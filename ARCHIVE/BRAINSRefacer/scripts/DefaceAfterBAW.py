import os

##https://www.hindawi.com/journals/joph/2014/503645/tab1/
GLB_EYE_DIAMETER = 19.0  ## CHANGE TO EYE_RADIUS with extra margins
GLB_MAX_SIZE = 3000

from support import *


class DefaceAfterBAW:
    """
    A class to obscure recognizable facial features from a subject by using
    BAW output information.  The FCSV landmarks file, brain mask segmentations
    are both used to identify regions of the image that must be preseved.
    """

    def __init__(self, acpc_aligned_face_image):
        self.acpc_aligned_face_image_fn = acpc_aligned_face_image
        self.experiment_dir = os.path.dirname(
            os.path.dirname(self.acpc_aligned_face_image_fn)
        )
        out_base_fn = os.path.basename(self.acpc_aligned_face_image_fn)

        self.result_dir = os.path.join(self.experiment_dir, "DEFACE")
        ## -- DEBUG: self.result_dir = "/tmp"

        self.out_preserve_mask_fn = os.path.join(
            self.result_dir, out_base_fn.replace(".nii.gz", "_deface_mask.nii.gz")
        )
        self.out_deface_image_fn = os.path.join(
            self.result_dir, out_base_fn.replace(".nii.gz", "_deface.nii.gz")
        )
        if os.path.exists(self.out_deface_image_fn):
            print(("FILE EXISTS SO QUITTING: {0}".format(self.out_deface_image_fn)))
            sys.exit(-1)
        self.out_re_png_fn = self.out_deface_image_fn.replace(
            "_deface.nii.gz", "_right_eye.png"
        )
        self.out_le_png_fn = self.out_deface_image_fn.replace(
            "_deface.nii.gz", "_left_eye.png"
        )

        self.brain_label_fn = os.path.join(
            self.experiment_dir,
            "JointFusion/JointFusion_HDAtlas20_2015_lobe_label.nii.gz",
        )

        lmks_fn = os.path.join(self.experiment_dir, "ACPCAlign/BCD_ACPC_Landmarks.fcsv")
        self.lndmk_pts = read_fcsv(lmks_fn)
        del lmks_fn

        self.IDTXFM = sitk.Transform()
        self.reference_image = sitk.ReadImage(self.acpc_aligned_face_image_fn)
        self.reference_image_storage_type = self.reference_image.GetPixelID()

        _tmp = make_identity_image()
        self.idimg = sitk.Cast(
            sitk.Resample(self.reference_image, _tmp, self.IDTXFM, sitk.sitkLinear),
            self.reference_image_storage_type,
        )
        mask = sitk.Image(self.idimg) * 0.0
        del _tmp

        _tmp = sitk.Resample(
            sitk.ReadImage(self.brain_label_fn),
            self.idimg,
            self.IDTXFM,
            sitk.sitkNearestNeighbor,
        )
        self.head_mask = quick_dilate(_tmp, 4, 5)
        del _tmp

    def do_defacing(self):
        # UPPER_CORNER=(-MAX_SIZE,self.lndmk_pts["RE"][1]+EYE_DIAMETER, self.lndmk_pts["RE"][2]-EYE_DIAMETER)
        # LOWER_CORNER=(+MAX_SIZE,+MAX_SIZE,self.lndmk_pts["RE"][2]-MAX_SIZE)
        # out_mask = make_mask_from_bb_corners(LOWER_CORNER,UPPER_CORNER,self.idimg)

        UPPER_CORNER = (
            -GLB_MAX_SIZE,
            self.lndmk_pts["RE"][1] + GLB_EYE_DIAMETER,
            self.lndmk_pts["RE"][2] + 20.0,
        )
        # -- LOWER_CORNER=(+MAX_SIZE,+MAX_SIZE, self.lndmk_pts["RE"][2]+MAX_SIZE)
        LOWER_CORNER = (+GLB_MAX_SIZE, +GLB_MAX_SIZE, +GLB_MAX_SIZE)
        top_of_head_mask = make_mask_from_bb_corners(LOWER_CORNER, UPPER_CORNER, self.idimg)
        UPPER_CORNER = (
            -GLB_MAX_SIZE,
            self.lndmk_pts["AC"][1],
            self.lndmk_pts["AC"][2] - 55,
        )
        LOWER_CORNER = (
            +GLB_MAX_SIZE,
            self.lndmk_pts["RE"][1] + GLB_MAX_SIZE,
            +GLB_MAX_SIZE,
        )
        behind_acpnt = make_mask_from_bb_corners(LOWER_CORNER, UPPER_CORNER, self.idimg)

        left_eye = draw_eye(self.lndmk_pts["LE"], self.head_mask)
        right_eye = draw_eye(self.lndmk_pts["RE"], self.head_mask)

        force_keep_InIDIMG = (
            (left_eye > 0)
            + (right_eye > 0)
            + (top_of_head_mask > 0)
            + (behind_acpnt > 0)
            + (self.head_mask > 0)
        ) > 0

        force_keep = sitk.Resample(
            force_keep_InIDIMG, self.reference_image, self.IDTXFM
        )
        self.in_mask = sitk.Cast(force_keep, self.reference_image_storage_type)
        del force_keep, force_keep_InIDIMG

        not_in_mask = sitk.Cast((1 - self.in_mask), self.reference_image_storage_type)

        self.out_img = sitk.Image(self.reference_image)  # Copy image

        smoothable_image_pixel_type = sitk.Cast(self.reference_image, sitk.sitkInt32)
        build_up_img = sitk.Cast(
            sitk.SmoothingRecursiveGaussian(smoothable_image_pixel_type, 7),
            self.reference_image_storage_type,
        )
        del smoothable_image_pixel_type

        mif = sitk.MedianImageFilter()
        mif.SetRadius(7)
        smoothed_img = mif.Execute(build_up_img)
        smoothed_img = mif.Execute(smoothed_img)

        self.out_img = sitk.Cast(
            (self.in_mask * self.reference_image + not_in_mask * smoothed_img),
            self.reference_image_storage_type,
        )

        re_2dimg, le_2dimg = defacer.get_eye_images()
        sitk.WriteImage(re_2dimg, self.out_re_png_fn)
        sitk.WriteImage(le_2dimg, self.out_le_png_fn)

    def get_eye_images(self):
        id_space_out_img = sitk.Resample(self.out_img, self.idimg, self.IDTXFM)
        right_eye = get_eye_image(
            self.lndmk_pts["RE"], id_space_out_img, GLB_EYE_DIAMETER * 3
        )
        left_eye = get_eye_image(
            self.lndmk_pts["LE"], id_space_out_img, GLB_EYE_DIAMETER * 3
        )
        return right_eye, left_eye

    def write_outputs(self):
        sitk.WriteImage(self.in_mask, self.out_preserve_mask_fn)
        # sitk.WriteImage(self.idimg, os.path.join(self.result_dir, "resamp.nii.gz"))
        sitk.WriteImage(self.out_img, self.out_deface_image_fn)


if __name__ == "__main__":
    import sys

    print((sys.argv[1]))
    if len(sys.argv) != 2:
        print(
            (
                """USAGE: {0} <FullPathToUnDefacedImage.nii.gz>\n
        """.format(
                    sys.argv[0]
                )
            )
        )
        sys.exit(-1)
    print(("Defacing: {0}".format(sys.argv[1])))
    ref_img = sys.argv[1]
    defacer = DefaceAfterBAW(ref_img)
    defacer.do_defacing()
    defacer.write_outputs()

    print("done")
