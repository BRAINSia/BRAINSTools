import math

import SimpleITK as sitk


def read_fcsv(lmks_fn):
    """
      lmks_fn: A slicer complant fcsv fiducial file
      This function returns a map of named landmark points.
    """
    from collections import (
        OrderedDict,
    )  # Need OrderedDict internally to ensure consistent ordering

    lmks_map = OrderedDict()
    with open(lmks_fn, "r") as fid:
        lines = fid.readlines()
        for line in lines:
            if line[0] == "#":
                continue
            string_list = line.split(",")
            lmkName = string_list[0]
            lmk_pnt = [float(v) for v in (string_list[1:4])]
            lmk_pnt[0] *= -1  # R -> L
            lmk_pnt[1] *= -1  # A -> P
            lmk_pnt[2] *= +1  # P -> S
            lmks_map[lmkName] = lmk_pnt

    return lmks_map


def make_identity_image(img_size=320):
    """
    make an image with idenity space, and origin at center of space
    """
    origin_pnt = -(img_size - 1) / 2.0
    idimg = sitk.Image(img_size, img_size, img_size, sitk.sitkInt32)
    idimg.SetSpacing([1.0, 1.0, 1.0])
    idimg.SetOrigin([origin_pnt, origin_pnt, origin_pnt])
    return idimg


def zero_bnd_idx(pnt, img):
    """
    Get the index from a point, but use edge of image if out of bounds
    :param pnt: The pysical point
    :param img: The image with bounds
    :return:
    """
    idx = img.TransformPhysicalPointToIndex(pnt)
    idx = [max(0, x) for x in idx]
    last_idx = img.GetSize()
    idx = [min(last_idx[i] - 1, idx[i]) for i in (0, 1, 2)]
    return idx


def fix_index_ranges(lower_corner, upper_corner):
    """
    Some images have index increases opposite physical space increases
    :return:
    """
    lc = [min(lower_corner[i], upper_corner[i]) for i in (0, 1, 2)]
    uc = [max(lower_corner[i], upper_corner[i]) for i in (0, 1, 2)]
    return lc, uc


def get_eye_image(one_eye_pnt, orig_img_any, fov):
    orig_img = sitk.Resample(orig_img_any)
    re_lower_idx = zero_bnd_idx([ii - fov for ii in one_eye_pnt], orig_img)
    re_upper_idx = zero_bnd_idx([ii + fov for ii in one_eye_pnt], orig_img)
    re_center_idx = orig_img.TransformPhysicalPointToIndex(one_eye_pnt)
    LOWER_CORNER = (re_center_idx[0], re_lower_idx[1], re_lower_idx[2])
    UPPER_CORNER = (re_center_idx[0], re_upper_idx[1], re_upper_idx[2])
    LOWER_CORNER, UPPER_CORNER = fix_index_ranges(LOWER_CORNER, UPPER_CORNER)
    re_2dimg = orig_img[
        LOWER_CORNER[0],
        LOWER_CORNER[1] : UPPER_CORNER[1],
        LOWER_CORNER[2] : UPPER_CORNER[2],
    ]
    ## Ensure that the image range of the image fits into 16 bits
    ThirtyTwoBitImage = sitk.RescaleIntensity(
        sitk.Cast(re_2dimg, sitk.sitkInt32), 0, (2 ** 16 - 1)
    )
    re_2dimg = sitk.Cast(ThirtyTwoBitImage, sitk.sitkUInt16)  # Fit into PNG img
    re_2dimg = sitk.Flip(re_2dimg, [False, True])  ## Flip for easy viewing
    return re_2dimg


def make_mask_from_bb_corners(LOWER_CORNER_pnt, UPPER_CORNER_pnt, ref_img):
    """
    Given two physical points, create a mask that is the bounding box between those two points

    """
    IDTXFM = sitk.Transform()

    UPPER_index_tmp = zero_bnd_idx(UPPER_CORNER_pnt, ref_img)
    LOWER_index_tmp = zero_bnd_idx(LOWER_CORNER_pnt, ref_img)

    def arrangeUpperLower(idx1, idx2):
        lower = [min(idx1[i], idx2[i]) for i in (0, 1, 2)]
        upper = [max(idx1[i], idx2[i]) for i in (0, 1, 2)]
        return lower, upper

    LOWER_index, UPPER_index = arrangeUpperLower(LOWER_index_tmp, UPPER_index_tmp)

    extractregion = ref_img[
        LOWER_index[0] : UPPER_index[0],
        LOWER_index[1] : UPPER_index[1],
        LOWER_index[2] : UPPER_index[2],
    ]
    extractregion = (extractregion * 0) + 1
    out_mask = sitk.Cast(sitk.Resample(extractregion, ref_img, IDTXFM), sitk.sitkInt32)
    return out_mask


def quick_dilate(img, count, iterdilate):
    """
    iteratively calling smaller dilations is faster than doing large dilations
    """
    outlbl = sitk.Image(img)
    for it in range(0, count):
        outlbl = sitk.DilateObjectMorphology((outlbl > 0), iterdilate)
    return outlbl


def draw_eye(pnt, img):
    """
    Eye diameters are less than 30mm (typically average about 24mm)
    """
    index = img.TransformPhysicalPointToIndex(pnt)
    spacing = img.GetSpacing()
    EYE_MAX_RADIUS = 15
    indexstart = [int(math.floor(-EYE_MAX_RADIUS * e)) for e in spacing]
    indexstop = [int(math.ceil(+EYE_MAX_RADIUS * e)) for e in spacing]
    for x in range(indexstart[0], indexstop[0]):
        for y in range(indexstart[1], indexstop[1]):
            for z in range(indexstart[2], indexstop[2]):
                if x ** 2 + y ** 2 + z ** 2 < EYE_MAX_RADIUS ** 2:
                    offset = (x, y, z)
                    setidx = [index[i] + offset[i] for i in (0, 1, 2)]
                    img[setidx] = 1
    return img
