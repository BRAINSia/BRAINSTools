from nipype.interfaces.base import BaseInterface, BaseInterfaceInputSpec, traits, TraitedSpec
from nipype.interfaces.freesurfer.base import FSCommand, FSTraitedSpec
import SimpleITK as sitk
import numpy as np
import pandas as pd
import os


def create_ones_image(in_volume, out_file, value=1):
    """Creates a volume filled with a scalar (1 by default)"""
    import nibabel as nb
    import os
    image = nb.load(in_volume)
    array = image.get_data()
    array.fill(value)
    image.to_filename(out_file)
    out_file = os.path.abspath(out_file)
    return out_file


def recode_labelmap(in_file, out_file, recode_file):
    """This function has been adapted from BRAINSTools and serves
    as a means to recode a label map based upon an input csv
    file."""
    # rewritten by Chris Markiewicz
    import sys
    import os
    import numpy as np
    import nibabel as nib
    import json

    if os.path.splitext(recode_file)[-1] == ".json":
        # Extract label-label map from a JSON file
        with open(recode_file, "rb") as json_file:
            abc_recode_data = json.load(json_file)
            recode_table = dict()
            for abc_label in abc_recode_data:
                fs_labels = abc_recode_data[abc_label]
                for fs_label in fs_labels:
                    recode_table[int(fs_label)] = int(abc_label)
    else:
        # Extract label-label map from a CSV file
        recode_data = np.loadtxt(recode_file, delimiter=',', skiprows=1,
                                 dtype='S20')
        # Permit (and ignore) label names
        if recode_data.shape[1] == 4:
            recode_data = recode_data[:, (0, 2)]
        if recode_data.shape[1] != 2:
            print("ERROR: input csv file for label recoding does meet "
                  "requirements")
            sys.exit()
        recode_table = dict(recode_data.astype(np.uint64))

    mapper = lambda x: recode_table[x] if x in recode_table else x

    img = nib.load(in_file)
    klass = img.__class__
    if ".nii" in out_file:
        klass = nib.Nifti1Image # allows for changing to NIFTI file type
    out_file = os.path.abspath(out_file)

    # Cast non-integer labels as unsigned, 32-bit integers
    dtype = img.get_data_dtype()
    if not np.issubdtype(img.get_data_dtype(), np.integer):
        dtype = np.uint32
    labels = np.asanyarray(img.get_data(), dtype=dtype)

    # Choose smallest integer type to contain all outputted values
    recode = np.vectorize(mapper, otypes=[np.uint64])
    max_val = recode(np.unique(labels)).max()
    for dtype in (np.uint8, np.uint16, np.uint32, np.uint64):
        if max_val <= np.iinfo(dtype).max:
            break
    recode = np.vectorize(mapper, otypes=[dtype])

    new_img = klass(recode(labels), img.affine, img.header, extra=img.extra)
    new_img.set_data_dtype(dtype)  # Output type defined in header
    new_img.to_filename(out_file)

    return out_file


# As advised on the ITK mailing list, label dilation can be implemented via
# distance transforms and watershed transforms. This algorithm is illustrated
# in SimpleITK python code below (courtesy of Bradely Lowekamp)
def multilabel_dilation(in_file, out_file, radius=1, kernel=None):
    import SimpleITK as sitk
    import os
    img = sitk.ReadImage(in_file)
    if not kernel:
        kernel = sitk.BinaryDilateImageFilter.Ball
    dilatImg = sitk.BinaryDilate(img != 0, radius, kernel)
    wsImg = create_label_watershed(img)
    sitk.WriteImage(sitk.Cast(dilatImg, wsImg.GetPixelID())*wsImg, out_file)
    out_file = os.path.abspath(out_file)
    return out_file


def create_label_watershed(labels_image, markWatershedLine=False):
    import SimpleITK as sitk
    distImg = sitk.SignedMaurerDistanceMap(labels_image != 0,
                                           insideIsPositive=False,
                                           squaredDistance=False,
                                           useImageSpacing=False)
    wsImg = sitk.MorphologicalWatershedFromMarkers(distImg, labels_image, markWatershedLine=markWatershedLine)
    return wsImg


class MultiLabelDilationInputSpec(BaseInterfaceInputSpec):
    in_file = traits.File(exists=True, mandatory=True)
    out_file = traits.File(mandatory=True)
    radius = traits.Int(1, use_default=True)


class MultiLabelDilationOutputSpec(TraitedSpec):
    out_file = traits.File()


class MultiLabelDilation(BaseInterface):
    input_spec = MultiLabelDilationInputSpec
    output_spec = MultiLabelDilationOutputSpec

    def _run_interface(self, runtime):
        self.output_spec.out_file = multilabel_dilation(self.inputs.in_file, self.inputs.out_file, self.inputs.radius)
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = self.output_spec.out_file
        return outputs


def create_image_like(array, image):
    """Takes an array and creates it into an image like the one given"""
    import SimpleITK as sitk
    image_array = sitk.GetArrayFromImage(image)
    ndims = len(array.shape)
    if ndims == 1:
        array = array.reshape(image_array.shape)
    new_image = sitk.GetImageFromArray(array)
    new_image.SetOrigin(image.GetOrigin())
    new_image.SetSpacing(image.GetSpacing())
    new_image.SetDirection(image.GetDirection())
    return new_image


def split_labels(labels_file, lut_file, out_file, left_label=1, right_label=2):
    """create a a hemisphere label map"""
    # read in the LUT
    df = pd.read_csv(lut_file, index_col=0)
    left_labels = df.index[df.Left == 1].values.tolist()
    right_labels = df.index[df.Right == 1].values.tolist()
    # Load the image
    image = sitk.ReadImage(labels_file)
    image_array = sitk.GetArrayFromImage(image)
    # create empty arrays
    left_array = np.zeros_like(image_array)
    right_array = np.zeros_like(image_array)
    labels = np.unique(image_array)
    for label in labels:
        if label in left_labels:
            index = image_array == label
            left_array[index] = 1
        elif label in right_labels:
            index = image_array == label
            right_array[index] = 1

    left_image = create_image_like(left_array, image)
    clean_left = sitk.BinaryMorphologicalClosing(left_image, 1)
    right_image = create_image_like(right_array, image)
    clean_right = sitk.BinaryMorphologicalClosing(right_image, 1)

    hemi_labels = (clean_left * (clean_right == 0)) * left_label + (clean_right * (clean_left == 0)) * right_label

    wsImg = create_label_watershed(hemi_labels)
    sitk.WriteImage(wsImg, out_file)
    out_file = os.path.abspath(out_file)
    return out_file


def apply_label_split(image_file, hemi_file, hemi, out_file, left_label=1, right_label=2):
    import SimpleITK as sitk
    import os
    image = sitk.ReadImage(image_file)
    hemis = sitk.ReadImage(hemi_file, image.GetPixelID())
    if hemi == 'lh':
        out_image = (hemis == left_label) * image
    elif hemi == 'rh':
        out_image = (hemis == right_label) * image
    else:
        print("ERROR: Invalid hemisphere name {0}".format(hemi))
        return
    sitk.WriteImage(out_image, out_file)
    out_file = os.path.abspath(out_file)
    return out_file


class SplitLabelsInputSpec(BaseInterfaceInputSpec):
    in_file = traits.File(exists=True, mandatory=True)
    labels_file = traits.File(exists=True, mandatory=True)
    lookup_table = traits.File(exists=True, mandatory=True)
    hemi = traits.Enum('lh', 'rh', mandatory=True)
    out_file = traits.File(mandatory=True)


class SplitLabelsOutputSpec(TraitedSpec):
    out_file = traits.File(exists=True)


class SplitLabels(BaseInterface):
    input_spec = SplitLabelsInputSpec
    output_spec = SplitLabelsOutputSpec

    def _run_interface(self, runtime):
        hemispheres_image_file = split_labels(self.inputs.labels_file, self.inputs.lookup_table, "hemispheres.nii.gz")
        self.output_spec.out_file = os.path.abspath(apply_label_split(self.inputs.in_file, hemispheres_image_file,
                                                                      self.inputs.hemi, self.inputs.out_file))
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = self.output_spec.out_file
        return outputs


class SurfaceMaskInputSpec(FSTraitedSpec):
    in_volume = traits.File(argstr="%s", position=-3, exist=True,
                            desc="Input volume to which mask is applied.")
    in_surface = traits.File(argstr="%s", position=-2, exist=True,
                             desc="Input surface. Values inside surface are set to the values of in_volume.")
    out_file = traits.File(argstr="%s", position=-1, exist=True,
                           desc="Output masked volume.")


class SurfaceMaskOutputSpec(TraitedSpec):
    out_file = traits.File(desc="Output masked volume.")


class SurfaceMask(FSCommand):
    """Purpose: Produce a new volume where all pixels outside the surface are set to zero.
    """

    _cmd = 'mri_surfacemask'
    input_spec = SurfaceMaskInputSpec
    output_spec = SurfaceMaskOutputSpec

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs["out_file"] = os.path.abspath(self.inputs.out_file)
        return outputs
