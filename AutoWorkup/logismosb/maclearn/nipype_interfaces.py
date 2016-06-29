from nipype.interfaces.base import BaseInterface, traits, BaseInterfaceInputSpec, TraitedSpec
from nipype.interfaces.semtools import BRAINSResample
from nipype.interfaces.freesurfer import MRIConvert
import os
from training import image_data
from sklearn.externals import joblib
from predict import image_file_from_array_with_reference_image_file
import SimpleITK as sitk


def run_resample(in_file, ref_file, transform, out_file, interpolation_mode='Linear', pixel_type='float',
                 inverse_transform=True):
    resample = BRAINSResample()
    resample.inputs.inputVolume = in_file
    resample.inputs.warpTransform = transform
    resample.inputs.pixelType = pixel_type
    resample.inputs.interpolationMode = interpolation_mode
    resample.inputs.outputVolume = os.path.abspath(out_file)
    resample.inputs.referenceVolume = ref_file
    resample.inputs.inverseTransform = inverse_transform
    result = resample.run()
    return result.outputs.outputVolume


class CollectFeatureFilesInputSpec(BaseInterfaceInputSpec):
    rho = traits.File(exists=True)
    phi = traits.File(exists=True)
    theta = traits.File(exists=True)
    posterior_files = traits.Dict(trait=traits.File(exists=True), desc="Dictionary of posterior probabilities")
    reference_file = traits.File(exists=True)
    transform_file = traits.File(exists=True)
    inverse_transform = traits.Bool(True, desc="if true, inverse transform will be used", use_default=True)


class CollectFeatureFilesOutputSpec(TraitedSpec):
    feature_files = traits.Dict(trait=traits.File(exists=True), desc="Output dictionary of feature files")


class CollectFeatureFiles(BaseInterface):
    input_spec = CollectFeatureFilesInputSpec
    output_spec = CollectFeatureFilesOutputSpec

    def _run_interface(self, runtime):
        list_of_feature_files = self.combine_files()
        self.resample_images_for_features(list_of_feature_files, self.inputs.reference_file, self.inputs.transform_file)
        return runtime

    @staticmethod
    def _get_dict_name(filename):
        return os.path.basename(filename).split(".")[0]

    def combine_files(self):
        feature_files = [self.inputs.rho, self.inputs.phi, self.inputs.theta] + self.inputs.posterior_files.values()
        return feature_files

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs["feature_files"] = self._list_resampled_feature_files(self.combine_files())
        return outputs

    def _list_resampled_feature_files(self, list_of_feature_files):
        resampled_feature_files_dict = dict()
        for _file in list_of_feature_files:
            basename = os.path.basename(_file)
            name = self._get_dict_name(_file)
            resampled_feature_files_dict[name] = os.path.abspath(basename)
        return resampled_feature_files_dict

    def resample_images_for_features(self, list_of_feature_files, ref_file, transform):
        resampled_feature_files = self._list_resampled_feature_files(list_of_feature_files)
        for _file in list_of_feature_files:
            name = self._get_dict_name(_file)
            resampled_feature_files[name] = run_resample(_file, ref_file, transform, resampled_feature_files[name],
                                                         "Linear", "float",
                                                         inverse_transform=self.inputs.inverse_transform)


def create_white_edge_cost_image(t1_file, t2_file, gm_proba_file, out_file):
    import SimpleITK as sitk
    import os
    gm_proba = sitk.ReadImage(gm_proba_file)
    negative_gm_proba = 1 - gm_proba
    t1 = sitk.ReadImage(t1_file)
    t2 = sitk.ReadImage(t2_file)
    t1_gradient = sitk.GradientMagnitude(t1)
    t2_gradient = sitk.GradientMagnitude(t2)
    multi_modal_gradient = sitk.Cast((t1_gradient + t2_gradient), negative_gm_proba.GetPixelID())
    cost_image = multi_modal_gradient * negative_gm_proba
    out_file = os.path.abspath(out_file)
    sitk.WriteImage(cost_image, out_file)
    return out_file


class PredictEdgeProbabilityInputSpec(BaseInterfaceInputSpec):
    t1_file = traits.File(exists=True)
    additional_files = traits.Dict(trait=traits.File(exists=True))
    gm_classifier_file = traits.File(exists=True, desc="Classifier file for predicting edge probability")
    wm_classifier_file = traits.File(exists=True, desc="Classifier file for predicting edge probability")


class PredictEdgeProbabilityOutputSpec(TraitedSpec):
    gm_edge_probability = traits.File()
    wm_edge_probability = traits.File()


class PredictEdgeProbability(BaseInterface):
    input_spec = PredictEdgeProbabilityInputSpec
    output_spec = PredictEdgeProbabilityOutputSpec

    def _run_interface(self, runtime):
        feature_data = image_data(self.inputs.t1_file, "T1", additional_images=self.inputs.additional_files)
        gm_classifier = joblib.load(self.inputs.gm_classifier_file)
        gm_probability_array = gm_classifier.predict_proba(feature_data.values)[:, 1]
        del gm_classifier
        gm_probability_image = image_file_from_array_with_reference_image_file(
            gm_probability_array, self.inputs.t1_file, self._list_outputs()["gm_edge_probability"])
        del gm_probability_array, gm_probability_image
        wm_classifier = joblib.load(self.inputs.wm_classifier_file)
        wm_probability_array = wm_classifier.predict_proba(feature_data.values)[:, 1]
        del wm_classifier
        wm_probability_image = image_file_from_array_with_reference_image_file(
            wm_probability_array, self.inputs.t1_file, self._list_outputs()["wm_edge_probability"])
        del wm_probability_array, wm_probability_image
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        for matter in ["gm", "wm"]:
            name = "{0}_edge_probability".format(matter)
            outputs[name] = os.path.abspath(name + ".nii.gz")
        return outputs


def create_identity_transform():
    dimension = 3
    offset = (0, 0, 0)
    transform = sitk.TranslationTransform(dimension , offset)
    transform.SetIdentity() # just to be safe
    return transform


def change_orientation(image_file, out_file, orientation="LPS"):
    convert = MRIConvert()
    convert.inputs.in_file = image_file
    convert.inputs.out_file = os.path.abspath(out_file)
    convert.inputs.out_orientation = orientation
    result = convert.run()
    return result.outputs.out_file


def sample_image(image, size, spacing=(1, 1, 1)):
    resample = sitk.ResampleImageFilter()
    resample.SetInterpolator(sitk.sitkLinear)
    resample.SetOutputSpacing(spacing)
    resample.SetTransform(create_identity_transform())
    resample.SetOutputDirection(image.GetDirection())
    resample.SetSize(size)
    resample.SetOutputOrigin(image.GetOrigin())
    return resample.Execute(image)


class CreateReferenceImageInputSpec(BaseInterfaceInputSpec):
    baw_t1 = traits.File(exists=True)
    orig_t1 = traits.File(exists=True)


class CreateReferenceImageOutputSpec(TraitedSpec):
    reference_file = traits.File(desc="reference file for resampling")


class CreateReferenceImage(BaseInterface):
    input_spec = CreateReferenceImageInputSpec
    output_spec = CreateReferenceImageOutputSpec

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs["reference_file"] = os.path.abspath("reference_file.nii.gz")
        return outputs

    def _run_interface(self, runtime):
        baw_image = sitk.ReadImage(str(self.inputs.baw_t1))
        spacing = baw_image.GetSpacing()
        size = baw_image.GetSize()
        output_image = sample_image(sitk.ReadImage(str(self.inputs.orig_t1)), size, spacing)
        out_file = self._list_outputs()["reference_file"]
        sitk.WriteImage(output_image, out_file)
        change_orientation(out_file, out_file)
        return runtime


def scale_image(image, scale):
    return image * scale


def increase_penalty(image, penalty, minimum):
    background = sitk.Cast(image < minimum, sitk.sitkFloat64)
    return image - (background * penalty)


class LOGISMOSBPreprocessingInputSpec(BaseInterfaceInputSpec):
    white_mask = traits.File(exists=True, mandatory=True)
    erode_mask = traits.Int(default_value=1, usedefault=True)
    gm_proba = traits.File(exists=True, mandatory=True)
    wm_proba = traits.File(exists=True, mandatory=True)
    background_penalty = traits.Int(default_value=100, usedefault=True, desc="Penalty for a zero probability.")
    proba_scale = traits.Int(default_value=50, usedefault=True, desc="Scale the probabilities.")
    min_probability = traits.Float(default_value=0.01, usedefault=True)


class LOGISMOSBPreprocessingOutputSpec(TraitedSpec):
    wm_proba = traits.File()
    gm_proba = traits.File()
    white_mask = traits.File()


class LOGISMOSBPreprocessing(BaseInterface):
    """Interface for playing with the inputs so that LOGISMOS-B is optimized for the probability maps."""
    input_spec = LOGISMOSBPreprocessingInputSpec
    output_spec = LOGISMOSBPreprocessingOutputSpec

    def _list_outputs(self):
        outputs = self._outputs().get()
        inputs = self.inputs.get()
        for output in outputs:
            outputs[output] = os.path.abspath(os.path.basename(inputs[output]))
        return outputs

    def _run_interface(self, runtime):
        inputs = self.inputs.get()
        for matter in ("gm", "wm"):
            proba_name = "{0}_proba".format(matter)
            proba_file = inputs[proba_name]
            proba = sitk.ReadImage(str(proba_file))
            if self.inputs.proba_scale:
                proba = scale_image(proba, self.inputs.proba_scale)
            if self.inputs.background_penalty:
                proba = increase_penalty(proba, self.inputs.background_penalty, self.inputs.min_probability)
            proba_out_file = self._list_outputs()[proba_name]
            sitk.WriteImage(proba, proba_out_file)
        white_mask = sitk.ReadImage(str(self.inputs.white_mask))
        white_out_file = self._list_outputs()["white_mask"]
        if self.inputs.erode_mask:
            eroded_white_mask = sitk.BinaryErode(sitk.Cast(white_mask, sitk.sitkUInt8), int(self.inputs.erode_mask))
            sitk.WriteImage(eroded_white_mask, white_out_file)
        else:
            sitk.WriteImage(white_mask, white_out_file)
        return runtime
