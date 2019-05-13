"""
BRAINSABCext.py
===========================
Description:

Author:

Usage:

"""
import os
from builtins import object
from builtins import range
from xml.etree import ElementTree as et

from nipype.interfaces.base import (
    CommandLine,
    CommandLineInputSpec,
    TraitedSpec,
    File,
    Directory,
    traits,
    isdefined,
    InputMultiPath,
    OutputMultiPath,
)
from nipype.interfaces.semtools.segmentation.specialized import (
    BRAINSABCOutputSpec,
    BRAINSABCInputSpec,
    BRAINSABC,
)


class GetPosteriorsFromAtlasXML(object):
    """
    This class represents a...
    """

    def __init__(self, xmlFile):
        self.xmlFile = xmlFile
        self.xmlString = self.get_xml_string(self.xmlFile)
        self.priorTypeNameList = self.get_prior_type_name_list(self.xmlString)

    def main(self):
        self.get_posterior_file_name_list(priorTypeNameList)

    def get_xml_string(self, xmlFile):
        """
        This function...

        :param xmlFile:
        :return:
        """
        Handle = open(xmlFile)
        _xmlString = Handle.read()
        Handle.close()
        return _xmlString

    def get_prior_type_name_list(self, xmlString):
        """
        This function...

        :param xmlString:
        :return:
        """
        myelem = et.fromstring(xmlString)
        elementsList = list(myelem.getiterator())
        iterator = list(range(len(elementsList)))
        priorTypeNameList = list()
        for i in iterator:
            ## the Prior type is the next item in elementsList after a Prior tag:
            if elementsList[i].tag == "Prior" and elementsList[i + 1].tag == "type":
                priorTypeNameList.append(elementsList[i + 1].text)
        return priorTypeNameList

    def get_posterior_file_name_list(self, posteriorTemplate):
        posteriorFileNameList = list()
        for priorType in self.priorTypeNameList:
            posteriorFileNameList.append(
                "POSTERIOR_{priorT}.nii.gz".format(priorT=priorType)
            )
            ## HACK:  The following is correct from the command line posteriorTemplate arguments
            # posteriorFileNameList.append(posteriorTemplate % priorType)
        return posteriorFileNameList


"""
class BRAINSABCextInputSpec(BRAINSABCInputSpec):
    outputAverageImages = traits.Either(traits.Bool(True,desc="The automatically generated average images"), InputMultiPath(File(),), hash_files = False,argstr = "")
    posteriorImages = traits.Either(traits.Bool(True,desc="The automatically generated posterior images"), InputMultiPath(File(),), hash_files = False,argstr = "")
"""


class BRAINSABCextOutputSpec(BRAINSABCOutputSpec):
    """
    This class represents a...
    """

    # Not convenient outputAverageImages = OutputMultiPath(File(exists=True), exists = True)
    outputT1AverageImage = traits.Either(File(exists=True), [None])
    outputT2AverageImage = traits.Either(File(exists=True), [None])
    outputPDAverageImage = traits.Either(File(exists=True), [None])
    outputFLAverageImage = traits.Either(File(exists=True), [None])
    posteriorImages = OutputMultiPath(File(exists=True), exists=True)
    atlasToSubjectInverseTransform = traits.Either(File(exists=True), [None])


class BRAINSABCext(BRAINSABC):
    """
    This class represents a...
    """

    # input_spec= BRAINSABCextInputSpec
    output_spec = BRAINSABCextOutputSpec

    def _list_outputs(self):
        """
        This function...

        :return:
        """
        from collections import (
            OrderedDict,
        )  # Need OrderedDict internally to ensure consistent ordering

        custom_implied_outputs_with_no_inputs = [
            "posteriorImages",
            "outputT1AverageImage",
            "outputT2AverageImage",
            "outputPDAverageImage",
            "outputFLAverageImage",
            "atlasToSubjectInverseTransform",
        ]
        full_outputs = self.output_spec().get()
        pruned_outputs = OrderedDict()
        for key, value in list(full_outputs.items()):
            if key not in custom_implied_outputs_with_no_inputs:
                pruned_outputs[key] = value
        outputs = super(BRAINSABCext, self)._outputs_from_inputs(pruned_outputs)
        input_check = OrderedDict(
            {
                "T1": ("outputT1AverageImage", "t1_average_BRAINSABC.nii.gz"),
                "T2": ("outputT2AverageImage", "t2_average_BRAINSABC.nii.gz"),
                "PD": ("outputPDAverageImage", "pd_average_BRAINSABC.nii.gz"),
                "FL": ("outputFLAverageImage", "fl_average_BRAINSABC.nii.gz"),
            }
        )
        for key, values in list(input_check.items()):
            if key in self.inputs.inputVolumeTypes:
                outputs[values[0]] = os.path.abspath(values[1])
            else:
                outputs[values[0]] = None

        PosteriorOutputs = GetPosteriorsFromAtlasXML(self.inputs.atlasDefinition)
        PosteriorPaths = PosteriorOutputs.get_posterior_file_name_list(
            self.inputs.posteriorTemplate
        )
        outputs["posteriorImages"] = [
            os.path.abspath(postPath) for postPath in PosteriorPaths
        ]

        fixed_inverse_name = os.path.abspath(
            outputs["atlasToSubjectTransform"].replace(".h5", "_Inverse.h5")
        )
        if os.path.exists(fixed_inverse_name):
            outputs["atlasToSubjectInverseTransform"] = fixed_inverse_name
        else:
            outputs["atlasToSubjectInverseTransform"] = None
        return outputs
