#! /usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Date:          2012-09-05
Author:        hans-johnson@uiowa.edu
Purpose:       Wrap a convenience function for the BRAINSCut program in Nipype

Requirements:  <<< Interface specifications >>>

"""
from nipype.interfaces.base import ( File, TraitedSpec, Interface, CommandLineInputSpec, CommandLine, traits, isdefined)
import sys
import os
import warnings

### CommandLine
class RF12BRAINSCutWrapperCLInputSpec(CommandLineInputSpec):
    ### subject specific
    inputSubjectT1Filename = File( desc="Subject T1 Volume", exists=True, mandatory=True, argstr="--inputSubjectT1Filename %s")
    inputSubjectT2Filename = File( desc="Subject T2 Volume", exists=True, mandatory=True, argstr="--inputSubjectT2Filename %s")
    inputSubjectGadSGFilename = File( desc="Subject SG Volume", exists=True, mandatory=True, argstr="--inputSubjectGadSGFilename %s")
    vectorNormalization = traits.Enum("IQR","Linear","Sigmoid_Q01","Sigmoid_Q05","ZScore","NONE",
        desc="The type of intensity normalization to use",exists=True,mandatory=True,argstr="--vectorNormalization %s")

    ### model specific
    modelFilename = File( desc="modelFilename", exists=True, mandatory=True, argstr="--modelFilename %s")
    trainingVectorFilename = File( desc="training vectof file name", exists=False, mandatory=False, argstr="--trainingVectorFilename %s")
    inputTemplateT1 = File( desc="Atlas Template T1 image", exists=False, mandatory=False, argstr="--inputTemplateT1 %s")
    inputTemplateRhoFilename = File( desc="Atlas Template rho image", exists=False, mandatory=False, argstr="--inputTemplateRhoFilename %s")
    inputTemplatePhiFilename = File( desc="Atlas Template phi image", exists=False, mandatory=False, argstr="--inputTemplatePhiFilename %s")
    inputTemplateThetaFilename = File( desc="Atlas Template theta image", exists=False, mandatory=False, argstr="--inputTemplateThetaFilename %s")
    deformationFromTemplateToSubject = File( desc="Atlas To subject Deformation", exists=False, mandatory=False, argstr="--deformationFromTemplateToSubject %s")

    ### probability maps
    probabilityMapsLeftAccumben = File( desc="Spatial probability map of left accumben", exists=True, mandatory=True, argstr="--probabilityMapsLeftAccumben %s")
    probabilityMapsRightAccumben = File( desc="Spatial probability map of right accumben", exists=True, mandatory=True, argstr="--probabilityMapsRightAccumben %s")

    probabilityMapsLeftCaudate = File( desc="Spatial probability map of left caudate", exists=True, mandatory=True, argstr="--probabilityMapsLeftCaudate %s")
    probabilityMapsRightCaudate = File( desc="Spatial probability map of right caudate", exists=True, mandatory=True, argstr="--probabilityMapsRightCaudate %s")

    probabilityMapsLeftGlobus = File( desc="Spatial probability map of left globus", exists=True, mandatory=True, argstr="--probabilityMapsLeftGlobus %s")
    probabilityMapsRightGlobus = File( desc="Spatial probability map of right globus", exists=True, mandatory=True, argstr="--probabilityMapsRightGlobus %s")

    probabilityMapsLeftHippocampus = File( desc="Spatial probability map of left hippocampus", exists=True, mandatory=True, argstr="--probabilityMapsLeftHippocampus %s")
    probabilityMapsRightHippocampus = File( desc="Spatial probability map of right hippocampus", exists=True, mandatory=True, argstr="--probabilityMapsRightHippocampus %s")

    probabilityMapsLeftPutamen = File( desc="Spatial probability map of left putamen", exists=True, mandatory=True, argstr="--probabilityMapsLeftPutamen %s")
    probabilityMapsRightPutamen = File( desc="Spatial probability map of right putamen", exists=True, mandatory=True, argstr="--probabilityMapsRightPutamen %s")

    probabilityMapsLeftThalamus = File( desc="Spatial probability map of left thalamus", exists=True, mandatory=True, argstr="--probabilityMapsLeftThalamus %s")
    probabilityMapsRightThalamus = File( desc="Spatial probability map of right thalamus", exists=True, mandatory=True, argstr="--probabilityMapsRightThalamus %s")

    xmlFilename = File( desc = "Net configuration xml file", exists = False, mandatory = False, argstr="--xmlFilename %s")

    outputBinaryLeftAccumben = File( desc = "Output binary file of left accumben", exists = False, mandatory = True, argstr="--outputBinaryLeftAccumben %s")
    outputBinaryRightAccumben = File( desc = "Output binary file of right accumben", exists = False, mandatory = True, argstr="--outputBinaryRightAccumben %s")

    outputBinaryLeftCaudate = File( desc = "Output binary file of left caudate", exists = False, mandatory = True, argstr="--outputBinaryLeftCaudate %s")
    outputBinaryRightCaudate = File( desc = "Output binary file of right caudate", exists = False, mandatory = True, argstr="--outputBinaryRightCaudate %s")

    outputBinaryLeftGlobus = File( desc = "Output binary file of left globus", exists = False, mandatory = True, argstr="--outputBinaryLeftGlobus %s")
    outputBinaryRightGlobus = File( desc = "Output binary file of right globus", exists = False, mandatory = True, argstr="--outputBinaryRightGlobus %s")

    outputBinaryLeftHippocampus = File( desc = "Output binary file of left hippocampus", exists = False, mandatory = True, argstr="--outputBinaryLeftHippocampus %s")
    outputBinaryRightHippocampus = File( desc = "Output binary file of right hippocampus", exists = False, mandatory = True, argstr="--outputBinaryRightHippocampus %s")

    outputBinaryLeftPutamen = File( desc = "Output binary file of left putamen", exists = False, mandatory = True, argstr="--outputBinaryLeftPutamen %s")
    outputBinaryRightPutamen = File( desc = "Output binary file of right putamen", exists = False, mandatory = True, argstr="--outputBinaryRightPutamen %s")

    outputBinaryLeftThalamus = File( desc = "Output binary file of left thalamus", exists = False, mandatory = True, argstr="--outputBinaryLeftThalamus %s")
    outputBinaryRightThalamus = File( desc = "Output binary file of right thalamus", exists = False, mandatory = True, argstr="--outputBinaryRightThalamus %s")

class RF12BRAINSCutWrapperCLOutputSpec(TraitedSpec):
    xmlFilename = File( desc = "Net configuration xml file", exists = False, mandatory = False)

    outputBinaryLeftAccumben = File( desc = "Output binary file of left accumben", exists = True, mandatory = True)
    outputBinaryRightAccumben = File( desc = "Output binary file of right accumben", exists = True, mandatory = True)

    outputBinaryLeftCaudate = File( desc = "Output binary file of left caudate", exists = True, mandatory = True)
    outputBinaryRightCaudate = File( desc = "Output binary file of right caudate", exists = True, mandatory = True)

    outputBinaryLeftGlobus = File( desc = "Output binary file of left globus", exists = True, mandatory = True)
    outputBinaryRightGlobus = File( desc = "Output binary file of right globus", exists = True, mandatory = True)

    outputBinaryLeftHippocampus = File( desc = "Output binary file of left hippocampus", exists = True, mandatory = True)
    outputBinaryRightHippocampus = File( desc = "Output binary file of right hippocampus", exists = True, mandatory = True)

    outputBinaryLeftPutamen = File( desc = "Output binary file of left putamen", exists = True, mandatory = True)
    outputBinaryRightPutamen = File( desc = "Output binary file of right putamen", exists = True, mandatory = True)

    outputBinaryLeftThalamus = File( desc = "Output binary file of left thalamus", exists = True, mandatory = True)
    outputBinaryRightThalamus = File( desc = "Output binary file of right thalamus", exists = True, mandatory = True)

class RF12BRAINSCutWrapper(CommandLine):
    """
    A script to wrap the complexity of BRAINSCut into a single script.
    """
    _cmd = 'BRAINSCutCMD.py'
    input_spec = RF12BRAINSCutWrapperCLInputSpec
    output_spec = RF12BRAINSCutWrapperCLOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        return self._outputs_from_inputs(outputs)

    def _outputs_from_inputs(self, outputs):
        for name in outputs.keys():
            coresponding_input = getattr(self.inputs, name)
            if isdefined(coresponding_input):
                if isinstance(coresponding_input, bool) and coresponding_input == True:
                    outputs[name] = os.path.abspath(self._outputs_filenames[name])
                else:
                    if isinstance(coresponding_input, list):
                        outputs[name] = [os.path.abspath(inp) for inp in coresponding_input]
                    else:
                        outputs[name] = os.path.abspath(coresponding_input)
        return outputs

#if __name__ == '__main__':
#    RF12Test = RF12BRAINSCutWrapper(sys.argv)
#    RF12Test.run()
