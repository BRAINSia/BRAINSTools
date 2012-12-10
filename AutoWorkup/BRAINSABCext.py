from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory, traits, isdefined, InputMultiPath, OutputMultiPath
import os
from SEMTools.segmentation.specialized import  BRAINSABCOutputSpec, BRAINSABCInputSpec, BRAINSABC
#from SEMTools import BRAINSABCInputSpec,BRAINSABCOutputSpec,BRAINSABC

from xml.etree import ElementTree as et
class GetPosteriorsFromAtlasXML():

    def __init__(self, xmlFile):
        self.xmlFile = xmlFile
        self.xmlString = self.getXMLstring(self.xmlFile)
        self.priorTypeNameList = self.getPriorTypeNameList(self.xmlString)

    def main(self):
        self.getPosteriorFileNameList(priorTypeNameList)

    def getXMLstring(self, xmlFile):
        Handle = open(xmlFile)
        _xmlString = Handle.read()
        Handle.close()
        return _xmlString

    def getPriorTypeNameList(self, xmlString):
        myelem = et.fromstring(xmlString)
        elementsList = myelem.getiterator()
        iterator = range(len(elementsList))
        priorTypeNameList = list()
        for i in iterator:
            ## the Prior type is the next item in elementsList after a Prior tag:
            if elementsList[i].tag == "Prior" and elementsList[i + 1].tag == "type":
                priorTypeNameList.append(elementsList[i + 1].text)
        return priorTypeNameList

    def getPosteriorFileNameList(self, posteriorTemplate):
        posteriorFileNameList = list()
        for priorType in self.priorTypeNameList:
            posteriorFileNameList.append("POSTERIOR_{priorT}.nii.gz".format(priorT=priorType))
            ## HACK:  The following is correct from the command line posteriorTemplate arguments
            #posteriorFileNameList.append(posteriorTemplate % priorType)
        return posteriorFileNameList

"""
class BRAINSABCextInputSpec(BRAINSABCInputSpec):
    outputAverageImages = traits.Either(traits.Bool(True,desc="The automatically generated average images"), InputMultiPath(File(),), hash_files = False,argstr = "")
    posteriorImages = traits.Either(traits.Bool(True,desc="The automatically generated posterior images"), InputMultiPath(File(),), hash_files = False,argstr = "")
"""

class BRAINSABCextOutputSpec(BRAINSABCOutputSpec):
    # Not convenient outputAverageImages = OutputMultiPath(File(exists=True), exists = True)
    outputT1AverageImage = traits.Either( File(exists=True), None )
    outputT2AverageImage = traits.Either( File(exists=True), None )
    outputPDAverageImage = traits.Either( File(exists=True), None )
    outputFLAverageImage = traits.Either( File(exists=True), None )
    posteriorImages = OutputMultiPath(File(exists=True), exists=True)
    atlasToSubjectInverseTransform  = traits.Either( File(exists=True), None )

class BRAINSABCext(BRAINSABC):
    #input_spec= BRAINSABCextInputSpec
    output_spec = BRAINSABCextOutputSpec

    def _list_outputs(self):
        ## HACK: This function is not being called properly.
        ## -- The assert should be removed, and teh fixed_inverse_name
        ## should be properly created.
        assert False, "KILL HERE"
        custom_implied_outputs_with_no_inputs = ['posteriorImages',
                                                 'outputT1AverageImage',
                                                 'outputT2AverageImage',
                                                 'outputPDAverageImage',
                                                 'outputFLAverageImage',
                                                 'atlasToSubjectInverseTransform'
                                                  ]
        full_outputs = self.output_spec().get()
        pruned_outputs = dict()
        for key, value in full_outputs.iteritems():
            if key not in custom_implied_outputs_with_no_inputs:
                pruned_outputs[key] = value
        outputs = super(BRAINSABCext,self)._outputs_from_inputs( pruned_outputs )
        input_check = {'T1':('outputT1AverageImage', 't1_average_BRAINSABC.nii.gz'),
                       'T2':('outputT2AverageImage', 't2_average_BRAINSABC.nii.gz'),
                       'PD':('outputPDAverageImage', 'pd_average_BRAINSABC.nii.gz'),
                       'FL':('outputFLAverageImage', 'fl_average_BRAINSABC.nii.gz')}
        for key, values in input_check.iteritems():
            if key in self.inputs.inputVolumeTypes:
                outputs[values[0]] = os.path.abspath(values[1])
            else:
                outputs[values[0]] = None

        PosteriorOutputs = GetPosteriorsFromAtlasXML(self.inputs.atlasDefinition)
        PosteriorPaths = PosteriorOutputs.getPosteriorFileNameList(self.inputs.posteriorTemplate)
        outputs['posteriorImages'] = [os.path.abspath(postPath) for postPath in PosteriorPaths]

        fixed_inverse_name=self.inputs.atlasToSubjectTransform.replace(".h5","_Inverse.h5")
        print "\n"*20
        print "LOOKHERE"
        fixed_inverse_name=self.inputs.atlasToSubjectTransform.replace(".h5","_Inverse.h5")
        print  fixed_inverse_name
        print "\n"*20

        outputs['atlasToSubjectInverseTransform'] = fixed_inverse_name

        return outputs

