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
    outputT1AverageImage = traits.Either( File(exits=True), None )
    outputT2AverageImage = traits.Either( File(exits=True), None )
    outputPDAverageImage = traits.Either( File(exits=True), None )
    outputFLAverageImage = traits.Either( File(exits=True), None )
    posteriorImages = OutputMultiPath(File(exists=True), exists = True)

class BRAINSABCext(BRAINSABC):
    #input_spec= BRAINSABCextInputSpec
    output_spec = BRAINSABCextOutputSpec

    def _list_outputs(self):
        custom_implied_outputs_with_no_inputs = [ 'posteriorImages','outputT1AverageImage','outputT2AverageImage','outputPDAverageImage','outputFLAverageImage' ]
        full_outputs = self.output_spec().get()
        pruned_outputs= dict((key,value) for key, value in full_outputs.iteritems() if key not in custom_implied_outputs_with_no_inputs )

        outputs=super(BRAINSABCext,self)._outputs_from_inputs( pruned_outputs )

        if 'T1' in self.inputs.inputVolumeTypes:
            outputs['outputT1AverageImage']=os.path.abspath('t1_average_BRAINSABC.nii.gz')
        else:
            outputs['outputT1AverageImage']=None
        if 'T2' in self.inputs.inputVolumeTypes:
            outputs['outputT2AverageImage']=os.path.abspath('t2_average_BRAINSABC.nii.gz')
        else:
            outputs['outputT2AverageImage']=None
        if 'PD' in self.inputs.inputVolumeTypes:
            outputs['outputPDAverageImage']=os.path.abspath('pd_average_BRAINSABC.nii.gz')
        else:
            outputs['outputPDAverageImage']=None
        if 'FL' in self.inputs.inputVolumeTypes:
            outputs['outputFLAverageImage']=os.path.abspath('fl_average_BRAINSABC.nii.gz')
        else:
            outputs['outputFLAverageImage']=None
        ## outputs['outputAverageImages'] = [ os.path.abspath(avgs) for avgs in templist ]
        ##  outputs['outputAverageImages'] = [ os.path.abspath('{imageType}_average.nii.gz'.format(imageType=test)) for test in set(self.inputs.inputVolumeTypes) ]

        PosteriorOutputs = GetPosteriorsFromAtlasXML(self.inputs.atlasDefinition)
        fullPosteriorPaths = [ os.path.abspath(post) for post in PosteriorOutputs.getPosteriorFileNameList(self.inputs.posteriorTemplate) ]
        outputs['posteriorImages']= fullPosteriorPaths
        return outputs
