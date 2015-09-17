from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory, traits, isdefined, InputMultiPath, OutputMultiPath
import os
from nipype.interfaces.semtools.diffusion.dtiprep import DTIPrepInputSpec, DTIPrepOutputSpec, DTIPrep


class DTIPrepextOutputSpec(DTIPrepOutputSpec):
	outputVolume =    traits.Either(File(exists=True), None)
	outputReportXML = traits.Either(File(exists=True), None)
	outputReportTxt = traits.Either(File(exists=True), None)
	# outputQCedBaseline = traits.Either(File(exists=True), None)
	# outputQCedDTI = traits.Either(File(exists=True), None)
	# outputQCedDTI_FA = traits.Either(File(exists=True), None)
	# outputQCedDTI_MD = traits.Either(File(exists=True), None)
	# outputQCedDTI_colorFA = traits.Either(File(exists=True), None)
	# outputQCedDTI_frobeniusnorm = traits.Either(File(exists=True), None)
	# outputQCedIDWI = traits.Either(File(exists=True), None)


class DTIPrepext(DTIPrep):
    #input_spec = DTIPrepextInputSpec
    output_spec = DTIPrepextOutputSpec

    def _list_outputs(self):
        custom_implied_outputs_with_no_inputs = ['outputVolume',
                                                 'outputReportXML',
                                                 'outputReportTxt'
                                                 ]
        full_outputs = self._outputs().get()
        pruned_outputs = dict()
        for key, value in list(full_outputs.items()):
            if key not in custom_implied_outputs_with_no_inputs:
                pruned_outputs[key] = value
        outputs = super(DTIPrepext, self)._outputs_from_inputs(pruned_outputs)
        inputDir, filename = os.path.split(self.inputs.DWINrrdFile)
        filenameList = filename.split(".")
        prefix = '.'.join(filenameList[:-1])
                                                                                             #ConcatenatedDWIFile.nrrd
                                                                                             #ConcatenatedDWIFile_QCed.nrrd
                                                                                             #ConcatenatedDWIFile_QCed.nrrd
        outputs['outputVolume'] =         os.path.abspath( os.path.join(self.inputs.outputFolder, prefix + "_QCed.nrrd") )
        outputs['outputReportXML'] =      os.path.abspath( os.path.join(self.inputs.outputFolder, prefix + "_XMLQCResult.xml") )
        outputs['outputReportTxt'] =      os.path.abspath( os.path.join(self.inputs.outputFolder, prefix + "_QCReport.txt") )
        # outputs['outputQCedBaseline'] = os.path.abspath( os.path.join(self.inputs.outputFolder, prefix + "_QCed_Baseline.nrrd") )
        # outputs['outputQCedDTI'] =      os.path.abspath( os.path.join(self.inputs.outputFolder, prefix + "_QCed_DTI.nrrd") )
        # outputs['outputQCedDTI_FA'] =   os.path.abspath( os.path.join(self.inputs.outputFolder, prefix + "_QCed_DTI_FA.nrrd") )
        # outputs['outputQCedDTI_MD'] =   os.path.abspath( os.path.join(self.inputs.outputFolder, prefix + "_QCed_DTI_MD.nrrd") )
        # outputs['outputQCedDTI_colorFA'] =       os.path.abspath( os.path.join(self.inputs.outputFolder, prefix + "_QCed_DTI_colorFA.nrrd") )
        # outputs['outputQCedDTI_frobeniusnorm'] = os.path.abspath( os.path.join(self.inputs.outputFolder, prefix + "_QCed_DTI_frobeniusnorm.nrrd") )
        # outputs['outputQCedIDWI'] =              os.path.abspath( os.path.join(self.inputs.outputFolder, prefix + "_QCed_IDWI.nrrd") )
        return outputs
