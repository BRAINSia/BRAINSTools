## \author Ali Ghayoor
##
## This workflow gets the output DWI scan of CompressSensing,
## and creates RISs and ukfTracts
##

import os
import SimpleITK as sitk
import nipype
from nipype.interfaces import ants
from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/oS
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.semtools import *

def CreateEstimationWorkflow(WFname):
    #### Utility function ####
    def RunDTIProcess(dti_image):
        import os
        from nipype.interfaces.semtools import dtiprocess
        DTIProcess = dtiprocess()
        DTIProcess.inputs.dti_image = dti_image
        DTIProcess.inputs.fa_output = 'FA.nrrd'
        DTIProcess.inputs.md_output = 'MD.nrrd'
        DTIProcess.inputs.RD_output = 'RD.nrrd'
        DTIProcess.inputs.frobenius_norm_output = 'frobenius_norm_output.nrrd'
        DTIProcess.inputs.lambda1_output = 'lambda1_output.nrrd'
        DTIProcess.inputs.lambda2_output = 'lambda2_output.nrrd'
        DTIProcess.inputs.lambda3_output = 'lambda3_output.nrrd'
        DTIProcess.inputs.correction = 'nearest'
        DTIProcess.inputs.scalar_float = True
        DTIProcess.inputs.ignore_exception = True
        DTIProcess.run()
        fa_output = os.path.join(os.getcwd(), 'FA.nrrd')
        md_output = os.path.join(os.getcwd(), 'MD.nrrd')
        RD_output = os.path.join(os.getcwd(), 'RD.nrrd')
        frobenius_norm_output = os.path.join(os.getcwd(), 'frobenius_norm_output.nrrd')
        lambda1_output = os.path.join(os.getcwd(), 'lambda1_output.nrrd')
        lambda2_output = os.path.join(os.getcwd(), 'lambda2_output.nrrd')
        lambda3_output = os.path.join(os.getcwd(), 'lambda3_output.nrrd')
        assert os.path.isfile(fa_output), "FA file is not found: %s" % fa_output
        assert os.path.isfile(md_output), "MD file is not found: %s" % md_output
        assert os.path.isfile(RD_output), "RD file is not found: %s" % RD_output
        assert os.path.isfile(frobenius_norm_output), "frobenius_norm file is not found: %s" % frobenius_norm_output
        assert os.path.isfile(lambda1_output), "lambda1 file is not found: %s" % lambda1_output
        assert os.path.isfile(lambda2_output), "lambda2 file is not found: %s" % lambda2_output
        assert os.path.isfile(lambda3_output), "lambda3 file is not found: %s" % lambda3_output
        return fa_output,md_output,RD_output,frobenius_norm_output,lambda1_output,lambda2_output,lambda3_output
    #########################

    EstimationWF = pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['inputDWIImage', 'DWIBrainMask']),
                         name='inputsSpec')

    outputsSpec = pe.Node(interface=IdentityInterface(fields=['tensor_image','idwi_image', # from dtiestim
                                                              'FAImage','MDImage','RDImage','FrobeniusNormImage', # from dtiprocess
                                                              'Lambda1Image','Lambda2Image','Lambda3Image']), # from dtiprocess
                          name='outputsSpec')

    # Step1: DTI estimation
    DTIEstim = pe.Node(interface=dtiestim(), name="DTIEstim")
    DTIEstim.inputs.method = 'wls'
    DTIEstim.inputs.threshold = 0
    DTIEstim.inputs.correctionType = 'nearest'
    DTIEstim.inputs.tensor_output = 'DTI_Output.nrrd'
    DTIEstim.inputs.idwi = 'IDWI_Output.nrrd'
    EstimationWF.connect(inputsSpec, 'inputDWIImage', DTIEstim, 'dwi_image')
    EstimationWF.connect(inputsSpec, 'DWIBrainMask', DTIEstim, 'brain_mask')
    EstimationWF.connect(DTIEstim, 'tensor_output', outputsSpec, 'tensor_image')
    EstimationWF.connect(DTIEstim, 'idwi', outputsSpec, 'idwi_image')

    # Step2: DTI process
    # HACK: In the linux, "dtiprocess" returns a segmentation fault after it finishes all its processing
    # at the moment of destroying all registered factories.
    # Therefore, we run "dtiprocess" through an utility function to ignore this probable exception.
    # However, some checkings are done to make sure all needed outputs of dtiprocess are written properly.
    DTIProcess = pe.Node(interface=Function(function = RunDTIProcess,
                                 input_names=['dti_image'],
                                 output_names=['fa_output','md_output','RD_output','frobenius_norm_output',
                                               'lambda1_output','lambda2_output','lambda3_output']),
                         name="DTIProcess")
    '''
    DTIProcess = pe.Node(interface=dtiprocess(), name='DTIProcess')
    DTIProcess.inputs.fa_output = 'FA.nrrd'
    DTIProcess.inputs.md_output = 'MD.nrrd'
    DTIProcess.inputs.RD_output = 'RD.nrrd'
    DTIProcess.inputs.frobenius_norm_output = 'frobenius_norm_output.nrrd'
    DTIProcess.inputs.lambda1_output = 'lambda1_output.nrrd'
    DTIProcess.inputs.lambda2_output = 'lambda2_output.nrrd'
    DTIProcess.inputs.lambda3_output = 'lambda3_output.nrrd'
    DTIProcess.inputs.scalar_float = True
    '''
    EstimationWF.connect(DTIEstim, 'tensor_output', DTIProcess, 'dti_image')
    EstimationWF.connect(DTIProcess, 'fa_output', outputsSpec, 'FAImage')
    EstimationWF.connect(DTIProcess, 'md_output', outputsSpec, 'MDImage')
    EstimationWF.connect(DTIProcess, 'RD_output', outputsSpec, 'RDImage')
    EstimationWF.connect(DTIProcess, 'frobenius_norm_output', outputsSpec, 'FrobeniusNormImage')
    EstimationWF.connect(DTIProcess, 'lambda1_output', outputsSpec, 'Lambda1Image')
    EstimationWF.connect(DTIProcess, 'lambda2_output', outputsSpec, 'Lambda2Image')
    EstimationWF.connect(DTIProcess, 'lambda3_output', outputsSpec, 'Lambda3Image')

    return EstimationWF
