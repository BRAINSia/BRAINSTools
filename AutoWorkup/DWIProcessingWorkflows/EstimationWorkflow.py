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
from SEMTools import *

def CreateEstimationWorkflow(WFname):

    EstimationWF = pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['DWI_Corrected_Aligned_CS', 'DWIBrainMask']),
                         name='inputsSpec')

    outputsSpec = pe.Node(interface=IdentityInterface(fields=['tensor_image','FAImage','MDImage',
                                                              'RDImage','FrobeniusNormImage',
                                                              'Lambda1Image','Lambda2Image','Lambda3Image',
                                                              'ukfTracks']),
                          name='outputsSpec')

    # Step1: DTI estimation
    DTIEstim = pe.Node(interface=dtiestim(), name="DTIEstim")
    DTIEstim.inputs.method = 'wls'
    DTIEstim.inputs.tensor_output = 'DTI_Output.nrrd'
    EstimationWF.connect(inputsSpec, 'DWI_Corrected_Aligned_CS', DTIEstim, 'dwi_image')
    EstimationWF.connect(inputsSpec, 'DWIBrainMask', DTIEstim, 'brain_mask')
    EstimationWF.connect(DTIEstim, 'tensor_output', outputsSpec, 'tensor_image')

    # Step2: DTI process
    DTIProcess = pe.Node(interface=dtiprocess(), name='DTIProcess')
    DTIProcess.inputs.fa_output = 'FA.nrrd'
    DTIProcess.inputs.md_output = 'MD.nrrd'
    DTIProcess.inputs.RD_output = 'RD.nrrd'
    DTIProcess.inputs.frobenius_norm_output = 'frobenius_norm_output.nrrd'
    DTIProcess.inputs.lambda1_output = 'lambda1_output.nrrd'
    DTIProcess.inputs.lambda2_output = 'lambda2_output.nrrd'
    DTIProcess.inputs.lambda3_output = 'lambda3_output.nrrd'
    DTIProcess.inputs.scalar_float = True

    EstimationWF.connect(DTIEstim, 'tensor_output', DTIProcess, 'dti_image')
    EstimationWF.connect(DTIProcess, 'fa_output', outputsSpec, 'FAImage')
    EstimationWF.connect(DTIProcess, 'md_output', outputsSpec, 'MDImage')
    EstimationWF.connect(DTIProcess, 'RD_output', outputsSpec, 'RDImage')
    EstimationWF.connect(DTIProcess, 'frobenius_norm_output', outputsSpec, 'FrobeniusNormImage')
    EstimationWF.connect(DTIProcess, 'lambda1_output', outputsSpec, 'Lambda1Image')
    EstimationWF.connect(DTIProcess, 'lambda2_output', outputsSpec, 'Lambda2Image')
    EstimationWF.connect(DTIProcess, 'lambda3_output', outputsSpec, 'Lambda3Image')

    # Step3: UKF Processing
    UKFNode = pe.Node(interface=UKFTractography(), name= "UKFRunRecordStates")
    UKFNode.inputs.tracts = "ukfTracts.vtk"
    #UKFNode.inputs.tractsWithSecondTensor = "ukfSecondTensorTracks.vtk"
    UKFNode.inputs.numTensor = '2'
    UKFNode.inputs.freeWater = True ## default False
    UKFNode.inputs.recordFA = True ## default False
    UKFNode.inputs.recordTensors = True ## default False
    #UKFNode.inputs.recordCovariance = True ## default False
    #UKFNode.inputs.recordState = True ## default False
    #UKFNode.inputs.recordFreeWater = True ## default False
    #UKFNode.inputs.recordTrace = True ## default False
    #UKFNode.inputs.recordNMSE = True ## default False

    EstimationWF.connect(inputsSpec, 'DWI_Corrected_Aligned_CS', UKFNode, 'dwiFile')
    EstimationWF.connect(inputsSpec, 'DWIBrainMask', UKFNode, 'maskFile')
    EstimationWF.connect(UKFNode,'tracts',outputsSpec,'ukfTracks')
    #DWIWorkflow.connect(UKFNode,'tractsWithSecondTensor',outputsSpec,'ukf2ndTracks')

    return EstimationWF
