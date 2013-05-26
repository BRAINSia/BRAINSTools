# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Purpose
# =======
#
# The purpose of this pipeline is to complete all the pre-processing steps needed to turn diffusion-weighted images
# into FA images that will be used to build a template diffusion tensor atlas for fiber tracking.
#
# Inputs
# ======
# The input to this pipeline is a list of subject IDs that is used to generate lists of the corresponding DWIs
# processed with automated quality control, T2s, and brain label images that are treated as brain masks.
#
# Pipeline Steps for CreateDWIWorkflow
# ====================================
# 1. A rigid transform from the b0 of the DWI to the T2 is first derived with BRAINSFit. This rigid transform is then
#  used to resample the DWI in place into the physical space of the T2 (with gtractResampleDWIInPlace) while preserving
#  the voxel lattice of the DWI.
#
# 1. The b0 from the DWI resampled in place is extracted with extractNrrdVectorIndex. A BSpline transform from the T2
#  to the b0 of the DWI resampled in place is then derived with BRAINSFit and used to resample the brain mask into the
#  the space of the DWI resampled in place with BRAINSResample.
#
# 1. A masked tensor image is estimated with dtiprocess using the DWI resampled in place and resampled brain mask.
# dtiprocess is used again to compute FA, MD, RD, Frobenius norm, lambda1 (AD), lambda2, and lambda3 images with
# the masked tensor image.

# <codecell>

import os
import glob
import sys

#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
#####################################################################################
#     Prepend the shell environment search paths
PROGRAM_PATHS = '/Users/johnsonhj/src/BT-build/bin:/Users/johnsonhj/src/ANTs-clang/bin:/Users/johnsonhj/src/DTIPrep-build/bin:/usr/local/bin'
PROGRAM_PATHS = PROGRAM_PATHS.split(':')
PROGRAM_PATHS.extend(os.environ['PATH'].split(':'))
os.environ['PATH'] = ':'.join(PROGRAM_PATHS)

CUSTOM_ENVIRONMENT = dict()

CLUSTER_QUEUE_LONG = '-q OSX'
CLUSTER_QUEUE = '-q OSX'


# Platform specific information
#     Prepend the python search paths
PYTHON_AUX_PATHS = '/raid0/homes/johnsonhj/src/BRAINSTools/AutoWorkup:/raid0/homes/johnsonhj/src/BSA-clang31/SimpleITK-build/XXXWrapping/:/raid0/homes/johnsonhj/src/BSA-clang31/NIPYPE'
PYTHON_AUX_PATHS = PYTHON_AUX_PATHS.split(':')
PYTHON_AUX_PATHS.extend(sys.path)
sys.path = PYTHON_AUX_PATHS

import SimpleITK as sitk
import nipype
from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/oS
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.freesurfer import ReconAll
from SEMTools import *


def get_global_sge_script(pythonPathsList, binPathsList, customEnvironment={}):
    """This is a wrapper script for running commands on an SGE cluster
so that all the python modules and commands are pathed properly"""

    custEnvString = ""
    for key, value in customEnvironment.items():
        custEnvString += "export " + key + "=" + value + "\n"

    PYTHONPATH = ":".join(pythonPathsList)
    BASE_BUILDS = ":".join(binPathsList)
    GLOBAL_SGE_SCRIPT = """#!/bin/bash
echo "STARTED at: $(date +'%F-%T')"
echo "Ran on: $(hostname)"
export PATH={BINPATH}
export PYTHONPATH={PYTHONPATH}

echo "========= CUSTOM ENVIORNMENT SETTINGS =========="
echo "export PYTHONPATH={PYTHONPATH}"
echo "export PATH={BINPATH}"
echo "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"

echo "With custom environment:"
echo {CUSTENV}
{CUSTENV}
## NOTE:  nipype inserts the actual commands that need running below this section.
""".format(PYTHONPATH=PYTHONPATH, BINPATH=BASE_BUILDS, CUSTENV=custEnvString)
    return GLOBAL_SGE_SCRIPT

## Create the shell wrapper script for ensuring that all jobs running on remote hosts from SGE
#  have the same environment as the job submission host.
JOB_SCRIPT = get_global_sge_script(sys.path, PROGRAM_PATHS, CUSTOM_ENVIRONMENT)
SGE_JOB_SCRIPT = JOB_SCRIPT


def GetDWIReferenceImagesFromSessionID(SESSION_TUPLE, BASE_STRUCT, BASE_DWI):
    """A function to extract file names from base parameters"""
    import os
    import glob
    PROJ_ID = SESSION_TUPLE[0]
    SUBJ_ID = SESSION_TUPLE[1]
    SESSION_ID = SESSION_TUPLE[2]
    FixImageList = glob.glob("{BASE_STRUCT}/{PROJ_ID}/{SUBJ_ID}/{SESSION_ID}/TissueClassify/t2_average_BRAINSABC.nii.gz".format(BASE_STRUCT=BASE_STRUCT, PROJ_ID=PROJ_ID, SUBJ_ID=SUBJ_ID, SESSION_ID=SESSION_ID))
    FixMaskImageList = glob.glob("{BASE_STRUCT}/{PROJ_ID}/{SUBJ_ID}/{SESSION_ID}/TissueClassify/fixed_brainlabels_seg.nii.gz".format(BASE_STRUCT=BASE_STRUCT, PROJ_ID=PROJ_ID, SUBJ_ID=SUBJ_ID, SESSION_ID=SESSION_ID))
    MovingGlobPattern = "{BASE_DWI}/{PROJ_ID}/{SUBJ_ID}/{SESSION_ID}/*/*_[Cc][oO][nN][Cc][aA][tT]_QCed.nrrd".format(
        BASE_DWI=BASE_DWI, PROJ_ID=PROJ_ID, SUBJ_ID=SUBJ_ID, SESSION_ID=SESSION_ID)

    print(MovingGlobPattern)

    MovingDWIList = glob.glob(MovingGlobPattern)

    #### Need better methods for finding files here

    ## Should check that each list has 1 element

    print "^" * 80
    print SESSION_TUPLE
    print BASE_STRUCT
    print BASE_DWI
    print "^" * 80
    print FixImageList
    print FixMaskImageList
    print MovingDWIList
    print "^" * 80
    FixImage = FixImageList[0]
    FixMaskImage = FixMaskImageList[0]
    MovingDWI = MovingDWIList[0]

    print "=" * 80
    print FixImage
    print FixMaskImage
    print MovingDWI
    print "=" * 80

    return PROJ_ID, SUBJ_ID, SESSION_ID, FixImage, FixMaskImage, MovingDWI


def MergeByExtendListElements(FAImageList):
    ## Initial list with empty dictionaries
    ListOfImagesDictionaries = list()
    for FAImage in FAImageList:
        ListOfImagesDictionaries.append({'FA': FAImage, 'DUMMY': FAImage})
    ## HACK:  Need to make it so that AVG_AIR.nii.gz is has a background value of 1
    registrationImageTypes = ['FA']  # ['T1','T2'] someday.
    """
    ***********
    HACK: Here's the deal - interpolationMapping is supposed to define the type of interpolation done in
          the call to antsRegistration(), NOT define the function to call in the registration node!
    ***********
    """
    # DefaultContinuousInterpolationType='LanczosWindowedSinc' ## Could also be Linear for speed.
    DefaultContinuousInterpolationType = 'Linear'
    DTIinterpolationType = 'ResampleDTIlogEuclidean'
    interpolationMapping = {'T1': DefaultContinuousInterpolationType,
                            'T2': DefaultContinuousInterpolationType,
                            'PD': DefaultContinuousInterpolationType,
                            'FL': DefaultContinuousInterpolationType,
                            'FA': DefaultContinuousInterpolationType,
                            'DUMMY': DefaultContinuousInterpolationType,
                            'BRAINMASK': 'MultiLabel',
                            'DTI': DTIinterpolationType}
    return ListOfImagesDictionaries, registrationImageTypes, interpolationMapping


def CreateDWIWorkFlow(SESSION_TUPLE, BASE_STRUCT, BASE_DWI):
    WFname = "DWIWorkflow_" + str(SESSION_TUPLE[0]) + "_" + str(SESSION_TUPLE[1]) + "_" + str(SESSION_TUPLE[2])
    DWIWorkflow = pe.Workflow(name=WFname)
    # DWIWorkflow.base_dir = os.path.join(CACHE_BASE,session_name)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['SESSION_TUPLE']),
                         name='inputspec')
    inputsSpec.inputs.SESSION_TUPLE = SESSION_TUPLE

    GetFileNamesNode = pe.Node(interface=Function(function=GetDWIReferenceImagesFromSessionID,
                               input_names=['SESSION_TUPLE', 'BASE_STRUCT', 'BASE_DWI'],
                               output_names=['PROJ_ID', 'SUBJ_ID', 'SESSION_ID', 'FixImage', 'FixMaskImage', 'MovingDWI']),
                               run_without_submitting=True, name="99_GetDWIReferenceImagesFromSessionID")
    GetFileNamesNode.inputs.BASE_STRUCT = BASE_STRUCT
    GetFileNamesNode.inputs.BASE_DWI = BASE_DWI

    DWIWorkflow.connect(inputsSpec, 'SESSION_TUPLE', GetFileNamesNode, 'SESSION_TUPLE')

    outputsSpec = pe.Node(interface=IdentityInterface(fields=['FAImage', 'MDImage', 'RDImage', 'FrobeniusNormImage',
                                                              'Lambda1Image', 'Lambda2Image', 'Lambda3Image', 'tensor_image']),
                          name='outputspec')

    BFitB0_T2 = pe.Node(interface=BRAINSFit(), name="B0ToT2_Rigid")
    # BF_cpu_sge_options_dictionary = {'qsub_args': '-S /bin/bash -pe smp1 2-12 -l h_vmem=14G,mem_free=4G -o /dev/null -e /dev/null ' + CLUSTER_QUEUE, 'overwrite': True}

    # BFitB0_T2.plugin_args = BF_cpu_sge_options_dictionary
    BFitB0_T2.inputs.costMetric = "MMI"
    BFitB0_T2.inputs.numberOfSamples = 100000
    BFitB0_T2.inputs.numberOfIterations = [1500]
    BFitB0_T2.inputs.numberOfHistogramBins = 50
    BFitB0_T2.inputs.maximumStepLength = 0.2
    BFitB0_T2.inputs.minimumStepLength = [0.00005]
    BFitB0_T2.inputs.useRigid = True
    BFitB0_T2.inputs.useAffine = True  # Use Affine/but extract Rigid. Using initial transform from BRAINSABC
    BFitB0_T2.inputs.maskInferiorCutOffFromCenter = 65
    BFitB0_T2.inputs.maskProcessingMode = "ROIAUTO"
    BFitB0_T2.inputs.ROIAutoDilateSize = 13
    BFitB0_T2.inputs.backgroundFillValue = 0.0
    BFitB0_T2.inputs.initializeTransformMode = 'useCenterOfHeadAlign'

    BFitB0_T2.inputs.outputTransform = "B0ToT2_AffineTransform.h5"
    BFitB0_T2.inputs.strippedOutputTransform = "B0ToT2_RigidTransform.h5"
    BFitB0_T2.inputs.outputVolume = "B0_in_T2Space_Output.nii.gz"

    ### DWIWorkflow.connect(inputsSpec, 'T2Volume', BFitB0_T2, 'fixedVolume')
    DWIWorkflow.connect(GetFileNamesNode, 'FixImage', BFitB0_T2, 'fixedVolume')
    ### DWIWorkflow.connect(inputsSpec, 'DWIVolume', BFitB0_T2, 'movingVolume')
    DWIWorkflow.connect(GetFileNamesNode, 'MovingDWI', BFitB0_T2, 'movingVolume')

    DWIRIP = pe.Node(interface=gtractResampleDWIInPlace(), name="DWIRIP_B0ToT2")
    DWIRIP.inputs.outputVolume = 'ACPC_DWI.nrrd'
    DWIRIP.inputs.outputResampledB0 = 'ACPC_B0.nrrd'
    # DWIRIP.inputs.imageOutputSize = [164,164,100]

    DWIWorkflow.connect(BFitB0_T2, 'strippedOutputTransform', DWIRIP, 'inputTransform')
    DWIWorkflow.connect(GetFileNamesNode, 'FixImage', DWIRIP, 'referenceVolume')
    DWIWorkflow.connect(GetFileNamesNode, 'MovingDWI', DWIRIP, 'inputVolume')

    ## TODO:  Replace with ANTS
    BSPLINE_T2_TO_RIPB0 = pe.Node(interface=BRAINSFit(), name="BSPLINE_T2_TO_RIPB0")
    # BSPLINE_T2_TO_RIPB0.plugin_args = BF_cpu_sge_options_dictionary
    BSPLINE_T2_TO_RIPB0.inputs.costMetric = "MMI"
    BSPLINE_T2_TO_RIPB0.inputs.numberOfSamples = 100000
    BSPLINE_T2_TO_RIPB0.inputs.numberOfIterations = [1500]
    BSPLINE_T2_TO_RIPB0.inputs.numberOfHistogramBins = 50
    BSPLINE_T2_TO_RIPB0.inputs.maximumStepLength = 0.2
    BSPLINE_T2_TO_RIPB0.inputs.minimumStepLength = [0.00025, 0.00025, 0.00025, 0.00025, 0.00025]
    BSPLINE_T2_TO_RIPB0.inputs.useRigid = True
    BSPLINE_T2_TO_RIPB0.inputs.useScaleVersor3D = True
    BSPLINE_T2_TO_RIPB0.inputs.useScaleSkewVersor3D = True
    BSPLINE_T2_TO_RIPB0.inputs.useAffine = True  # Using initial transform from BRAINSABC

    BSPLINE_T2_TO_RIPB0.inputs.useBSpline = True
    # BSPLINE_T2_TO_RIPB0.inputs.useROIBSpline = True

    ##  This needs to be debugged, it should work. BSPLINE_T2_TO_RIPB0.inputs.useROIBSpline = True
    BSPLINE_T2_TO_RIPB0.inputs.useExplicitPDFDerivativesMode = "AUTO"
    BSPLINE_T2_TO_RIPB0.inputs.useCachingOfBSplineWeightsMode = "ON"
    BSPLINE_T2_TO_RIPB0.inputs.maxBSplineDisplacement = 24
    BSPLINE_T2_TO_RIPB0.inputs.splineGridSize = [14, 10, 12]

    BSPLINE_T2_TO_RIPB0.inputs.maskInferiorCutOffFromCenter = 65
    BSPLINE_T2_TO_RIPB0.inputs.maskProcessingMode = "ROIAUTO"
    BSPLINE_T2_TO_RIPB0.inputs.ROIAutoDilateSize = 13
    BSPLINE_T2_TO_RIPB0.inputs.backgroundFillValue = 0.0
    BSPLINE_T2_TO_RIPB0.inputs.initializeTransformMode = 'useCenterOfHeadAlign'

    BSPLINE_T2_TO_RIPB0.inputs.bsplineTransform = "T2ToRIPB0_BSplineTransform.h5"
    BSPLINE_T2_TO_RIPB0.inputs.outputVolume = "T2ToRIPB0_Output.nii.gz"

    DWIWorkflow.connect(DWIRIP, 'outputVolume', BSPLINE_T2_TO_RIPB0, 'fixedVolume')
    ### DWIWorkflow.connect(inputsSpec, 'T2Volume', BSPLINE_T2_TO_RIPB0, 'movingVolume')
    DWIWorkflow.connect(GetFileNamesNode, 'FixImage', BSPLINE_T2_TO_RIPB0, 'movingVolume')

    RESAMPLE_BRAINMASK = pe.Node(interface=BRAINSResample(), name="RESAMPLE_BRAINMASK")
    RESAMPLE_BRAINMASK.inputs.interpolationMode = 'NearestNeighbor'  # This needs to be debugged'Binary'
    RESAMPLE_BRAINMASK.inputs.outputVolume = 'DeformedBrainMaskDWIRIP.nrrd'
    RESAMPLE_BRAINMASK.inputs.pixelType = 'uchar'

    DWIWorkflow.connect(BSPLINE_T2_TO_RIPB0, 'bsplineTransform', RESAMPLE_BRAINMASK, 'warpTransform')
    ### DWIWorkflow.connect(inputsSpec, 'BrainMask',RESAMPLE_BRAINMASK,'inputVolume')
    DWIWorkflow.connect(GetFileNamesNode, 'FixMaskImage', RESAMPLE_BRAINMASK, 'inputVolume')
    DWIWorkflow.connect(DWIRIP, 'outputResampledB0', RESAMPLE_BRAINMASK, 'referenceVolume')

    DTIEstim = pe.Node(interface=dtiestim(), name="DTIEstim_Process")
    DTIEstim.inputs.method = "wls"
    DTIEstim.inputs.tensor_output = 'DTI_Output.nrrd'
    DWIWorkflow.connect(DWIRIP, 'outputVolume', DTIEstim, 'dwi_image')
    DWIWorkflow.connect(RESAMPLE_BRAINMASK, 'outputVolume', DTIEstim, 'brain_mask')

    DTIProcess = pe.Node(interface=dtiprocess(), name="DTIProcess")
    DTIProcess.inputs.fa_output = "FA.nrrd"
    DTIProcess.inputs.md_output = "MD.nrrd"
    DTIProcess.inputs.RD_output = "RD.nrrd"
    DTIProcess.inputs.frobenius_norm_output = "frobenius_norm_output.nrrd"
    DTIProcess.inputs.lambda1_output = "lambda1_output.nrrd"
    DTIProcess.inputs.lambda2_output = "lambda2_output.nrrd"
    DTIProcess.inputs.lambda3_output = "lambda3_output.nrrd"
    DTIProcess.inputs.scalar_float = True

    DWIWorkflow.connect(DTIEstim, 'tensor_output', DTIProcess, 'dti_image')
    DWIWorkflow.connect(DTIEstim, 'tensor_output', outputsSpec, 'tensor_image')
    DWIWorkflow.connect(DTIProcess, 'fa_output', outputsSpec, 'FAImage')
    DWIWorkflow.connect(DTIProcess, 'md_output', outputsSpec, 'MDImage')
    DWIWorkflow.connect(DTIProcess, 'RD_output', outputsSpec, 'RDImage')
    DWIWorkflow.connect(DTIProcess, 'frobenius_norm_output', outputsSpec, 'FrobeniusNormImage')
    DWIWorkflow.connect(DTIProcess, 'lambda1_output', outputsSpec, 'Lambda1Image')
    DWIWorkflow.connect(DTIProcess, 'lambda2_output', outputsSpec, 'Lambda2Image')
    DWIWorkflow.connect(DTIProcess, 'lambda3_output', outputsSpec, 'Lambda3Image')

    DWIDataSink = pe.Node(interface=nio.DataSink(), name='DWIDataSink')
    DWIDataSink.inputs.base_directory = '/scratch/20130214_DWIPROCESSING_NIPYPE/DWIPrototype_Results'
    # DWIDataSink.inputs.regex_substitutions = [('/Output/*/','/')]

    def sinkContainer(_tuple):
        import os.path
        return os.path.join(*_tuple)

    DWIWorkflow.connect([(inputsSpec, DWIDataSink, [(('SESSION_TUPLE', sinkContainer), 'container')])])

    DWIWorkflow.connect(outputsSpec, 'FAImage', DWIDataSink, 'Output.@FAImage')
    DWIWorkflow.connect(outputsSpec, 'MDImage', DWIDataSink, 'Output.@MDImage')
    DWIWorkflow.connect(outputsSpec, 'RDImage', DWIDataSink, 'Output.@RDImage')
    DWIWorkflow.connect(outputsSpec, 'FrobeniusNormImage', DWIDataSink, 'Output.@FrobeniusNormImage')
    DWIWorkflow.connect(outputsSpec, 'Lambda1Image', DWIDataSink, 'Output.@Lambda1Image')
    DWIWorkflow.connect(outputsSpec, 'Lambda2Image', DWIDataSink, 'Output.@Lambda2Image')
    DWIWorkflow.connect(outputsSpec, 'Lambda3Image', DWIDataSink, 'Output.@Lambda3Image')
    DWIWorkflow.connect(outputsSpec, 'tensor_image', DWIDataSink, 'Output.@tensor_image')

    return DWIWorkflow

# <codecell>

SESSIONS_TO_PROCESS = [('HDNI_001', '068044003', '068044003_20120522_30')]
#\"\"\"                 ,('PHD_024','0029','84091'),
#                     ('PHD_024','0091','60387'),('PHD_024','0091','78867'),('PHD_024','0093','50120'),
#                     ('PHD_024','0093','60307'),('PHD_024','0093','88775'),('PHD_024','0122','42742'),
#                     ('PHD_024','0122','63892'),('PHD_024','0131','76658'),('PHD_024','0131','90863'),
#                     ('PHD_024','0132','38235'),('PHD_024','0132','43991'),('PHD_024','0132','74443'),
#      "                     ('PHD_024','0133','63793'),('PHD_024','0133','81826'),('PHD_024','0137','11834'),
#      "                     ('PHD_024','0138','84460'),('PHD_024','0140','31352')]\"\"\"
"""
SLICER_REFERENCE_DIR='/scratch/DWI_DATA'
SLICER_RESULTS_DIR='/scratch/DWI_DATA/SlicerResults'
    DWIQCed='20121019_DTIPrep/HDNI_001/068044003/068044003_20120522_30/DTIPrepOutput/068044003_068044003_20120522_30_DWI_CONCAT_QCed.nrrd',
    ACPCT1='20130109_TrackOn_Results/HDNI_001/068044003/068044003_20120522_30/TissueClassify/t1_average_BRAINSABC.nii.gz',
    ACPCT2='20130109_TrackOn_Results/HDNI_001/068044003/068044003_20120522_30/TissueClassify/t2_average_BRAINSABC.nii.gz',
    FST1='20130109_TrackOn_Results/FREESURFER52_SUBJECTS/068044003_068044003_20120522_30/mri/T1.mgz',
    FSWMParc='20130109_TrackOn_Results/FREESURFER52_SUBJECTS/068044003_068044003_20120522_30/mri/wmparc.mgz'
"""
BASE_STRUCT = '/scratch/DWI_DATA/20130109_TrackOn_Results'
BASE_DWI = '/scratch/DWI_DATA/20121019_DTIPrep'

# <markdowncell>

# The desired behavior would be to use a MapNode here
MasterWFname = "ManySubjectDWIPrototype"
MasterDWIWorkflow = pe.Workflow(name=MasterWFname)
#      "\"\"\"
#      ***********
#      "Changed the cache directory"
#      "***********
#      "\"\"\"
BASE_DIR = os.path.join('/scratch/DWI_DATA/HANS_test', "HansMasterWFname")
MasterDWIWorkflow.base_dir = BASE_DIR

MasterDWIWorkflow.config['execution'] = {
    'plugin': 'Linear',
    #'stop_on_first_crash':'true',
    #'stop_on_first_rerun': 'true',
    'stop_on_first_crash': 'false',
    'stop_on_first_rerun': 'false',  # This stops at first attempt to rerun, before running, and before deleting previous results.
    'hash_method': 'timestamp',
    'single_thread_matlab': 'true',  # Multi-core 2011a  multi-core for matrix multiplication.
    'remove_unnecessary_outputs': 'false',
    'use_relative_paths': 'false',  # relative paths should be on, require hash update when changed.
    'remove_node_directories': 'false',  # Experimental
    'local_hash_check': 'true',
    'job_finished_timeout': 45
}
MasterDWIWorkflow.config['logging'] = {
    'workflow_level': 'DEBUG',
    'filemanip_level': 'DEBUG',
    'interface_level': 'DEBUG',
    'log_directory': BASE_DIR
}


if True:
    sessionNode = pe.Node(interface=IdentityInterface(fields=['SESSIONS_TO_PROCESS', 'BASE_STRUCT', 'BASE_DWI']), name='sessionNode')
    sessionNode.inputs.SESSIONS_TO_PROCESS = SESSIONS_TO_PROCESS
    sessionNode.inputs.BASE_STRUCT = BASE_STRUCT
    sessionNode.inputs.BASE_DWI = BASE_DWI
    myDWIWorkflow = pe.MapNode(interface=Function(function=CreateDWIWorkFlow,
                                                  input_names=['SESSIONS_TO_PROCESS', 'BASE_STRUCT', 'BASE_DWI'],
                                                  output_names=['FAImage', 'MDImage', 'RDImage', 'FrobeniusNormImage',
                                                                'Lambda1Image', 'Lambda2Image', 'Lambda3Image', 'tensor_output']),
                               iterfield=['SESSIONS_TO_PROCESS'],
                               name='myDWIWorkflow')
    MasterDWIWorkflow.connect(sessionNode, 'SESSIONS_TO_PROCESS', myDWIWorkflow, 'SESSIONS_TO_PROCESS')
    MasterDWIWorkflow.connect(sessionNode, 'BASE_STRUCT', myDWIWorkflow, 'BASE_STRUCT')
    MasterDWIWorkflow.connect(sessionNode, 'BASE_DWI', myDWIWorkflow, 'BASE_DWI')

    ## Now do template building with FA's
    import nipype.interfaces.ants as ants

    initAvg = pe.Node(interface=ants.AverageImages(), name='initAvg')
    initAvg.inputs.dimension = 3
    initAvg.inputs.normalize = True

    MergeByExtendListElementsNode = pe.Node(interface=Function(function=MergeByExtendListElements, input_names=['FAImageList'], output_names=['ListOfImagesDictionaries', 'registrationImageTypes', 'interpolationMapping']), run_without_submitting=True, name="99_FAMergeByExtendListElements")

    from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory, traits, isdefined, InputMultiPath, OutputMultiPath
    import os

    class ResampleDTIlogEuclideanInputSpec(CommandLineInputSpec):
        in_file = File(argstr='%s', exists=True, mandatory=True, position=-2, desc="The input file for Resampledtilogeuclidean")
        out_file = File(argstr='%s', exists=False, mandatory=True, position=-1, desc="The output file for Resampledtilogeuclidean")
        transformationFile = File(argstr='--transformationFile %s', exists=True, mandatory=True, desc="The affine transformation file for Resampledtilogeuclidean")
        Reference = File(argstr='--Reference %s', exists=True, mandatory=True, desc="The reference file for ResampleDTIlogEuclidean output")
        defField = File(argstr='--defField %s', exists=True, mandatory=True, desc="The deformation file for ResampleDTIlogEuclidean output")
        hfieldtype = traits.Enum('hfield', 'displacement', argstr="--hfieldtype %s", desc="Set if the deformation field is an -Field")
        interpolation = traits.Enum('linear', 'nn', 'ws', 'bs', argstr="--interpolation %s", desc="Sampling algorithm")
        transform_tensor_method = traits.Enum('PPD', 'FS', argstr="--transform_tensor_method %s", desc="Chooses between 2 methods to transform the tensors: Finite Strain (FS), faster but less accurate, or Preservation of the Principal Direction (PPD)")

    class ResampleDTIlogEuclideanOutputSpec(TraitedSpec):
        out_file = traits.File(exists=False, desc="The output file for Resampledtilogeuclidean")

    class ResampleDTIlogEuclidean(CommandLine):
        _cmd = "ResampleDTIlogEuclidean"
        input_spec = ResampleDTIlogEuclideanInputSpec
        output_spec = ResampleDTIlogEuclideanOutputSpec

        def _list_outputs(self):
            outputs = self.output_spec().get()
            outputs['out_file'] = self.inputs.out_file
            return outputs

    def getElement(inputList, index):
        """ Custom function to get out affine and deformation fields from BeginANTS output list """
        return inputList[int(index)]

    ResampleDTI = pe.MapNode(interface=ResampleDTIlogEuclidean(), iterfield=['in_file'], name='ResampleDTI')
    ResampleDTI.inputs.hfieldtype = 'displacement'
    ResampleDTI.inputs.interpolation = 'linear'
    ResampleDTI.inputs.transform_tensor_method = 'PPD'
    """
    ***********
    """
    mergeFA = pe.Node(interface=Merge(len(SESSIONS_TO_PROCESS)), name='99_mergeFA')
    mergeTensor = pe.Node(interface=Merge(len(SESSIONS_TO_PROCESS)), name='99_mergeTensor')
    count = 1
    for SESSION_TUPLE in SESSIONS_TO_PROCESS:
        myDWIWorkflow = CreateDWIWorkFlow(SESSION_TUPLE, BASE_STRUCT, BASE_DWI)
        MasterDWIWorkflow.connect(myDWIWorkflow, "outputspec.FAImage", mergeFA, 'in' + str(count))
        MasterDWIWorkflow.connect(myDWIWorkflow, 'outputspec.tensor_image', mergeTensor, 'in' + str(count))
        count += 1

    MasterDWIWorkflow.connect(mergeFA, 'out', MergeByExtendListElementsNode, 'FAImageList')
    MasterDWIWorkflow.connect(mergeFA, 'out', initAvg, "images")

    # resampleIDNode = pe.Node(interface=IdentityInterface(fields=['in_file', 'SESSIONS_TO_PROCESS']), name='99_resampleIdentity')
    # MasterDWIWorkflow.connect(mergeTensor, 'out', resampleIDNode, 'in_file')
    # MasterDWIWorkflow.connect(resampleIDNode, 'in_file', ResampleDTI, 'in_file')
    MasterDWIWorkflow.connect(mergeTensor, 'out', ResampleDTI, 'in_file')

    ### USE ANTS REGISTRATION
    from BAWantsRegistrationBuildTemplate import BAWantsRegistrationTemplateBuildSingleIterationWF
    buildTemplateIteration1 = BAWantsRegistrationTemplateBuildSingleIterationWF('iteration01')

    ## TODO:  Change these parameters
    BeginANTS_iter1 = buildTemplateIteration1.get_node("BeginANTS")
    BeginANTS_iter1.plugin_args = {'template': SGE_JOB_SCRIPT, 'qsub_args': '-S /bin/bash -cwd -pe smp1 4-8 -l mem_free=9000M -o /dev/null -e /dev/null {QUEUE_OPTIONS}'.format(QUEUE_OPTIONS=CLUSTER_QUEUE_LONG), 'overwrite': True}

    wimtdeformed_iter1 = buildTemplateIteration1.get_node("wimtdeformed")
    wimtdeformed_iter1.plugin_args = {'template': SGE_JOB_SCRIPT, 'qsub_args': '-S /bin/bash -cwd -pe smp1 1-2 -l mem_free=2000M -o /dev/null -e /dev/null {QUEUE_OPTIONS}'.format(QUEUE_OPTIONS=CLUSTER_QUEUE), 'overwrite': True}

    AvgAffineTransform_iter1 = buildTemplateIteration1.get_node("AvgAffineTransform")
    AvgAffineTransform_iter1.plugin_args = {'template': SGE_JOB_SCRIPT, 'qsub_args': '-S /bin/bash -cwd -pe smp1 1 -l mem_free=2000M -o /dev/null -e /dev/null {QUEUE_OPTIONS}'.format(QUEUE_OPTIONS=CLUSTER_QUEUE), 'overwrite': True}

    wimtPassivedeformed_iter1 = buildTemplateIteration1.get_node("wimtPassivedeformed")
    wimtPassivedeformed_iter1.plugin_args = {'template': SGE_JOB_SCRIPT, 'qsub_args': '-S /bin/bash -cwd -pe smp1 1-2 -l mem_free=2000M -o /dev/null -e /dev/null {QUEUE_OPTIONS}'.format(QUEUE_OPTIONS=CLUSTER_QUEUE), 'overwrite': True}

    MasterDWIWorkflow.connect(initAvg, 'output_average_image', buildTemplateIteration1, 'inputspec.fixed_image')
    MasterDWIWorkflow.connect(MergeByExtendListElementsNode, 'ListOfImagesDictionaries', buildTemplateIteration1, 'inputspec.ListOfImagesDictionaries')
    MasterDWIWorkflow.connect(MergeByExtendListElementsNode, 'registrationImageTypes', buildTemplateIteration1, 'inputspec.registrationImageTypes')
    MasterDWIWorkflow.connect(MergeByExtendListElementsNode, 'interpolationMapping', buildTemplateIteration1, 'inputspec.interpolationMapping')

    buildTemplateIteration2 = buildTemplateIteration1.clone(name='buildTemplateIteration2')
    buildTemplateIteration2 = BAWantsRegistrationTemplateBuildSingleIterationWF('Iteration02')
    MasterDWIWorkflow.connect([(buildTemplateIteration2, ResampleDTI, [(('BeginANTS.forward_transforms', getElement, 0), 'transformationFile')])])
    MasterDWIWorkflow.connect([(buildTemplateIteration2, ResampleDTI, [(('BeginANTS.forward_transforms', getElement, 1), 'defField')])])
    MasterDWIWorkflow.connect(buildTemplateIteration2, 'outputspec.template', ResampleDTI, 'Reference')
    ## TODO:  Change these parameters
    BeginANTS_iter2 = buildTemplateIteration2.get_node("BeginANTS")
    BeginANTS_iter2.plugin_args = {'template': SGE_JOB_SCRIPT, 'qsub_args': '-S /bin/bash -cwd -pe smp1 4-8 -l mem_free=9000M -o /dev/null -e /dev/null {QUEUE_OPTIONS}'.format(QUEUE_OPTIONS=CLUSTER_QUEUE_LONG), 'overwrite': True}

    wimtdeformed_iter2 = buildTemplateIteration2.get_node("wimtdeformed")
    wimtdeformed_iter2.plugin_args = {'template': SGE_JOB_SCRIPT, 'qsub_args': '-S /bin/bash -cwd -pe smp1 1-2 -l mem_free=2000M -o /dev/null -e /dev/null {QUEUE_OPTIONS}'.format(QUEUE_OPTIONS=CLUSTER_QUEUE), 'overwrite': True}

    AvgAffineTransform_iter2 = buildTemplateIteration2.get_node("AvgAffineTransform")
    AvgAffineTransform_iter2.plugin_args = {'template': SGE_JOB_SCRIPT, 'qsub_args': '-S /bin/bash -cwd -pe smp1 1 -l mem_free=2000M -o /dev/null -e /dev/null {QUEUE_OPTIONS}'.format(QUEUE_OPTIONS=CLUSTER_QUEUE), 'overwrite': True}

    wimtPassivedeformed_iter2 = buildTemplateIteration2.get_node("wimtPassivedeformed")
    wimtPassivedeformed_iter2.plugin_args = {'template': SGE_JOB_SCRIPT, 'qsub_args': '-S /bin/bash -cwd -pe smp1 1-2 -l mem_free=2000M -o /dev/null -e /dev/null {QUEUE_OPTIONS}'.format(QUEUE_OPTIONS=CLUSTER_QUEUE), 'overwrite': True}

    MasterDWIWorkflow.connect(buildTemplateIteration1, 'outputspec.template', buildTemplateIteration2, 'inputspec.fixed_image')
    MasterDWIWorkflow.connect(MergeByExtendListElementsNode, 'ListOfImagesDictionaries', buildTemplateIteration2, 'inputspec.ListOfImagesDictionaries')
    MasterDWIWorkflow.connect(MergeByExtendListElementsNode, 'registrationImageTypes', buildTemplateIteration2, 'inputspec.registrationImageTypes')
    MasterDWIWorkflow.connect(MergeByExtendListElementsNode, 'interpolationMapping', buildTemplateIteration2, 'inputspec.interpolationMapping')

    DWIAverageSink = pe.Node(interface=nio.DataSink(), name="DWIAverageSink")
    #""" 'Average DTIEstim_Process and average deformed DTI image' """
    DWIAverageSink = pe.Node(interface=nio.DataSink(), name='DWIAverageSink')
    DWIAverageSink.inputs.base_directory = '/scratch/20130214_DWIPROCESSING_NIPYPE'
    DWIAverageSink.inputs.container = os.path.join("DWIPrototype_Results")
    # DWIAverageSink.inputs.regexp_substitutions = [('Average/*/', 'Average/')]

    MasterDWIWorkflow.connect(buildTemplateIteration2, 'AvgAffineTransform.affine_transform',
                              DWIAverageSink, 'Average.@transforms')
    MasterDWIWorkflow.connect(buildTemplateIteration2, 'AvgWarpImages.output_average_image',
                              DWIAverageSink, 'Average.@passive_transforms')
    MasterDWIWorkflow.connect(buildTemplateIteration2, 'outputspec.template',
                              DWIAverageSink, 'Average.@template')
    MasterDWIWorkflow.connect(buildTemplateIteration2, 'outputspec.passive_deformed_templates',
                              DWIAverageSink, 'Average.@passive')
    MasterDWIWorkflow.connect(initAvg, 'output_average_image', DWIAverageSink, 'Average.@image')

    DTIDataSink = pe.MapNode(interface=nio.DataSink(infields=['Output.@out_file']),
                             iterfield=['container', 'Output.@out_file'],
                             name='DTIDataSink')
    DTIDataSink.inputs.base_directory = '/scratch/20130214_DWIPROCESSING_NIPYPE/DWIPrototype_Results'
    # DWIDataSink.inputs.regex_substitutions = [('/Output/*/','/')]

    def sinkContainer(args):
        import os.path
        retval = []
        for tup in args:
            retval.append(os.path.join(*tup))
        return retval

"""
    MasterDWIWorkflow.connect([(resampleIDNode, DTIDataSink,[(('SESSIONS_TO_PROCESS', sinkContainer), 'container')])])
    MasterDWIWorkflow.connect(ResampleDTI, 'out_file', DTIDataSink, 'Output.@out_file')
"""


# <markdowncell>

# Now Start Processing
# ====================

# <codecell>
import multiprocessing
total_CPUS = multiprocessing.cpu_count()
NUMPARALLEL = 1
os.environ['NSLOTS'] = "{0}".format(total_CPUS / NUMPARALLEL)

MasterDWIWorkflow.write_graph()

SGEFlavor = 'SGE'

if True:
    MasterDWIWorkflow.write_graph()
    MasterDWIWorkflow.run()
else:
    MasterDWIWorkflow.run(plugin=SGEFlavor,
                          plugin_args=dict(template=JOB_SCRIPT,
                                           qsub_args="-S /bin/bash -cwd -pe smp1 1-12 -l h_vmem=19G,mem_free=2G -o /dev/null -e /dev/null " + "-q OSX"))


print sys.argv
print sys.api_version
