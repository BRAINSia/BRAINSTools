#!/usr/bin/env python

from __future__ import print_function
from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine

import os

from nipype.interfaces.semtools import *
from utilities.misc import CommonANTsRegistrationSettings

"""
    from WorkupT1T2ANTS import CreateANTSRegistrationWorkflow
    myLocalAntsWF = CreateANTSRegistrationWorkflow("ANTSRegistration",CLUSTER_QUEUE)
    ANTSWF.connect( SplitAvgBABC,'avgBABCT1',myLocalAntsWF,"inputspec.fixedVolumesList")
    ANTSWF.connect( BAtlas,'template_t1_denoised_gaussian',    myLocalAntsWF,"inputspec.movingVolumesList")
    ANTSWF.connect(myLocalLMIWF,'outputspec.atlasToSubjectTransform',myLocalAntsWF,'inputspec.initial_moving_transform')
"""


def CreateANTSRegistrationWorkflow(WFname, CLUSTER_QUEUE, CLUSTER_QUEUE_LONG, NumberOfThreads=-1):
    ANTSWF = pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['fixedVolumesList', 'movingVolumesList', 'initial_moving_transform',
                                                             'fixedBinaryVolume', 'movingBinaryVolume', 'warpFixedVolumesList'
                                                             ]), name='inputspec')

    print("""Run ANTS Registration""")

    BFitAtlasToSubject = pe.Node(interface=BRAINSFit(), name="bfA2S")
    BF_cpu_sge_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,8,4,16), 'overwrite': True}
    BFitAtlasToSubject.plugin_args = BF_cpu_sge_options_dictionary
    BFitAtlasToSubject.inputs.costMetric = "MMI"
    BFitAtlasToSubject.inputs.numberOfSamples = 1000000
    BFitAtlasToSubject.inputs.numberOfIterations = [1500]
    BFitAtlasToSubject.inputs.numberOfHistogramBins = 50
    BFitAtlasToSubject.inputs.maximumStepLength = 0.2
    BFitAtlasToSubject.inputs.minimumStepLength = [0.000005]
    BFitAtlasToSubject.inputs.useAffine = True  # Using initial transform from BRAINSABC
    BFitAtlasToSubject.inputs.maskInferiorCutOffFromCenter = 65
    BFitAtlasToSubject.inputs.outputVolume = "Trial_Initializer_Output.nii.gz"
    # Bug in BRAINSFit PREDICTIMG-1379 BFitAtlasToSubject.inputs.outputFixedVolumeROI="FixedROI.nii.gz"
    # Bug in BRAINSFit PREDICTIMG-1379 BFitAtlasToSubject.inputs.outputMovingVolumeROI="MovingROI.nii.gz"
    BFitAtlasToSubject.inputs.outputTransform = "Trial_Initializer_Output.h5"
    BFitAtlasToSubject.inputs.maskProcessingMode = "ROIAUTO"
    BFitAtlasToSubject.inputs.ROIAutoDilateSize = 4
    BFitAtlasToSubject.inputs.writeOutputTransformInFloat = True
    # BFitAtlasToSubject.inputs.maskProcessingMode="ROI"
   # ANTSWF.connect(inputsSpec,'fixedBinaryVolume',BFitAtlasToSubject,'fixedBinaryVolume')
   # ANTSWF.connect(inputsSpec,'movingBinaryVolume',BFitAtlasToSubject,'movingBinaryVolume')
    ANTSWF.connect(inputsSpec, 'fixedVolumesList', BFitAtlasToSubject, 'fixedVolume')
    ANTSWF.connect(inputsSpec, 'movingVolumesList', BFitAtlasToSubject, 'movingVolume')
    ANTSWF.connect(inputsSpec, 'initial_moving_transform', BFitAtlasToSubject, 'initialTransform')

    ComputeAtlasToSubjectTransform = pe.Node(interface=antsRegistration(), name="antsA2S")
    many_cpu_sge_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,8,8,16), 'overwrite': True}
    ComputeAtlasToSubjectTransform.plugin_args = many_cpu_sge_options_dictionary
    CommonANTsRegistrationSettings(
            antsRegistrationNode=ComputeAtlasToSubjectTransform,
            registrationTypeDescription="FromWorkupT1T2ANTS.py",
            output_transform_prefix='antsRegPrefix_',
            output_warped_image='moving_to_fixed.nii.gz',
            output_inverse_warped_image='fixed_to_moving.nii.gz',
            save_state = None)


    ANTSWF.connect(inputsSpec, 'fixedVolumesList', ComputeAtlasToSubjectTransform, "fixed_image")
    ANTSWF.connect(inputsSpec, 'movingVolumesList', ComputeAtlasToSubjectTransform, "moving_image")
    ANTSWF.connect(BFitAtlasToSubject, 'outputTransform', ComputeAtlasToSubjectTransform, 'initial_moving_transform')

    if 1 == 1:
        mergeAffineWarp = pe.Node(interface=Merge(2), name="Merge_AffineWarp")
        ANTSWF.connect(ComputeAtlasToSubjectTransform, 'warp_transform', mergeAffineWarp, 'in1')
        ANTSWF.connect(BFitAtlasToSubject, 'outputTransform', mergeAffineWarp, 'in2')

        from nipype.interfaces.ants import WarpImageMultiTransform
        debugWarpTest = pe.Node(interface=WarpImageMultiTransform(), name="dbgWarpTest")
        # Not allowed as an input debugWarpTest.inputs.output_image = 'debugWarpedMovingToFixed.nii.gz'

        ANTSWF.connect(inputsSpec, 'fixedVolumesList', debugWarpTest, 'reference_image')
        ANTSWF.connect(inputsSpec, 'movingVolumesList', debugWarpTest, 'moving_image')
        ANTSWF.connect(mergeAffineWarp, 'out', debugWarpTest, 'transformation_series')

    #############
    outputsSpec = pe.Node(interface=IdentityInterface(fields=['warped_image', 'inverse_warped_image', 'warp_transform',
                                                              'inverse_warp_transform', 'affine_transform'
                                                              ]), name='outputspec')

    ANTSWF.connect(ComputeAtlasToSubjectTransform, 'warped_image', outputsSpec, 'warped_image')
    ANTSWF.connect(ComputeAtlasToSubjectTransform, 'inverse_warped_image', outputsSpec, 'inverse_warped_image')
    ANTSWF.connect(ComputeAtlasToSubjectTransform, 'warp_transform', outputsSpec, 'warp_transform')
    ANTSWF.connect(ComputeAtlasToSubjectTransform, 'inverse_warp_transform', outputsSpec, 'inverse_warp_transform')
    ANTSWF.connect(BFitAtlasToSubject, 'outputTransform', outputsSpec, 'affine_transform')

    return ANTSWF
