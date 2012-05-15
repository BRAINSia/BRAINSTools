#!/usr/bin/env python

from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine

import os

from BRAINSTools import *

from BRAINSTools.ants.antsRegistration import *

"""
    from WorkupT1T2ANTS import CreateANTSRegistrationWorkflow
    myLocalAntsWF = CreateANTSRegistrationWorkflow("31_ANTSRegistration",CLUSTER_QUEUE)
    ANTSWF.connect( SplitAvgBABC,'avgBABCT1',myLocalAntsWF,"InputSpec.fixedVolumesList")
    ANTSWF.connect( BAtlas,'template_t1',    myLocalAntsWF,"InputSpec.movingVolumesList")
    ANTSWF.connect(myLocalLMIWF,'OutputSpec.atlasToSubjectTransform',myLocalAntsWF,'InputSpec.initial_moving_transform')
"""
def CreateANTSRegistrationWorkflow(WFname,CLUSTER_QUEUE,NumberOfThreads=-1):
    ANTSWF= pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['fixedVolumesList','movingVolumesList','initial_moving_transform'
        ]), name='InputSpec' )

    many_cpu_sge_options_dictionary={'qsub_args': '-S /bin/bash -pe smp1 2-8 -l mem_free=5000M -o /dev/null -e /dev/null '+CLUSTER_QUEUE, 'overwrite': True}
    print("""Run ANTS Registration""")
    ComputeAtlasToSubjectTransform = pe.Node(interface=antsRegistration(), name="19_ComputeAtlasToSubjectTransform")
    ComputeAtlasToSubjectTransform.plugin_args=many_cpu_sge_options_dictionary
    
    ComputeAtlasToSubjectTransform.inputs.dimension=3
    ComputeAtlasToSubjectTransform.inputs.metric='CC'                  ## This is a family of interfaces, CC,MeanSquares,Demons,GC,MI,Mattes
    ComputeAtlasToSubjectTransform.inputs.transform='SyN[0.25,3.0,0.0]'
    ComputeAtlasToSubjectTransform.inputs.n_iterations=[1]
    ComputeAtlasToSubjectTransform.inputs.convergence_threshold=1e-6
    ComputeAtlasToSubjectTransform.inputs.smoothing_sigmas=[0]
    ComputeAtlasToSubjectTransform.inputs.shrink_factors=[1]
    ComputeAtlasToSubjectTransform.inputs.use_histogram_matching=True
    ComputeAtlasToSubjectTransform.inputs.invert_initial_moving_transform=False
    ComputeAtlasToSubjectTransform.inputs.output_transform_prefix='antsRegPrefix_'
    ComputeAtlasToSubjectTransform.inputs.output_warped_image='fixed_to_moving.nii.gz'
    ComputeAtlasToSubjectTransform.inputs.output_inverse_warped_image='moving_to_fixed.nii.gz'
    #if os.environ.has_key('NSLOTS'):
    #    ComputeAtlasToSubjectTransform.inputs.num_threads=int(os.environ.has_key('NSLOTS'))
    #else:
    #    ComputeAtlasToSubjectTransform.inputs.num_threads=NumberOfThreads
    # ComputeAtlasToSubjectTransform.inputs.fixedMask=SUBJ_A_small_T2_mask.nii.gz
    # ComputeAtlasToSubjectTransform.inputs.movingMask=SUBJ_B_small_T2_mask.nii.gz

    ANTSWF.connect( inputsSpec,'fixedVolumesList', ComputeAtlasToSubjectTransform,"fixed_image")
    ANTSWF.connect( inputsSpec,'movingVolumesList',ComputeAtlasToSubjectTransform,"moving_image")
    ANTSWF.connect( inputsSpec,'initial_moving_transform', ComputeAtlasToSubjectTransform,'initial_moving_transform')
    
    if 1 == 1:
        TestResampleMovingImage=pe.Node(interface=BRAINSResample(),name="99_TestAffineRegistration")
        TestResampleMovingImage.inputs.interpolationMode = "Linear"
        TestResampleMovingImage.inputs.outputVolume = "atlasToSubjectTest.nii.gz"
        ANTSWF.connect(inputsSpec,'initial_moving_transform',TestResampleMovingImage,'warpTransform')
        ANTSWF.connect(inputsSpec,'fixedVolumesList',TestResampleMovingImage,'referenceVolume')
        ANTSWF.connect(inputsSpec,'movingVolumesList',TestResampleMovingImage,'inputVolume')
    
    #############
    outputsSpec = pe.Node(interface=IdentityInterface(fields=['warped_image','inverse_warped_image','warped_transform',
            'inverse_warp_transform'
            ]), name='OutputSpec' )

    ANTSWF.connect(ComputeAtlasToSubjectTransform,'warped_image',          outputsSpec,'warped_image')
    ANTSWF.connect(ComputeAtlasToSubjectTransform,'inverse_warped_image',  outputsSpec,'inverse_warped_image')
    ANTSWF.connect(ComputeAtlasToSubjectTransform,'warped_transform',      outputsSpec,'warped_transform')
    ANTSWF.connect(ComputeAtlasToSubjectTransform,'inverse_warp_transform',outputsSpec,'inverse_warp_transform')

    return ANTSWF