#!/usr/bin/env python

from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine

import os

from BRAINSTools import *

from BRAINSTools.ants.antsRegistration import *

from BRAINSTools.ants.ants import *

"""
    from WorkupT1T2ANTS import CreateANTSRegistrationWorkflow
    myLocalAntsWF = CreateANTSRegistrationWorkflow("31_ANTSRegistration",CLUSTER_QUEUE)
    ANTSWF.connect( SplitAvgBABC,'avgBABCT1',myLocalAntsWF,"InputSpec.fixedVolumesList")
    ANTSWF.connect( BAtlas,'template_t1',    myLocalAntsWF,"InputSpec.movingVolumesList")
    ANTSWF.connect(myLocalLMIWF,'OutputSpec.atlasToSubjectTransform',myLocalAntsWF,'InputSpec.initial_moving_transform')
"""
def CreateANTSRegistrationWorkflow(WFname,CLUSTER_QUEUE,NumberOfThreads=-1):
    ANTSWF= pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['fixedVolumesList','movingVolumesList','initial_moving_transform',
                                                             'fixedBinaryVolume','movingBinaryVolume'
        ]), name='InputSpec' )

    many_cpu_sge_options_dictionary={'qsub_args': '-S /bin/bash -pe smp1 2-8 -l mem_free=5000M -o /dev/null -e /dev/null '+CLUSTER_QUEUE, 'overwrite': True}
    print("""Run ANTS Registration""")

    BFitAtlasToSubject = pe.Node(interface=BRAINSFit(),name="30_BFitAtlasToSubjectAffine")
    BFitAtlasToSubject.inputs.costMetric="MMI"
    BFitAtlasToSubject.inputs.numberOfSamples=1000000
    BFitAtlasToSubject.inputs.numberOfIterations=[1500]
    BFitAtlasToSubject.inputs.numberOfHistogramBins=50
    BFitAtlasToSubject.inputs.maximumStepLength=0.2
    BFitAtlasToSubject.inputs.minimumStepLength=[0.000005]
    BFitAtlasToSubject.inputs.transformType= ["Affine"]
    BFitAtlasToSubject.inputs.maskInferiorCutOffFromCenter=65
    BFitAtlasToSubject.inputs.outputVolume="Trial_Initializer_Output.nii.gz"
    BFitAtlasToSubject.inputs.outputTransform="Trial_Initializer_Output.mat"
    BFitAtlasToSubject.inputs.maskProcessingMode="ROIAUTO"
    BFitAtlasToSubject.inputs.ROIAutoDilateSize=4
    #BFitAtlasToSubject.inputs.maskProcessingMode="ROI"
   # ANTSWF.connect(inputsSpec,'fixedBinaryVolume',BFitAtlasToSubject,'fixedBinaryVolume')
   # ANTSWF.connect(inputsSpec,'movingBinaryVolume',BFitAtlasToSubject,'movingBinaryVolume')
    ANTSWF.connect(inputsSpec,'fixedVolumesList',BFitAtlasToSubject,'fixedVolume')
    ANTSWF.connect(inputsSpec,'movingVolumesList',BFitAtlasToSubject,'movingVolume')
    ANTSWF.connect(inputsSpec,'initial_moving_transform',BFitAtlasToSubject,'initialTransform')

    if 1 == 1:
        ComputeAtlasToSubjectTransform = pe.Node(interface=antsRegistration(), name="19_ComputeAtlasToSubjectTransform")
        ComputeAtlasToSubjectTransform.plugin_args=many_cpu_sge_options_dictionary

        ComputeAtlasToSubjectTransform.inputs.dimension=3
        ComputeAtlasToSubjectTransform.inputs.metric='CC'                  ## This is a family of interfaces, CC,MeanSquares,Demons,GC,MI,Mattes
        ComputeAtlasToSubjectTransform.inputs.transform='SyN[0.25,3.0,0.0]'
        ComputeAtlasToSubjectTransform.inputs.number_of_iterations=[100,70,20]
        ComputeAtlasToSubjectTransform.inputs.convergence_threshold=1e-6
        ComputeAtlasToSubjectTransform.inputs.smoothing_sigmas=[0,0,0]
        ComputeAtlasToSubjectTransform.inputs.shrink_factors=[3,2,1]
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
        ANTSWF.connect( BFitAtlasToSubject,'outputTransform', ComputeAtlasToSubjectTransform,'initial_moving_transform')

        #############
        outputsSpec = pe.Node(interface=IdentityInterface(fields=['warped_image','inverse_warped_image','warp_transform',
                'inverse_warp_transform','affine_transform'
                ]), name='OutputSpec' )

        ANTSWF.connect(ComputeAtlasToSubjectTransform,'warped_image',          outputsSpec,'warped_image')
        ANTSWF.connect(ComputeAtlasToSubjectTransform,'inverse_warped_image',  outputsSpec,'inverse_warped_image')
        ANTSWF.connect(ComputeAtlasToSubjectTransform,'warp_transform',       outputsSpec,'warp_transform')
        ANTSWF.connect(ComputeAtlasToSubjectTransform,'inverse_warp_transform',outputsSpec,'inverse_warp_transform')
        ANTSWF.connect(BFitAtlasToSubject,'outputTransform',      outputsSpec,'affine_transform')

    else:
        ANTS_AtlasToSubjectTransform = pe.Node(interface=ANTS(), name="ANTS_ANTS_AtlasToSubjectTransform")
        ANTS_AtlasToSubjectTransform.plugin_args=many_cpu_sge_options_dictionary

        ANTS_AtlasToSubjectTransform.inputs.dimension=3
        ANTS_AtlasToSubjectTransform.inputs.output_transform_prefix='antsRegPrefix_'
        ANTS_AtlasToSubjectTransform.inputs.metric=['CC']                  ## This is a family of interfaces, CC,MeanSquares,Demons,GC,MI,Mattes
        ANTS_AtlasToSubjectTransform.inputs.metric_weight= [1.0]
        ANTS_AtlasToSubjectTransform.inputs.radius = [5]

        ANTS_AtlasToSubjectTransform.inputs.radius = [5]
        ANTS_AtlasToSubjectTransform.inputs.affine_gradient_descent_option = [0.25,0.05,0.0001,0.0001]
        ANTS_AtlasToSubjectTransform.inputs.transformation_model = 'SyN'
        ANTS_AtlasToSubjectTransform.inputs.gradient_step_length = 0.25
        ANTS_AtlasToSubjectTransform.inputs.number_of_time_steps = 3.0
        ANTS_AtlasToSubjectTransform.inputs.delta_time = 0.0
        ANTS_AtlasToSubjectTransform.inputs.number_of_iterations = [100,35,10]
        ANTS_AtlasToSubjectTransform.inputs.subsampling_factors = [3,2,1]
        ANTS_AtlasToSubjectTransform.inputs.smoothing_sigmas = [0,0,0]
        ANTS_AtlasToSubjectTransform.inputs.use_histogram_matching = True

        #ANTS_AtlasToSubjectTransform.inputs.output_warped_image='fixed_to_moving.nii.gz'
        #ANTS_AtlasToSubjectTransform.inputs.output_inverse_warped_image='moving_to_fixed.nii.gz'
        #if os.environ.has_key('NSLOTS'):
        #    ANTS_AtlasToSubjectTransform.inputs.num_threads=int(os.environ.has_key('NSLOTS'))
        #else:
        #    ANTS_AtlasToSubjectTransform.inputs.num_threads=NumberOfThreads
        # ANTS_AtlasToSubjectTransform.inputs.fixedMask=SUBJ_A_small_T2_mask.nii.gz
        # ANTS_AtlasToSubjectTransform.inputs.movingMask=SUBJ_B_small_T2_mask.nii.gz

        ANTSWF.connect( inputsSpec,'fixedVolumesList', ANTS_AtlasToSubjectTransform,"fixed_image")
        ANTSWF.connect( inputsSpec,'movingVolumesList',ANTS_AtlasToSubjectTransform,"moving_image")
        #ANTSWF.connect( BFitAtlasToSubject,'outputTransform', ANTS_AtlasToSubjectTransform,'initial_moving_transform')

        #############
        outputsSpec = pe.Node(interface=IdentityInterface(fields=['warped_image','inverse_warped_image','warp_transform',
                'inverse_warp_transform','affine_transform'
                ]), name='OutputSpec' )

        #ANTSWF.connect(ANTS_AtlasToSubjectTransform,'warped_image',          outputsSpec,'warped_image')
        #ANTSWF.connect(ANTS_AtlasToSubjectTransform,'inverse_warped_image',  outputsSpec,'inverse_warped_image')
        ANTSWF.connect(ANTS_AtlasToSubjectTransform,'affine_transform',      outputsSpec,'affine_transform')
        ANTSWF.connect(ANTS_AtlasToSubjectTransform,'warp_transform',       outputsSpec,'warp_transform')
        ANTSWF.connect(ANTS_AtlasToSubjectTransform,'inverse_warp_transform',outputsSpec,'inverse_warp_transform')

    if 0 == 1:
        TestResampleMovingImage=pe.Node(interface=BRAINSResample(),name="99_TestAffineRegistration")
        TestResampleMovingImage.inputs.interpolationMode = "Linear"
        TestResampleMovingImage.inputs.outputVolume = "atlasToSubjectTest.nii.gz"
        ANTSWF.connect(inputsSpec,'initial_moving_transform',TestResampleMovingImage,'warpTransform')
        ANTSWF.connect(inputsSpec,'fixedVolumesList',TestResampleMovingImage,'referenceVolume')
        ANTSWF.connect(inputsSpec,'movingVolumesList',TestResampleMovingImage,'inputVolume')


    return ANTSWF
