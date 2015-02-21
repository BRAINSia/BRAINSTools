#!/usr/bin/env python

from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine

from utilities.misc import *
from utilities.distributed import modify_qsub_args
from SEMTools.utilities.brains import BRAINSLandmarkInitializer
from SEMTools import BRAINSSnapShotWriter

"""
    from WorkupT1T2MALF import CreateMALFWorkflow
    myLocalTCWF= CreateMALFWorkflow("MALF")
    MALFWF.connect( [ (uidSource, myLocalTCWF, [(('uid', getT1s, subjectDatabaseFile ), 'T1List')] ), ])
    MALFWF.connect( [ (uidSource, myLocalTCWF, [(('uid', getT2s, subjectDatabaseFile ), 'T2List')] ), ])
    MALFWF.connect( [ (uidSource, myLocalTCWF, [(('uid', getT1sLength, subjectDatabaseFile ), 'T1_count')] ), ])
    MALFWF.connect( BCD,    'outputResampledVolume', myLocalTCWF, 'PrimaryT1' )
    MALFWF.connect(BAtlas,'ExtendedAtlasDefinition.xml',myLocalTCWF,'atlasDefinition')
    MALFWF.connect(BLI,'outputTransformFilename',myLocalTCWF,'atlasToSubjectInitialTransform')
"""

def getListIndexOrNoneIfOutOfRange(imageList, index):
    if index < len(imageList):
        return imageList[index]
    else:
        return None


def CreateMALFWorkflow(WFname, master_config,good_subjects,BASE_DATA_GRABBER_DIR):
    from nipype.interfaces import ants

    CLUSTER_QUEUE=master_config['queue']
    CLUSTER_QUEUE_LONG=master_config['long_q']

    MALFWF = pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['subj_t1_image', #Desired image to create label map for
                                                             'subj_lmks', #The landmarks corresponding to t1_image
                                                             'subj_fixed_head_labels', #The fixed head labels from BABC
                                                             'subj_left_hemisphere', #The warped left hemisphere mask
                                                             'atlasWeightFilename'  #The static weights file name
                                                            ]),
                         run_without_submitting=True,
                         name='inputspec')
    outputsSpec = pe.Node(interface=IdentityInterface(fields=['MALF_neuro2012_labelmap',
                                                       'MALF_fswm_extended_neuro2012_labelmap',
                                                       'MALF_fswm_standard_neuro2012_labelmap',
                                                       'MALF_extended_snapshot']),
                          run_without_submitting=True,
                          name='outputspec')

    BLICreator = dict()
    MALF_DG = dict()
    A2SantsRegistrationPreABCRigid =dict()
    A2SantsRegistrationPreABCSyN = dict()
    fixedROIAuto = dict()
    movingROIAuto = dict()
    labelMapResample = dict()
    NewlabelMapResample = dict()

    warpedAtlasT1MergeNode = pe.Node(interface=Merge(len(good_subjects)),name="T1sMergeAtlas")
    warpedAtlasLblMergeNode = pe.Node(interface=Merge(len(good_subjects)),name="LblMergeAtlas")
    NewwarpedAtlasLblMergeNode = pe.Node(interface=Merge(len(good_subjects)),name="fswmLblMergeAtlas")
    malf_atlas_mergeindex = 1;
    for malf_atlas_subject in good_subjects:
        ## Need DataGrabber Here For the Atlas
        MALF_DG[malf_atlas_subject] = pe.Node(interface=nio.DataGrabber(infields=['subject'],
                                                        outfields=['malf_atlas_t1',
                                                                   'malf_atlas_lbls',
                                                                   'malf_fswm_atlas_lbls',
                                                                   'malf_atlas_lmks'
                                                        ]),
                              run_without_submitting=True,name='MALF_DG_'+malf_atlas_subject)
        #MALF_DG[malf_atlas_subject].inputs.base_directory = master_config['previousresult']
        MALF_DG[malf_atlas_subject].inputs.base_directory = BASE_DATA_GRABBER_DIR

        MALF_DG[malf_atlas_subject].inputs.subject = malf_atlas_subject
        MALF_DG[malf_atlas_subject].inputs.field_template = {
                                             'malf_atlas_t1': '%s/TissueClassify/t1_average_BRAINSABC.nii.gz',
                                             'malf_atlas_lbls': '%s/TissueClassify/neuro_lbls.nii.gz',
                                             'malf_fswm_atlas_lbls': '%s/TissueClassify/neuro_lbls_MALF_In_FSLabel.nii.gz',
                                             'malf_atlas_lmks': '%s/ACPCAlign/BCD_ACPC_Landmarks.fcsv',
        }
        MALF_DG[malf_atlas_subject].inputs.template_args = {
                                            'malf_atlas_t1':   [['subject']],
                                            'malf_atlas_lbls': [['subject']],
                                            'malf_fswm_atlas_lbls': [['subject']],
                                            'malf_atlas_lmks': [['subject']],
        }
        MALF_DG[malf_atlas_subject].inputs.template = '*'
        MALF_DG[malf_atlas_subject].inputs.sort_filelist = True
        MALF_DG[malf_atlas_subject].inputs.raise_on_empty = True

        ## Create BLI first
        ########################################################
        # Run BLI atlas_to_subject
        ########################################################
        BLICreator[malf_atlas_subject] = pe.Node(interface=BRAINSLandmarkInitializer(), name="BLI_"+malf_atlas_subject)
        BLICreator[malf_atlas_subject].inputs.outputTransformFilename = "landmarkInitializer_{0}_to_subject_transform.h5".format(malf_atlas_subject)

        MALFWF.connect(inputsSpec, 'atlasWeightFilename', BLICreator[malf_atlas_subject], 'inputWeightFilename')
        MALFWF.connect(MALF_DG[malf_atlas_subject], 'malf_atlas_lmks', BLICreator[malf_atlas_subject], 'inputMovingLandmarkFilename')
        MALFWF.connect(inputsSpec, 'subj_lmks', BLICreator[malf_atlas_subject], 'inputFixedLandmarkFilename')


        ##### Initialize with ANTS Transform For AffineComponentBABC
        currentAtlasToSubjectantsRigidRegistration = 'Rigid_AtlasToSubjectANTsPreABC_'+malf_atlas_subject
        A2SantsRegistrationPreABCRigid[malf_atlas_subject] = pe.Node(interface=ants.Registration(), name=currentAtlasToSubjectantsRigidRegistration)
        many_cpu_ANTsRigid_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,2,1,1), 'overwrite': True}
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].plugin_args = many_cpu_ANTsRigid_options_dictionary

        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.num_threads   = -1
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.dimension = 3
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.transforms = ["Affine",]
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.transform_parameters = [[0.1]]
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.metric = ['MI']
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.sampling_strategy = ['Regular']
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.sampling_percentage = [0.5]
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.metric_weight = [1.0]
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.radius_or_number_of_bins = [32]
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.number_of_iterations = [[1000,1000, 500, 100]]

        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.convergence_threshold = [1e-8]

        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.convergence_window_size = [10]
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.use_histogram_matching = [True]
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.shrink_factors = [[8, 4, 2, 1]]
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.smoothing_sigmas = [[3, 2, 1, 0]]
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.sigma_units = ["vox"]
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.use_estimate_learning_rate_once = [False]
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.write_composite_transform = True  # Required for initialize_transforms_per_stage
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.collapse_output_transforms = False # Mutually Exclusive with initialize_transforms_per_stage
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.initialize_transforms_per_stage = True
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.output_transform_prefix = 'AtlasToSubjectPreBABC_Rigid'
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.winsorize_lower_quantile = 0.01
        A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.winsorize_upper_quantile = 0.99
        ## NO NEED FOR THIS A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.output_warped_image = 'atlas2subjectRigid.nii.gz'
        ## NO NEED FOR THIS A2SantsRegistrationPreABCRigid[malf_atlas_subject].inputs.output_inverse_warped_image = 'subject2atlasRigid.nii.gz'

        MALFWF.connect(BLICreator[malf_atlas_subject], 'outputTransformFilename',A2SantsRegistrationPreABCRigid[malf_atlas_subject],'initial_moving_transform')
        MALFWF.connect(inputsSpec, 'subj_t1_image',A2SantsRegistrationPreABCRigid[malf_atlas_subject],'fixed_image')
        MALFWF.connect(MALF_DG[malf_atlas_subject], 'malf_atlas_t1',A2SantsRegistrationPreABCRigid[malf_atlas_subject],'moving_image')


        ##### Initialize with ANTS Transform For SyN component BABC
        currentAtlasToSubjectantsRegistration = 'SyN_AtlasToSubjectANTsPreABC_'+malf_atlas_subject
        A2SantsRegistrationPreABCSyN[malf_atlas_subject] = pe.Node(interface=ants.Registration(), name=currentAtlasToSubjectantsRegistration)
        many_cpu_ANTsSyN_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE_LONG,4,2,16), 'overwrite': True}
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].plugin_args = many_cpu_ANTsSyN_options_dictionary

        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.num_threads   = -1
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.dimension = 3
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.transforms = ["SyN","SyN"]
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.transform_parameters = [[0.1, 3, 0],[0.1, 3, 0] ]
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.metric = ['MI','MI']
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.sampling_strategy = [None,None]
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.sampling_percentage = [1.0,1.0]
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.metric_weight = [1.0,1.0]
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.radius_or_number_of_bins = [32,32]
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.number_of_iterations = [[500, 500, 500, 500 ], [70]]

        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.convergence_threshold = [1e-8,1e-4]

        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.convergence_window_size = [12]
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.use_histogram_matching = [True,True]
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.shrink_factors = [[8, 4, 3, 2], [1]]
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.smoothing_sigmas = [[3, 2, 2, 1], [0]]
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.sigma_units = ["vox","vox"]
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.use_estimate_learning_rate_once = [False,False]
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.write_composite_transform = True # Required for initialize_transforms_per_stage
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.collapse_output_transforms = False # Mutually Exclusive with initialize_transforms_per_stage
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.initialize_transforms_per_stage = True
        ## NO NEED FOR THIS A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.save_state = 'SavedInternalSyNState.h5'
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.output_transform_prefix = malf_atlas_subject+'_ToSubjectPreBABC_SyN'
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.winsorize_lower_quantile = 0.01
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.winsorize_upper_quantile = 0.99
        A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.output_warped_image = malf_atlas_subject + '_2subject.nii.gz'
        ## NO NEED FOR THIS A2SantsRegistrationPreABCSyN[malf_atlas_subject].inputs.output_inverse_warped_image = 'subject2atlas.nii.gz'

        ## if using Registration masking, then do ROIAuto on fixed and moving images and connect to registraitons
        UseRegistrationMasking = True
        if UseRegistrationMasking == True:
            from SEMTools.segmentation.specialized import BRAINSROIAuto

            fixedROIAuto[malf_atlas_subject] = pe.Node(interface=BRAINSROIAuto(), name="fixedROIAUTOMask_"+malf_atlas_subject)
            fixedROIAuto[malf_atlas_subject].inputs.ROIAutoDilateSize=10
            fixedROIAuto[malf_atlas_subject].inputs.outputROIMaskVolume = "fixedImageROIAutoMask.nii.gz"

            movingROIAuto[malf_atlas_subject] = pe.Node(interface=BRAINSROIAuto(), name="movingROIAUTOMask_"+malf_atlas_subject)
            fixedROIAuto[malf_atlas_subject].inputs.ROIAutoDilateSize=10
            movingROIAuto[malf_atlas_subject].inputs.outputROIMaskVolume = "movingImageROIAutoMask.nii.gz"

            MALFWF.connect(inputsSpec, 'subj_t1_image',fixedROIAuto[malf_atlas_subject],'inputVolume')
            MALFWF.connect(MALF_DG[malf_atlas_subject], 'malf_atlas_t1', movingROIAuto[malf_atlas_subject],'inputVolume')

            MALFWF.connect(fixedROIAuto[malf_atlas_subject], 'outputROIMaskVolume',A2SantsRegistrationPreABCRigid[malf_atlas_subject],'fixed_image_mask')
            MALFWF.connect(movingROIAuto[malf_atlas_subject], 'outputROIMaskVolume',A2SantsRegistrationPreABCRigid[malf_atlas_subject],'moving_image_mask')

            MALFWF.connect(fixedROIAuto[malf_atlas_subject], 'outputROIMaskVolume',A2SantsRegistrationPreABCSyN[malf_atlas_subject],'fixed_image_mask')
            MALFWF.connect(movingROIAuto[malf_atlas_subject], 'outputROIMaskVolume',A2SantsRegistrationPreABCSyN[malf_atlas_subject],'moving_image_mask')

        MALFWF.connect(A2SantsRegistrationPreABCRigid[malf_atlas_subject],
                                 ('composite_transform', getListIndexOrNoneIfOutOfRange, 0 ),
                                 A2SantsRegistrationPreABCSyN[malf_atlas_subject],'initial_moving_transform')
        MALFWF.connect(inputsSpec, 'subj_t1_image',A2SantsRegistrationPreABCSyN[malf_atlas_subject],'fixed_image')
        MALFWF.connect(MALF_DG[malf_atlas_subject], 'malf_atlas_t1',A2SantsRegistrationPreABCSyN[malf_atlas_subject],'moving_image')
        MALFWF.connect(A2SantsRegistrationPreABCSyN[malf_atlas_subject],'warped_image',warpedAtlasT1MergeNode,'in'+str(malf_atlas_mergeindex) )

        ### Original labelmap resampling
        labelMapResample[malf_atlas_subject] = pe.Node(interface=ants.ApplyTransforms(),name="WLABEL_"+malf_atlas_subject)
        many_cpu_labelMapResample_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,1,1,1), 'overwrite': True}
        labelMapResample[malf_atlas_subject].plugin_args = many_cpu_labelMapResample_options_dictionary
        labelMapResample[malf_atlas_subject].inputs.dimension=3
        labelMapResample[malf_atlas_subject].inputs.output_image=malf_atlas_subject+'_2_subj_lbl.nii.gz'
        labelMapResample[malf_atlas_subject].inputs.interpolation='MultiLabel'
        labelMapResample[malf_atlas_subject].inputs.default_value=0
        labelMapResample[malf_atlas_subject].inputs.invert_transform_flags=[False]

        MALFWF.connect( A2SantsRegistrationPreABCSyN[malf_atlas_subject],'composite_transform',
                        labelMapResample[malf_atlas_subject],'transforms')
        MALFWF.connect( inputsSpec, 'subj_t1_image',
                        labelMapResample[malf_atlas_subject],'reference_image')
        MALFWF.connect( MALF_DG[malf_atlas_subject], 'malf_atlas_lbls',
                        labelMapResample[malf_atlas_subject],'input_image')


        MALFWF.connect(labelMapResample[malf_atlas_subject],'output_image',warpedAtlasLblMergeNode,'in'+str(malf_atlas_mergeindex) )

        ### New labelmap resampling
        NewlabelMapResample[malf_atlas_subject] = pe.Node(interface=ants.ApplyTransforms(),name="FSWM_WLABEL_"+malf_atlas_subject)
        many_cpu_NewlabelMapResample_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,1,1,1), 'overwrite': True}
        NewlabelMapResample[malf_atlas_subject].plugin_args = many_cpu_NewlabelMapResample_options_dictionary
        NewlabelMapResample[malf_atlas_subject].inputs.dimension=3
        NewlabelMapResample[malf_atlas_subject].inputs.output_image=malf_atlas_subject+'fswm_2_subj_lbl.nii.gz'
        NewlabelMapResample[malf_atlas_subject].inputs.interpolation='MultiLabel'
        NewlabelMapResample[malf_atlas_subject].inputs.default_value=0
        NewlabelMapResample[malf_atlas_subject].inputs.invert_transform_flags=[False]

        MALFWF.connect( A2SantsRegistrationPreABCSyN[malf_atlas_subject],'composite_transform',
                        NewlabelMapResample[malf_atlas_subject],'transforms')
        MALFWF.connect( inputsSpec, 'subj_t1_image',
                        NewlabelMapResample[malf_atlas_subject],'reference_image')
        MALFWF.connect( MALF_DG[malf_atlas_subject], 'malf_fswm_atlas_lbls',
                        NewlabelMapResample[malf_atlas_subject],'input_image')


        MALFWF.connect(NewlabelMapResample[malf_atlas_subject],'output_image',NewwarpedAtlasLblMergeNode,'in'+str(malf_atlas_mergeindex) )

        malf_atlas_mergeindex += 1

    ## Now work on cleaning up the label maps
    from FixLabelMapsTools import FixLabelMapFromNeuromorphemetrics2012
    from FixLabelMapsTools import RecodeLabelMap

    ### Original NeuroMorphometrica merged fusion
    jointFusion = pe.Node(interface=ants.JointFusion(),name="JointFusion")
    many_cpu_JointFusion_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,8,4,4), 'overwrite': True}
    jointFusion.plugin_args = many_cpu_JointFusion_options_dictionary
    jointFusion.inputs.dimension=3
    jointFusion.inputs.num_modalities=1
    jointFusion.inputs.method='Joint[0.1,2]'
    jointFusion.inputs.output_label_image='fusion_neuro2012_20.nii.gz'

    MALFWF.connect(warpedAtlasT1MergeNode,'out',jointFusion,'warped_intensity_images')
    MALFWF.connect(warpedAtlasLblMergeNode,'out',jointFusion,'warped_label_images')
    MALFWF.connect(inputsSpec, 'subj_t1_image',jointFusion,'target_image')

    NEUROLABELS_DICT = { 'BRAINSTEM': 35, 'RH_CSF': 51, 'LH_CSF': 52, 'BLOOD': 230 , 'UNKNOWN': 255 ,
    'CONNECTED': [36,37,57,58,55,56,59,60,47,48,23,30]
    }
    fixFusionLabelMap = pe.Node(Function(function=FixLabelMapFromNeuromorphemetrics2012,
                                                   input_names=['fusionFN','FixedHeadFN','LeftHemisphereFN','outFN',
                                                   'OUT_DICT'],
                                                   output_names=['fixedFusionLabelFN']), name="FixedFusionLabelmap")
    fixFusionLabelMap.inputs.outFN = 'neuro2012_20fusion_merge_seg.nii.gz'
    fixFusionLabelMap.inputs.OUT_DICT = NEUROLABELS_DICT
    MALFWF.connect(jointFusion, 'output_label_image', fixFusionLabelMap, 'fusionFN')
    MALFWF.connect(inputsSpec, 'subj_fixed_head_labels', fixFusionLabelMap, 'FixedHeadFN')
    MALFWF.connect(inputsSpec, 'subj_left_hemisphere', fixFusionLabelMap, 'LeftHemisphereFN')

    MALFWF.connect(fixFusionLabelMap,'fixedFusionLabelFN',outputsSpec,'MALF_neuro2012_labelmap')

    ## 2014-02-19 Updated fs_wmparcelation_improved malf
    newJointFusion = pe.Node(interface=ants.JointFusion(),name="FSWM_JointFusion")
    many_cpu_JointFusion_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,8,4,4), 'overwrite': True}
    newJointFusion.plugin_args = many_cpu_JointFusion_options_dictionary
    newJointFusion.inputs.dimension=3
    newJointFusion.inputs.num_modalities=1
    newJointFusion.inputs.method='Joint[0.1,2]'
    newJointFusion.inputs.output_label_image='fswm_neuro2012_20.nii.gz'

    MALFWF.connect(warpedAtlasT1MergeNode,'out',newJointFusion,'warped_intensity_images')
    MALFWF.connect(NewwarpedAtlasLblMergeNode,'out',newJointFusion,'warped_label_images')
    MALFWF.connect(inputsSpec, 'subj_t1_image',newJointFusion,'target_image')

    FREESURFER_DICT = { 'BRAINSTEM': 16, 'RH_CSF':24, 'LH_CSF':24, 'BLOOD': 15000, 'UNKNOWN': 999,
                        'CONNECTED': [11,12,13,9,17,26,50,51,52,48,53,58]
    }
    injectSurfaceCSFandVBIntoLabelMap = pe.Node(Function(function=FixLabelMapFromNeuromorphemetrics2012,
                                                   input_names=['fusionFN','FixedHeadFN','LeftHemisphereFN','outFN',
                                                   'OUT_DICT'],
                                                   output_names=['fixedFusionLabelFN']), name="injectSurfaceCSFandVBIntoLabelMap")
    injectSurfaceCSFandVBIntoLabelMap.inputs.outFN = 'fswm_neuro2012_20_merge_seg.nii.gz'
    injectSurfaceCSFandVBIntoLabelMap.inputs.OUT_DICT = FREESURFER_DICT
    MALFWF.connect(newJointFusion, 'output_label_image', injectSurfaceCSFandVBIntoLabelMap, 'fusionFN')
    MALFWF.connect(inputsSpec, 'subj_fixed_head_labels', injectSurfaceCSFandVBIntoLabelMap, 'FixedHeadFN')
    MALFWF.connect(inputsSpec, 'subj_left_hemisphere', injectSurfaceCSFandVBIntoLabelMap, 'LeftHemisphereFN')

    ## We need to recode values to ensure that there are no conflicts in the future
    RECODE_LABELS_2_Extended_FSWM = [
                                (7071,15071),(7072,15072),(7073,15073),(7145,15145),(7157,15157),
                                (7161,15161),(7179,15179),(7141,15141),(7151,15151),(7163,15163),
                                (7165,15165),(7143,15143),(7191,15191),(7193,15193),(7185,15185),
                                (7201,15201),(7175,15175),(7195,15195),(7173,15173),(7144,15144),
                                (7156,15156),(7160,15160),(7178,15178),(7140,15140),(7150,15150),
                                (7162,15162),(7164,15164),(7142,15142),(7190,15190),(7192,15192),
                                (7184,15184),(7174,15174),(7194,15194),(7172,15172)]
    ## def RecodeLabelMap(InputFileName,OutputFileName,RECODE_TABLE):
    RecodeToExtended = pe.Node(Function(function=RecodeLabelMap,
                                                   input_names=['InputFileName','OutputFileName','RECODE_TABLE'],
                                                   output_names=['OutputFileName']),
                                                   name="RecodeToExteneded")
    RecodeToExtended.inputs.RECODE_TABLE = RECODE_LABELS_2_Extended_FSWM
    RecodeToExtended.inputs.OutputFileName = 'fswm_extended_neuro2012_20_merge_seg.nii.gz'
    MALFWF.connect(injectSurfaceCSFandVBIntoLabelMap, 'fixedFusionLabelFN',RecodeToExtended,'InputFileName')

    ## We need to recode values to ensure that the labels match FreeSurer as close as possible by merging
    ## some labels together to standard FreeSurfer confenventions (i.e. for WMQL)
    RECODE_LABELS_2_Standard_FSWM = [
                                (15071,47),(15072,47),(15073,47),(15145,1011),(15157,1011),(15161,1011),
                                (15179,1012),(15141,1014),(15151,1017),(15163,1018),(15165,1019),(15143,1027),
                                (15191,1028),(15193,1028),(15185,1030),(15201,1030),(15175,1031),(15195,1031),
                                (15173,1035),(15144,2011),(15156,2011),(15160,2011),(15178,2012),(15140,2014),
                                (15150,2017),(15162,2018),(15164,2019),(15142,2027),(15190,2028),(15192,2028),
                                (15184,2030),(15174,2031),(15194,2031),(15172,2035)]
    ## def RecodeLabelMap(InputFileName,OutputFileName,RECODE_TABLE):
    RecodeToStandardFSWM = pe.Node(Function(function=RecodeLabelMap,
                                                   input_names=['InputFileName','OutputFileName','RECODE_TABLE'],
                                                   output_names=['OutputFileName']),
                                                   name="RecodeToStandardFSWM")
    RecodeToStandardFSWM.inputs.RECODE_TABLE = RECODE_LABELS_2_Standard_FSWM
    RecodeToStandardFSWM.inputs.OutputFileName = 'fswm_standard_neuro2012_20_merge_seg.nii.gz'
    MALFWF.connect(RecodeToExtended, 'OutputFileName',RecodeToStandardFSWM,'InputFileName')


    MALFWF.connect(RecodeToExtended,'OutputFileName',outputsSpec,'MALF_fswm_extended_neuro2012_labelmap')
    MALFWF.connect(RecodeToStandardFSWM,'OutputFileName',outputsSpec,'MALF_fswm_standard_neuro2012_labelmap')

    ## MALF_SNAPSHOT_WRITER for Segmented result checking:
    MALF_SNAPSHOT_WRITERNodeName = "MALF_ExtendedMALF_SNAPSHOT_WRITER"
    MALF_SNAPSHOT_WRITER = pe.Node(interface=BRAINSSnapShotWriter(), name=MALF_SNAPSHOT_WRITERNodeName)

    MALF_SNAPSHOT_WRITER.inputs.outputFilename = 'fswm_extended_neuro2012_labelmap.png'  # output specification
    MALF_SNAPSHOT_WRITER.inputs.inputPlaneDirection = [2, 1, 1, 1, 1, 0, 0]
    MALF_SNAPSHOT_WRITER.inputs.inputSliceToExtractInPhysicalPoint = [-3, -7, -3, 5, 7, 22, -22]

    MALFWF.connect([(inputsSpec, MALF_SNAPSHOT_WRITER, [( 'subj_t1_image','inputVolumes')]),
                    (RecodeToExtended, MALF_SNAPSHOT_WRITER, [('OutputFileName', 'inputBinaryVolumes')])
                   ])
    MALFWF.connect(MALF_SNAPSHOT_WRITER,'outputFilename',outputsSpec,'MALF_extended_snapshot')

    return MALFWF
