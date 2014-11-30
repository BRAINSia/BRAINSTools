#!/usr/bin/env python

from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine

from utilities.misc import *
from utilities.distributed import modify_qsub_args
from SEMTools.utilities.brains import BRAINSLandmarkInitializer

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
    outputsSpec = pe.Node(interface=IdentityInterface(fields=['MALF_neuro2012_labelmap']),
                          run_without_submitting=True,
                          name='outputspec')

    BLICreator = dict()
    MALF_DG = dict()
    A2SantsRegistrationPreABCRigid =dict()
    A2SantsRegistrationPreABCSyN = dict()
    fixedROIAuto = dict()
    movingROIAuto = dict()
    labelMapResample = dict()

    warpedAtlasT1MergeNode = pe.Node(interface=Merge(len(good_subjects)),name="T1sMergeAtlas")
    warpedAtlasLblMergeNode = pe.Node(interface=Merge(len(good_subjects)),name="LblMergeAtlas")
    malf_atlas_mergeindex = 1;
    for malf_atlas_subject in good_subjects:
        ## Need DataGrabber Here For the Atlas
        MALF_DG[malf_atlas_subject] = pe.Node(interface=nio.DataGrabber(infields=['subject'],
                                                        outfields=['malf_atlas_t1',
                                                                   'malf_atlas_lbls',
                                                                   'malf_atlas_lmks'
                                                        ]),
                              run_without_submitting=True,name='MALF_DG_'+malf_atlas_subject)
        #MALF_DG[malf_atlas_subject].inputs.base_directory = master_config['previousresult']
        MALF_DG[malf_atlas_subject].inputs.base_directory = BASE_DATA_GRABBER_DIR

        MALF_DG[malf_atlas_subject].inputs.subject = malf_atlas_subject
        MALF_DG[malf_atlas_subject].inputs.field_template = {
                                             'malf_atlas_t1': '%s/TissueClassify/t1_average_BRAINSABC.nii.gz',
                                             'malf_atlas_lbls': '%s/TissueClassify/neuro_lbls.nii.gz',
                                             'malf_atlas_lmks': '%s/ACPCAlign/BCD_ACPC_Landmarks.fcsv',
        }
        MALF_DG[malf_atlas_subject].inputs.template_args = {
                                            'malf_atlas_t1':   [['subject']],
                                            'malf_atlas_lbls': [['subject']],
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
        many_cpu_ANTsSyN_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE_LONG,4,2,4), 'overwrite': True}
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


        MALFWF.connect(A2SantsRegistrationPreABCSyN[malf_atlas_subject],'warped_image',warpedAtlasT1MergeNode,'in'+str(malf_atlas_mergeindex) )
        MALFWF.connect(labelMapResample[malf_atlas_subject],'output_image',warpedAtlasLblMergeNode,'in'+str(malf_atlas_mergeindex) )
        malf_atlas_mergeindex += 1

    jointFusion = pe.Node(interface=ants.JointFusion(),name="JointFusion")
    many_cpu_JointFusion_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,8,4,4), 'overwrite': True}
    jointFusion.plugin_args = many_cpu_JointFusion_options_dictionary
    jointFusion.inputs.dimension=3
    jointFusion.inputs.num_modalities=1
    jointFusion.inputs.method='Joint[0.1,2]'
    jointFusion.inputs.output_label_image='fusion_neuro2012_20.nii.gz'

    def FixLabelMapFromNeuromorphemetrics2012(fusionFN,FixedHeadFN,LeftHemisphereFN,outFN):
        import SimpleITK as sitk
        import os

        def ForceMaskInsert(inlabels,newmask,newmaskvalue):
            newmask = sitk.Cast( (newmask>0) , sitk.sitkUInt8)
            outlabels=inlabels*sitk.Cast( (1-newmask), sitk.sitkUInt8)
            outlabels = outlabels + newmask*newmaskvalue
            return sitk.Cast(outlabels,sitk.sitkUInt8)
        ## TODO: GetLargestLabel is copied from elsewhere
        def GetLargestLabel(inputMask, UseErosionCleaning):
            LargestComponentCode = 1
            if UseErosionCleaning:
                erosionMask = sitk.ErodeObjectMorphology(inputMask, 1)
            else:
                erosionMask = inputMask
            CC = sitk.ConnectedComponent(erosionMask)
            Rlabel = sitk.RelabelComponent(CC)
            largestMask = ( Rlabel == LargestComponentCode)
            if UseErosionCleaning:
                dilateMask = sitk.DilateObjectMorphology(largestMask, 1)
            else:
                dilateMask = largestMask

            return (largestMask * dilateMask > 0)

        def RecodeNonLargest(outlabels,keepCode,UNKNOWN_LABEL_CODE):
            orig_mask = (outlabels ==  keepCode)
            connected_mask = GetLargestLabel(orig_mask,False)
            small_regions = ( orig_mask - connected_mask )
            outlabels = ForceMaskInsert(outlabels,connected_mask,keepCode)
            outlabels = ForceMaskInsert(outlabels,small_regions,UNKNOWN_LABEL_CODE)
            return outlabels

        fusionIm=sitk.Cast(sitk.ReadImage(fusionFN),sitk.sitkUInt8)
        FixedHead=sitk.Cast(sitk.ReadImage(FixedHeadFN),sitk.sitkUInt8)
        LeftHemisphereIm=sitk.Cast(sitk.ReadImage(LeftHemisphereFN),sitk.sitkUInt8)

        csf_labels=(FixedHead == 4)
        outlabels= ForceMaskInsert(fusionIm,csf_labels,51)
        blood_labels=(FixedHead == 5)
        BLOOD_CODE=230
        outlabels = ForceMaskInsert(outlabels,blood_labels,BLOOD_CODE)  ## Add blood as value 230
        left_hemi_pre = ( outlabels == 52 )
        outlabels = ForceMaskInsert(outlabels,left_hemi_pre,51)  ## Make all CSF Right hemisphere
        left_hemi_post =  (LeftHemisphereIm * sitk.Cast ( ( outlabels == 51 ),sitk.sitkUInt8) > 0 )
        outlabels = ForceMaskInsert(outlabels,left_hemi_post,52)  ## Make all CSF Right hemisphere
        ## Now extend brainstem lower
        brain_stem = (FixedHead == 30) * (outlabels == 0) ## Only extend to areas where there is not already a label
        outlabels = ForceMaskInsert(outlabels,brain_stem,35)  ## Make all CSF Right hemisphere
        BRAIN_MASK=sitk.Cast( (FixedHead > 0),sitk.sitkUInt8)
        outlabels = outlabels * BRAIN_MASK

        ## Caudate = 36 37
        ## Putamen = 57 58
        ## Pallidus = 55,56
        ## Thalamus = 59,60
        ## Hippocampus = 47,48
        ## Accumbens  = 23,30
        UNKNOWN_LABEL_CODE=255
        labels_to_ensure_connected = [36,37,57,58,55,56,59,60,47,48,23,30]
        for keepCode in labels_to_ensure_connected:
            outlabels = RecodeNonLargest(outlabels,keepCode,UNKNOWN_LABEL_CODE)

        ## FILL IN HOLES
        unkown_holes = ( BRAIN_MASK > 0 ) * ( outlabels == 0 )
        outlabels = ForceMaskInsert(outlabels,unkown_holes,UNKNOWN_LABEL_CODE)  ## Fill unkown regeions with unkown code

        fixedFusionLabelFN=os.path.realpath(outFN)
        sitk.WriteImage(outlabels,fixedFusionLabelFN)
        #print("\n\n\n\n\n\n{0}\n\n\n\nXXXXXXXX".format(fixedFusionLabelFN))
        return fixedFusionLabelFN

    fixFusionLabelMap = pe.Node(Function(function=FixLabelMapFromNeuromorphemetrics2012,
                                                   input_names=['fusionFN','FixedHeadFN','LeftHemisphereFN','outFN' ],
                                                   output_names=['fixedFusionLabelFN']), name="FixedFusionLabelmap")
    fixFusionLabelMap.inputs.outFN = 'neuro2012_20fusion_merge_seg.nii.gz'
    MALFWF.connect(jointFusion, 'output_label_image', fixFusionLabelMap, 'fusionFN')
    MALFWF.connect(inputsSpec, 'subj_fixed_head_labels', fixFusionLabelMap, 'FixedHeadFN')
    MALFWF.connect(inputsSpec, 'subj_left_hemisphere', fixFusionLabelMap, 'LeftHemisphereFN')



    MALFWF.connect(warpedAtlasT1MergeNode,'out',jointFusion,'warped_intensity_images')
    MALFWF.connect(warpedAtlasLblMergeNode,'out',jointFusion,'warped_label_images')
    MALFWF.connect(inputsSpec, 'subj_t1_image',jointFusion,'target_image')
    MALFWF.connect(fixFusionLabelMap,'fixedFusionLabelFN',outputsSpec,'MALF_neuro2012_labelmap')

    return MALFWF
