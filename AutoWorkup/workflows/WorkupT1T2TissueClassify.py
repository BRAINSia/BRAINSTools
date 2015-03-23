#!/usr/bin/env python

from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine

from BRAINSABCext import *

from utilities.misc import *
from utilities.distributed import modify_qsub_args


"""
    from WorkupT1T2TissueClassify import CreateTissueClassifyWorkflow
    myLocalTCWF= CreateTissueClassifyWorkflow("TissueClassify")
    tissueClassifyWF.connect( [ (uidSource, myLocalTCWF, [(('uid', getT1s, subjectDatabaseFile ), 'T1List')] ), ])
    tissueClassifyWF.connect( [ (uidSource, myLocalTCWF, [(('uid', getT2s, subjectDatabaseFile ), 'T2List')] ), ])
    tissueClassifyWF.connect( [ (uidSource, myLocalTCWF, [(('uid', getT1sLength, subjectDatabaseFile ), 'T1_count')] ), ])
    tissueClassifyWF.connect( BCD,    'outputResampledVolume', myLocalTCWF, 'PrimaryT1' )
    tissueClassifyWF.connect(BAtlas,'ExtendedAtlasDefinition.xml',myLocalTCWF,'atlasDefinition')
    tissueClassifyWF.connect(BLI,'outputTransformFilename',myLocalTCWF,'atlasToSubjectInitialTransform')
"""



def getListIndexOrNoneIfOutOfRange(imageList, index):
    if index < len(imageList):
        return imageList[index]
    else:
        return None


def MakePosteriorDictionaryFunc(posteriorImages):
    from PipeLineFunctionHelpers import POSTERIORS
    if len(POSTERIORS) != len(posteriorImages):
        print "ERROR: ", posteriorNames
        print "ERROR: ", POSTERIORS
        return -1
    temp_dictionary = dict(zip(POSTERIORS, posteriorImages))
    return temp_dictionary


def CreateTissueClassifyWorkflow(WFname, master_config, InterpolationMode,UseRegistrationMasking):
    from nipype.interfaces import ants

    CLUSTER_QUEUE=master_config['queue']
    CLUSTER_QUEUE_LONG=master_config['long_q']

    tissueClassifyWF = pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['T1List', 'T2List', 'PDList', 'FLList',
                                                             'OtherList', 'T1_count', 'PrimaryT1',
                                                             'atlasDefinition',
                                                             'atlasToSubjectInitialTransform','atlasVolume'
                                                            ]),
                         run_without_submitting=True,
                         name='inputspec')
    outputsSpec = pe.Node(interface=IdentityInterface(fields=['atlasToSubjectTransform',
                                                              'atlasToSubjectInverseTransform',
                                                              'atlasToSubjectRegistrationState',
                                                              'outputLabels',
                                                              'outputHeadLabels',  # ???
                                                              #'t1_corrected', 't2_corrected',
                                                              't1_average',
                                                              't2_average',
                                                              'pd_average',
                                                              'fl_average',
                                                              'posteriorImages',
                                                              ]),
                          run_without_submitting=True,
                          name='outputspec')


    ########################################################
    # Run BABCext on Multi-modal images
    ########################################################
    makeOutImageList = pe.Node(Function(function=MakeOutFileList,
                                        input_names=['T1List', 'T2List', 'PDList', 'FLList',
                                                     'OtherList','postfix','PrimaryT1'],
                                        output_names=['inImageList','outImageList','imageTypeList']),
                                        run_without_submitting=True, name="99_makeOutImageList")
    tissueClassifyWF.connect(inputsSpec, 'T1List', makeOutImageList, 'T1List')
    tissueClassifyWF.connect(inputsSpec, 'T2List', makeOutImageList, 'T2List')
    tissueClassifyWF.connect(inputsSpec, 'PDList', makeOutImageList, 'PDList')
    tissueClassifyWF.connect(inputsSpec, 'PrimaryT1', makeOutImageList, 'PrimaryT1')
    makeOutImageList.inputs.FLList = []  # an emptyList HACK
    makeOutImageList.inputs.postfix = "_corrected.nii.gz"
    # HACK tissueClassifyWF.connect( inputsSpec, 'FLList', makeOutImageList, 'FLList' )
    tissueClassifyWF.connect(inputsSpec, 'OtherList', makeOutImageList, 'OtherList')

    ##### Initialize with ANTS Transform For AffineComponentBABC
    currentAtlasToSubjectantsRigidRegistration = 'AtlasToSubjectANTsPreABC_Rigid'
    A2SantsRegistrationPreABCRigid = pe.Node(interface=ants.Registration(), name=currentAtlasToSubjectantsRigidRegistration)
    many_cpu_ANTsRigid_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,4,2,8), 'overwrite': True}
    A2SantsRegistrationPreABCRigid.plugin_args = many_cpu_ANTsRigid_options_dictionary

    A2SantsRegistrationPreABCRigid.inputs.num_threads   = -1
    A2SantsRegistrationPreABCRigid.inputs.dimension = 3
    A2SantsRegistrationPreABCRigid.inputs.transforms = ["Affine",]
    A2SantsRegistrationPreABCRigid.inputs.transform_parameters = [[0.1]]
    A2SantsRegistrationPreABCRigid.inputs.metric = ['MI']
    A2SantsRegistrationPreABCRigid.inputs.sampling_strategy = ['Regular']
    A2SantsRegistrationPreABCRigid.inputs.sampling_percentage = [0.5]
    A2SantsRegistrationPreABCRigid.inputs.metric_weight = [1.0]
    A2SantsRegistrationPreABCRigid.inputs.radius_or_number_of_bins = [32]
    A2SantsRegistrationPreABCRigid.inputs.number_of_iterations = [[1000,1000, 500, 100]]

    A2SantsRegistrationPreABCRigid.inputs.convergence_threshold = [1e-8]

    A2SantsRegistrationPreABCRigid.inputs.convergence_window_size = [10]
    A2SantsRegistrationPreABCRigid.inputs.use_histogram_matching = [True]
    A2SantsRegistrationPreABCRigid.inputs.shrink_factors = [[8, 4, 2, 1]]
    A2SantsRegistrationPreABCRigid.inputs.smoothing_sigmas = [[3, 2, 1, 0]]
    A2SantsRegistrationPreABCRigid.inputs.sigma_units = ["vox"]
    A2SantsRegistrationPreABCRigid.inputs.use_estimate_learning_rate_once = [False]
    A2SantsRegistrationPreABCRigid.inputs.write_composite_transform = True  # Required for initialize_transforms_per_stage
    A2SantsRegistrationPreABCRigid.inputs.collapse_output_transforms = False # Mutually Exclusive with initialize_transforms_per_stage
    A2SantsRegistrationPreABCRigid.inputs.initialize_transforms_per_stage = True
    A2SantsRegistrationPreABCRigid.inputs.output_transform_prefix = 'AtlasToSubjectPreBABC_Rigid'
    A2SantsRegistrationPreABCRigid.inputs.winsorize_lower_quantile = 0.01
    A2SantsRegistrationPreABCRigid.inputs.winsorize_upper_quantile = 0.99
    A2SantsRegistrationPreABCRigid.inputs.output_warped_image = 'atlas2subjectRigid.nii.gz'
    A2SantsRegistrationPreABCRigid.inputs.output_inverse_warped_image = 'subject2atlasRigid.nii.gz'
    A2SantsRegistrationPreABCRigid.inputs.float = True

    tissueClassifyWF.connect(inputsSpec, 'atlasToSubjectInitialTransform',A2SantsRegistrationPreABCRigid,'initial_moving_transform')
    tissueClassifyWF.connect(inputsSpec, 'PrimaryT1',A2SantsRegistrationPreABCRigid,'fixed_image')
    tissueClassifyWF.connect(inputsSpec, 'atlasVolume',A2SantsRegistrationPreABCRigid,'moving_image')


    ##### Initialize with ANTS Transform For SyN component BABC
    currentAtlasToSubjectantsRegistration = 'AtlasToSubjectANTsPreABC_SyN'
    A2SantsRegistrationPreABCSyN = pe.Node(interface=ants.Registration(), name=currentAtlasToSubjectantsRegistration)
    many_cpu_ANTsSyN_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE_LONG,8,8,12), 'overwrite': True}
    A2SantsRegistrationPreABCSyN.plugin_args = many_cpu_ANTsSyN_options_dictionary

    A2SantsRegistrationPreABCSyN.inputs.num_threads   = -1
    A2SantsRegistrationPreABCSyN.inputs.dimension = 3
    A2SantsRegistrationPreABCSyN.inputs.transforms = ["SyN","SyN"]
    A2SantsRegistrationPreABCSyN.inputs.transform_parameters = [[0.1, 3, 0],[0.1, 3, 0]]
    A2SantsRegistrationPreABCSyN.inputs.metric = ['CC','CC']
    A2SantsRegistrationPreABCSyN.inputs.sampling_strategy = [None,None]
    A2SantsRegistrationPreABCSyN.inputs.sampling_percentage = [1.0,1.0]
    A2SantsRegistrationPreABCSyN.inputs.metric_weight = [1.0,1.0]
    A2SantsRegistrationPreABCSyN.inputs.radius_or_number_of_bins = [4,4]
    A2SantsRegistrationPreABCSyN.inputs.number_of_iterations = [[500, 500], [500, 70]]

    A2SantsRegistrationPreABCSyN.inputs.convergence_threshold = [1e-8,1e-6]

    A2SantsRegistrationPreABCSyN.inputs.convergence_window_size = [12]
    A2SantsRegistrationPreABCSyN.inputs.use_histogram_matching = [True,True]
    A2SantsRegistrationPreABCSyN.inputs.shrink_factors = [[8, 4], [2, 1]]
    A2SantsRegistrationPreABCSyN.inputs.smoothing_sigmas = [[3, 2], [1, 0]]
    A2SantsRegistrationPreABCSyN.inputs.sigma_units = ["vox","vox"]
    A2SantsRegistrationPreABCSyN.inputs.use_estimate_learning_rate_once = [False,False]
    A2SantsRegistrationPreABCSyN.inputs.write_composite_transform = True # Required for initialize_transforms_per_stage
    A2SantsRegistrationPreABCSyN.inputs.collapse_output_transforms = False # Mutually Exclusive with initialize_transforms_per_stage
    A2SantsRegistrationPreABCSyN.inputs.initialize_transforms_per_stage = True
    A2SantsRegistrationPreABCSyN.inputs.save_state = 'SavedInternalSyNState.h5'
    A2SantsRegistrationPreABCSyN.inputs.output_transform_prefix = 'AtlasToSubjectPreBABC_SyN'
    A2SantsRegistrationPreABCSyN.inputs.winsorize_lower_quantile = 0.01
    A2SantsRegistrationPreABCSyN.inputs.winsorize_upper_quantile = 0.99
    A2SantsRegistrationPreABCSyN.inputs.output_warped_image = 'atlas2subject.nii.gz'
    A2SantsRegistrationPreABCSyN.inputs.output_inverse_warped_image = 'subject2atlas.nii.gz'
    A2SantsRegistrationPreABCSyN.inputs.float = True

    ## if using Registration masking, then do ROIAuto on fixed and moving images and connect to registraitons
    if UseRegistrationMasking == True:
        from SEMTools.segmentation.specialized import BRAINSROIAuto

        fixedROIAuto = pe.Node(interface=BRAINSROIAuto(), name="fixedImageROIAUTOMask")
        fixedROIAuto.inputs.ROIAutoDilateSize=10
        fixedROIAuto.inputs.outputROIMaskVolume = "fixedImageROIAutoMask.nii.gz"

        movingROIAuto = pe.Node(interface=BRAINSROIAuto(), name="movingImageROIAUTOMask")
        fixedROIAuto.inputs.ROIAutoDilateSize=10
        movingROIAuto.inputs.outputROIMaskVolume = "movingImageROIAutoMask.nii.gz"

        tissueClassifyWF.connect(inputsSpec, 'PrimaryT1',fixedROIAuto,'inputVolume')
        tissueClassifyWF.connect(inputsSpec, 'atlasVolume',movingROIAuto,'inputVolume')

        tissueClassifyWF.connect(fixedROIAuto, 'outputROIMaskVolume',A2SantsRegistrationPreABCRigid,'fixed_image_mask')
        tissueClassifyWF.connect(movingROIAuto, 'outputROIMaskVolume',A2SantsRegistrationPreABCRigid,'moving_image_mask')

        tissueClassifyWF.connect(fixedROIAuto, 'outputROIMaskVolume',A2SantsRegistrationPreABCSyN,'fixed_image_mask')
        tissueClassifyWF.connect(movingROIAuto, 'outputROIMaskVolume',A2SantsRegistrationPreABCSyN,'moving_image_mask')

    tissueClassifyWF.connect(A2SantsRegistrationPreABCRigid,
                             ('composite_transform', getListIndexOrNoneIfOutOfRange, 0 ),
                             A2SantsRegistrationPreABCSyN,'initial_moving_transform')
    tissueClassifyWF.connect(inputsSpec, 'PrimaryT1',A2SantsRegistrationPreABCSyN,'fixed_image')
    tissueClassifyWF.connect(inputsSpec, 'atlasVolume',A2SantsRegistrationPreABCSyN,'moving_image')

    BABCext = pe.Node(interface=BRAINSABCext(), name="BABC")
    many_cpu_BABC_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,8,2,4), 'overwrite': True}
    BABCext.plugin_args = many_cpu_BABC_options_dictionary
    tissueClassifyWF.connect(makeOutImageList, 'inImageList', BABCext, 'inputVolumes')
    tissueClassifyWF.connect(makeOutImageList, 'imageTypeList', BABCext, 'inputVolumeTypes')
    tissueClassifyWF.connect(makeOutImageList, 'outImageList', BABCext, 'outputVolumes')
    BABCext.inputs.debuglevel = 0
    BABCext.inputs.useKNN = True
    BABCext.inputs.maxIterations = 3
    BABCext.inputs.maxBiasDegree = 4
    BABCext.inputs.filterIteration = 3
    BABCext.inputs.filterMethod = 'GradientAnisotropicDiffusion'
    BABCext.inputs.atlasToSubjectTransformType = 'SyN'
    BABCext.inputs.gridSize = [10, 10, 10]
    BABCext.inputs.outputFormat = "NIFTI"
    BABCext.inputs.outputLabels = "brain_label_seg.nii.gz"
    BABCext.inputs.outputDirtyLabels = "volume_label_seg.nii.gz"
    BABCext.inputs.posteriorTemplate = "POSTERIOR_%s.nii.gz"
    BABCext.inputs.atlasToSubjectTransform = "atlas_to_subject.h5"
    # BABCext.inputs.implicitOutputs = ['t1_average_BRAINSABC.nii.gz', 't2_average_BRAINSABC.nii.gz']
    BABCext.inputs.interpolationMode = InterpolationMode
    BABCext.inputs.outputDir = './'
    BABCext.inputs.saveState = 'SavedBABCInternalSyNState.h5'

    tissueClassifyWF.connect(inputsSpec, 'atlasDefinition', BABCext, 'atlasDefinition')
    # NOTE: MUTUALLY EXCLUSIVE with restoreState
    #tissueClassifyWF.connect(A2SantsRegistrationPreABCSyN,
    #                         ( 'composite_transform', getListIndexOrNoneIfOutOfRange, 0 ),
    #                         BABCext, 'atlasToSubjectInitialTransform')
    tissueClassifyWF.connect(A2SantsRegistrationPreABCSyN,'save_state',
                             BABCext, 'restoreState')

    """
    Get the first T1 and T2 corrected images from BABCext
    """

    """ HACK:  THIS IS NOT NEEDED!  We should use the averged t1 and averaged t2 images instead!
    def get_first_T1_and_T2(in_files,T1_count):
        '''
        Returns the first T1 and T2 file in in_files, based on offset in T1_count.
        '''
        return in_files[0],in_files[T1_count]
    bfc_files = pe.Node(Function(input_names=['in_files','T1_count'],
                               output_names=['t1_corrected','t2_corrected'],
                               function=get_first_T1_and_T2), run_without_submitting=True, name='99_bfc_files' )
    tissueClassifyWF.connect( inputsSpec, 'T1_count', bfc_files, 'T1_count')
    tissueClassifyWF.connect(BABCext,'outputVolumes',bfc_files, 'in_files')


    tissueClassifyWF.connect(bfc_files,'t1_corrected',outputsSpec,'t1_corrected')
    tissueClassifyWF.connect(bfc_files,'t2_corrected',outputsSpec,'t2_corrected')
    #tissueClassifyWF.connect(bfc_files,'pd_corrected',outputsSpec,'pd_corrected')
    #tissueClassifyWF.connect(bfc_files,'fl_corrected',outputsSpec,'fl_corrected')

    """

    #############
    tissueClassifyWF.connect(BABCext, 'saveState', outputsSpec, 'atlasToSubjectRegistrationState')

    tissueClassifyWF.connect(BABCext, 'atlasToSubjectTransform', outputsSpec, 'atlasToSubjectTransform')


    def MakeInverseTransformFileName(TransformFileName):
        """### HACK:  This function is to work around a deficiency in BRAINSABCext where the inverse transform name is not being computed properly
          in the list outputs"""
        fixed_inverse_name = TransformFileName.replace(".h5", "_Inverse.h5")
        return [fixed_inverse_name]

    tissueClassifyWF.connect([(BABCext, outputsSpec, [(('atlasToSubjectTransform', MakeInverseTransformFileName), "atlasToSubjectInverseTransform")]), ])
    tissueClassifyWF.connect(BABCext, 'outputLabels', outputsSpec, 'outputLabels')
    tissueClassifyWF.connect(BABCext, 'outputDirtyLabels', outputsSpec, 'outputHeadLabels')

    tissueClassifyWF.connect(BABCext, 'outputT1AverageImage', outputsSpec, 't1_average')
    tissueClassifyWF.connect(BABCext, 'outputT2AverageImage', outputsSpec, 't2_average')
    tissueClassifyWF.connect(BABCext, 'outputPDAverageImage', outputsSpec, 'pd_average')
    tissueClassifyWF.connect(BABCext, 'outputFLAverageImage', outputsSpec, 'fl_average')
    ##  remove tissueClassifyWF.connect( [ ( BABCext, outputsSpec, [ (( 'outputAverageImages', getListIndexOrNoneIfOutOfRange, 0 ), "t1_average")] ), ] )
    ##  remove tissueClassifyWF.connect( [ ( BABCext, outputsSpec, [ (( 'outputAverageImages', getListIndexOrNoneIfOutOfRange, 1 ), "t2_average")] ), ] )
    ##  remove tissueClassifyWF.connect( [ ( BABCext, outputsSpec, [ (( 'outputAverageImages', getListIndexOrNoneIfOutOfRange, 2 ), "pd_average")] ), ] )

    MakePosteriorDictionaryNode = pe.Node(Function(function=MakePosteriorDictionaryFunc,
                                                   input_names=['posteriorImages'],
                                                   output_names=['posteriorDictionary']), run_without_submitting=True, name="99_makePosteriorDictionary")
    tissueClassifyWF.connect(BABCext, 'posteriorImages', MakePosteriorDictionaryNode, 'posteriorImages')

    tissueClassifyWF.connect(MakePosteriorDictionaryNode, 'posteriorDictionary', outputsSpec, 'posteriorImages')

    return tissueClassifyWF
