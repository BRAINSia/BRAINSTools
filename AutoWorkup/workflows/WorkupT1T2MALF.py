#!/usr/bin/env python

from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import range
from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine

from utilities.misc import *
from utilities.distributed import modify_qsub_args
from nipype.interfaces.semtools.utilities.brains import BRAINSLandmarkInitializer
# HACK Remove due to bugs from nipype.interfaces.semtools import BRAINSSnapShotWriter

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

def MakeVector(inFN1, inFN2=None):
    #print("inFN1: {0}".format(inFN1))
    #print("inFN2: {0}".format(inFN2))
    if inFN2 == None:
        return inFN1
    else:
        return [inFN1, inFN2]

def readRecodingList( recodeLabelFilename ):
    recodeLabelPairList = []
    import csv
    with open( recodeLabelFilename, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for line in reader:
           if line[0].startswith("#"):
               pass
           else:
               origLbl, origName, targetLbl, targetName= line
               origLbl = int(origLbl)
               targetLbl = int(targetLbl)
               recodeLabelPairList.append( (origLbl, targetLbl) )
               #print(" relabel: {0:30} ({1:5d}) --> {2:30} ({3:5d})".format(origName, origLbl, targetName, targetLbl))
    return recodeLabelPairList

def readMalfAtlasDbBase( dictionaryFilename ):
    malfAtlasDict = {}
    #scanID, ['atlasID', 't1', 't2' ,'label', 'lmks']
    import ast
    with open(dictionaryFilename, 'r') as f:
        for line in f.readlines():
            scanID, scanDict=ast.literal_eval(line)
            malfAtlasDict[scanID]=scanDict

    return malfAtlasDict

def getListIndexOrNoneIfOutOfRange(imageList, index):
    if index < len(imageList):
        return imageList[index]
    else:
        return None


def CreateMALFWorkflow(WFname, onlyT1, master_config,BASE_DATA_GRABBER_DIR=None, runFixFusionLabelMap=True):
    from nipype.interfaces import ants

    if onlyT1:
      n_modality = 1
    else:
      n_modality = 2
    CLUSTER_QUEUE=master_config['queue']
    CLUSTER_QUEUE_LONG=master_config['long_q']

    MALFWF = pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['subj_t1_image', #Desired image to create label map for
                                                             'subj_t2_image', #Desired image to create label map for
                                                             'subj_lmks', #The landmarks corresponding to t1_image
                                                             'subj_fixed_head_labels', #The fixed head labels from BABC
                                                             'subj_left_hemisphere', #The warped left hemisphere mask
                                                             'atlasWeightFilename',  #The static weights file name
                                                             'labelBaseFilename' #Atlas label base name ex) neuro_lbls.nii.gz
                                                            ]),
                         run_without_submitting=True,
                         name='inputspec')
    outputsSpec = pe.Node(interface=IdentityInterface(fields=['MALF_HDAtlas20_2015_label',
                                                       'MALF_HDAtlas20_2015_CSFVBInjected_label',
                                                       'MALF_HDAtlas20_2015_fs_standard_label',
                                                       'MALF_HDAtlas20_2015_lobar_label',
                                                       'MALF_extended_snapshot']),
                          run_without_submitting=True,
                          name='outputspec')

    BLICreator = dict()
    A2SantsRegistrationPreMALF_SyN = dict()
    fixedROIAuto = dict()
    movingROIAuto = dict()
    labelMapResample = dict()
    NewlabelMapResample = dict()

    malf_atlas_mergeindex = 0
    merge_input_offset = 1 #Merge nodes are indexed from 1, not zero!
    """
    multimodal ants registration if t2 exists
    """
    sessionMakeMultimodalInput = pe.Node(Function(function=MakeVector,
                                                                      input_names=['inFN1', 'inFN2'],
                                                                      output_names=['outFNs']),
                                run_without_submitting=True, name="sessionMakeMultimodalInput")
    MALFWF.connect(inputsSpec, 'subj_t1_image', sessionMakeMultimodalInput, 'inFN1')
    if not onlyT1:
        MALFWF.connect(inputsSpec, 'subj_t2_image', sessionMakeMultimodalInput, 'inFN2')
    else:
        pass


    #print('malf_atlas_db_base')
    #print(master_config['malf_atlas_db_base'])
    malfAtlasDict = readMalfAtlasDbBase( master_config['malf_atlas_db_base'] )
    number_of_atlas_sources = len(malfAtlasDict)
    malfAtlases = dict()
    atlasMakeMultimodalInput = dict()
    t2Resample = dict()
    warpedAtlasLblMergeNode = pe.Node(interface=Merge(number_of_atlas_sources),name="LblMergeAtlas")
    NewwarpedAtlasLblMergeNode = pe.Node(interface=Merge(number_of_atlas_sources),name="fswmLblMergeAtlas")
    warpedAtlasesMergeNode = pe.Node(interface=Merge(number_of_atlas_sources*n_modality),name="MergeAtlases")

    for malf_atlas_subject in list(malfAtlasDict.keys()):
        ## Need DataGrabber Here For the Atlas
        malfAtlases[malf_atlas_subject] = pe.Node(interface = IdentityInterface(
                                                                  fields=['t1', 't2', 'label', 'lmks']),
                                                                  name='malfAtlasInput'+malf_atlas_subject)
        malfAtlases[malf_atlas_subject].inputs.t1 = malfAtlasDict[malf_atlas_subject]['t1']
        malfAtlases[malf_atlas_subject].inputs.t2 = malfAtlasDict[malf_atlas_subject]['t2']
        malfAtlases[malf_atlas_subject].inputs.label = malfAtlasDict[malf_atlas_subject]['label']
        malfAtlases[malf_atlas_subject].inputs.lmks = malfAtlasDict[malf_atlas_subject]['lmks']
        ## Create BLI first
        ########################################################
        # Run BLI atlas_to_subject
        ########################################################
        BLICreator[malf_atlas_subject] = pe.Node(interface=BRAINSLandmarkInitializer(), name="BLI_"+malf_atlas_subject)
        BLICreator[malf_atlas_subject].inputs.outputTransformFilename = "landmarkInitializer_{0}_to_subject_transform.h5".format(malf_atlas_subject)

        MALFWF.connect(inputsSpec, 'atlasWeightFilename', BLICreator[malf_atlas_subject], 'inputWeightFilename')
        MALFWF.connect(malfAtlases[malf_atlas_subject], 'lmks', BLICreator[malf_atlas_subject], 'inputMovingLandmarkFilename')
        MALFWF.connect(inputsSpec, 'subj_lmks', BLICreator[malf_atlas_subject], 'inputFixedLandmarkFilename')

        ##### Initialize with ANTS Transform For SyN
        currentAtlasToSubjectantsRegistration = 'SyN_AtlasToSubjectANTsPreMALF_'+malf_atlas_subject
        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject] = pe.Node(interface=ants.Registration(), name=currentAtlasToSubjectantsRegistration)
        many_cpu_ANTsSyN_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE_LONG,4,2,16), 'overwrite': True}
        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].plugin_args = many_cpu_ANTsSyN_options_dictionary

        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.num_threads   = -1
        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.dimension = 3
        #### DEBUGGIN
        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.transforms = ["Affine","Affine","SyN","SyN"]
        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.transform_parameters = [[0.1],[0.1],[0.1, 3, 0],[0.1, 3, 0]]
        if onlyT1:
            A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.metric = ['MI','MI','CC','CC']
            A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.metric_weight = [1.0,1.0,1.0,1.0]
            A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.sampling_percentage = [.5,.5,1.0,1.0]
            A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.radius_or_number_of_bins = [32,32,4,4]
            A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.sampling_strategy = ['Regular','Regular',None,None]
        else:
            A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.metric = ['MI',['MI','MI'],'CC',['CC','CC']]
            A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.metric_weight = [1.0,[1.0,1.0],1.0,[1.0,1.0]]
            A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.sampling_percentage = [.5,[.5,0.5],1.0,[1.0,1.0]]
            A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.radius_or_number_of_bins = [32,[32,32],4,[4,4]]
            A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.sampling_strategy = ['Regular',['Regular','Regular'],None,[None,None]]


        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.number_of_iterations = [[1000,1000,500],[500,500],[500,500],[500,70]]

        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.convergence_threshold = [1e-8,1e-6,1e-8,1e-6]

        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.convergence_window_size = [12]
        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.use_histogram_matching = [True,True,True,True]
        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.shrink_factors = [[8, 4, 2],[2, 1],[8, 4],[2, 1]]
        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.smoothing_sigmas = [[3, 2, 1],[1, 0],[3, 2],[1, 0]]
        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.sigma_units = ["vox","vox","vox","vox"]
        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.use_estimate_learning_rate_once = [False,False,False,False]
        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.write_composite_transform = True # Required for initialize_transforms_per_stage
        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.collapse_output_transforms = False # Mutually Exclusive with initialize_transforms_per_stage
        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.initialize_transforms_per_stage = True
        ## NO NEED FOR THIS A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.save_state = 'SavedInternalSyNState.h5'
        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.output_transform_prefix = malf_atlas_subject+'_ToSubjectPreMALF_SyN'
        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.winsorize_lower_quantile = 0.01
        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.winsorize_upper_quantile = 0.99
        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.output_warped_image = malf_atlas_subject + '_2subject.nii.gz'
        ## NO NEED FOR THIS A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.output_inverse_warped_image = 'subject2atlas.nii.gz'
        A2SantsRegistrationPreMALF_SyN[malf_atlas_subject].inputs.float = True

        ## if using Registration masking, then do ROIAuto on fixed and moving images and connect to registraitons
        UseRegistrationMasking = True
        if UseRegistrationMasking == True:
            from nipype.interfaces.semtools.segmentation.specialized import BRAINSROIAuto

            fixedROIAuto[malf_atlas_subject] = pe.Node(interface=BRAINSROIAuto(), name="fixedROIAUTOMask_"+malf_atlas_subject)
            fixedROIAuto[malf_atlas_subject].inputs.ROIAutoDilateSize=10
            fixedROIAuto[malf_atlas_subject].inputs.outputROIMaskVolume = "fixedImageROIAutoMask.nii.gz"

            movingROIAuto[malf_atlas_subject] = pe.Node(interface=BRAINSROIAuto(), name="movingROIAUTOMask_"+malf_atlas_subject)
            fixedROIAuto[malf_atlas_subject].inputs.ROIAutoDilateSize=10
            movingROIAuto[malf_atlas_subject].inputs.outputROIMaskVolume = "movingImageROIAutoMask.nii.gz"

            MALFWF.connect(inputsSpec, 'subj_t1_image',fixedROIAuto[malf_atlas_subject],'inputVolume')
            MALFWF.connect(malfAtlases[malf_atlas_subject], 't1', movingROIAuto[malf_atlas_subject],'inputVolume')

            MALFWF.connect(fixedROIAuto[malf_atlas_subject], 'outputROIMaskVolume',A2SantsRegistrationPreMALF_SyN[malf_atlas_subject],'fixed_image_mask')
            MALFWF.connect(movingROIAuto[malf_atlas_subject], 'outputROIMaskVolume',A2SantsRegistrationPreMALF_SyN[malf_atlas_subject],'moving_image_mask')

        MALFWF.connect(BLICreator[malf_atlas_subject],'outputTransformFilename',
                       A2SantsRegistrationPreMALF_SyN[malf_atlas_subject],'initial_moving_transform')

        """
        make multimodal input for atlases
        """
        atlasMakeMultimodalInput[malf_atlas_subject] = pe.Node(Function(function=MakeVector, input_names=['inFN1', 'inFN2'], output_names=['outFNs']),
                                  run_without_submitting=True, name="atlasMakeMultimodalInput"+malf_atlas_subject)
        MALFWF.connect(malfAtlases[malf_atlas_subject], 't1', atlasMakeMultimodalInput[malf_atlas_subject], 'inFN1')
        if not onlyT1:
            MALFWF.connect(malfAtlases[malf_atlas_subject], 't2', atlasMakeMultimodalInput[malf_atlas_subject], 'inFN2')
        else:
            pass

        MALFWF.connect(sessionMakeMultimodalInput, 'outFNs',
                       A2SantsRegistrationPreMALF_SyN[malf_atlas_subject],'fixed_image')
        MALFWF.connect(atlasMakeMultimodalInput[malf_atlas_subject], 'outFNs',
                       A2SantsRegistrationPreMALF_SyN[malf_atlas_subject],'moving_image')
        MALFWF.connect(A2SantsRegistrationPreMALF_SyN[malf_atlas_subject],'warped_image',
                       warpedAtlasesMergeNode,'in'+str(merge_input_offset + malf_atlas_mergeindex*n_modality) )

        """
        Original t2 resampling
        """
        for modality_index in range(1,n_modality):
            t2Resample[malf_atlas_subject] = pe.Node(interface=ants.ApplyTransforms(),name="resampledT2"+malf_atlas_subject)
            many_cpu_t2Resample_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,1,1,1), 'overwrite': True}
            t2Resample[malf_atlas_subject].plugin_args = many_cpu_t2Resample_options_dictionary
            t2Resample[malf_atlas_subject].inputs.dimension=3
            t2Resample[malf_atlas_subject].inputs.output_image=malf_atlas_subject+'_t2.nii.gz'
            t2Resample[malf_atlas_subject].inputs.interpolation='BSpline'
            t2Resample[malf_atlas_subject].inputs.default_value=0
            t2Resample[malf_atlas_subject].inputs.invert_transform_flags=[False]

            MALFWF.connect( A2SantsRegistrationPreMALF_SyN[malf_atlas_subject],'composite_transform',
                            t2Resample[malf_atlas_subject],'transforms')
            MALFWF.connect( inputsSpec, 'subj_t1_image',
                            t2Resample[malf_atlas_subject],'reference_image')
            MALFWF.connect( malfAtlases[malf_atlas_subject], 't2',
                            t2Resample[malf_atlas_subject],'input_image')
            MALFWF.connect(t2Resample[malf_atlas_subject],'output_image',
                           warpedAtlasesMergeNode,'in'+str(merge_input_offset + malf_atlas_mergeindex*n_modality+modality_index) )

        """
        Original labelmap resampling
        """
        labelMapResample[malf_atlas_subject] = pe.Node(interface=ants.ApplyTransforms(),name="resampledLabel"+malf_atlas_subject)
        many_cpu_labelMapResample_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,1,1,1), 'overwrite': True}
        labelMapResample[malf_atlas_subject].plugin_args = many_cpu_labelMapResample_options_dictionary
        labelMapResample[malf_atlas_subject].inputs.dimension=3
        labelMapResample[malf_atlas_subject].inputs.output_image=malf_atlas_subject+'_2_subj_lbl.nii.gz'
        labelMapResample[malf_atlas_subject].inputs.interpolation='MultiLabel'
        labelMapResample[malf_atlas_subject].inputs.default_value=0
        labelMapResample[malf_atlas_subject].inputs.invert_transform_flags=[False]

        MALFWF.connect( A2SantsRegistrationPreMALF_SyN[malf_atlas_subject],'composite_transform',
                        labelMapResample[malf_atlas_subject],'transforms')
        MALFWF.connect( inputsSpec, 'subj_t1_image',
                        labelMapResample[malf_atlas_subject],'reference_image')
        MALFWF.connect( malfAtlases[malf_atlas_subject], 'label',
                        labelMapResample[malf_atlas_subject],'input_image')


        MALFWF.connect(labelMapResample[malf_atlas_subject],'output_image',warpedAtlasLblMergeNode,'in'+str(merge_input_offset + malf_atlas_mergeindex) )

        ### New labelmap resampling
        NewlabelMapResample[malf_atlas_subject] = pe.Node(interface=ants.ApplyTransforms(),name="FSWM_WLABEL_"+malf_atlas_subject)
        many_cpu_NewlabelMapResample_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,1,1,1), 'overwrite': True}
        NewlabelMapResample[malf_atlas_subject].plugin_args = many_cpu_NewlabelMapResample_options_dictionary
        NewlabelMapResample[malf_atlas_subject].inputs.dimension=3
        NewlabelMapResample[malf_atlas_subject].inputs.output_image=malf_atlas_subject+'fswm_2_subj_lbl.nii.gz'
        NewlabelMapResample[malf_atlas_subject].inputs.interpolation='MultiLabel'
        NewlabelMapResample[malf_atlas_subject].inputs.default_value=0
        NewlabelMapResample[malf_atlas_subject].inputs.invert_transform_flags=[False]

        MALFWF.connect( A2SantsRegistrationPreMALF_SyN[malf_atlas_subject],'composite_transform',
                        NewlabelMapResample[malf_atlas_subject],'transforms')
        MALFWF.connect( inputsSpec, 'subj_t1_image',
                        NewlabelMapResample[malf_atlas_subject],'reference_image')
        MALFWF.connect( malfAtlases[malf_atlas_subject], 'label',
                        NewlabelMapResample[malf_atlas_subject],'input_image')


        MALFWF.connect(NewlabelMapResample[malf_atlas_subject],'output_image',NewwarpedAtlasLblMergeNode,'in'+str(merge_input_offset + malf_atlas_mergeindex) )

        malf_atlas_mergeindex += 1


    ## Now work on cleaning up the label maps
    from .FixLabelMapsTools import FixLabelMapFromNeuromorphemetrics2012
    from .FixLabelMapsTools import RecodeLabelMap

    ### Original NeuroMorphometrica merged fusion
    jointFusion = pe.Node(interface=ants.JointFusion(),name="JointFusion")
    many_cpu_JointFusion_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,8,4,4), 'overwrite': True}
    jointFusion.plugin_args = many_cpu_JointFusion_options_dictionary
    jointFusion.inputs.dimension=3
    jointFusion.inputs.method='Joint[0.1,2]'
    jointFusion.inputs.output_label_image='MALF_HDAtlas20_2015_label.nii.gz'

    MALFWF.connect(warpedAtlasesMergeNode,'out',jointFusion,'warped_intensity_images')
    MALFWF.connect(warpedAtlasLblMergeNode,'out',jointFusion,'warped_label_images')
    #MALFWF.connect(inputsSpec, 'subj_t1_image',jointFusion,'target_image')
    MALFWF.connect(sessionMakeMultimodalInput, 'outFNs',jointFusion,'target_image')
    MALFWF.connect(jointFusion, 'output_label_image', outputsSpec,'MALF_HDAtlas20_2015_label')

    if onlyT1:
        jointFusion.inputs.modalities=1
    else:
        jointFusion.inputs.modalities=2


    ## We need to recode values to ensure that the labels match FreeSurer as close as possible by merging
    ## some labels together to standard FreeSurfer confenventions (i.e. for WMQL)
    RECODE_LABELS_2_Standard_FSWM = [
                                (15071,47),(15072,47),(15073,47),(15145,1011),(15157,1011),(15161,1011),
                                (15179,1012),(15141,1014),(15151,1017),(15163,1018),(15165,1019),(15143,1027),
                                (15191,1028),(15193,1028),(15185,1030),(15201,1030),(15175,1031),(15195,1031),
                                (15173,1035),(15144,2011),(15156,2011),(15160,2011),(15178,2012),(15140,2014),
                                (15150,2017),(15162,2018),(15164,2019),(15142,2027),(15190,2028),(15192,2028),
                                (15184,2030),(15174,2031),(15194,2031),(15172,2035),(15200,2030)]
    ## def RecodeLabelMap(InputFileName,OutputFileName,RECODE_TABLE):
    RecodeToStandardFSWM = pe.Node(Function(function=RecodeLabelMap,
                                                   input_names=['InputFileName','OutputFileName','RECODE_TABLE'],
                                                   output_names=['OutputFileName']),
                                                   name="RecodeToStandardFSWM")
    RecodeToStandardFSWM.inputs.RECODE_TABLE = RECODE_LABELS_2_Standard_FSWM
    RecodeToStandardFSWM.inputs.OutputFileName = 'MALF_HDAtlas20_2015_fs_standard_label.nii.gz'

    MALFWF.connect(RecodeToStandardFSWM,'OutputFileName',outputsSpec,'MALF_HDAtlas20_2015_fs_standard_label')

    ## MALF_SNAPSHOT_WRITER for Segmented result checking:
#    MALF_SNAPSHOT_WRITERNodeName = "MALF_ExtendedMALF_SNAPSHOT_WRITER"
#    MALF_SNAPSHOT_WRITER = pe.Node(interface=BRAINSSnapShotWriter(), name=MALF_SNAPSHOT_WRITERNodeName)

#    MALF_SNAPSHOT_WRITER.inputs.outputFilename = 'MALF_HDAtlas20_2015_CSFVBInjected_label.png'  # output specification
#    MALF_SNAPSHOT_WRITER.inputs.inputPlaneDirection = [2, 1, 1, 1, 1, 0, 0]
#    MALF_SNAPSHOT_WRITER.inputs.inputSliceToExtractInPhysicalPoint = [-3, -7, -3, 5, 7, 22, -22]

#    MALFWF.connect(MALF_SNAPSHOT_WRITER,'outputFilename',outputsSpec,'MALF_extended_snapshot')

    if runFixFusionLabelMap:
        ## post processing of jointfusion
        injectSurfaceCSFandVBIntoLabelMap = pe.Node(Function(function=FixLabelMapFromNeuromorphemetrics2012,
                                                      input_names=['fusionFN',
                                                        'FixedHeadFN',
                                                        'LeftHemisphereFN',
                                                        'outFN',
                                                        'OUT_DICT'],
                                                      output_names=['fixedFusionLabelFN']),
                                               name="injectSurfaceCSFandVBIntoLabelMap")
        injectSurfaceCSFandVBIntoLabelMap.inputs.outFN = 'MALF_HDAtlas20_2015_CSFVBInjected_label.nii.gz'
        FREESURFER_DICT = { 'BRAINSTEM': 16, 'RH_CSF':24, 'LH_CSF':24, 'BLOOD': 15000, 'UNKNOWN': 999,
                            'CONNECTED': [11,12,13,9,17,26,50,51,52,48,53,58]
                          }
        injectSurfaceCSFandVBIntoLabelMap.inputs.OUT_DICT = FREESURFER_DICT
        MALFWF.connect(jointFusion, 'output_label_image', injectSurfaceCSFandVBIntoLabelMap, 'fusionFN')
        MALFWF.connect(inputsSpec, 'subj_fixed_head_labels', injectSurfaceCSFandVBIntoLabelMap, 'FixedHeadFN')
        MALFWF.connect(inputsSpec, 'subj_left_hemisphere', injectSurfaceCSFandVBIntoLabelMap, 'LeftHemisphereFN')

        MALFWF.connect(injectSurfaceCSFandVBIntoLabelMap, 'fixedFusionLabelFN',
                       RecodeToStandardFSWM,'InputFileName')

        MALFWF.connect(injectSurfaceCSFandVBIntoLabelMap,'fixedFusionLabelFN',
                       outputsSpec,'MALF_HDAtlas20_2015_CSFVBInjected_label')
#        MALFWF.connect([(inputsSpec, MALF_SNAPSHOT_WRITER, [( 'subj_t1_image','inputVolumes')]),
#                    (injectSurfaceCSFandVBIntoLabelMap, MALF_SNAPSHOT_WRITER,
#                      [('fixedFusionLabelFN', 'inputBinaryVolumes')])
#                   ])
    else:
        MALFWF.connect(jointFusion, 'output_label_image',
                       RecodeToStandardFSWM,'InputFileName')
        MALFWF.connect(jointFusion, 'output_label_image',
                       outputsSpec,'MALF_HDAtlas20_2015_CSFVBInjected_label')
#        MALFWF.connect([(inputsSpec, MALF_SNAPSHOT_WRITER, [( 'subj_t1_image','inputVolumes')]),
#                    (jointFusion, MALF_SNAPSHOT_WRITER,
#                      [('output_label_image', 'inputBinaryVolumes')])
#                   ])

    ## Lobar Pacellation by recoding
    if master_config['relabel2lobes_filename'] != None:
        #print("Generate relabeled version based on {0}".format(master_config['relabel2lobes_filename']))

        RECODE_LABELS_2_LobarPacellation = readRecodingList( master_config['relabel2lobes_filename'] )
        RecordToFSLobes = pe.Node(Function(function=RecodeLabelMap,
                                                    input_names=['InputFileName','OutputFileName','RECODE_TABLE'],
                                                    output_names=['OutputFileName']),
                                                    name="RecordToFSLobes")
        RecordToFSLobes.inputs.RECODE_TABLE = RECODE_LABELS_2_LobarPacellation
        RecordToFSLobes.inputs.OutputFileName = 'MALF_HDAtlas20_2015_lobar_label.nii.gz'
        MALFWF.connect(RecodeToStandardFSWM, 'OutputFileName',RecordToFSLobes,'InputFileName')
        MALFWF.connect(RecordToFSLobes,'OutputFileName',outputsSpec,'MALF_HDAtlas20_2015_lobar_label')

    return MALFWF
