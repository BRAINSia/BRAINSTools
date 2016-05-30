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

from nipype.interfaces.semtools.registration.brainsresample import BRAINSResample

from utilities.misc import *
from utilities.distributed import modify_qsub_args
from nipype.interfaces.semtools.utilities.brains import BRAINSLandmarkInitializer
from .WorkupAtlasDustCleanup import CreateDustCleanupWorkflow

from utilities.misc import CommonANTsRegistrationSettings
from .WorkupComputeLabelVolume import *

# HACK Remove due to bugs from nipype.interfaces.semtools import BRAINSSnapShotWriter

"""
    from WorkupJointFusion import CreateJointFusionWorkflow
    myLocalTCWF= CreateJointFusionWorkflow("JointFusion")
    JointFusionWF.connect( [ (uidSource, myLocalTCWF, [(('uid', getT1s, subjectDatabaseFile ), 'T1List')] ), ])
    JointFusionWF.connect( [ (uidSource, myLocalTCWF, [(('uid', getT2s, subjectDatabaseFile ), 'T2List')] ), ])
    JointFusionWF.connect( [ (uidSource, myLocalTCWF, [(('uid', getT1sLength, subjectDatabaseFile ), 'T1_count')] ), ])
    JointFusionWF.connect( BCD,    'outputResampledVolume', myLocalTCWF, 'PrimaryT1' )
    JointFusionWF.connect(BAtlas,'ExtendedAtlasDefinition.xml',myLocalTCWF,'atlasDefinition')
    JointFusionWF.connect(BLI,'outputTransformFilename',myLocalTCWF,'atlasToSubjectInitialTransform')
"""
def MakeVector(inFN1, inFN2=None, jointFusion =False):
    #print("inFN1: {0}".format(inFN1))
    #print("inFN2: {0}".format(inFN2))
    if inFN2 == None:
        returnVector = [ inFN1 ]
    else:
        returnVector = [inFN1, inFN2]

    if jointFusion:
        returnVector = [returnVector]

    print ("jointFusion: ")
    print (str(jointFusion))
    print (returnVector)
    print ("============================================")
    return returnVector


def adjustMergeList(allList, n_modality):
    def yieldList(inList, n):
        for i in xrange(0, len(inList), n):
            yield inList[i:i+n]
    return list(yieldList(allList, n_modality))

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
    jointFusionAtlasDict = {}
    #scanID, ['atlasID', 't1', 't2' ,'label', 'lmks']
    import ast
    with open(dictionaryFilename, 'r') as f:
        for line in f.readlines():
            scanID, scanDict=ast.literal_eval(line)
            jointFusionAtlasDict[scanID]=scanDict
            print( line )

    return jointFusionAtlasDict

def getListIndexOrNoneIfOutOfRange(imageList, index):
    if index < len(imageList):
        return imageList[index]
    else:
        return None


def CreateJointFusionWorkflow(WFname, onlyT1, master_config, runFixFusionLabelMap=True):
    from nipype.interfaces import ants

    if onlyT1:
      n_modality = 1
    else:
      n_modality = 2
    CLUSTER_QUEUE=master_config['queue']
    CLUSTER_QUEUE_LONG=master_config['long_q']

    JointFusionWF = pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['subj_t1_image', #Desired image to create label map for
                                                             'subj_t2_image', #Desired image to create label map for
                                                             'subj_lmks', #The landmarks corresponding to t1_image
                                                             'subj_fixed_head_labels', #The fixed head labels from BABC
                                                             'subj_posteriors', #The BABC posteriors
                                                             'subj_left_hemisphere', #The warped left hemisphere mask
                                                             'atlasWeightFilename',  #The static weights file name
                                                             'labelBaseFilename' #Atlas label base name ex) neuro_lbls.nii.gz
                                                            ]),
                         run_without_submitting=True,
                         name='inputspec')
    outputsSpec = pe.Node(interface=IdentityInterface(fields=['JointFusion_HDAtlas20_2015_label',
                                                       'JointFusion_HDAtlas20_2015_CSFVBInjected_label',
                                                       'JointFusion_HDAtlas20_2015_fs_standard_label',
                                                       'JointFusion_HDAtlas20_2015_lobe_label',
                                                       'JointFusion_extended_snapshot',
                                                       'JointFusion_HDAtlas20_2015_dustCleaned_label',
                                                       'JointFusion_volumes_csv',
                                                       'JointFusion_volumes_json',
                                                       'JointFusion_lobe_volumes_csv',
                                                       'JointFusion_lobe_volumes_json']),
                          run_without_submitting=True,
                          name='outputspec')

    BLICreator = dict()
    A2SantsRegistrationPreJointFusion_SyN = dict()
    movingROIAuto = dict()
    labelMapResample = dict()
    NewlabelMapResample = dict()

    jointFusion_atlas_mergeindex = 0
    merge_input_offset = 1 #Merge nodes are indexed from 1, not zero!
    """
    multimodal ants registration if t2 exists
    """
    sessionMakeMultimodalInput = pe.Node(Function(function=MakeVector,
                                         input_names=['inFN1', 'inFN2', 'jointFusion'],
                                         output_names=['outFNs']),
                                run_without_submitting=True, name="sessionMakeMultimodalInput")
    sessionMakeMultimodalInput.inputs.jointFusion = False
    JointFusionWF.connect(inputsSpec, 'subj_t1_image', sessionMakeMultimodalInput, 'inFN1')
    """
    T2 resample to T1 average image
    :: BRAINSABC changed its behavior to retain image's original spacing & origin
    :: Since antsJointFusion only works for the identical origin images for targets,
    :: Resampling is placed at this stage
    """
    subjectT2Resample = pe.Node(interface=BRAINSResample(), name="BRAINSResample_T2_forAntsJointFusion")
    if not onlyT1:
        subjectT2Resample.plugin_args = {'qsub_args': modify_qsub_args(master_config['queue'], 1, 1, 1),
                                  'overwrite': True}
        subjectT2Resample.inputs.pixelType = 'short'
        subjectT2Resample.inputs.interpolationMode = 'Linear'
        subjectT2Resample.inputs.outputVolume = "t2_resampled_in_t1.nii.gz"
        #subjectT2Resample.inputs.warpTransform= "Identity" # Default is "Identity"

        JointFusionWF.connect(inputsSpec, 'subj_t1_image', subjectT2Resample, 'referenceVolume')
        JointFusionWF.connect(inputsSpec, 'subj_t2_image', subjectT2Resample, 'inputVolume')

        JointFusionWF.connect(subjectT2Resample, 'outputVolume', sessionMakeMultimodalInput, 'inFN2')
    else:
        pass


    #print('jointFusion_atlas_db_base')
    print("master_config")
    print(master_config)
    print("master_config['jointfusion_atlas_db_base']")
    print(master_config['jointfusion_atlas_db_base'])
    jointFusionAtlasDict = readMalfAtlasDbBase( master_config['jointfusion_atlas_db_base'] )
    number_of_atlas_sources = len(jointFusionAtlasDict)
    jointFusionAtlases = dict()
    atlasMakeMultimodalInput = dict()
    t2Resample = dict()
    warpedAtlasLblMergeNode = pe.Node(interface=Merge(number_of_atlas_sources),name="LblMergeAtlas")
    NewwarpedAtlasLblMergeNode = pe.Node(interface=Merge(number_of_atlas_sources),name="fswmLblMergeAtlas")
    # "HACK NOT to use T2 for JointFusion only"
    #warpedAtlasesMergeNode = pe.Node(interface=Merge(number_of_atlas_sources*n_modality),name="MergeAtlases")
    warpedAtlasesMergeNode = pe.Node(interface=Merge(number_of_atlas_sources*1),name="MergeAtlases")

    ## if using Registration masking, then do ROIAuto on fixed and moving images and connect to registraitons
    UseRegistrationMasking = True
    if UseRegistrationMasking == True:
       from nipype.interfaces.semtools.segmentation.specialized import BRAINSROIAuto

       fixedROIAuto = pe.Node(interface=BRAINSROIAuto(), name="fixedROIAUTOMask")
       fixedROIAuto.inputs.ROIAutoDilateSize=10
       fixedROIAuto.inputs.outputROIMaskVolume = "fixedImageROIAutoMask.nii.gz"
       JointFusionWF.connect(inputsSpec, 'subj_t1_image',fixedROIAuto,'inputVolume')


    for jointFusion_atlas_subject in list(jointFusionAtlasDict.keys()):
        ## Need DataGrabber Here For the Atlas
        jointFusionAtlases[jointFusion_atlas_subject] = pe.Node(interface = IdentityInterface(
                                                                  fields=['t1','t2','label','lmks','regisration_mask']),
                                                                  name='jointFusionAtlasInput'+jointFusion_atlas_subject)
        jointFusionAtlases[jointFusion_atlas_subject].inputs.t1 = jointFusionAtlasDict[jointFusion_atlas_subject]['t1']
        jointFusionAtlases[jointFusion_atlas_subject].inputs.t2 = jointFusionAtlasDict[jointFusion_atlas_subject]['t2']
        jointFusionAtlases[jointFusion_atlas_subject].inputs.label = jointFusionAtlasDict[jointFusion_atlas_subject]['label']
        jointFusionAtlases[jointFusion_atlas_subject].inputs.lmks = jointFusionAtlasDict[jointFusion_atlas_subject]['lmks']
        jointFusionAtlases[jointFusion_atlas_subject].inputs.regisration_mask = jointFusionAtlasDict[jointFusion_atlas_subject]['regisration_mask']
        ## Create BLI first
        ########################################################
        # Run BLI atlas_to_subject
        ########################################################
        BLICreator[jointFusion_atlas_subject] = pe.Node(interface=BRAINSLandmarkInitializer(), name="BLI_"+jointFusion_atlas_subject)
        BLICreator[jointFusion_atlas_subject].inputs.outputTransformFilename = "landmarkInitializer_{0}_to_subject_transform.h5".format(jointFusion_atlas_subject)

        JointFusionWF.connect(inputsSpec, 'atlasWeightFilename', BLICreator[jointFusion_atlas_subject], 'inputWeightFilename')
        JointFusionWF.connect(jointFusionAtlases[jointFusion_atlas_subject], 'lmks', BLICreator[jointFusion_atlas_subject], 'inputMovingLandmarkFilename')
        JointFusionWF.connect(inputsSpec, 'subj_lmks', BLICreator[jointFusion_atlas_subject], 'inputFixedLandmarkFilename')

        ##### Initialize with ANTS Transform For SyN
        currentAtlasToSubjectantsRegistration = 'SyN_AtlasToSubjectANTsPreJointFusion_'+jointFusion_atlas_subject
        A2SantsRegistrationPreJointFusion_SyN[jointFusion_atlas_subject] = pe.Node(interface=ants.Registration(), name=currentAtlasToSubjectantsRegistration)
        many_cpu_ANTsSyN_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE_LONG,4,2,16), 'overwrite': True}
        A2SantsRegistrationPreJointFusion_SyN[jointFusion_atlas_subject].plugin_args = many_cpu_ANTsSyN_options_dictionary
        if onlyT1:
            JFregistrationTypeDescription="FiveStageAntsRegistrationT1Only"
        else:
            JFregistrationTypeDescription="FiveStageAntsRegistrationMultiModal"
        CommonANTsRegistrationSettings(
                      antsRegistrationNode=A2SantsRegistrationPreJointFusion_SyN[jointFusion_atlas_subject],
                      registrationTypeDescription=JFregistrationTypeDescription,
                      output_transform_prefix=jointFusion_atlas_subject+'_ToSubjectPreJointFusion_SyN',
                      output_warped_image=jointFusion_atlas_subject + '_2subject.nii.gz',
                      output_inverse_warped_image=None, #NO NEED FOR THIS
                      save_state=None,                  #NO NEED FOR THIS
                      invert_initial_moving_transform=False)

        ## if using Registration masking, then do ROIAuto on fixed and moving images and connect to registraitons
        if UseRegistrationMasking == True:
            from nipype.interfaces.semtools.segmentation.specialized import BRAINSROIAuto
            JointFusionWF.connect(fixedROIAuto, 'outputROIMaskVolume',A2SantsRegistrationPreJointFusion_SyN[jointFusion_atlas_subject],'fixed_image_mask')
            # JointFusionWF.connect(inputsSpec, 'subj_fixed_head_labels',
            #                       A2SantsRegistrationPreJointFusion_SyN[jointFusion_atlas_subject],'fixed_image_mask')

            # NOTE: Moving image mask can be taken from Atlas directly so that it does not need to be read in
            #movingROIAuto[jointFusion_atlas_subject] = pe.Node(interface=BRAINSROIAuto(), name="movingROIAUTOMask_"+jointFusion_atlas_subject)
            #movingROIAuto.inputs.ROIAutoDilateSize=10
            #movingROIAuto[jointFusion_atlas_subject].inputs.outputROIMaskVolume = "movingImageROIAutoMask.nii.gz"
            #JointFusionWF.connect(jointFusionAtlases[jointFusion_atlas_subject], 't1', movingROIAuto[jointFusion_atlas_subject],'inputVolume')
            #JointFusionWF.connect(movingROIAuto[jointFusion_atlas_subject], 'outputROIMaskVolume',A2SantsRegistrationPreJointFusion_SyN[jointFusion_atlas_subject],'moving_image_mask')
            JointFusionWF.connect(jointFusionAtlases[jointFusion_atlas_subject], 'regisration_mask',
                                  A2SantsRegistrationPreJointFusion_SyN[jointFusion_atlas_subject],'moving_image_mask')

        JointFusionWF.connect(BLICreator[jointFusion_atlas_subject],'outputTransformFilename',
                       A2SantsRegistrationPreJointFusion_SyN[jointFusion_atlas_subject],'initial_moving_transform')

        """
        make multimodal input for atlases
        """
        atlasMakeMultimodalInput[jointFusion_atlas_subject] = pe.Node(Function(function=MakeVector,
                                                                               input_names=['inFN1', 'inFN2','jointFusion'],
                                                                               output_names=['outFNs']),
                                  run_without_submitting=True, name="atlasMakeMultimodalInput"+jointFusion_atlas_subject)
        atlasMakeMultimodalInput[jointFusion_atlas_subject].inputs.jointFusion = False
        JointFusionWF.connect(jointFusionAtlases[jointFusion_atlas_subject], 't1',
                              atlasMakeMultimodalInput[jointFusion_atlas_subject], 'inFN1')

        if not onlyT1:
            JointFusionWF.connect(jointFusionAtlases[jointFusion_atlas_subject], 't2', atlasMakeMultimodalInput[jointFusion_atlas_subject], 'inFN2')
        else:
            pass

        JointFusionWF.connect(sessionMakeMultimodalInput, 'outFNs',
                       A2SantsRegistrationPreJointFusion_SyN[jointFusion_atlas_subject],'fixed_image')
        JointFusionWF.connect(atlasMakeMultimodalInput[jointFusion_atlas_subject], 'outFNs',
                       A2SantsRegistrationPreJointFusion_SyN[jointFusion_atlas_subject],'moving_image')
        "HACK NOT to use T2 for JointFusion"
        #JointFusionWF.connect(A2SantsRegistrationPreJointFusion_SyN[jointFusion_atlas_subject],'warped_image',
        #               warpedAtlasesMergeNode,'in'+str(merge_input_offset + jointFusion_atlas_mergeindex*n_modality) )
        JointFusionWF.connect(A2SantsRegistrationPreJointFusion_SyN[jointFusion_atlas_subject],'warped_image',
                       warpedAtlasesMergeNode,'in'+str(merge_input_offset + jointFusion_atlas_mergeindex*1) )

        """
        Original t2 resampling
        """
        for modality_index in range(1,n_modality):
            t2Resample[jointFusion_atlas_subject] = pe.Node(interface=ants.ApplyTransforms(),name="resampledT2"+jointFusion_atlas_subject)
            many_cpu_t2Resample_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,1,1,1), 'overwrite': True}
            t2Resample[jointFusion_atlas_subject].plugin_args = many_cpu_t2Resample_options_dictionary
            t2Resample[jointFusion_atlas_subject].inputs.num_threads=-1
            t2Resample[jointFusion_atlas_subject].inputs.dimension=3
            t2Resample[jointFusion_atlas_subject].inputs.output_image=jointFusion_atlas_subject+'_t2.nii.gz'
            t2Resample[jointFusion_atlas_subject].inputs.interpolation='BSpline'
            t2Resample[jointFusion_atlas_subject].inputs.default_value=0
            t2Resample[jointFusion_atlas_subject].inputs.invert_transform_flags=[False]

            JointFusionWF.connect( A2SantsRegistrationPreJointFusion_SyN[jointFusion_atlas_subject],'composite_transform',
                            t2Resample[jointFusion_atlas_subject],'transforms')
            JointFusionWF.connect( inputsSpec, 'subj_t1_image',
                            t2Resample[jointFusion_atlas_subject],'reference_image')
            JointFusionWF.connect( jointFusionAtlases[jointFusion_atlas_subject], 't2',
                            t2Resample[jointFusion_atlas_subject],'input_image')
            "HACK NOT to use T2 for JointFusion only"
            #JointFusionWF.connect(t2Resample[jointFusion_atlas_subject],'output_image',
            #               warpedAtlasesMergeNode,'in'+str(merge_input_offset + jointFusion_atlas_mergeindex*n_modality+modality_index) )

        """
        Original labelmap resampling
        """
        labelMapResample[jointFusion_atlas_subject] = pe.Node(interface=ants.ApplyTransforms(),name="resampledLabel"+jointFusion_atlas_subject)
        many_cpu_labelMapResample_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,1,1,1), 'overwrite': True}
        labelMapResample[jointFusion_atlas_subject].plugin_args = many_cpu_labelMapResample_options_dictionary
        labelMapResample[jointFusion_atlas_subject].inputs.num_threads=-1
        labelMapResample[jointFusion_atlas_subject].inputs.dimension=3
        labelMapResample[jointFusion_atlas_subject].inputs.output_image=jointFusion_atlas_subject+'_2_subj_lbl.nii.gz'
        labelMapResample[jointFusion_atlas_subject].inputs.interpolation='MultiLabel'
        labelMapResample[jointFusion_atlas_subject].inputs.default_value=0
        labelMapResample[jointFusion_atlas_subject].inputs.invert_transform_flags=[False]

        JointFusionWF.connect( A2SantsRegistrationPreJointFusion_SyN[jointFusion_atlas_subject],'composite_transform',
                        labelMapResample[jointFusion_atlas_subject],'transforms')
        JointFusionWF.connect( inputsSpec, 'subj_t1_image',
                        labelMapResample[jointFusion_atlas_subject],'reference_image')
        JointFusionWF.connect( jointFusionAtlases[jointFusion_atlas_subject], 'label',
                        labelMapResample[jointFusion_atlas_subject],'input_image')


        JointFusionWF.connect(labelMapResample[jointFusion_atlas_subject],'output_image',warpedAtlasLblMergeNode,'in'+str(merge_input_offset + jointFusion_atlas_mergeindex) )

        ### New labelmap resampling
        NewlabelMapResample[jointFusion_atlas_subject] = pe.Node(interface=ants.ApplyTransforms(),name="FSWM_WLABEL_"+jointFusion_atlas_subject)
        many_cpu_NewlabelMapResample_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,1,1,1), 'overwrite': True}
        NewlabelMapResample[jointFusion_atlas_subject].plugin_args = many_cpu_NewlabelMapResample_options_dictionary
        NewlabelMapResample[jointFusion_atlas_subject].inputs.num_threads=-1
        NewlabelMapResample[jointFusion_atlas_subject].inputs.dimension=3
        NewlabelMapResample[jointFusion_atlas_subject].inputs.output_image=jointFusion_atlas_subject+'fswm_2_subj_lbl.nii.gz'
        NewlabelMapResample[jointFusion_atlas_subject].inputs.interpolation='MultiLabel'
        NewlabelMapResample[jointFusion_atlas_subject].inputs.default_value=0
        NewlabelMapResample[jointFusion_atlas_subject].inputs.invert_transform_flags=[False]

        JointFusionWF.connect( A2SantsRegistrationPreJointFusion_SyN[jointFusion_atlas_subject],'composite_transform',
                        NewlabelMapResample[jointFusion_atlas_subject],'transforms')
        JointFusionWF.connect( inputsSpec, 'subj_t1_image',
                        NewlabelMapResample[jointFusion_atlas_subject],'reference_image')
        JointFusionWF.connect( jointFusionAtlases[jointFusion_atlas_subject], 'label',
                        NewlabelMapResample[jointFusion_atlas_subject],'input_image')


        JointFusionWF.connect(NewlabelMapResample[jointFusion_atlas_subject],'output_image',NewwarpedAtlasLblMergeNode,'in'+str(merge_input_offset + jointFusion_atlas_mergeindex) )

        jointFusion_atlas_mergeindex += 1


    ## Now work on cleaning up the label maps
    from .FixLabelMapsTools import FixLabelMapFromNeuromorphemetrics2012
    from .FixLabelMapsTools import RecodeLabelMap

    ### Original NeuroMorphometrica merged fusion
    jointFusion = pe.Node(interface=ants.AntsJointFusion(),name="AntsJointFusion")
    many_cpu_JointFusion_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE,10,8,16), 'overwrite': True}
    jointFusion.plugin_args = many_cpu_JointFusion_options_dictionary
    jointFusion.inputs.num_threads = -1
    jointFusion.inputs.dimension=3
    jointFusion.inputs.search_radius=[3]
    #jointFusion.inputs.method='Joint[0.1,2]'
    jointFusion.inputs.out_label_fusion='JointFusion_HDAtlas20_2015_label.nii.gz'
    #JointFusionWF.connect(inputsSpec, 'subj_fixed_head_labels', jointFusion, 'mask_image')
    JointFusionWF.connect(fixedROIAuto, 'outputROIMaskVolume', jointFusion, 'mask_image')

    JointFusionWF.connect(warpedAtlasLblMergeNode,'out',
                          jointFusion,'atlas_segmentation_image')

    AdjustMergeListNode = pe.Node(Function(function=adjustMergeList,
                                                   input_names=['allList','n_modality'],
                                                   output_names=['out']),
                                                   name="AdjustMergeListNode")
    "*** HACK JointFusion only uses T1"
    #AdjustMergeListNode.inputs.n_modality = n_modality
    AdjustMergeListNode.inputs.n_modality = 1

    JointFusionWF.connect(warpedAtlasesMergeNode,'out',AdjustMergeListNode,'allList')
    JointFusionWF.connect(AdjustMergeListNode,'out',jointFusion,'atlas_image')

    AdjustTargetImageListNode = pe.Node(Function(function=adjustMergeList,
                                                   input_names=['allList','n_modality'],
                                                   output_names=['out']),
                                                   name="AdjustTargetImageListNode")
    AdjustTargetImageListNode.inputs.n_modality = n_modality

    "*** HACK JointFusion only uses T1"
    """ Once JointFusion works with T2 properly,
        delete sessionMakeListSingleModalInput and use sessionMakeMultimodalInput instead
    """
    sessionMakeListSingleModalInput = pe.Node(Function(function=MakeVector,
                                         input_names=['inFN1', 'inFN2', 'jointFusion'],
                                         output_names=['outFNs']),
                                run_without_submitting=True, name="sessionMakeListSingleModalInput")
    sessionMakeListSingleModalInput.inputs.jointFusion = False
    JointFusionWF.connect(inputsSpec, 'subj_t1_image', sessionMakeListSingleModalInput, 'inFN1')
    JointFusionWF.connect(sessionMakeListSingleModalInput, 'outFNs', jointFusion,'target_image')

    JointFusionWF.connect(jointFusion, 'out_label_fusion', outputsSpec,'JointFusion_HDAtlas20_2015_label')

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
    RecodeToStandardFSWM.inputs.OutputFileName = 'JointFusion_HDAtlas20_2015_fs_standard_label.nii.gz'

    JointFusionWF.connect(RecodeToStandardFSWM,'OutputFileName',outputsSpec,'JointFusion_HDAtlas20_2015_fs_standard_label')

    ## JointFusion_SNAPSHOT_WRITER for Segmented result checking:
#    JointFusion_SNAPSHOT_WRITERNodeName = "JointFusion_ExtendedJointFusion_SNAPSHOT_WRITER"
#    JointFusion_SNAPSHOT_WRITER = pe.Node(interface=BRAINSSnapShotWriter(), name=JointFusion_SNAPSHOT_WRITERNodeName)

#    JointFusion_SNAPSHOT_WRITER.inputs.outputFilename = 'JointFusion_HDAtlas20_2015_CSFVBInjected_label.png'  # output specification
#    JointFusion_SNAPSHOT_WRITER.inputs.inputPlaneDirection = [2, 1, 1, 1, 1, 0, 0]
#    JointFusion_SNAPSHOT_WRITER.inputs.inputSliceToExtractInPhysicalPoint = [-3, -7, -3, 5, 7, 22, -22]

#    JointFusionWF.connect(JointFusion_SNAPSHOT_WRITER,'outputFilename',outputsSpec,'JointFusion_extended_snapshot')

    myLocalDustCleanup = CreateDustCleanupWorkflow("DUST_CLEANUP", onlyT1, master_config)
    JointFusionWF.connect(inputsSpec, 'subj_t1_image', myLocalDustCleanup, 'inputspec.subj_t1_image')
    if not onlyT1:
        JointFusionWF.connect(subjectT2Resample, 'outputVolume', myLocalDustCleanup, 'inputspec.subj_t2_image')
    if runFixFusionLabelMap:
        ## post processing of jointfusion
        injectSurfaceCSFandVBIntoLabelMap = pe.Node(Function(function=FixLabelMapFromNeuromorphemetrics2012,
                                                      input_names=['fusionFN',
                                                        'FixedHeadFN',
                                                        'posterior_dict',
                                                        'LeftHemisphereFN',
                                                        'outFN',
                                                        'OUT_DICT'],
                                                      output_names=['fixedFusionLabelFN']),
                                               name="injectSurfaceCSFandVBIntoLabelMap")
        injectSurfaceCSFandVBIntoLabelMap.inputs.outFN = 'JointFusion_HDAtlas20_2015_CSFVBInjected_label.nii.gz'
        FREESURFER_DICT = { 'BRAINSTEM': 16, 'RH_CSF':24, 'LH_CSF':24, 'BLOOD': 15000, 'UNKNOWN': 999,
                            'CONNECTED': [11,12,13,9,17,26,50,51,52,48,53,58]
                          }
        injectSurfaceCSFandVBIntoLabelMap.inputs.OUT_DICT = FREESURFER_DICT
        JointFusionWF.connect(jointFusion, 'out_label_fusion', injectSurfaceCSFandVBIntoLabelMap, 'fusionFN')
        JointFusionWF.connect(inputsSpec, 'subj_fixed_head_labels', injectSurfaceCSFandVBIntoLabelMap, 'FixedHeadFN')
        JointFusionWF.connect(inputsSpec, 'subj_posteriors', injectSurfaceCSFandVBIntoLabelMap, 'posterior_dict')
        JointFusionWF.connect(inputsSpec, 'subj_left_hemisphere', injectSurfaceCSFandVBIntoLabelMap, 'LeftHemisphereFN')

        JointFusionWF.connect(injectSurfaceCSFandVBIntoLabelMap, 'fixedFusionLabelFN',
                       myLocalDustCleanup, 'inputspec.subj_label_atlas')

        JointFusionWF.connect(injectSurfaceCSFandVBIntoLabelMap,'fixedFusionLabelFN',
                       outputsSpec,'JointFusion_HDAtlas20_2015_CSFVBInjected_label')

        JointFusionWF.connect(myLocalDustCleanup, 'outputspec.JointFusion_HDAtlas20_2015_dustCleaned_label',
                       RecodeToStandardFSWM,'InputFileName')

        JointFusionWF.connect(myLocalDustCleanup, 'outputspec.JointFusion_HDAtlas20_2015_dustCleaned_label',
                       outputsSpec, 'JointFusion_HDAtlas20_2015_dustCleaned_label')

#        JointFusionWF.connect([(inputsSpec, JointFusion_SNAPSHOT_WRITER, [( 'subj_t1_image','inputVolumes')]),
#                    (injectSurfaceCSFandVBIntoLabelMap, JointFusion_SNAPSHOT_WRITER,
#                      [('fixedFusionLabelFN', 'inputBinaryVolumes')])
#                   ])
    else:
        JointFusionWF.connect(jointFusion, 'output_label_image',
                       myLocalDustCleanup, 'inputspec.subj_label_atlas')

        JointFusionWF.connect(jointFusion, 'output_label_image',
                       outputsSpec,'JointFusion_HDAtlas20_2015_CSFVBInjected_label')

        JointFusionWF.connect(myLocalDustCleanup, 'outputspec.JointFusion_HDAtlas20_2015_dustCleaned_label',
                       RecodeToStandardFSWM,'InputFileName')

        JointFusionWF.connect(myLocalDustCleanup, 'outputspec.JointFusion_HDAtlas20_2015_dustCleaned_label',
                       outputsSpec, 'JointFusion_HDAtlas20_2015_dustCleaned_label')

#        JointFusionWF.connect([(inputsSpec, JointFusion_SNAPSHOT_WRITER, [( 'subj_t1_image','inputVolumes')]),
#                    (jointFusion, JointFusion_SNAPSHOT_WRITER,
#                      [('output_label_image', 'inputBinaryVolumes')])
#                   ])

    """
    Compute label volumes
    """
    computeLabelVolumes = CreateVolumeMeasureWorkflow("LabelVolume", master_config)
    JointFusionWF.connect( inputsSpec, 'subj_t1_image',
                           computeLabelVolumes, 'inputspec.subj_t1_image')
    JointFusionWF.connect( myLocalDustCleanup, 'outputspec.JointFusion_HDAtlas20_2015_dustCleaned_label',
                           computeLabelVolumes, 'inputspec.subj_label_image')
    JointFusionWF.connect( computeLabelVolumes, 'outputspec.csvFilename',
                           outputsSpec, 'JointFusion_volumes_csv')
    JointFusionWF.connect( computeLabelVolumes, 'outputspec.jsonFilename',
                           outputsSpec, 'JointFusion_volumes_json')

    ## Lobe Pacellation by recoding
    if master_config['relabel2lobes_filename'] != None:
        #print("Generate relabeled version based on {0}".format(master_config['relabel2lobes_filename']))

        RECODE_LABELS_2_LobePacellation = readRecodingList( master_config['relabel2lobes_filename'] )
        RecordToFSLobes = pe.Node(Function(function=RecodeLabelMap,
                                                    input_names=['InputFileName','OutputFileName','RECODE_TABLE'],
                                                    output_names=['OutputFileName']),
                                                    name="RecordToFSLobes")
        RecordToFSLobes.inputs.RECODE_TABLE = RECODE_LABELS_2_LobePacellation
        RecordToFSLobes.inputs.OutputFileName = 'JointFusion_HDAtlas20_2015_lobe_label.nii.gz'
        JointFusionWF.connect(RecodeToStandardFSWM, 'OutputFileName',RecordToFSLobes,'InputFileName')
        JointFusionWF.connect(RecordToFSLobes,'OutputFileName',outputsSpec,'JointFusion_HDAtlas20_2015_lobe_label')

        """
        Compute lobe volumes
        """
        computeLobeVolumes = CreateVolumeMeasureWorkflow("LobeVolume", master_config)
        JointFusionWF.connect( inputsSpec, 'subj_t1_image',
                               computeLobeVolumes, 'inputspec.subj_t1_image')
        JointFusionWF.connect( RecordToFSLobes, 'OutputFileName',
                               computeLobeVolumes, 'inputspec.subj_label_image')
        JointFusionWF.connect( computeLobeVolumes, 'outputspec.csvFilename',
                               outputsSpec, 'JointFusion_lobe_volumes_csv')
        JointFusionWF.connect( computeLobeVolumes, 'outputspec.jsonFilename',
                               outputsSpec, 'JointFusion_lobe_volumes_json')

    return JointFusionWF
