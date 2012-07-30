#!/usr/bin/python
#################################################################################
## Program:   BRAINS (Brain Research: Analysis of Images, Networks, and Systems)
## Language:  Python
##
## Author:  Hans J. Johnson
##
##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.
##
#################################################################################

import os
import csv
import sys
import string
import argparse
#"""Import necessary modules from nipype."""
#from nipype.utils.config import config
#config.set('logging', 'log_to_file', 'false')
#config.set_log_dir(os.getcwd())
#--config.set('logging', 'workflow_level', 'DEBUG')
#--config.set('logging', 'interface_level', 'DEBUG')
#--config.set('execution','remove_unnecessary_outputs','false')

from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.freesurfer import ReconAll

from nipype.utils.misc import package_check
#package_check('nipype', '5.4', 'tutorial1') ## HACK: Check nipype version
package_check('numpy', '1.3', 'tutorial1')
package_check('scipy', '0.7', 'tutorial1')
package_check('networkx', '1.0', 'tutorial1')
package_check('IPython', '0.10', 'tutorial1')

from BRAINSTools import *
from BRAINSTools.ANTSWrapper import *
from BRAINSTools.WarpAllAtlas import *
from BRAINSTools.ants.normalize import WarpImageMultiTransform

from WorkupT1T2AtlasNode import MakeAtlasNode


#############################################################################
#############################################################################
## Utility functions for the pipeline
#############################################################################
#############################################################################
def get_first_T1_and_T2(in_files,T1_count):
    '''
    Returns the first T1 and T2 file in in_files, based on offset in T1_count.
    '''
    return in_files[0],in_files[T1_count]

def GetExtensionlessBaseName(filename):
    '''
    Get the filename without the extension.  Works for .ext and .ext.gz
    '''
    import os
    currBaseName = os.path.basename(filename)
    currExt = os.path.splitext(currBaseName)[1]
    currBaseName = os.path.splitext(currBaseName)[0]
    if currExt == ".gz":
        currBaseName = os.path.splitext(currBaseName)[0]
        currExt = os.path.splitext(currBaseName)[1]
    return currBaseName

def get_list_element( nestedList, index ):
    return nestedList[index]
def MakeList(firstElement,secondElement):
    return [firstElement, secondElement]

def GenerateOutputPattern(subjectDatabaseFile,DefaultNodeName,uidIsFirst=True):
    """ This function generates output path substitutions for workflows and nodes that conform to a common standard.
    """
    patternList=[]
    """
    import SessionDB
    dbObject=SessionDB.SessionDB(subjectDatabaseFile)
    for session in dbObject.getAllSessions():
        #print session,"HACK"
        if uidIsFirst == True:
            find_pat=os.path.join(DefaultNodeName,'_uid_'+session)
        else:
            find_pat=os.path.join('_uid_'+session,DefaultNodeName)
        currProj=dbObject.getProjFromSession(session)
        currSubj=dbObject.getSubjFromSession(session)
        replace_pat=os.path.join(currProj,currSubj,session,DefaultNodeName)
        patternList.append( (find_pat,replace_pat) )
    """
    return patternList

###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
## WorkupT1T2 is the main workflow to be run
###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################
def MakeOneSubWorkFlow(sessionid, BAtlas, WORKFLOW_COMPONENTS, BCD_model_path, InterpolationMode, CLUSTER_QUEUE, ExperimentBaseDirectoryResults, subjectDatabaseFile):
    """
    Run autoworkup on a single Subject

    This is the main function to call when processing a data set with T1 & T2
    data.  ExperimentBaseDirectoryPrefix is the base of the directory to place results, T1Images & T2Images
    are the lists of images to be used in the auto-workup. atlas_fname_wpath is
    the path and filename of the atlas to use.
    """

    print "Building Pipeline for ",sessionid
    ########### PIPELINE INITIALIZATION #############
    subWorkflowName='WF_'+str(sessionid)
    T1T2WorkupSingle = pe.Workflow(name=subWorkflowName)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=
                        ['sessionid','subjectid','projectid',
                         'ReferenceT1'
                         ]),
                         run_without_submitting=True,
                         name='InputSpec' )

    outputsSpec = pe.Node(interface=IdentityInterface(fields=['BCD_ACPC_T1'
            ]),
            run_without_submitting=True,
            name='OutputSpec' )


    """
    ### Get the first T1
    def helperGetFirstT1(databaseName,sessionid):
        import SessionDB
        dbObject=SessionDB.SessionDB(databaseName)
        ReferenceT1 = dbObject.getFirstT1(sessionid)
        return ReferenceT1
    ReferenceT1Source = pe.Node(interface=Function(
                               input_names=['databaseName','sessionid'],
                               output_names=['ReferenceT1'],
                               function=helperGetFirstT1),
                      run_without_submitting=True,
                      name='99_ReferenceT1')
    ReferenceT1Source.inputs.databaseName = subjectDatabaseFile
    T1T2WorkupSingle.connect(inputsSpec,'sessionid',ReferenceT1Source,'sessionid')
    """

    ### Get all the T1's
    def helperGetFilenamesByScantype(databaseName, scantype, sessionid ):
        import SessionDB
        dbObject=SessionDB.SessionDB(databaseName)
        allScans = dbObject.getFilenamesByScantype(sessionid,scantype)
        return allScans
    allT1Source = pe.Node(interface=Function(
                               input_names=['databaseName','scantype','sessionid'],
                               output_names=['allT1s'],
                               function=helperGetFilenamesByScantype),
                      run_without_submitting=True,
                      name='99_allT1s')
    allT1Source.inputs.databaseName = subjectDatabaseFile
    allT1Source.inputs.scantype = 'T1-30'
    T1T2WorkupSingle.connect(inputsSpec,'sessionid',allT1Source,'sessionid')

    ### Get all the T2's
    allT2Source = pe.Node(interface=Function(
                               input_names=['databaseName','scantype','sessionid'],
                               output_names=['allT2s'],
                               function=helperGetFilenamesByScantype),
                      run_without_submitting=True,
                      name='99_allT2s')
    allT2Source.inputs.databaseName = subjectDatabaseFile
    allT2Source.inputs.scantype = 'T2-30'
    T1T2WorkupSingle.connect(inputsSpec,'sessionid',allT2Source,'sessionid')

    def getAllT1sLength(allT1s):
        return len(allT1s)

    if 'BASIC' in WORKFLOW_COMPONENTS:
        from WorkupT1T2LandmarkInitialization import CreateLandmarkInitializeWorkflow
        DoReverseMapping = False   # Set to true for debugging outputs
        if 'AUXLMK' in WORKFLOW_COMPONENTS:
            DoReverseMapping = True
        myLocalLMIWF= CreateLandmarkInitializeWorkflow("LandmarkInitialize", BCD_model_path, InterpolationMode,DoReverseMapping)
        T1T2WorkupSingle.connect( inputsSpec, 'ReferenceT1', myLocalLMIWF, 'InputSpec.inputVolume')
        T1T2WorkupSingle.connect( BAtlas, 'template_landmarks_31_fcsv', myLocalLMIWF,'InputSpec.atlasLandmarkFilename')
        T1T2WorkupSingle.connect( BAtlas, 'template_landmark_weights_31_csv', myLocalLMIWF,'InputSpec.atlasWeightFilename')
        if 'AUXLMK' in WORKFLOW_COMPONENTS:
            T1T2WorkupSingle.connect(BAtlas,'template_t1',myLocalLMIWF,'InputSpec.atlasVolume')

        ### Now define where the final organized outputs should go.
        BASIC_DataSink=pe.Node(nio.DataSink(),name="BASIC_DS")
        BASIC_DataSink.inputs.base_directory=ExperimentBaseDirectoryResults
        BASIC_DataSink.inputs.regexp_substitutions = GenerateOutputPattern(subjectDatabaseFile,'ACPCAlign')

        T1T2WorkupSingle.connect(myLocalLMIWF,'OutputSpec.outputLandmarksInACPCAlignedSpace',BASIC_DataSink,'ACPCAlign.@outputLandmarksInACPCAlignedSpace')
        T1T2WorkupSingle.connect(myLocalLMIWF,'OutputSpec.outputResampledVolume',BASIC_DataSink,'ACPCAlign.@outputResampledVolume')
        T1T2WorkupSingle.connect(myLocalLMIWF,'OutputSpec.outputLandmarksInInputSpace',BASIC_DataSink,'ACPCAlign.@outputLandmarksInInputSpace')
        T1T2WorkupSingle.connect(myLocalLMIWF,'OutputSpec.outputTransform',BASIC_DataSink,'ACPCAlign.@outputTransform')
        T1T2WorkupSingle.connect(myLocalLMIWF,'OutputSpec.atlasToSubjectTransform',BASIC_DataSink,'ACPCAlign.@atlasToSubjectTransform')
        ### Now connect OutputSpec
        T1T2WorkupSingle.connect( myLocalLMIWF, 'OutputSpec.outputResampledVolume', outputsSpec, 'BCD_ACPC_T1' )

    if 'TISSUE_CLASSIFY' in WORKFLOW_COMPONENTS:
        from WorkupT1T2TissueClassifiy import CreateTissueClassifyWorkflow
        myLocalTCWF= CreateTissueClassifyWorkflow("TissueClassify",CLUSTER_QUEUE,InterpolationMode)
        T1T2WorkupSingle.connect( allT1Source, 'allT1s', myLocalTCWF, 'InputSpec.T1List')
        T1T2WorkupSingle.connect( allT2Source, 'allT2s', myLocalTCWF, 'InputSpec.T2List')
        T1T2WorkupSingle.connect( [ (allT1Source, myLocalTCWF, [(('allT1s', getAllT1sLength), 'InputSpec.T1_count')] ), ])
        T1T2WorkupSingle.connect( BAtlas,'AtlasPVDefinition_xml',myLocalTCWF,'InputSpec.atlasDefinition')
        T1T2WorkupSingle.connect( myLocalLMIWF, 'OutputSpec.outputResampledVolume', myLocalTCWF, 'InputSpec.PrimaryT1' )
        T1T2WorkupSingle.connect( myLocalLMIWF,'OutputSpec.atlasToSubjectTransform',myLocalTCWF,'InputSpec.atlasToSubjectInitialTransform')

        ### Now define where the final organized outputs should go.
        TC_DataSink=pe.Node(nio.DataSink(),name="TISSUE_CLASSIFY_DS")
        TC_DataSink.inputs.base_directory=ExperimentBaseDirectoryResults
        TC_DataSink.inputs.regexp_substitutions = GenerateOutputPattern(subjectDatabaseFile,'TissueClassify')
        T1T2WorkupSingle.connect(myLocalTCWF, 'OutputSpec.TissueClassifyOutputDir', TC_DataSink,'TissueClassify.@TissueClassifyOutputDir')

    ## Make deformed Atlas image space
    if 'ANTS' in WORKFLOW_COMPONENTS:
        from WorkupT1T2ANTS import CreateANTSRegistrationWorkflow
        myLocalAntsWF = CreateANTSRegistrationWorkflow("ANTSRegistration",CLUSTER_QUEUE,-1)
        T1T2WorkupSingle.connect( myLocalTCWF,'OutputSpec.t1_corrected',myLocalAntsWF,"InputSpec.fixedVolumesList")
        T1T2WorkupSingle.connect( BAtlas,'template_t1',    myLocalAntsWF,"InputSpec.movingVolumesList")
        T1T2WorkupSingle.connect(myLocalLMIWF,'OutputSpec.atlasToSubjectTransform',myLocalAntsWF,'InputSpec.initial_moving_transform')
        # Must register the entire head, not just the brain!
        T1T2WorkupSingle.connect(myLocalTCWF,'OutputSpec.outputHeadLabels',myLocalAntsWF,'InputSpec.fixedBinaryVolume')
        T1T2WorkupSingle.connect(BAtlas,'template_headregion',myLocalAntsWF,'InputSpec.movingBinaryVolume')

        ### Now define where the final organized outputs should go.
        ANTS_DataSink=pe.Node(nio.DataSink(),name="ANTSRegistration_DS")
        ANTS_DataSink.inputs.base_directory=ExperimentBaseDirectoryResults
        ANTS_DataSink.inputs.regexp_substitutions = GenerateOutputPattern(subjectDatabaseFile,'ANTSRegistration')
        T1T2WorkupSingle.connect(myLocalAntsWF, 'OutputSpec.warped_image', ANTS_DataSink,'ANTSRegistration.@warped_image')
        T1T2WorkupSingle.connect(myLocalAntsWF, 'OutputSpec.inverse_warped_image', ANTS_DataSink,'ANTSRegistration.@inverse_warped_image')
        T1T2WorkupSingle.connect(myLocalAntsWF, 'OutputSpec.affine_transform', ANTS_DataSink,'ANTSRegistration.@affine_transform')
        T1T2WorkupSingle.connect(myLocalAntsWF, 'OutputSpec.warp_transform', ANTS_DataSink,'ANTSRegistration.@warp_transform')
        T1T2WorkupSingle.connect(myLocalAntsWF, 'OutputSpec.inverse_warp_transform', ANTS_DataSink,'ANTSRegistration.@inverse_warp_transform')

    if 'SEGMENTATION' in WORKFLOW_COMPONENTS:
        def getListIndex( imageList, index):
            return imageList[index]
        from WorkupT1T2BRAINSCut import CreateBRAINSCutWorkflow
        ## TODO:  Remove BAtlas From Here as well!
        myLocalSegWF = CreateBRAINSCutWorkflow("Segmentation",CLUSTER_QUEUE,BAtlas) ##Note:  Passing in the entire BAtlas Object here!
        T1T2WorkupSingle.connect( [ ( myLocalTCWF, myLocalSegWF, [ (( 'OutputSpec.outputAverageImages', getListIndex, 0 ), "InputSpec.T1Volume")] ), ] )
        T1T2WorkupSingle.connect( [ ( myLocalTCWF, myLocalSegWF, [ (( 'OutputSpec.outputAverageImages', getListIndex, 1 ), "InputSpec.T2Volume")] ), ] )
        T1T2WorkupSingle.connect( myLocalTCWF,'OutputSpec.atlasToSubjectTransform',myLocalSegWF,'InputSpec.atlasToSubjectTransform')

        ### Now define where the final organized outputs should go.
        SEGMENTATION_DataSink=pe.Node(nio.DataSink(),name="SEGMENTATION_DS")
        SEGMENTATION_DataSink.inputs.base_directory=ExperimentBaseDirectoryResults
        SEGMENTATION_DataSink.inputs.regexp_substitutions = GenerateOutputPattern(subjectDatabaseFile,'BRAINSCut')
        T1T2WorkupSingle.connect(myLocalSegWF, 'OutputSpec.outputBinaryLeftAccumben',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryLeftAccumben')
        T1T2WorkupSingle.connect(myLocalSegWF, 'OutputSpec.outputBinaryRightAccumben',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryRightAccumben')
        T1T2WorkupSingle.connect(myLocalSegWF, 'OutputSpec.outputBinaryLeftCaudate',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryLeftCaudate')
        T1T2WorkupSingle.connect(myLocalSegWF, 'OutputSpec.outputBinaryRightCaudate',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryRightCaudate')
        T1T2WorkupSingle.connect(myLocalSegWF, 'OutputSpec.outputBinaryLeftGlobus',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryLeftGlobus')
        T1T2WorkupSingle.connect(myLocalSegWF, 'OutputSpec.outputBinaryRightGlobus',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryRightGlobus')
        T1T2WorkupSingle.connect(myLocalSegWF, 'OutputSpec.outputBinaryLeftHippocampus',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryLeftHippocampus')
        T1T2WorkupSingle.connect(myLocalSegWF, 'OutputSpec.outputBinaryRightHippocampus',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryRightHippocampus')
        T1T2WorkupSingle.connect(myLocalSegWF, 'OutputSpec.outputBinaryLeftPutamen',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryLeftPutamen')
        T1T2WorkupSingle.connect(myLocalSegWF, 'OutputSpec.outputBinaryRightPutamen',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryRightPutamen')
        T1T2WorkupSingle.connect(myLocalSegWF, 'OutputSpec.outputBinaryLeftThalamus',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryLeftThalamus')
        T1T2WorkupSingle.connect(myLocalSegWF, 'OutputSpec.outputBinaryRightThalamus',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryRightThalamus')
        T1T2WorkupSingle.connect(myLocalSegWF, 'OutputSpec.outputLabelImageName', SEGMENTATION_DataSink,'BRAINSCut.@outputLabelImageName')
        T1T2WorkupSingle.connect(myLocalSegWF, 'OutputSpec.outputCSVFileName', SEGMENTATION_DataSink,'BRAINSCut.@outputCSVFileName')

    if 'FREESURFER' in WORKFLOW_COMPONENTS:
        RunAllFSComponents=True ## A hack to avoid 26 hour run of freesurfer
        from WorkupT1T2FreeSurfer import CreateFreeSurferWorkflow
        myLocalFSWF= CreateFreeSurferWorkflow("Level1_FSTest",CLUSTER_QUEUE,RunAllFSComponents)
        T1T2WorkupSingle.connect(inputsSpec,'sessionid',myLocalFSWF,'InputSpec.subject_id')
        T1T2WorkupSingle.connect(myLocalTCWF,'OutputSpec.t1_corrected',myLocalFSWF,'InputSpec.T1_files')
        T1T2WorkupSingle.connect(myLocalTCWF,'OutputSpec.t2_corrected',myLocalFSWF,'InputSpec.T2_files')
        T1T2WorkupSingle.connect(myLocalTCWF,'OutputSpec.outputLabels',myLocalFSWF,'InputSpec.label_file')
        #T1T2WorkupSingle.connect(myLocalTCWF,'OutputSpec.outputLabels',myLocalFSWF,'InputSpec.mask_file') #Yes, the same file as label_file!

        ### Now define where the final organized outputs should go.
        if RunAllFSComponents == True:
            T1T2WorkupSingleDataSink=pe.Node(nio.DataSink(),name="FREESURFER_DS")
            T1T2WorkupSingleDataSink.inputs.base_directory=ExperimentBaseDirectoryResults
            T1T2WorkupSingleDataSink.inputs.regexp_substitutions = [
                ('/_uid_(?P<myuid>[^/]*)',r'/\g<myuid>')
                ]
            T1T2WorkupSingle.connect(myLocalFSWF, 'OutputSpec.FreesurferOutputDirectory', T1T2WorkupSingleDataSink,'FREESURFER_SUBJ.@FreesurferOutputDirectory')
        ### Now define where the final organized outputs should go.
        FSPREP_DataSink=pe.Node(nio.DataSink(),name="FREESURFER_PREP")
        FSPREP_DataSink.inputs.base_directory=ExperimentBaseDirectoryResults
        FREESURFER_PREP_PATTERNS = GenerateOutputPattern(subjectDatabaseFile,'FREESURFER_PREP')
        FSPREP_DataSink.inputs.regexp_substitutions = FREESURFER_PREP_PATTERNS
        print "========================="
        print "========================="
        print "========================="
        print FREESURFER_PREP_PATTERNS
        print "========================="
        print "========================="
        print "========================="
        T1T2WorkupSingle.connect(myLocalFSWF, 'OutputSpec.cnr_optimal_image', FSPREP_DataSink,'FREESURFER_PREP.@cnr_optimal_image')

    else:
        print "Skipping freesurfer"

    try:
        T1T2WorkupSingle.write_graph()
    except:
        pass

    return T1T2WorkupSingle

#############
## The following are just notes, and not really part of this script.
##
        #T1T2WorkupSingle.connect(myLocalLMIWF, 'OutputSpec.outputLandmarksInACPCAlignedSpace', T1T2WorkupSingleDataSink,'foo.@outputLandmarksInACPCAlignedSpace')
        #T1T2WorkupSingle.connect(myLocalLMIWF, 'OutputSpec.outputResampledVolume', T1T2WorkupSingleDataSink,'foo.@outputResampledVolume')
        #T1T2WorkupSingle.connect(myLocalLMIWF, 'OutputSpec.outputLandmarksInInputSpace', T1T2WorkupSingleDataSink,'foo.@outputLandmarksInInputSpace')
        #T1T2WorkupSingle.connect(myLocalLMIWF, 'OutputSpec.outputTransform', T1T2WorkupSingleDataSink,'foo.@outputTransform')
        #T1T2WorkupSingle.connect(myLocalLMIWF, 'OutputSpec.outputMRML', T1T2WorkupSingleDataSink,'foo.@outputMRML')
"""
    subs=r'test/\g<project>/\g<subject>/\g<session>'
pe.sub(subs,test)
pat=r'foo/_uid_(?P<project>PHD_[0-9][0-9][0-9])_(?P<subject>[0-9][0-9][0-9][0-9])_(?P<session>[0-9][0-9][0-9][0-9][0-9])'
pe=re.compile(pat)
pe.sub(subs,test)
test
test='foo/_uid_PHD_024_0003_12345'
pe.sub(subs,test)
pat=r'(?P<modulename>[^/]*)/_uid_(?P<project>PHD_[0-9][0-9][0-9])_(?P<subject>[0-9][0-9][0-9][0-9])_(?P<session>[0-9][0-9][0-9][0-9][0-9])'
subs=r'test/\g<project>/\g<subject>/\g<session>/\g<modulename>'
pe.sub(subs,test)
pe=re.compile(pat)
pe.sub(subs,test)

"""
