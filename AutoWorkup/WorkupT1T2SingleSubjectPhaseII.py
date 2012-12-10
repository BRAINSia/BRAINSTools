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

from SEMTools import *
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
def getAllT1sLength(allT1s):
    return len(allT1s)
def get_list_element( nestedList, index ):
    return nestedList[index]
def MakeList(firstElement,secondElement):
    return [firstElement, secondElement]

def GenerateWFName(projectid, subjectid, sessionid):
    return 'WF_'+str(subjectid)+"_"+str(sessionid)+"_"+str(projectid)

def GenerateOutputPattern(projectid, subjectid, sessionid, DefaultNodeName, uidIsFirst):
    """ This function generates output path substitutions for workflows and nodes that conform to a common standard.
    """
    WFName=""#GenerateWFName(projectid,subjectid,sessionid)
    patternList = []
    if uidIsFirst == True:
        find_pat = os.path.join(DefaultNodeName, WFName)
    else:
        find_pat = os.path.join(WFName, DefaultNodeName)
    replace_pat = os.path.join(projectid, subjectid, sessionid, DefaultNodeName)
    patternList.append((find_pat,replace_pat))
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
def MakeOneSubWorkFlow(projectid, subjectid, sessionid, BAtlas, WORKFLOW_COMPONENTS, BCD_model_path, InterpolationMode, CLUSTER_QUEUE, ExperimentBaseDirectoryResults):
    """
    Run autoworkup on a single Subject

    This is the main function to call when processing a data set with T1 & T2
    data.  ExperimentBaseDirectoryPrefix is the base of the directory to place results, T1Images & T2Images
    are the lists of images to be used in the auto-workup. atlas_fname_wpath is
    the path and filename of the atlas to use.
    """

    # print "Building Pipeline for ",sessionid
    ########### PIPELINE INITIALIZATION #############
    T1T2WorkupSingle = pe.Workflow(name=GenerateWFName(projectid, subjectid, sessionid))

    inputsSpec = pe.Node(interface=IdentityInterface(fields=
                        ['sessionid','subjectid','projectid',
                         'allT1s',
                         'allT2s',
                         'allPDs',
                         'allOthers'
                         ]),
                         run_without_submitting=True,
                         name='inputspec' )

    outputsSpec = pe.Node(interface=IdentityInterface(fields=['BCD_ACPC_T1',
            't1_average','t2_average'
            ]),
            run_without_submitting=True,
            name='outputspec' )

    if 'BASIC' in WORKFLOW_COMPONENTS:
        from WorkupT1T2LandmarkInitialization import CreateLandmarkInitializeWorkflow
        DoReverseMapping = False   # Set to true for debugging outputs
        if 'AUXLMK' in WORKFLOW_COMPONENTS:
            DoReverseMapping = True
        myLocalLMIWF= CreateLandmarkInitializeWorkflow("LandmarkInitialize", BCD_model_path, InterpolationMode,DoReverseMapping)

        T1T2WorkupSingle.connect( [ ( inputsSpec, myLocalLMIWF, [ ( ( 'allT1s', get_list_element, 0 ), 'inputspec.inputVolume') ] ), ] )
        T1T2WorkupSingle.connect( BAtlas, 'template_landmarks_31_fcsv', myLocalLMIWF,'inputspec.atlasLandmarkFilename')
        T1T2WorkupSingle.connect( BAtlas, 'template_landmark_weights_31_csv', myLocalLMIWF,'inputspec.atlasWeightFilename')
        if 'AUXLMK' in WORKFLOW_COMPONENTS:
            T1T2WorkupSingle.connect(BAtlas,'template_t1',myLocalLMIWF,'inputspec.atlasVolume')

        ### Now define where the final organized outputs should go.
        BASIC_DataSink=pe.Node(nio.DataSink(),name="BASIC_DS")
        BASIC_DataSink.inputs.base_directory=ExperimentBaseDirectoryResults
        BASIC_DataSink.inputs.regexp_substitutions = GenerateOutputPattern(projectid, subjectid, sessionid,'ACPCAlign',False)

        T1T2WorkupSingle.connect(myLocalLMIWF,'outputspec.outputLandmarksInACPCAlignedSpace',BASIC_DataSink,'ACPCAlign.@outputLandmarksInACPCAlignedSpace')
        T1T2WorkupSingle.connect(myLocalLMIWF,'outputspec.outputResampledVolume',BASIC_DataSink,'ACPCAlign.@outputResampledVolume')
        T1T2WorkupSingle.connect(myLocalLMIWF,'outputspec.outputLandmarksInInputSpace',BASIC_DataSink,'ACPCAlign.@outputLandmarksInInputSpace')
        T1T2WorkupSingle.connect(myLocalLMIWF,'outputspec.outputTransform',BASIC_DataSink,'ACPCAlign.@outputTransform')
        T1T2WorkupSingle.connect(myLocalLMIWF,'outputspec.atlasToSubjectTransform',BASIC_DataSink,'ACPCAlign.@atlasToSubjectTransform')
        ### Now connect outputspec
        T1T2WorkupSingle.connect( myLocalLMIWF, 'outputspec.outputResampledVolume', outputsSpec, 'BCD_ACPC_T1' )

    if 'TISSUE_CLASSIFY' in WORKFLOW_COMPONENTS:
        from WorkupT1T2TissueClassifiy import CreateTissueClassifyWorkflow
        myLocalTCWF = CreateTissueClassifyWorkflow("TissueClassify", CLUSTER_QUEUE, InterpolationMode)
        T1T2WorkupSingle.connect( inputsSpec, 'allT1s', myLocalTCWF, 'inputspec.T1List')
        T1T2WorkupSingle.connect( inputsSpec, 'allT2s', myLocalTCWF, 'inputspec.T2List')
        T1T2WorkupSingle.connect( inputsSpec, 'allPDs', myLocalTCWF, 'inputspec.PDList')
        T1T2WorkupSingle.connect( inputsSpec, 'allOthers', myLocalTCWF, 'inputspec.OtherList')
        T1T2WorkupSingle.connect( [ (inputsSpec, myLocalTCWF, [(('allT1s', getAllT1sLength), 'inputspec.T1_count')] ), ])
        T1T2WorkupSingle.connect( BAtlas,'ExtendedAtlasDefinition.xml',myLocalTCWF,'inputspec.atlasDefinition')
        T1T2WorkupSingle.connect( myLocalLMIWF, 'outputspec.outputResampledVolume', myLocalTCWF, 'inputspec.PrimaryT1' )
        T1T2WorkupSingle.connect( myLocalLMIWF,'outputspec.atlasToSubjectTransform',myLocalTCWF,'inputspec.atlasToSubjectInitialTransform')

        ### Now define where the final organized outputs should go.
        ###     For posterior probability files, we need to use a MapNode for the keys from the
        ### myLocalTCWF.outputspec.posteriorImages dictionary.  To use a MapNode, we must know
        ### the list BEFORE we run the pipeline...
        from PipeLineFunctionHelpers import POSTERIORS
        TC_DataSink = pe.MapNode(nio.DataSink(infields=['TissueClassify.@posteriors']), name="TISSUE_CLASSIFY_DS")
        TC_DataSink.inputs.base_directory = ExperimentBaseDirectoryResults
        TC_DataSink.inputs.regexp_substitutions = GenerateOutputPattern(projectid, subjectid,
                                                                        sessionid, 'TissueClassify',
                                                                        False)
        setattr(TC_DataSink.inputs, 'TissueClassify.@posteriors',
                ['POSTERIOR_%s.nii.gz' % label for label in POSTERIORS])
        T1T2WorkupSingle.connect(myLocalTCWF, 'outputspec.t1_average', TC_DataSink, 'TissueClassify.@t1_average')
        T1T2WorkupSingle.connect(myLocalTCWF, 'outputspec.t2_average', TC_DataSink, 'TissueClassify.@t2_average')
        T1T2WorkupSingle.connect(myLocalTCWF, 'outputspec.posteriorImages', TC_DataSink, 'TissueClassify.@posteriorImages')
        ### Now connect outputspec
        T1T2WorkupSingle.connect(TC_DataSink, 'TissueClassify.@t1_average', outputsSpec,'t1_average')
        T1T2WorkupSingle.connect(TC_DataSink, 'TissueClassify.@t2_average', outputsSpec,'t2_average')

    ## Make deformed Atlas image space
    if 'ANTS' in WORKFLOW_COMPONENTS:
        from WorkupT1T2ANTS import CreateANTSRegistrationWorkflow
        myLocalAntsWF = CreateANTSRegistrationWorkflow("ANTSRegistration",CLUSTER_QUEUE,-1)
        T1T2WorkupSingle.connect( myLocalTCWF,'outputspec.t1_average',myLocalAntsWF,"inputspec.fixedVolumesList")
        T1T2WorkupSingle.connect( BAtlas,'template_t1',    myLocalAntsWF,"inputspec.movingVolumesList")
        T1T2WorkupSingle.connect(myLocalLMIWF,'outputspec.atlasToSubjectTransform',myLocalAntsWF,'inputspec.initial_moving_transform')
        # Must register the entire head, not just the brain!
        T1T2WorkupSingle.connect(myLocalTCWF,'outputspec.outputHeadLabels',myLocalAntsWF,'inputspec.fixedBinaryVolume')
        T1T2WorkupSingle.connect(BAtlas,'template_headregion',myLocalAntsWF,'inputspec.movingBinaryVolume')

        ### Now define where the final organized outputs should go.
        ANTS_DataSink=pe.Node(nio.DataSink(),name="ANTSRegistration_DS")
        ANTS_DataSink.inputs.base_directory=ExperimentBaseDirectoryResults
        ANTS_DataSink.inputs.regexp_substitutions = GenerateOutputPattern(projectid, subjectid, sessionid,'ANTSRegistration',False)
        T1T2WorkupSingle.connect(myLocalAntsWF, 'outputspec.warped_image', ANTS_DataSink,'ANTSRegistration.@warped_image')
        T1T2WorkupSingle.connect(myLocalAntsWF, 'outputspec.inverse_warped_image', ANTS_DataSink,'ANTSRegistration.@inverse_warped_image')
        T1T2WorkupSingle.connect(myLocalAntsWF, 'outputspec.affine_transform', ANTS_DataSink,'ANTSRegistration.@affine_transform')
        T1T2WorkupSingle.connect(myLocalAntsWF, 'outputspec.warp_transform', ANTS_DataSink,'ANTSRegistration.@warp_transform')
        T1T2WorkupSingle.connect(myLocalAntsWF, 'outputspec.inverse_warp_transform', ANTS_DataSink,'ANTSRegistration.@inverse_warp_transform')

    if 'SEGMENTATION' in WORKFLOW_COMPONENTS:
        from WorkupT1T2BRAINSCut import CreateBRAINSCutWorkflow
        ## TODO:  Remove BAtlas From Here as well!
        myLocalSegWF = CreateBRAINSCutWorkflow("Segmentation",CLUSTER_QUEUE,BAtlas) ##Note:  Passing in the entire BAtlas Object here!
        T1T2WorkupSingle.connect( myLocalTCWF,'outputspec.t1_average',myLocalSegWF,'inputspec.T1Volume')
        T1T2WorkupSingle.connect( myLocalTCWF,'outputspec.t2_average',myLocalSegWF,'inputspec.T2Volume')
        T1T2WorkupSingle.connect( myLocalTCWF,'outputspec.atlasToSubjectTransform',myLocalSegWF,'inputspec.atlasToSubjectTransform')

        ### Now define where the final organized outputs should go.
        SEGMENTATION_DataSink=pe.Node(nio.DataSink(),name="SEGMENTATION_DS")
        SEGMENTATION_DataSink.inputs.base_directory=ExperimentBaseDirectoryResults
        SEGMENTATION_DataSink.inputs.regexp_substitutions = GenerateOutputPattern(projectid, subjectid, sessionid,'BRAINSCut',False)
        T1T2WorkupSingle.connect(myLocalSegWF, 'outputspec.outputBinaryLeftAccumben',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryLeftAccumben')
        T1T2WorkupSingle.connect(myLocalSegWF, 'outputspec.outputBinaryRightAccumben',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryRightAccumben')
        T1T2WorkupSingle.connect(myLocalSegWF, 'outputspec.outputBinaryLeftCaudate',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryLeftCaudate')
        T1T2WorkupSingle.connect(myLocalSegWF, 'outputspec.outputBinaryRightCaudate',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryRightCaudate')
        T1T2WorkupSingle.connect(myLocalSegWF, 'outputspec.outputBinaryLeftGlobus',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryLeftGlobus')
        T1T2WorkupSingle.connect(myLocalSegWF, 'outputspec.outputBinaryRightGlobus',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryRightGlobus')
        T1T2WorkupSingle.connect(myLocalSegWF, 'outputspec.outputBinaryLeftHippocampus',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryLeftHippocampus')
        T1T2WorkupSingle.connect(myLocalSegWF, 'outputspec.outputBinaryRightHippocampus',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryRightHippocampus')
        T1T2WorkupSingle.connect(myLocalSegWF, 'outputspec.outputBinaryLeftPutamen',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryLeftPutamen')
        T1T2WorkupSingle.connect(myLocalSegWF, 'outputspec.outputBinaryRightPutamen',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryRightPutamen')
        T1T2WorkupSingle.connect(myLocalSegWF, 'outputspec.outputBinaryLeftThalamus',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryLeftThalamus')
        T1T2WorkupSingle.connect(myLocalSegWF, 'outputspec.outputBinaryRightThalamus',SEGMENTATION_DataSink, 'BRAINSCut.@outputBinaryRightThalamus')
        T1T2WorkupSingle.connect(myLocalSegWF, 'outputspec.outputLabelImageName', SEGMENTATION_DataSink,'BRAINSCut.@outputLabelImageName')
        T1T2WorkupSingle.connect(myLocalSegWF, 'outputspec.outputCSVFileName', SEGMENTATION_DataSink,'BRAINSCut.@outputCSVFileName')

    if 'FREESURFER' in WORKFLOW_COMPONENTS:
        RunAllFSComponents=True ## A hack to avoid 26 hour run of freesurfer
        from WorkupT1T2FreeSurfer import CreateFreeSurferWorkflow
        myLocalFSWF= CreateFreeSurferWorkflow("Level1_FSTest",CLUSTER_QUEUE,RunAllFSComponents)
        T1T2WorkupSingle.connect(inputsSpec,'sessionid',myLocalFSWF,'inputspec.subject_id')
        T1T2WorkupSingle.connect(myLocalTCWF,'outputspec.t1_average',myLocalFSWF,'inputspec.T1_files')
        T1T2WorkupSingle.connect(myLocalTCWF,'outputspec.t2_average',myLocalFSWF,'inputspec.T2_files')
        T1T2WorkupSingle.connect(myLocalTCWF,'outputspec.outputLabels',myLocalFSWF,'inputspec.label_file')
        #T1T2WorkupSingle.connect(myLocalTCWF,'outputspec.outputLabels',myLocalFSWF,'inputspec.mask_file') #Yes, the same file as label_file!

        ### Now define where the final organized outputs should go.
        if RunAllFSComponents == True:
            T1T2WorkupSingleDataSink=pe.Node(nio.DataSink(),name="FREESURFER_DS")
            T1T2WorkupSingleDataSink.inputs.base_directory=ExperimentBaseDirectoryResults
            T1T2WorkupSingleDataSink.inputs.regexp_substitutions = [
                ('/_uid_(?P<myuid>[^/]*)',r'/\g<myuid>')
                ]
            T1T2WorkupSingle.connect(myLocalFSWF, 'outputspec.FreesurferOutputDirectory', T1T2WorkupSingleDataSink,'FREESURFER_SUBJ.@FreesurferOutputDirectory')
        ### Now define where the final organized outputs should go.
        FSPREP_DataSink=pe.Node(nio.DataSink(),name="FREESURFER_PREP")
        FSPREP_DataSink.inputs.base_directory=ExperimentBaseDirectoryResults
        FREESURFER_PREP_PATTERNS = GenerateOutputPattern(projectid, subjectid, sessionid,'FREESURFER_PREP',False)
        FSPREP_DataSink.inputs.regexp_substitutions = FREESURFER_PREP_PATTERNS
        print "========================="
        print "========================="
        print "========================="
        print FREESURFER_PREP_PATTERNS
        print "========================="
        print "========================="
        print "========================="
        T1T2WorkupSingle.connect(myLocalFSWF, 'outputspec.cnr_optimal_image', FSPREP_DataSink,'FREESURFER_PREP.@cnr_optimal_image')
    else:
        pass
        #print "Skipping freesurfer"
    """
    try:
        T1T2WorkupSingle.write_graph()
    except:
        pass
    """

    return T1T2WorkupSingle

#############
## The following are just notes, and not really part of this script.
##
        #T1T2WorkupSingle.connect(myLocalLMIWF, 'outputspec.outputLandmarksInACPCAlignedSpace', T1T2WorkupSingleDataSink,'foo.@outputLandmarksInACPCAlignedSpace')
        #T1T2WorkupSingle.connect(myLocalLMIWF, 'outputspec.outputResampledVolume', T1T2WorkupSingleDataSink,'foo.@outputResampledVolume')
        #T1T2WorkupSingle.connect(myLocalLMIWF, 'outputspec.outputLandmarksInInputSpace', T1T2WorkupSingleDataSink,'foo.@outputLandmarksInInputSpace')
        #T1T2WorkupSingle.connect(myLocalLMIWF, 'outputspec.outputTransform', T1T2WorkupSingleDataSink,'foo.@outputTransform')
        #T1T2WorkupSingle.connect(myLocalLMIWF, 'outputspec.outputMRML', T1T2WorkupSingleDataSink,'foo.@outputMRML')
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
