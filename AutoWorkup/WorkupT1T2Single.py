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

def GenerateWFName(projectid, subjectid, sessionid,processing_phase):
    return 'WF_'+str(subjectid)+"_"+str(sessionid)+"_"+str(projectid)+processing_phase

def GenerateOutputPattern(projectid, subjectid, sessionid,DefaultNodeName,uidIsFirst):
    """ This function generates output path substitutions for workflows and nodes that conform to a common standard.
    """
    WFName=""#GenerateWFName(projectid,subjectid,sessionid)
    patternList=[]
    if uidIsFirst == True:
        find_pat=os.path.join(DefaultNodeName,WFName)
    else:
        find_pat=os.path.join(WFName,DefaultNodeName)
    replace_pat=os.path.join(projectid,subjectid,sessionid,DefaultNodeName)
    patternList.append( (find_pat,replace_pat) )
    print "HACK: ", patternList
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
def MakeOneSubWorkFlow(projectid, subjectid, sessionid,processing_phase, WORKFLOW_COMPONENTS, BCD_model_path, InterpolationMode, CLUSTER_QUEUE, ExperimentBaseDirectoryResults):
    """
    Run autoworkup on a single Subject

    This is the main function to call when processing a data set with T1 & T2
    data.  ExperimentBaseDirectoryPrefix is the base of the directory to place results, T1Images & T2Images
    are the lists of images to be used in the auto-workup. atlas_fname_wpath is
    the path and filename of the atlas to use.
    """

    # print "Building Pipeline for ",sessionid
    ########### PIPELINE INITIALIZATION #############
    T1T2WorkupSingle = pe.Workflow(name=GenerateWFName(projectid, subjectid, sessionid,processing_phase))

    inputsSpec = pe.Node(interface=IdentityInterface(fields=
                        ['sessionid','subjectid','projectid',
                         'allT1s',
                         'allT2s',
                         'allPDs',
                         'allOthers',
                         'template_landmarks_31_fcsv',
                         'template_landmark_weights_31_csv',
                         'template_t1',
                         'AtlasPVDefinition_xml'
                         ]),
                         run_without_submitting=True,
                         name='InputSpec' )

    outputsSpec = pe.Node(interface=IdentityInterface(fields=['BCD_ACPC_T1',
            't1_average','t2_average',
            'posteriorImages','outputLabels'
            ]),
            run_without_submitting=True,
            name='OutputSpec' )

    if 'BASIC' in WORKFLOW_COMPONENTS:
        from WorkupT1T2LandmarkInitialization import CreateLandmarkInitializeWorkflow
        DoReverseMapping = False   # Set to true for debugging outputs
        if 'AUXLMK' in WORKFLOW_COMPONENTS:
            DoReverseMapping = True
        myLocalLMIWF= CreateLandmarkInitializeWorkflow("LandmarkInitialize", BCD_model_path, InterpolationMode,DoReverseMapping)

        def get_list_element( nestedList, index ):
            return nestedList[index]
        T1T2WorkupSingle.connect( [ ( inputsSpec, myLocalLMIWF, [ ( ( 'allT1s', get_list_element, 0 ), 'InputSpec.inputVolume') ] ), ] )
        T1T2WorkupSingle.connect( inputsSpec, 'template_landmarks_31_fcsv', myLocalLMIWF,'InputSpec.atlasLandmarkFilename')
        T1T2WorkupSingle.connect( inputsSpec, 'template_landmark_weights_31_csv', myLocalLMIWF,'InputSpec.atlasWeightFilename')
        if 'AUXLMK' in WORKFLOW_COMPONENTS:
            T1T2WorkupSingle.connect(inputsSpec,'template_t1',myLocalLMIWF,'InputSpec.atlasVolume')

        ### Now define where the final organized outputs should go.
        BASIC_DataSink=pe.Node(nio.DataSink(),name="BASIC_DS")
        BASIC_DataSink.inputs.base_directory=ExperimentBaseDirectoryResults
        BASIC_DataSink.inputs.regexp_substitutions = GenerateOutputPattern(projectid, subjectid, sessionid,'ACPCAlign',False)

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
        T1T2WorkupSingle.connect( inputsSpec, 'allT1s', myLocalTCWF, 'InputSpec.T1List')
        T1T2WorkupSingle.connect( inputsSpec, 'allT2s', myLocalTCWF, 'InputSpec.T2List')
        T1T2WorkupSingle.connect( inputsSpec, 'allPDs', myLocalTCWF, 'InputSpec.PDList')
        T1T2WorkupSingle.connect( inputsSpec, 'allOthers', myLocalTCWF, 'InputSpec.OtherList')
        def getAllT1sLength(allT1s):
            return len(allT1s)
        T1T2WorkupSingle.connect( [ (inputsSpec, myLocalTCWF, [(('allT1s', getAllT1sLength), 'InputSpec.T1_count')] ), ])
        T1T2WorkupSingle.connect( inputsSpec,'AtlasPVDefinition_xml',myLocalTCWF,'InputSpec.atlasDefinition')
        T1T2WorkupSingle.connect( myLocalLMIWF, 'OutputSpec.outputResampledVolume', myLocalTCWF, 'InputSpec.PrimaryT1' )
        T1T2WorkupSingle.connect( myLocalLMIWF,'OutputSpec.atlasToSubjectTransform',myLocalTCWF,'InputSpec.atlasToSubjectInitialTransform')

        ### Now define where the final organized outputs should go.
        TC_DataSink=pe.Node(nio.DataSink(),name="TISSUE_CLASSIFY_DS")
        TC_DataSink.inputs.base_directory=ExperimentBaseDirectoryResults
        TC_DataSink.inputs.regexp_substitutions = GenerateOutputPattern(projectid, subjectid, sessionid,'TissueClassify',False)
        T1T2WorkupSingle.connect(myLocalTCWF, 'OutputSpec.TissueClassifyOutputDir', TC_DataSink,'TissueClassify.@TissueClassifyOutputDir')

        ### Now connect OutputSpec
        T1T2WorkupSingle.connect(myLocalTCWF, 'OutputSpec.t1_average', outputsSpec,'t1_average')
        T1T2WorkupSingle.connect(myLocalTCWF, 'OutputSpec.t2_average', outputsSpec,'t2_average')
        T1T2WorkupSingle.connect(myLocalTCWF, 'OutputSpec.posteriorImages', outputsSpec,'posteriorImages')
        T1T2WorkupSingle.connect(myLocalTCWF, 'OutputSpec.outputLabels', outputsSpec,'outputLabels')

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
