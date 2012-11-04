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

def GenerateWFName(projectid, subjectid, sessionid,processing_phase):
    return 'WF_'+str(subjectid)+"_"+str(sessionid)+"_"+str(projectid)+"_"+processing_phase

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
def MakeOneSubWorkFlow(projectid, subjectid, sessionid,processing_phase, WORKFLOW_COMPONENTS, BCD_model_path, InterpolationMode, CLUSTER_QUEUE):
    """
    Run autoworkup on a single Subject

    This is the main function to call when processing a data set with T1 & T2
    data.  ExperimentBaseDirectoryPrefix is the base of the directory to place results, T1Images & T2Images
    are the lists of images to be used in the auto-workup. atlas_fname_wpath is
    the path and filename of the atlas to use.
    """

    print "Building Pipeline for ",sessionid
    ########### PIPELINE INITIALIZATION #############
    T1T2WorkupSingle = pe.Workflow(name=GenerateWFName(projectid, subjectid, sessionid,processing_phase))

    inputsSpec = pe.Node(interface=IdentityInterface(fields=
                        ['sessionid','subjectid','projectid',
                         'allT1s',
                         'allT2s',
                         'allPDs',
                         'allFLs',
                         'allOthers',
                         'template_landmarks_31_fcsv',
                         'template_landmark_weights_31_csv',
                         'template_t1',
                         'atlasDefinition'
                         ]),
                         run_without_submitting=True,
                         name='inputspec' )

    outputsSpec = pe.Node(interface=IdentityInterface(fields=[
            't1_average','t2_average',
            'pd_average','fl_average',
            'posteriorImages','outputLabels',
            'TissueClassifyOutputDir',
            'TissueClassifyatlasToSubjectTransform',

#            'BCD_ACPC_T1',
            'BCD_ACPC_T1_CROPPED',
            'outputLandmarksInACPCAlignedSpace',
            'outputLandmarksInInputSpace',
            'outputTransform','LMIatlasToSubjectTransform'
            ]),
            run_without_submitting=True,
            name='outputspec' )

    if True: #'BASIC' in WORKFLOW_COMPONENTS:
        from WorkupT1T2LandmarkInitialization import CreateLandmarkInitializeWorkflow
        DoReverseMapping = False   # Set to true for debugging outputs
        if 'AUXLMK' in WORKFLOW_COMPONENTS:
            DoReverseMapping = True
        myLocalLMIWF= CreateLandmarkInitializeWorkflow("LandmarkInitialize", BCD_model_path, InterpolationMode,DoReverseMapping)

        T1T2WorkupSingle.connect( [ ( inputsSpec, myLocalLMIWF, [ ( ( 'allT1s', get_list_element, 0 ), 'inputspec.inputVolume') ] ), ] )
        T1T2WorkupSingle.connect( inputsSpec, 'template_landmarks_31_fcsv', myLocalLMIWF,'inputspec.atlasLandmarkFilename')
        T1T2WorkupSingle.connect( inputsSpec, 'template_landmark_weights_31_csv', myLocalLMIWF,'inputspec.atlasWeightFilename')
        if 'AUXLMK' in WORKFLOW_COMPONENTS:
            T1T2WorkupSingle.connect(inputsSpec,'template_t1',myLocalLMIWF,'inputspec.atlasVolume')
        ### Now connect outputspec
#        T1T2WorkupSingle.connect(myLocalLMIWF,'outputspec.outputResampledVolume', outputsSpec, 'BCD_ACPC_T1' )
        T1T2WorkupSingle.connect(myLocalLMIWF,'outputspec.outputResampledCroppedVolume', outputsSpec, 'BCD_ACPC_T1_CROPPED' )
        T1T2WorkupSingle.connect(myLocalLMIWF,'outputspec.outputLandmarksInACPCAlignedSpace',outputsSpec,'outputLandmarksInACPCAlignedSpace')
        T1T2WorkupSingle.connect(myLocalLMIWF,'outputspec.outputLandmarksInInputSpace',outputsSpec,'outputLandmarksInInputSpace')
        T1T2WorkupSingle.connect(myLocalLMIWF,'outputspec.outputTransform',outputsSpec,'outputTransform')
        T1T2WorkupSingle.connect(myLocalLMIWF,'outputspec.atlasToSubjectTransform',outputsSpec,'LMIatlasToSubjectTransform')

    if 'TISSUE_CLASSIFY' in WORKFLOW_COMPONENTS:
        from WorkupT1T2TissueClassifiy import CreateTissueClassifyWorkflow
        myLocalTCWF= CreateTissueClassifyWorkflow("TissueClassify",CLUSTER_QUEUE,InterpolationMode)
        T1T2WorkupSingle.connect( inputsSpec, 'allT1s', myLocalTCWF, 'inputspec.T1List')
        T1T2WorkupSingle.connect( inputsSpec, 'allT2s', myLocalTCWF, 'inputspec.T2List')
        T1T2WorkupSingle.connect( inputsSpec, 'allPDs', myLocalTCWF, 'inputspec.PDList')
        T1T2WorkupSingle.connect( inputsSpec, 'allFLs', myLocalTCWF, 'inputspec.FLList')
        T1T2WorkupSingle.connect( inputsSpec, 'allOthers', myLocalTCWF, 'inputspec.OtherList')
        T1T2WorkupSingle.connect( [ (inputsSpec, myLocalTCWF, [(('allT1s', getAllT1sLength), 'inputspec.T1_count')] ), ])
        T1T2WorkupSingle.connect( inputsSpec,'atlasDefinition',myLocalTCWF,'inputspec.atlasDefinition')
        T1T2WorkupSingle.connect( myLocalLMIWF, 'outputspec.outputResampledCroppedVolume', myLocalTCWF, 'inputspec.PrimaryT1' )
        T1T2WorkupSingle.connect( myLocalLMIWF,'outputspec.atlasToSubjectTransform',myLocalTCWF,'inputspec.atlasToSubjectInitialTransform')

        ### Now connect outputspec
        T1T2WorkupSingle.connect(myLocalTCWF, 'outputspec.t1_average', outputsSpec,'t1_average')
        T1T2WorkupSingle.connect(myLocalTCWF, 'outputspec.t2_average', outputsSpec,'t2_average')
        T1T2WorkupSingle.connect(myLocalTCWF, 'outputspec.pd_average', outputsSpec,'pd_average')
        T1T2WorkupSingle.connect(myLocalTCWF, 'outputspec.fl_average', outputsSpec,'fl_average')
        T1T2WorkupSingle.connect(myLocalTCWF, 'outputspec.posteriorImages', outputsSpec,'posteriorImages')
        T1T2WorkupSingle.connect(myLocalTCWF, 'outputspec.outputLabels', outputsSpec,'outputLabels')
        T1T2WorkupSingle.connect(myLocalTCWF, 'outputspec.TissueClassifyOutputDir', outputsSpec,'TissueClassifyOutputDir')
        T1T2WorkupSingle.connect(myLocalTCWF, 'outputspec.atlasToSubjectTransform', outputsSpec,'TissueClassifyatlasToSubjectTransform')

    return T1T2WorkupSingle
