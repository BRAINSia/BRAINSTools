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
# from nipype.utils.config import config
# config.set('logging', 'log_to_file', 'false')
# config.set_log_dir(os.getcwd())
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
# package_check('nipype', '5.4', 'tutorial1') ## HACK: Check nipype version
package_check('numpy', '1.3', 'tutorial1')
package_check('scipy', '0.7', 'tutorial1')
package_check('networkx', '1.0', 'tutorial1')
package_check('IPython', '0.10', 'tutorial1')

from SEMTools import *

from WorkupT1T2AtlasNode import MakeAtlasNode
from PipeLineFunctionHelpers import getListIndex

# HACK:  [('buildTemplateIteration2', 'SUBJECT_TEMPLATES/0249/buildTemplateIteration2')]


def GenerateSubjectOutputPattern(subjectid):
    """ This function generates output path substitutions for workflows and nodes that conform to a common standard.
HACK:  [('ANTSTemplate/Iteration02_Reshaped.nii.gz', 'SUBJECT_TEMPLATES/0668/T1_RESHAPED.nii.gz'),
('ANTSTemplate/_ReshapeAveragePassiveImageWithShapeUpdate[0-9]*', 'SUBJECT_TEMPLATES/0668')]
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
    patternList = []

    find_pat = os.path.join('ANTSTemplate', 'Iteration02_Reshaped.nii.gz')
    replace_pat = os.path.join('SUBJECT_TEMPLATES', subjectid, r'AVG_T1.nii.gz')
    patternList.append((find_pat, replace_pat))

    find_pat = os.path.join('ANTSTemplate', r'_ReshapeAveragePassiveImageWithShapeUpdate[0-9]*/AVG_[A-Z0-9]*WARP_(?P<structure>AVG_[A-Z0-9]*.nii.gz)')
    replace_pat = os.path.join('SUBJECT_TEMPLATES', subjectid, r'\g<structure>')
    patternList.append((find_pat, replace_pat))

    find_pat = os.path.join('ANTSTemplate', r'CLIPPED_AVG_[A-Z]*WARP_(?P<structure>AVG_[A-Z]*.nii.gz)')
    replace_pat = os.path.join('SUBJECT_TEMPLATES', subjectid, r'\g<structure>')
    patternList.append((find_pat, replace_pat))

    print "HACK: ", patternList
    return patternList


def GenerateOutputPattern(projectid, subjectid, sessionid, DefaultNodeName):
    """ This function generates output path substitutions for workflows and nodes that conform to a common standard.
    """
    patternList = []
    find_pat = os.path.join(DefaultNodeName)
    replace_pat = os.path.join(projectid, subjectid, sessionid, DefaultNodeName)
    patternList.append((find_pat, replace_pat))
    print "HACK: ", patternList
    return patternList


def GenerateAccumulatorImagesOutputPattern(projectid, subjectid, sessionid):
    """ This function generates output path substitutions for workflows and nodes that conform to a common standard.
    """
    patternList = []
    find_pat = "POSTERIOR_"
    replace_pat = os.path.join(projectid, subjectid, sessionid) + "/POSTERIOR_"
    patternList.append((find_pat, replace_pat))
    print "HACK: ", patternList
    return patternList

## This takes several lists and merges them, but it also removes all empty values from the lists


def MergeByExtendListElements(t2_averageList, pd_averageList, fl_averageList, outputLabels_averageList, ListOfPosteriorImagesDictionary):
    for t2_index in range(0, len(t2_averageList)):
        if t2_averageList[t2_index] is not None:
            ListOfPosteriorImagesDictionary[t2_index]['T2'] = t2_averageList[t2_index]
        if pd_averageList[t2_index] is not None:
            ListOfPosteriorImagesDictionary[t2_index]['PD'] = pd_averageList[t2_index]
        if fl_averageList[t2_index] is not None:
            ListOfPosteriorImagesDictionary[t2_index]['FL'] = fl_averageList[t2_index]
        if outputLabels_averageList[t2_index] is not None:
            ListOfPosteriorImagesDictionary[t2_index]['BRAINMASK'] = outputLabels_averageList[t2_index]
    return ListOfPosteriorImagesDictionary


def MakeNewAtlasTemplate(t1_image, deformed_list,
                         AtlasTemplate, outDefinition):
    import os
    import sys
    import SimpleITK as sitk
    patternDict = {
        'AVG_AIRWARP_AVG_AIR.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_AIR.nii.gz',
        'AVG_BGMWARP_AVG_BGM.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_BASALTISSUE.nii.gz',
        'AVG_CRBLGMWARP_AVG_CRBLGM.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_CRBLGM.nii.gz',
        'AVG_CRBLWMWARP_AVG_CRBLWM.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_CRBLWM.nii.gz',
        'AVG_CSFWARP_AVG_CSF.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_CSF.nii.gz',
        'AVG_NOTCSFWARP_AVG_NOTCSF.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_NOTCSF.nii.gz',
        'AVG_NOTGMWARP_AVG_NOTGM.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_NOTGM.nii.gz',
        'AVG_NOTVBWARP_AVG_NOTVB.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_NOTVB.nii.gz',
        'AVG_NOTWMWARP_AVG_NOTWM.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_NOTWM.nii.gz',
        'AVG_SURFGMWARP_AVG_SURFGM.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_SURFGM.nii.gz',
        'AVG_VBWARP_AVG_VB.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_VB.nii.gz',
        'AVG_WMWARP_AVG_WM.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_WM.nii.gz',
        'AVG_ACCUMBENWARP_AVG_ACCUMBEN.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_ACCUMBEN.nii.gz',
        'AVG_CAUDATEWARP_AVG_CAUDATE.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_CAUDATE.nii.gz',
        'AVG_PUTAMENWARP_AVG_PUTAMEN.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_PUTAMEN.nii.gz',
        'AVG_GLOBUSWARP_AVG_GLOBUS.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_GLOBUS.nii.gz',
        'AVG_THALAMUSWARP_AVG_THALAMUS.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_THALAMUS.nii.gz',
        'AVG_HIPPOCAMPUSWARP_AVG_HIPPOCAMPUS.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_HIPPOCAMPUS.nii.gz',
        'AVG_T2WARP_AVG_T2.nii.gz': '@ATLAS_DIRECTORY@/template_t2.nii.gz',
        'AVG_BRAINMASKWARP_AVG_BRAINMASK.nii.gz': '@ATLAS_DIRECTORY@/template_brain.nii.gz',
        'T1_RESHAPED.nii.gz': '@ATLAS_DIRECTORY@/template_t1.nii.gz'
    }
    templateFile = open(AtlasTemplate, 'r')
    content = templateFile.read()              # read entire file into memory
    templateFile.close()

    ## Now clean up the posteriors based on anatomical knowlege.
    ## sometimes the posteriors are not relevant for priors
    ## due to anomolies around the edges.
    load_images_list = dict()
    for full_pathname in deformed_list:
        base_name = os.path.basename(full_pathname)
        if base_name in patternDict.keys():
            load_images_list[base_name] = sitk.ReadImage(full_pathname)
    ## Make binary dilated mask
    binmask = sitk.BinaryThreshold(load_images_list['AVG_BRAINMASKWARP_AVG_BRAINMASK.nii.gz'], 1, 1000000)
    dilated5 = sitk.DilateObjectMorphology(binmask, 5)
    dilated5 = sitk.Cast(dilated5, sitk.sitkFloat32)  # Convert to Float32 for multiply
    ## Now clip the interior brain mask with dilated5
    interiorPriors = [
        'AVG_BGMWARP_AVG_BGM.nii.gz',
        'AVG_CRBLGMWARP_AVG_CRBLGM.nii.gz',
        'AVG_CRBLWMWARP_AVG_CRBLWM.nii.gz',
        'AVG_CSFWARP_AVG_CSF.nii.gz',
        'AVG_SURFGMWARP_AVG_SURFGM.nii.gz',
        'AVG_VBWARP_AVG_VB.nii.gz',
        'AVG_WMWARP_AVG_WM.nii.gz',
        'AVG_ACCUMBENWARP_AVG_ACCUMBEN.nii.gz',
        'AVG_CAUDATEWARP_AVG_CAUDATE.nii.gz',
        'AVG_PUTAMENWARP_AVG_PUTAMEN.nii.gz',
        'AVG_GLOBUSWARP_AVG_GLOBUS.nii.gz',
        'AVG_THALAMUSWARP_AVG_THALAMUS.nii.gz',
        'AVG_HIPPOCAMPUSWARP_AVG_HIPPOCAMPUS.nii.gz',
    ]
    clean_deformed_list = deformed_list
    for index in range(0, len(deformed_list)):
        full_pathname = deformed_list[index]
        base_name = os.path.basename(full_pathname)
        if base_name == 'AVG_BRAINMASKWARP_AVG_BRAINMASK.nii.gz':
            ### Make Brain Mask Binary
            clipped_name = 'CLIPPED_' + base_name
            patternDict[clipped_name] = patternDict[base_name]
            sitk.WriteImage(binmask, clipped_name)
            clean_deformed_list[index] = os.path.realpath(clipped_name)
        if base_name in interiorPriors:
            ### Make clipped posteriors for brain regions
            curr = sitk.Cast(sitk.ReadImage(full_pathname), sitk.sitkFloat32)
            curr = curr * dilated5
            clipped_name = 'CLIPPED_' + base_name
            patternDict[clipped_name] = patternDict[base_name]
            sitk.WriteImage(curr, clipped_name)
            clean_deformed_list[index] = os.path.realpath(clipped_name)
            print "HACK: ", clean_deformed_list[index]
            curr = None
    binmask = None
    dilated5 = None

    for full_pathname in clean_deformed_list:
        base_name = os.path.basename(full_pathname)
        if base_name in patternDict.keys():
            content = content.replace(patternDict[base_name], full_pathname)
    content = content.replace('@ATLAS_DIRECTORY@/template_t1.nii.gz', t1_image)
    ## NOTE:  HEAD REGION CAN JUST BE T1 image.
    content = content.replace('@ATLAS_DIRECTORY@/template_headregion.nii.gz', t1_image)
    ## NOTE:  BRAIN REGION CAN JUST BE the label images.
    outAtlasFullPath = os.path.realpath(outDefinition)
    newFile = open(outAtlasFullPath, 'w')
    newFile.write(content)  # write the file with the text substitution
    newFile.close()
    return outAtlasFullPath, clean_deformed_list


def AccumulateLikeTissuePosteriors(posteriorImages):
    import os
    import sys
    import SimpleITK as sitk
    ## Now clean up the posteriors based on anatomical knowlege.
    ## sometimes the posteriors are not relevant for priors
    ## due to anomolies around the edges.
    load_images_list = dict()
    for full_pathname in posteriorImages.values():
        base_name = os.path.basename(full_pathname)
        load_images_list[base_name] = sitk.ReadImage(full_pathname)
    GM_ACCUM = [
        'POSTERIOR_ACCUMBEN.nii.gz',
        'POSTERIOR_CAUDATE.nii.gz',
        'POSTERIOR_CRBLGM.nii.gz',
        'POSTERIOR_HIPPOCAMPUS.nii.gz',
        'POSTERIOR_PUTAMEN.nii.gz',
        'POSTERIOR_THALAMUS.nii.gz',
        'POSTERIOR_SURFGM.nii.gz',
    ]
    WM_ACCUM = [
        'POSTERIOR_CRBLWM.nii.gz',
        'POSTERIOR_WM.nii.gz'
    ]
    CSF_ACCUM = [
        'POSTERIOR_CSF.nii.gz',
    ]
    VB_ACCUM = [
        'POSTERIOR_VB.nii.gz',
    ]
    GLOBUS_ACCUM = [
        'POSTERIOR_GLOBUS.nii.gz',
    ]
    BACKGROUND_ACCUM = [
        'POSTERIOR_AIR.nii.gz',
        'POSTERIOR_NOTCSF.nii.gz',
        'POSTERIOR_NOTGM.nii.gz',
        'POSTERIOR_NOTVB.nii.gz',
        'POSTERIOR_NOTWM.nii.gz',
    ]
    ## The next 2 items MUST be syncronized
    AccumulatePriorsNames = ['POSTERIOR_GM_TOTAL.nii.gz', 'POSTERIOR_WM_TOTAL.nii.gz',
                             'POSTERIOR_CSF_TOTAL.nii.gz', 'POSTERIOR_VB_TOTAL.nii.gz',
                             'POSTERIOR_GLOBUS_TOTAL.nii.gz', 'POSTERIOR_BACKGROUND_TOTAL.nii.gz']
    ForcedOrderingLists = [GM_ACCUM, WM_ACCUM, CSF_ACCUM, VB_ACCUM, GLOBUS_ACCUM, BACKGROUND_ACCUM]
    AccumulatePriorsList = list()
    for index in range(0, len(ForcedOrderingLists)):
        outname = AccumulatePriorsNames[index]
        inlist = ForcedOrderingLists[index]
        accum_image = load_images_list[inlist[0]]  # copy first image
        for curr_image in range(1, len(inlist)):
            accum_image = accum_image + load_images_list[inlist[curr_image]]
        sitk.WriteImage(accum_image, outname)
        AccumulatePriorsList.append(os.path.realpath(outname))
    print "HACK \n\n\n\n\n\n\n HACK \n\n\n: {APL}\n".format(APL=AccumulatePriorsList)
    print ": {APN}\n".format(APN=AccumulatePriorsNames)
    return AccumulatePriorsList, AccumulatePriorsNames
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


def ShortWorkupT1T2(subjectid, mountPrefix, ExperimentBaseDirectoryCache, ExperimentBaseDirectoryResults, ExperimentDatabase, atlas_fname_wpath, BCD_model_path,
                    GLOBAL_DATA_SINK_REWRITE,
                    InterpolationMode="Linear", Mode=10, DwiList=[], WORKFLOW_COMPONENTS=[], CLUSTER_QUEUE='', CLUSTER_QUEUE_LONG=''):
    """
    Run autoworkup on all subjects data defined in the ExperimentDatabase

    This is the main function to call when processing a data set with T1 & T2
    data.  ExperimentBaseDirectoryPrefix is the base of the directory to place results, T1Images & T2Images
    are the lists of images to be used in the auto-workup. atlas_fname_wpath is
    the path and filename of the atlas to use.
    """

    print "Building Pipeline"
    ########### PIPELINE INITIALIZATION #############
    baw200 = pe.Workflow(name="BAW_20120813")
    baw200.config['execution'] = {
        'plugin': 'Linear',
        #'stop_on_first_crash':'true',
        #'stop_on_first_rerun': 'true',
        'stop_on_first_crash': 'false',
        'stop_on_first_rerun': 'false',  # This stops at first attempt to rerun, before running, and before deleting previous results.
        'hash_method': 'timestamp',
        'single_thread_matlab': 'true',  # Multi-core 2011a  multi-core for matrix multiplication.
        'remove_unnecessary_outputs': 'false',
        'use_relative_paths': 'false',  # relative paths should be on, require hash update when changed.
        'remove_node_directories': 'false',  # Experimental
        'local_hash_check': 'true',
        'job_finished_timeout': 15
    }
    baw200.config['logging'] = {
        'workflow_level': 'DEBUG',
        'filemanip_level': 'DEBUG',
        'interface_level': 'DEBUG',
        'log_directory': ExperimentBaseDirectoryCache
    }
    baw200.base_dir = ExperimentBaseDirectoryCache

    import WorkupT1T2Single
    MergeT1s = dict()
    MergeT2s = dict()
    MergePDs = dict()
    MergeFLs = dict()
    MergeOutputLabels = dict()
    MergePosteriors = dict()
    BAtlas = dict()
    if True:
        print("===================== SUBJECT: {0} ===========================".format(subjectid))
        PHASE_1_oneSubjWorkflow = dict()
        PHASE_1_subjInfoNode = dict()
        allSessions = ExperimentDatabase.getSessionsFromSubject(subjectid)
        print("Running sessions: {ses} for subject {sub}".format(ses=allSessions, sub=subjectid))
        BAtlas[subjectid] = MakeAtlasNode(atlas_fname_wpath, "BAtlas_" + str(subjectid))  # Call function to create node

        for sessionid in allSessions:
            global_AllT1s = ExperimentDatabase.getFilenamesByScantype(sessionid, ['T1-30', 'T1-15'])
            global_AllT2s = ExperimentDatabase.getFilenamesByScantype(sessionid, ['T2-30', 'T2-15'])
            global_AllPDs = ExperimentDatabase.getFilenamesByScantype(sessionid, ['PD-30', 'PD-15'])
            global_AllFLs = ExperimentDatabase.getFilenamesByScantype(sessionid, ['FL-30', 'FL-15'])
            global_AllOthers = ExperimentDatabase.getFilenamesByScantype(sessionid, ['OTHER-30', 'OTHER-15'])
            print("HACK:  all T1s: {0} {1}".format(global_AllT1s, len(global_AllT1s)))
            print("HACK:  all T2s: {0} {1}".format(global_AllT2s, len(global_AllT2s)))
            print("HACK:  all PDs: {0} {1}".format(global_AllPDs, len(global_AllPDs)))
            print("HACK:  all FLs: {0} {1}".format(global_AllFLs, len(global_AllFLs)))
            print("HACK:  all Others: {0} {1}".format(global_AllOthers, len(global_AllOthers)))

            projectid = ExperimentDatabase.getProjFromSession(sessionid)
            print("PROJECT: {0} SUBJECT: {1} SESSION: {2}".format(projectid, subjectid, sessionid))
            PHASE_1_subjInfoNode[sessionid] = pe.Node(interface=IdentityInterface(fields=
                                                                                  ['sessionid', 'subjectid', 'projectid',
                                                                                   'allT1s',
                                                                                   'allT2s',
                                                                                   'allPDs',
                                                                                   'allFLs',
                                                                                   'allOthers']),
                                                      run_without_submitting=True,
                                                      name='99_PHASE_1_SubjInfoNode_' + str(subjectid) + "_" + str(sessionid))
            PHASE_1_subjInfoNode[sessionid].inputs.projectid = projectid
            PHASE_1_subjInfoNode[sessionid].inputs.subjectid = subjectid
            PHASE_1_subjInfoNode[sessionid].inputs.sessionid = sessionid
            PHASE_1_subjInfoNode[sessionid].inputs.allT1s = global_AllT1s
            PHASE_1_subjInfoNode[sessionid].inputs.allT2s = global_AllT2s
            PHASE_1_subjInfoNode[sessionid].inputs.allPDs = global_AllPDs
            PHASE_1_subjInfoNode[sessionid].inputs.allFLs = global_AllFLs
            PHASE_1_subjInfoNode[sessionid].inputs.allOthers = global_AllOthers

            PROCESSING_PHASE = 'PHASE_1'
            PHASE_1_WORKFLOW_COMPONENTS = ['BASIC', 'TISSUE_CLASSIFY']
            PHASE_1_oneSubjWorkflow[sessionid] = WorkupT1T2Single.MakeOneSubWorkFlow(
                projectid, subjectid, sessionid, PROCESSING_PHASE,
                PHASE_1_WORKFLOW_COMPONENTS,
                BCD_model_path, InterpolationMode, CLUSTER_QUEUE, CLUSTER_QUEUE_LONG)
            baw200.connect(PHASE_1_subjInfoNode[sessionid], 'projectid', PHASE_1_oneSubjWorkflow[sessionid], 'inputspec.projectid')
            baw200.connect(PHASE_1_subjInfoNode[sessionid], 'subjectid', PHASE_1_oneSubjWorkflow[sessionid], 'inputspec.subjectid')
            baw200.connect(PHASE_1_subjInfoNode[sessionid], 'sessionid', PHASE_1_oneSubjWorkflow[sessionid], 'inputspec.sessionid')
            baw200.connect(PHASE_1_subjInfoNode[sessionid], 'allT1s', PHASE_1_oneSubjWorkflow[sessionid], 'inputspec.allT1s')
            baw200.connect(PHASE_1_subjInfoNode[sessionid], 'allT2s', PHASE_1_oneSubjWorkflow[sessionid], 'inputspec.allT2s')
            baw200.connect(PHASE_1_subjInfoNode[sessionid], 'allPDs', PHASE_1_oneSubjWorkflow[sessionid], 'inputspec.allPDs')
            baw200.connect(PHASE_1_subjInfoNode[sessionid], 'allFLs', PHASE_1_oneSubjWorkflow[sessionid], 'inputspec.allFLs')
            baw200.connect(PHASE_1_subjInfoNode[sessionid], 'allOthers', PHASE_1_oneSubjWorkflow[sessionid], 'inputspec.allOthers')

            baw200.connect(BAtlas[subjectid], 'template_landmarks_31_fcsv', PHASE_1_oneSubjWorkflow[sessionid], 'inputspec.template_landmarks_31_fcsv')
            baw200.connect(BAtlas[subjectid], 'template_landmark_weights_31_csv', PHASE_1_oneSubjWorkflow[sessionid], 'inputspec.template_landmark_weights_31_csv')
            baw200.connect(BAtlas[subjectid], 'template_t1', PHASE_1_oneSubjWorkflow[sessionid], 'inputspec.template_t1')
            baw200.connect(BAtlas[subjectid], 'ExtendedAtlasDefinition_xml', PHASE_1_oneSubjWorkflow[sessionid], 'inputspec.atlasDefinition')

            BASIC_DataSink = dict()
            TC_DataSink = dict()
            AccumulateLikeTissuePosteriorsNode = dict()
            AddLikeTissueSink = dict()
            if True:
                ### Now define where the final organized outputs should go.
                BASIC_DataSink[sessionid] = pe.Node(nio.DataSink(), name="BASIC_DS_" + str(subjectid) + "_" + str(sessionid))
                BASIC_DataSink[sessionid].inputs.base_directory = ExperimentBaseDirectoryResults
                BASIC_DataSink[sessionid].inputs.regexp_substitutions = GenerateOutputPattern(projectid, subjectid, sessionid, 'ACPCAlign')

                baw200.connect(PHASE_1_oneSubjWorkflow[sessionid], 'outputspec.outputLandmarksInACPCAlignedSpace', BASIC_DataSink[sessionid], 'ACPCAlign.@outputLandmarksInACPCAlignedSpace')
                baw200.connect(PHASE_1_oneSubjWorkflow[sessionid], 'outputspec.BCD_ACPC_T1', BASIC_DataSink[sessionid], 'ACPCAlign.@BCD_ACPC_T1')
                baw200.connect(PHASE_1_oneSubjWorkflow[sessionid], 'outputspec.outputLandmarksInInputSpace', BASIC_DataSink[sessionid], 'ACPCAlign.@outputLandmarksInInputSpace')
                baw200.connect(PHASE_1_oneSubjWorkflow[sessionid], 'outputspec.outputTransform', BASIC_DataSink[sessionid], 'ACPCAlign.@outputTransform')
                baw200.connect(PHASE_1_oneSubjWorkflow[sessionid], 'outputspec.atlasToSubjectTransform', BASIC_DataSink[sessionid], 'ACPCAlign.@atlasToSubjectTransform')

                ### Now define where the final organized outputs should go.
                TC_DataSink[sessionid] = pe.Node(nio.DataSink(), name="TISSUE_CLASSIFY_DS_" + str(subjectid) + "_" + str(sessionid))
                TC_DataSink[sessionid].inputs.base_directory = ExperimentBaseDirectoryResults
                TC_DataSink[sessionid].inputs.regexp_substitutions = GenerateOutputPattern(projectid, subjectid, sessionid, 'TissueClassify')
                baw200.connect(PHASE_1_oneSubjWorkflow[sessionid], 'outputspec.TissueClassifyOutputDir', TC_DataSink[sessionid], 'TissueClassify.@TissueClassifyOutputDir')

                ### Now clean up by adding together many of the items PHASE_1_oneSubjWorkflow
                currentAccumulateLikeTissuePosteriorsName = 'AccumulateLikeTissuePosteriors_' + str(subjectid) + "_" + str(sessionid)
                AccumulateLikeTissuePosteriorsNode[sessionid] = pe.Node(interface=Function(function=AccumulateLikeTissuePosteriors,
                                                                                           input_names=['posteriorImages'],
                                                                                           output_names=['AccumulatePriorsList', 'AccumulatePriorsNames']),
                                                                        name=currentAccumulateLikeTissuePosteriorsName)
                baw200.connect(PHASE_1_oneSubjWorkflow[sessionid], 'outputspec.posteriorImages',
                               AccumulateLikeTissuePosteriorsNode[sessionid], 'posteriorImages')

                ### Now define where the final organized outputs should go.
                AddLikeTissueSink[sessionid] = pe.Node(nio.DataSink(), name="ACCUMULATED_POSTERIORS_" + str(subjectid) + "_" + str(sessionid))
                AddLikeTissueSink[sessionid].inputs.base_directory = ExperimentBaseDirectoryResults
                # AddLikeTissueSink[sessionid].inputs.regexp_substitutions = GenerateAccumulatorImagesOutputPattern(projectid, subjectid, sessionid)
                AddLikeTissueSink[sessionid].inputs.regexp_substitutions = GenerateOutputPattern(projectid, subjectid, sessionid, 'ACCUMULATED_POSTERIORS')
                baw200.connect(AccumulateLikeTissuePosteriorsNode[sessionid], 'AccumulatePriorsList', AddLikeTissueSink[sessionid], 'ACCUMULATED_POSTERIORS.@AccumulateLikeTissuePosteriorsOutputDir')

    return baw200
