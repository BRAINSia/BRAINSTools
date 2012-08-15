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
from BRAINSTools.ants import antsAverageImages
#BRAINSTools/ants/buildtemplateparallel.py:import antsAverageImages

from WorkupT1T2AtlasNode import MakeAtlasNode

#HACK:  [('buildTemplateIteration2', 'SUBJECT_TEMPLATES/0249/buildTemplateIteration2')]
def GenerateSubjectOutputPattern(subjectid):
    """ This function generates output path substitutions for workflows and nodes that conform to a common standard.
    """
    patternList=[]
    find_pat=os.path.join('ANTSTemplate','Iteration02_Reshaped.nii.gz')
    replace_pat=os.path.join('SUBJECT_TEMPLATES',subjectid,'T1_RESHAPED.nii.gz')
    patternList.append( (find_pat,replace_pat) )
    find_pat=os.path.join('ANTSTemplate','_ReshapeAveragePassiveImageWithShapeUpdate[0-9]*')
    replace_pat=os.path.join('SUBJECT_TEMPLATES',subjectid)
    patternList.append( (find_pat,replace_pat) )
    print "HACK: ", patternList
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
def WorkupT1T2(subjectid,mountPrefix,ExperimentBaseDirectoryCache, ExperimentBaseDirectoryResults, ExperimentDatabase, atlas_fname_wpath, BCD_model_path,
               InterpolationMode="Linear", Mode=10,DwiList=[],WORKFLOW_COMPONENTS=[],CLUSTER_QUEUE=''):
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
                                     'plugin':'Linear',
                                     #'stop_on_first_crash':'true',
                                     #'stop_on_first_rerun': 'true',
                                     'stop_on_first_crash':'false',
                                     'stop_on_first_rerun': 'false',      ## This stops at first attempt to rerun, before running, and before deleting previous results.
                                     'hash_method': 'timestamp',
                                     'single_thread_matlab':'true',       ## Multi-core 2011a  multi-core for matrix multiplication.
                                     'remove_unnecessary_outputs':'false',
                                     'use_relative_paths':'false',         ## relative paths should be on, require hash update when changed.
                                     'remove_node_directories':'false',   ## Experimental
                                     'local_hash_check':'true',           ##
                                     'job_finished_timeout':15            ##
                                     }
    baw200.config['logging'] = {
      'workflow_level':'DEBUG',
      'filemanip_level':'DEBUG',
      'interface_level':'DEBUG',
      'log_directory': ExperimentBaseDirectoryCache
    }
    baw200.base_dir = ExperimentBaseDirectoryCache


    import WorkupT1T2Single
    MergeT1s=dict()
    MergeT2s=dict()
    MergeOutputLabels=dict()
    MergePosteriors=dict()
    BAtlas=dict()
    if True:
        print("===================== SUBJECT: {0} ===========================".format(subjectid))
        PHASE_1_oneSubjWorkflow=dict()
        PHASE_1_subjInfoNode=dict()
        allSessions = ExperimentDatabase.getSessionsFromSubject(subjectid)
        BAtlas[subjectid] = MakeAtlasNode(atlas_fname_wpath,"BAtlas_"+str(subjectid)) ## Call function to create node
        for sessionid in allSessions:
            projectid = ExperimentDatabase.getProjFromSession(sessionid)
            print("PROJECT: {0} SUBJECT: {1} SESSION: {2}".format(projectid,subjectid,sessionid))
            PHASE_1_subjInfoNode[sessionid] = pe.Node(interface=IdentityInterface(fields=
                    ['sessionid','subjectid','projectid',
                     'allT1s',
                     'allT2s',
                     'allPDs',
                     'allOthers']),
                    run_without_submitting=True,
                    name='99_PHASE_1_SubjInfoNode_'+str(subjectid)+"_"+str(sessionid) )
            PHASE_1_subjInfoNode[sessionid].inputs.projectid=projectid
            PHASE_1_subjInfoNode[sessionid].inputs.subjectid=subjectid
            PHASE_1_subjInfoNode[sessionid].inputs.sessionid=sessionid
            PHASE_1_subjInfoNode[sessionid].inputs.allT1s=ExperimentDatabase.getFilenamesByScantype(sessionid,['T1-30','T1-15'])
            PHASE_1_subjInfoNode[sessionid].inputs.allT2s=ExperimentDatabase.getFilenamesByScantype(sessionid,['T2-30','T2-15'])
            PHASE_1_subjInfoNode[sessionid].inputs.allPDs=ExperimentDatabase.getFilenamesByScantype(sessionid,['PD-30','PD-15'])
            PHASE_1_subjInfoNode[sessionid].inputs.allOthers=ExperimentDatabase.getFilenamesByScantype(sessionid,['OTHER-30','OTHER-15'])

            PROCESSING_PHASE='PHASE_1'
            PHASE_1_WORKFLOW_COMPONENTS =  ['BASIC','TISSUE_CLASSIFY']
            PHASE_1_oneSubjWorkflow[sessionid]=WorkupT1T2Single.MakeOneSubWorkFlow(
                              projectid, subjectid, sessionid,PROCESSING_PHASE,
                              PHASE_1_WORKFLOW_COMPONENTS,
                              BCD_model_path, InterpolationMode, CLUSTER_QUEUE,
                              ExperimentBaseDirectoryResults)
            baw200.connect(PHASE_1_subjInfoNode[sessionid],'projectid',PHASE_1_oneSubjWorkflow[sessionid],'InputSpec.projectid')
            baw200.connect(PHASE_1_subjInfoNode[sessionid],'subjectid',PHASE_1_oneSubjWorkflow[sessionid],'InputSpec.subjectid')
            baw200.connect(PHASE_1_subjInfoNode[sessionid],'sessionid',PHASE_1_oneSubjWorkflow[sessionid],'InputSpec.sessionid')
            baw200.connect(PHASE_1_subjInfoNode[sessionid],'allT1s',PHASE_1_oneSubjWorkflow[sessionid],'InputSpec.allT1s')
            baw200.connect(PHASE_1_subjInfoNode[sessionid],'allT2s',PHASE_1_oneSubjWorkflow[sessionid],'InputSpec.allT2s')
            baw200.connect(PHASE_1_subjInfoNode[sessionid],'allPDs',PHASE_1_oneSubjWorkflow[sessionid],'InputSpec.allPDs')
            baw200.connect(PHASE_1_subjInfoNode[sessionid],'allOthers',PHASE_1_oneSubjWorkflow[sessionid],'InputSpec.allOthers')

            baw200.connect(BAtlas[subjectid],'template_landmarks_31_fcsv', PHASE_1_oneSubjWorkflow[sessionid],'InputSpec.template_landmarks_31_fcsv')
            baw200.connect(BAtlas[subjectid],'template_landmark_weights_31_csv', PHASE_1_oneSubjWorkflow[sessionid],'InputSpec.template_landmark_weights_31_csv')
            baw200.connect(BAtlas[subjectid],'template_t1', PHASE_1_oneSubjWorkflow[sessionid],'InputSpec.template_t1')
            baw200.connect(BAtlas[subjectid],'ExtendedAtlasDefinition_xml', PHASE_1_oneSubjWorkflow[sessionid],'InputSpec.atlasDefinition')

        numSessions=len(allSessions)
        if numSessions > 1: ## Merge all BCD_Results into a global average
            mergeSubjectSessionNamesT1="99_MergeAllSessions_T1"+str(subjectid)
            MergeT1s[subjectid] = pe.Node(interface=Merge(numSessions),
                                          run_without_submitting=True,
                                          name=mergeSubjectSessionNamesT1)
            mergeSubjectSessionNamesT2="99_MergeAllSessions_T2"+str(subjectid)
            MergeT2s[subjectid] = pe.Node(interface=Merge(numSessions),
                                          run_without_submitting=True,
                                          name=mergeSubjectSessionNamesT2)
            mergeSubjectSessionNamesoutputLabels="99_MergeAllSessions_outputLabels"+str(subjectid)
            MergeOutputLabels[subjectid] = pe.Node(interface=Merge(numSessions),
                                          run_without_submitting=True,
                                          name=mergeSubjectSessionNamesoutputLabels)
            mergeSubjectSessionNamesPosteriors="99_MergeAllSessions_Posteriors"+str(subjectid)
            MergePosteriors[subjectid] = pe.Node(interface=Merge(numSessions),
                                          run_without_submitting=True,
                                          name=mergeSubjectSessionNamesPosteriors)
            index=1
            #print("HACK: HACK: HACK:  {0}".format(allSessions))
            for sessionid in allSessions:
                index_name='in'+str(index)
                index+=1
                baw200.connect(PHASE_1_oneSubjWorkflow[sessionid],'OutputSpec.t1_average',MergeT1s[subjectid],index_name)
                baw200.connect(PHASE_1_oneSubjWorkflow[sessionid],'OutputSpec.t2_average',MergeT2s[subjectid],index_name)
                baw200.connect(PHASE_1_oneSubjWorkflow[sessionid],'OutputSpec.outputLabels',MergeOutputLabels[subjectid],index_name)
                baw200.connect(PHASE_1_oneSubjWorkflow[sessionid],'OutputSpec.posteriorImages',MergePosteriors[subjectid],index_name)

            def MergeByExtendListElements(t2_averageList,outputLabels_averageList,ListOfPosteriorImagesDictionary):
                for t2_index in range(0,len(t2_averageList)):
                    ListOfPosteriorImagesDictionary[t2_index]['T2']=t2_averageList[t2_index]
                    ListOfPosteriorImagesDictionary[t2_index]['BRAINMASK']=outputLabels_averageList[t2_index]
                return ListOfPosteriorImagesDictionary
            MergeByExtendListElementsNode = pe.Node( Function(function=MergeByExtendListElements,
                                          input_names = ['t2_averageList','outputLabels_averageList','ListOfPosteriorImagesDictionary'],
                                          output_names = ['ListOfExtendedPassiveImages']),
                                          run_without_submitting=True, name="99_MergeByExtendListElements")
            #MergeByExtendListElementsNode.inputs.preserve_nested_lists = True
            baw200.connect( MergeT2s[subjectid],'out', MergeByExtendListElementsNode, 't2_averageList' )
            baw200.connect( MergeOutputLabels[subjectid],'out', MergeByExtendListElementsNode, 'outputLabels_averageList' )
            baw200.connect( MergePosteriors[subjectid],'out', MergeByExtendListElementsNode, 'ListOfPosteriorImagesDictionary' )

            """ Now part of ants directly
            InitAvgImages=pe.Node(interface=antsAverageImages.AntsAverageImages(), name ='InitBCDAvgImages_'+str(subjectid))
            InitAvgImages.inputs.dimension = 3
            InitAvgImages.inputs.output_average_image = str(subjectid)+"_TissueClassAVG_T1.nii.gz"
            InitAvgImages.inputs.normalize = 1
            baw200.connect(MergeT1s[subjectid], 'out', InitAvgImages, 'images')
            """
            from BRAINSTools.ants.antsSimpleAverageWF import antsSimpleAverageWF
            ### USE ANTS
            from BRAINSTools.ants.buildtemplateparallel import ANTSTemplateBuildSingleIterationWF
            ### USE ANTS REGISTRATION
            #from BRAINSTools.ants.buildtemplateparallel_antsRegistration import antsRegistrationTemplateBuildSingleIterationWF

            myInitAvgWF = antsSimpleAverageWF()
            baw200.connect(MergeT1s[subjectid], 'out', myInitAvgWF, 'InputSpec.images')

            buildTemplateIteration1 = ANTSTemplateBuildSingleIterationWF('Iteration01',CLUSTER_QUEUE)
            baw200.connect(myInitAvgWF, 'OutputSpec.average_image', buildTemplateIteration1, 'InputSpec.fixed_image')
            baw200.connect(MergeT1s[subjectid], 'out', buildTemplateIteration1, 'InputSpec.images')
            baw200.connect(MergeByExtendListElementsNode, 'ListOfExtendedPassiveImages', buildTemplateIteration1, 'InputSpec.ListOfPassiveImagesDictionararies')

            buildTemplateIteration2 = buildTemplateIteration1.clone(name='buildTemplateIteration2')
            buildTemplateIteration2 = ANTSTemplateBuildSingleIterationWF('Iteration02',CLUSTER_QUEUE)
            baw200.connect(buildTemplateIteration1, 'OutputSpec.template', buildTemplateIteration2, 'InputSpec.fixed_image')
            baw200.connect(MergeT1s[subjectid], 'out', buildTemplateIteration2, 'InputSpec.images')
            baw200.connect(MergeByExtendListElementsNode, 'ListOfExtendedPassiveImages', buildTemplateIteration2, 'InputSpec.ListOfPassiveImagesDictionararies')

            #baw200.connect(InitAvgImages, 'average_image', outputSpec, 'average_image')

            ### Now define where the final organized outputs should go.
            SubjectTemplate_DataSink=pe.Node(nio.DataSink(),name="SubjectTemplate_DS")
            SubjectTemplate_DataSink.inputs.base_directory=ExperimentBaseDirectoryResults
            SubjectTemplate_DataSink.inputs.regexp_substitutions = GenerateSubjectOutputPattern(subjectid)
            baw200.connect(buildTemplateIteration2,'OutputSpec.template',SubjectTemplate_DataSink,'ANTSTemplate.@template')
            baw200.connect(buildTemplateIteration2,'OutputSpec.passive_deformed_templates',SubjectTemplate_DataSink,'ANTSTemplate.@passive_deformed_templates')

            def MakeNewAtlasTemplate(t1_image,deformed_list,
                        AtlasTemplate,outDefinition):
                import os
                import sys
                patternDict= {
                    'AVG_AIRWARP_AVG_AIR.nii.gz':'@ATLAS_DIRECTORY@/EXTENDED_AIR.nii.gz',
                    'AVG_BGMWARP_AVG_BGM.nii.gz':'@ATLAS_DIRECTORY@/EXTENDED_BASALTISSUE.nii.gz',
                    'AVG_CRBLGMWARP_AVG_CRBLGM.nii.gz':'@ATLAS_DIRECTORY@/EXTENDED_CRBLGM.nii.gz',
                    'AVG_CRBLWMWARP_AVG_CRBLWM.nii.gz':'@ATLAS_DIRECTORY@/EXTENDED_CRBLWM.nii.gz',
                    'AVG_CSFWARP_AVG_CSF.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_CSF.nii.gz',
                    'AVG_NOTCSFWARP_AVG_NOTCSF.nii.gz' :'@ATLAS_DIRECTORY@/EXTENDED_NOTCSF.nii.gz',
                    'AVG_NOTGMWARP_AVG_NOTGM.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_NOTGM.nii.gz',
                    'AVG_NOTVBWARP_AVG_NOTVB.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_NOTVB.nii.gz',
                    'AVG_NOTWMWARP_AVG_NOTWM.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_NOTWM.nii.gz',
                    'AVG_SURFGMWARP_AVG_SURFGM.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_SURFGM.nii.gz',
                    'AVG_VBWARP_AVG_VB.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_VB.nii.gz',
                    'AVG_WMWARP_AVG_WM.nii.gz':'@ATLAS_DIRECTORY@/EXTENDED_WM.nii.gz',
                    'AVG_ACCUMBENWARP_AVG_ACCUMBEN.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_ACCUMBEN.nii.gz',
                    'AVG_CAUDATEWARP_AVG_CAUDATE.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_CAUDATE.nii.gz',
                    'AVG_PUTAMENWARP_AVG_PUTAMEN.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_PUTAMEN.nii.gz',
                    'AVG_GLOBUSWARP_AVG_GLOBUS.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_GLOBUS.nii.gz',
                    'AVG_THALAMUS_AVG_THALAMUS.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_THALAMUS.nii.gz',
                    'AVG_HIPPOCAMPUS_AVG_HIPPOCAMPUS.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_HIPPOCAMPUS.nii.gz',
                    'AVG_T2WARP_AVG_T2.nii.gz':'@ATLAS_DIRECTORY@/template_t2.nii.gz',
                    'AVG_BRAINMASKWARP_AVG_BRAINMASK.nii.gz':'@ATLAS_DIRECTORY@/template_brain.nii.gz',
                    'T1_RESHAPED.nii.gz':'@ATLAS_DIRECTORY@/template_t1.nii.gz'
                    }
                templateFile = open(AtlasTemplate,'r')
                content = templateFile.read()              # read entire file into memory
                templateFile.close()
                for full_pathname in deformed_list:
                        base_name=os.path.basename(full_pathname)
                        if base_name in patternDict.keys():
                                content=content.replace(patternDict[base_name],full_pathname)
                content=content.replace('@ATLAS_DIRECTORY@/template_t1.nii.gz',t1_image)
                ## NOTE:  HEAD REGION CAN JUST BE T1 image.
                content=content.replace('@ATLAS_DIRECTORY@/template_headregion.nii.gz',t1_image)
                ## NOTE:  BRAIN REGION CAN JUST BE the label images.
                outAtlasFullPath=os.path.realpath(outDefinition)
                newFile = open(outAtlasFullPath, 'w')
                newFile.write(content)  # write the file with the text substitution
                newFile.close()
                return outAtlasFullPath
            MakeNewAtlasTemplateNode = pe.Node(interface=Function(function=MakeNewAtlasTemplate,
                    input_names=['t1_image', 'deformed_list','AtlasTemplate','outDefinition'],
                    output_names=['outAtlasFullPath']),
                    run_without_submitting=True,
                    name='99_MakeNewAtlasTemplate')
            MakeNewAtlasTemplateNode.inputs.outDefinition='AtlasDefinition_'+subjectid+'.xml'
            baw200.connect(BAtlas[subjectid],'ExtendedAtlasDefinition_xml_in',MakeNewAtlasTemplateNode,'AtlasTemplate')
            baw200.connect(buildTemplateIteration2,'OutputSpec.template',MakeNewAtlasTemplateNode,'t1_image')
            baw200.connect(buildTemplateIteration2,'OutputSpec.passive_deformed_templates',MakeNewAtlasTemplateNode,'deformed_list')

            ###### Starting Phase II
            PHASE_2_oneSubjWorkflow=dict()
            PHASE_2_subjInfoNode=dict()
            for sessionid in allSessions:
               projectid = ExperimentDatabase.getProjFromSession(sessionid)
               print("PHASE II PROJECT: {0} SUBJECT: {1} SESSION: {2}".format(projectid,subjectid,sessionid))
               PHASE_2_subjInfoNode[sessionid] = pe.Node(interface=IdentityInterface(fields=
                       ['sessionid','subjectid','projectid',
                        'allT1s',
                        'allT2s',
                        'allPDs',
                        'allOthers']),
                       run_without_submitting=True,
                       name='99_PHASE_2_SubjInfoNode_'+str(subjectid)+"_"+str(sessionid) )
               PHASE_2_subjInfoNode[sessionid].inputs.projectid=projectid
               PHASE_2_subjInfoNode[sessionid].inputs.subjectid=subjectid
               PHASE_2_subjInfoNode[sessionid].inputs.sessionid=sessionid
               PHASE_2_subjInfoNode[sessionid].inputs.allT1s=ExperimentDatabase.getFilenamesByScantype(sessionid,['T1-30','T1-15'])
               PHASE_2_subjInfoNode[sessionid].inputs.allT2s=ExperimentDatabase.getFilenamesByScantype(sessionid,['T2-30','T2-15'])
               PHASE_2_subjInfoNode[sessionid].inputs.allPDs=ExperimentDatabase.getFilenamesByScantype(sessionid,['PD-30','PD-15'])
               PHASE_2_subjInfoNode[sessionid].inputs.allOthers=ExperimentDatabase.getFilenamesByScantype(sessionid,['OTHER-30','OTHER-15'])

               PROCESSING_PHASE='PHASE_2'
               PHASE_2_oneSubjWorkflow[sessionid]=WorkupT1T2Single.MakeOneSubWorkFlow(
                                 projectid, subjectid, sessionid,PROCESSING_PHASE,
                                 WORKFLOW_COMPONENTS,
                                 BCD_model_path, InterpolationMode, CLUSTER_QUEUE,
                                 ExperimentBaseDirectoryResults)
               baw200.connect(PHASE_2_subjInfoNode[sessionid],'projectid',PHASE_2_oneSubjWorkflow[sessionid],'InputSpec.projectid')
               baw200.connect(PHASE_2_subjInfoNode[sessionid],'subjectid',PHASE_2_oneSubjWorkflow[sessionid],'InputSpec.subjectid')
               baw200.connect(PHASE_2_subjInfoNode[sessionid],'sessionid',PHASE_2_oneSubjWorkflow[sessionid],'InputSpec.sessionid')
               baw200.connect(PHASE_2_subjInfoNode[sessionid],'allT1s',PHASE_2_oneSubjWorkflow[sessionid],'InputSpec.allT1s')
               baw200.connect(PHASE_2_subjInfoNode[sessionid],'allT2s',PHASE_2_oneSubjWorkflow[sessionid],'InputSpec.allT2s')
               baw200.connect(PHASE_2_subjInfoNode[sessionid],'allPDs',PHASE_2_oneSubjWorkflow[sessionid],'InputSpec.allPDs')
               baw200.connect(PHASE_2_subjInfoNode[sessionid],'allOthers',PHASE_2_oneSubjWorkflow[sessionid],'InputSpec.allOthers')

               baw200.connect(BAtlas[subjectid],'template_landmarks_31_fcsv', PHASE_2_oneSubjWorkflow[sessionid],'InputSpec.template_landmarks_31_fcsv')
               baw200.connect(BAtlas[subjectid],'template_landmark_weights_31_csv', PHASE_2_oneSubjWorkflow[sessionid],'InputSpec.template_landmark_weights_31_csv')
               baw200.connect(buildTemplateIteration2,'OutputSpec.template', PHASE_2_oneSubjWorkflow[sessionid],'InputSpec.template_t1')
               baw200.connect(MakeNewAtlasTemplateNode,'outAtlasFullPath', PHASE_2_oneSubjWorkflow[sessionid],'InputSpec.atlasDefinition')
    return baw200
