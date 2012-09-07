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

def getListIndex( imageList, index):
    return imageList[index]

#HACK:  [('buildTemplateIteration2', 'SUBJECT_TEMPLATES/0249/buildTemplateIteration2')]
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
    patternList=[]

    find_pat=os.path.join('ANTSTemplate','Iteration02_Reshaped.nii.gz')
    replace_pat=os.path.join('SUBJECT_TEMPLATES',subjectid,r'AVG_T1.nii.gz')
    patternList.append( (find_pat,replace_pat) )

    find_pat=os.path.join('ANTSTemplate',r'_ReshapeAveragePassiveImageWithShapeUpdate[0-9]*/AVG_[A-Z0-9]*WARP_(?P<structure>AVG_[A-Z0-9]*.nii.gz)')
    replace_pat=os.path.join('SUBJECT_TEMPLATES',subjectid,r'\g<structure>')
    patternList.append( (find_pat,replace_pat) )

    find_pat=os.path.join('ANTSTemplate',r'CLIPPED_AVG_[A-Z]*WARP_(?P<structure>AVG_[A-Z]*.nii.gz)')
    replace_pat=os.path.join('SUBJECT_TEMPLATES',subjectid,r'\g<structure>')
    patternList.append( (find_pat,replace_pat) )

    print "HACK: ", patternList
    return patternList

def GenerateOutputPattern(projectid, subjectid, sessionid,DefaultNodeName):
    """ This function generates output path substitutions for workflows and nodes that conform to a common standard.
    """
    patternList=[]
    find_pat=os.path.join(DefaultNodeName)
    replace_pat=os.path.join(projectid,subjectid,sessionid,DefaultNodeName)
    patternList.append( (find_pat,replace_pat) )
    print "HACK: ", patternList
    return patternList

def GenerateAccumulatorImagesOutputPattern(projectid, subjectid, sessionid):
    """ This function generates output path substitutions for workflows and nodes that conform to a common standard.
    """
    patternList=[]
    find_pat="POSTERIOR_"
    replace_pat=os.path.join(projectid,subjectid,sessionid)+"/POSTERIOR_"
    patternList.append( (find_pat,replace_pat) )
    print "HACK: ", patternList
    return patternList

def MergeByExtendListElements(t2_averageList,outputLabels_averageList,ListOfPosteriorImagesDictionary):
    for t2_index in range(0,len(t2_averageList)):
        ListOfPosteriorImagesDictionary[t2_index]['T2']=t2_averageList[t2_index]
        ListOfPosteriorImagesDictionary[t2_index]['BRAINMASK']=outputLabels_averageList[t2_index]
    return ListOfPosteriorImagesDictionary

def MakeNewAtlasTemplate(t1_image,deformed_list,
            AtlasTemplate,outDefinition):
    import os
    import sys
    import SimpleITK as sitk
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
        'AVG_THALAMUSWARP_AVG_THALAMUS.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_THALAMUS.nii.gz',
        'AVG_HIPPOCAMPUSWARP_AVG_HIPPOCAMPUS.nii.gz': '@ATLAS_DIRECTORY@/EXTENDED_HIPPOCAMPUS.nii.gz',
        'AVG_T2WARP_AVG_T2.nii.gz':'@ATLAS_DIRECTORY@/template_t2.nii.gz',
        'AVG_BRAINMASKWARP_AVG_BRAINMASK.nii.gz':'@ATLAS_DIRECTORY@/template_brain.nii.gz',
        'T1_RESHAPED.nii.gz':'@ATLAS_DIRECTORY@/template_t1.nii.gz'
        }
    templateFile = open(AtlasTemplate,'r')
    content = templateFile.read()              # read entire file into memory
    templateFile.close()

    ## Now clean up the posteriors based on anatomical knowlege.
    ## sometimes the posteriors are not relevant for priors
    ## due to anomolies around the edges.
    load_images_list=dict()
    for full_pathname in deformed_list:
        base_name=os.path.basename(full_pathname)
        if base_name in patternDict.keys():
            load_images_list[base_name]=sitk.ReadImage(full_pathname)
    ## Make binary dilated mask
    binmask=sitk.BinaryThreshold(load_images_list['AVG_BRAINMASKWARP_AVG_BRAINMASK.nii.gz'],1,1000000)
    dilated5=sitk.DilateObjectMorphology(binmask,5)
    dilated5=sitk.Cast(dilated5,sitk.sitkFloat32) # Convert to Float32 for multiply
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
    clean_deformed_list=deformed_list
    for index in range(0,len(deformed_list)):
        full_pathname=deformed_list[index]
        base_name=os.path.basename(full_pathname)
        if base_name == 'AVG_BRAINMASKWARP_AVG_BRAINMASK.nii.gz':
            ### Make Brain Mask Binary
            clipped_name='CLIPPED_'+base_name
            patternDict[clipped_name]=patternDict[base_name]
            sitk.WriteImage(binmask,clipped_name)
            clean_deformed_list[index]=os.path.realpath(clipped_name)
        if base_name in interiorPriors:
            ### Make clipped posteriors for brain regions
            curr=sitk.Cast(sitk.ReadImage(full_pathname),sitk.sitkFloat32)
            curr=curr*dilated5
            clipped_name='CLIPPED_'+base_name
            patternDict[clipped_name]=patternDict[base_name]
            sitk.WriteImage(curr,clipped_name)
            clean_deformed_list[index]=os.path.realpath(clipped_name)
            print "HACK: ", clean_deformed_list[index]
            curr=None
    binmask=None
    dilated5=None

    for full_pathname in clean_deformed_list:
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
    return outAtlasFullPath,clean_deformed_list

def AccumulateLikeTissuePosteriors(posteriorImages):
    import os
    import sys
    import SimpleITK as sitk
    ## Now clean up the posteriors based on anatomical knowlege.
    ## sometimes the posteriors are not relevant for priors
    ## due to anomolies around the edges.
    load_images_list=dict()
    for full_pathname in posteriorImages.values():
        base_name=os.path.basename(full_pathname)
        load_images_list[base_name]=sitk.ReadImage(full_pathname)
    GM_ACCUM=[
              'POSTERIOR_ACCUMBEN.nii.gz',
              'POSTERIOR_CAUDATE.nii.gz',
              'POSTERIOR_CRBLGM.nii.gz',
              'POSTERIOR_HIPPOCAMPUS.nii.gz',
              'POSTERIOR_PUTAMEN.nii.gz',
              'POSTERIOR_THALAMUS.nii.gz',
              'POSTERIOR_SURFGM.nii.gz',
             ]
    WM_ACCUM=[
              'POSTERIOR_CRBLWM.nii.gz',
              'POSTERIOR_WM.nii.gz'
              ]
    CSF_ACCUM=[
              'POSTERIOR_CSF.nii.gz',
              ]
    VB_ACCUM=[
              'POSTERIOR_VB.nii.gz',
              ]
    GLOBUS_ACCUM=[
              'POSTERIOR_GLOBUS.nii.gz',
              ]
    BACKGROUND_ACCUM=[
              'POSTERIOR_AIR.nii.gz',
              'POSTERIOR_NOTCSF.nii.gz',
              'POSTERIOR_NOTGM.nii.gz',
              'POSTERIOR_NOTVB.nii.gz',
              'POSTERIOR_NOTWM.nii.gz',
              ]
    ## The next 2 items MUST be syncronized
    AccumulatePriorsNames=['POSTERIOR_GM_TOTAL.nii.gz','POSTERIOR_WM_TOTAL.nii.gz',
                        'POSTERIOR_CSF_TOTAL.nii.gz','POSTERIOR_VB_TOTAL.nii.gz',
                        'POSTERIOR_GLOBUS_TOTAL.nii.gz','POSTERIOR_BACKGROUND_TOTAL.nii.gz']
    ForcedOrderingLists=[GM_ACCUM,WM_ACCUM,CSF_ACCUM,VB_ACCUM,GLOBUS_ACCUM,BACKGROUND_ACCUM]
    AccumulatePriorsList=list()
    for index in range(0,len(ForcedOrderingLists)):
        outname=AccumulatePriorsNames[index]
        inlist=ForcedOrderingLists[index]
        accum_image= load_images_list[ inlist[0] ] # copy first image
        for curr_image in range(1,len(inlist)):
            accum_image=accum_image + load_images_list[ inlist[curr_image] ]
        sitk.WriteImage(accum_image,outname)
        AccumulatePriorsList.append(os.path.realpath(outname))
    print "HACK \n\n\n\n\n\n\n HACK \n\n\n: {APL}\n".format(APL=AccumulatePriorsList)
    print ": {APN}\n".format(APN=AccumulatePriorsNames)
    return AccumulatePriorsList,AccumulatePriorsNames
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
        print("Running sessions: {ses} for subject {sub}".format(ses=allSessions,sub=subjectid))
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
                              BCD_model_path, InterpolationMode, CLUSTER_QUEUE)
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
        if True or numSessions > 1: ## Merge all BCD_Results into a global average
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

            TEMPLATE_BUILD_RUN_MODE='MULTI_IMAGE'
            if numSessions == 1:
                TEMPLATE_BUILD_RUN_MODE='SINGLE_IMAGE'

            buildTemplateIteration1 = ANTSTemplateBuildSingleIterationWF('Iteration01',CLUSTER_QUEUE,TEMPLATE_BUILD_RUN_MODE)
            baw200.connect(myInitAvgWF, 'OutputSpec.average_image', buildTemplateIteration1, 'InputSpec.fixed_image')
            baw200.connect(MergeT1s[subjectid], 'out', buildTemplateIteration1, 'InputSpec.images')
            baw200.connect(MergeByExtendListElementsNode, 'ListOfExtendedPassiveImages', buildTemplateIteration1, 'InputSpec.ListOfPassiveImagesDictionararies')

            buildTemplateIteration2 = buildTemplateIteration1.clone(name='buildTemplateIteration2')
            buildTemplateIteration2 = ANTSTemplateBuildSingleIterationWF('Iteration02',CLUSTER_QUEUE,TEMPLATE_BUILD_RUN_MODE)
            baw200.connect(buildTemplateIteration1, 'OutputSpec.template', buildTemplateIteration2, 'InputSpec.fixed_image')
            baw200.connect(MergeT1s[subjectid], 'out', buildTemplateIteration2, 'InputSpec.images')
            baw200.connect(MergeByExtendListElementsNode, 'ListOfExtendedPassiveImages', buildTemplateIteration2, 'InputSpec.ListOfPassiveImagesDictionararies')

            #baw200.connect(InitAvgImages, 'average_image', outputSpec, 'average_image')

            ### Now define where the final organized outputs should go.
            SubjectTemplate_DataSink=pe.Node(nio.DataSink(),name="SubjectTemplate_DS")
            SubjectTemplate_DataSink.inputs.base_directory=ExperimentBaseDirectoryResults
            SubjectTemplate_DataSink.inputs.regexp_substitutions = GenerateSubjectOutputPattern(subjectid)
            baw200.connect(buildTemplateIteration2,'OutputSpec.template',SubjectTemplate_DataSink,'ANTSTemplate.@template')

            MakeNewAtlasTemplateNode = pe.Node(interface=Function(function=MakeNewAtlasTemplate,
                    input_names=['t1_image', 'deformed_list','AtlasTemplate','outDefinition'],
                    output_names=['outAtlasFullPath','clean_deformed_list']),
                    # This is a lot of work, so submit it run_without_submitting=True,
                    name='99_MakeNewAtlasTemplate')
            MakeNewAtlasTemplateNode.inputs.outDefinition='AtlasDefinition_'+subjectid+'.xml'
            baw200.connect(BAtlas[subjectid],'ExtendedAtlasDefinition_xml_in',MakeNewAtlasTemplateNode,'AtlasTemplate')
            baw200.connect(buildTemplateIteration2,'OutputSpec.template',MakeNewAtlasTemplateNode,'t1_image')
            baw200.connect(buildTemplateIteration2,'OutputSpec.passive_deformed_templates',MakeNewAtlasTemplateNode,'deformed_list')
            baw200.connect(MakeNewAtlasTemplateNode,'clean_deformed_list',SubjectTemplate_DataSink,'ANTSTemplate.@passive_deformed_templates')

            ###### Starting Phase II
            PHASE_2_oneSubjWorkflow=dict()
            PHASE_2_subjInfoNode=dict()
            BASIC_DataSink=dict()
            TC_DataSink=dict()
            AddLikeTissueSink=dict()
            AccumulateLikeTissuePosteriorsNode=dict()
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
                                  BCD_model_path, InterpolationMode, CLUSTER_QUEUE)
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

                ### Now define where the final organized outputs should go.
                BASIC_DataSink[sessionid]=pe.Node(nio.DataSink(),name="BASIC_DS_"+str(subjectid)+"_"+str(sessionid))
                BASIC_DataSink[sessionid].inputs.base_directory=ExperimentBaseDirectoryResults
                BASIC_DataSink[sessionid].inputs.regexp_substitutions = GenerateOutputPattern(projectid, subjectid, sessionid,'ACPCAlign')

                baw200.connect(PHASE_2_oneSubjWorkflow[sessionid],'OutputSpec.outputLandmarksInACPCAlignedSpace',BASIC_DataSink[sessionid],'ACPCAlign.@outputLandmarksInACPCAlignedSpace')
                baw200.connect(PHASE_2_oneSubjWorkflow[sessionid],'OutputSpec.BCD_ACPC_T1',BASIC_DataSink[sessionid],'ACPCAlign.@BCD_ACPC_T1')
                baw200.connect(PHASE_2_oneSubjWorkflow[sessionid],'OutputSpec.outputLandmarksInInputSpace',BASIC_DataSink[sessionid],'ACPCAlign.@outputLandmarksInInputSpace')
                baw200.connect(PHASE_2_oneSubjWorkflow[sessionid],'OutputSpec.outputTransform',BASIC_DataSink[sessionid],'ACPCAlign.@outputTransform')
                baw200.connect(PHASE_2_oneSubjWorkflow[sessionid],'OutputSpec.atlasToSubjectTransform',BASIC_DataSink[sessionid],'ACPCAlign.@atlasToSubjectTransform')

                ### Now define where the final organized outputs should go.
                TC_DataSink[sessionid]=pe.Node(nio.DataSink(),name="TISSUE_CLASSIFY_DS_"+str(subjectid)+"_"+str(sessionid))
                TC_DataSink[sessionid].inputs.base_directory=ExperimentBaseDirectoryResults
                TC_DataSink[sessionid].inputs.regexp_substitutions = GenerateOutputPattern(projectid, subjectid, sessionid,'TissueClassify')
                baw200.connect(PHASE_2_oneSubjWorkflow[sessionid], 'OutputSpec.TissueClassifyOutputDir', TC_DataSink[sessionid],'TissueClassify.@TissueClassifyOutputDir')

                ### Now clean up by adding together many of the items PHASE_2_oneSubjWorkflow
                currentAccumulateLikeTissuePosteriorsName='AccumulateLikeTissuePosteriors_'+str(subjectid)+"_"+str(sessionid)
                AccumulateLikeTissuePosteriorsNode[sessionid] = pe.Node(interface=Function(function=AccumulateLikeTissuePosteriors,
                     input_names=['posteriorImages'],
                     output_names=['AccumulatePriorsList','AccumulatePriorsNames']),
                     name=currentAccumulateLikeTissuePosteriorsName)
                baw200.connect(PHASE_2_oneSubjWorkflow[sessionid],'OutputSpec.posteriorImages',
                               AccumulateLikeTissuePosteriorsNode[sessionid],'posteriorImages')

                ### Now define where the final organized outputs should go.
                AddLikeTissueSink[sessionid]=pe.Node(nio.DataSink(),name="ACCUMULATED_POSTERIORS"+str(subjectid)+"_"+str(sessionid))
                AddLikeTissueSink[sessionid].inputs.base_directory=ExperimentBaseDirectoryResults
                #AddLikeTissueSink[sessionid].inputs.regexp_substitutions = GenerateAccumulatorImagesOutputPattern(projectid, subjectid, sessionid)
                AddLikeTissueSink[sessionid].inputs.regexp_substitutions = GenerateOutputPattern(projectid, subjectid, sessionid,'ACCUMULATED_POSTERIORS')
                baw200.connect(AccumulateLikeTissuePosteriorsNode[sessionid], 'AccumulatePriorsList', AddLikeTissueSink[sessionid],'ACCUMULATED_POSTERIORS.@AccumulateLikeTissuePosteriorsOutputDir')

                ClipT1ImageWithBrainMaskNode=dict()
                AtlasToSubjectantsRegistration=dict()
                myLocalSegWF=dict()
                SEGMENTATION_DataSink=dict()
                myLocalFSWF=dict()
                FSPREP_DataSink=dict()
                FS_DS=dict()

                if True: ## Run the ANTS Registration from Atlas to Subject for BCut spatial priors propagation.
                    import PipeLineFunctionHelpers
                    import BRAINSTools.ants.antsRegistration
                    ## Second clip to brain tissue region
                    ### Now clean up by adding together many of the items PHASE_2_oneSubjWorkflow
                    currentClipT1ImageWithBrainMaskName='ClipT1ImageWithBrainMask_'+str(subjectid)+"_"+str(sessionid)
                    ClipT1ImageWithBrainMaskNode[sessionid] = pe.Node(interface=Function(function=PipeLineFunctionHelpers.ClipT1ImageWithBrainMask,
                          input_names=['t1_image','brain_labels','clipped_file_name'],
                          output_names=['clipped_file']),
                          name=currentClipT1ImageWithBrainMaskName)
                    ClipT1ImageWithBrainMaskNode[sessionid].inputs.clipped_file_name = 'clipped_t1.nii.gz'
                    baw200.connect(PHASE_2_oneSubjWorkflow[sessionid],'OutputSpec.t1_average',ClipT1ImageWithBrainMaskNode[sessionid],'t1_image')
                    baw200.connect(BAtlas[subjectid],'template_t1_clipped',ClipT1ImageWithBrainMaskNode[sessionid],'brain_labels')


                    currentAtlasToSubjectantsRegistration='AtlasToSubjectantsRegistration_'+str(subjectid)+"_"+str(sessionid)
                    AtlasToSubjectantsRegistration[subjectid]=pe.Node(interface=BRAINSTools.ants.antsRegistration.antsRegistration(), name =
 currentAtlasToSubjectantsRegistration)
                    AtlasToSubjectantsRegistration[subjectid].inputs.dimension = 3
                    AtlasToSubjectantsRegistration[subjectid].inputs.output_transform_prefix = 'AtlasToSubject_'
                    AtlasToSubjectantsRegistration[subjectid].inputs.metric = ['Mattes','Mattes']
                    AtlasToSubjectantsRegistration[subjectid].inputs.metric_weight = [1,1]
                    AtlasToSubjectantsRegistration[subjectid].inputs.radius_or_number_of_bins = [32,32] ## This is really number of bins
                    AtlasToSubjectantsRegistration[subjectid].inputs.transforms =           [ "Affine",  "SyN"              ]
                    AtlasToSubjectantsRegistration[subjectid].inputs.transform_parameters = [ [ 2 ],     [ 0.25, 3.0, 0.0 ] ]
                    AtlasToSubjectantsRegistration[subjectid].inputs.number_of_iterations = [[1500, 200], [100, 50, 30]     ]
                    AtlasToSubjectantsRegistration[subjectid].inputs.use_histogram_matching = [ True, True ]
                    AtlasToSubjectantsRegistration[subjectid].inputs.shrink_factors = [[2,1],[3,2,1]]
                    AtlasToSubjectantsRegistration[subjectid].inputs.smoothing_sigmas = [[1,0],[2,1,0]]
                    AtlasToSubjectantsRegistration[subjectid].inputs.use_estimate_learning_rate_once = [True,True]
                    baw200.connect(BAtlas[subjectid],'template_t1_clipped',AtlasToSubjectantsRegistration[subjectid], 'moving_image')
                    baw200.connect(ClipT1ImageWithBrainMaskNode[sessionid], 'clipped_file', AtlasToSubjectantsRegistration[subjectid], 'fixed_image')
                    baw200.connect(PHASE_2_oneSubjWorkflow[sessionid],'OutputSpec.atlasToSubjectTransform',AtlasToSubjectantsRegistration[subjectid],'initial_moving_transform')

                if 'SEGMENTATION' in WORKFLOW_COMPONENTS:
                    from WorkupT1T2BRAINSCut import CreateBRAINSCutWorkflow
                    myLocalSegWF[subjectid] = CreateBRAINSCutWorkflow(projectid, subjectid, sessionid,'Segmentation',CLUSTER_QUEUE,BAtlas[subjectid]) ##Note:  Passing in the entire BAtlas Object here!
                    baw200.connect( PHASE_2_oneSubjWorkflow[sessionid], 'OutputSpec.t1_average', myLocalSegWF[subjectid], "InputSpec.T1Volume" )
                    baw200.connect( PHASE_2_oneSubjWorkflow[sessionid], 'OutputSpec.t2_average', myLocalSegWF[subjectid], "InputSpec.T2Volume")
                    baw200.connect( PHASE_2_oneSubjWorkflow[sessionid], 'OutputSpec.outputLabels', myLocalSegWF[subjectid],"InputSpec.RegistrationROI")
                    ## NOTE: Element 0 of AccumulatePriorsList is the accumulated GM tissue
                    baw200.connect( [ ( AccumulateLikeTissuePosteriorsNode[sessionid], myLocalSegWF[subjectid],
                                      [ (( 'AccumulatePriorsList', getListIndex, 0 ), "InputSpec.TotalGM")] ),
] )

                    baw200.connect( PHASE_2_oneSubjWorkflow[sessionid],'OutputSpec.atlasToSubjectTransform',myLocalSegWF[subjectid],'InputSpec.atlasToSubjectTransform')

                    ### Now define where the final organized outputs should go.
                    SEGMENTATION_DataSink[subjectid]=pe.Node(nio.DataSink(),name="SEGMENTATION_DS_"+str(subjectid)+"_"+str(sessionid))
                    SEGMENTATION_DataSink[subjectid].inputs.base_directory=ExperimentBaseDirectoryResults
                    SEGMENTATION_DataSink[subjectid].inputs.regexp_substitutions = GenerateOutputPattern(projectid, subjectid, sessionid,'BRAINSCut')
                    baw200.connect(myLocalSegWF[subjectid], 'OutputSpec.outputBinaryLeftCaudate',SEGMENTATION_DataSink[subjectid], 'BRAINSCut.@outputBinaryLeftCaudate')
                    baw200.connect(myLocalSegWF[subjectid], 'OutputSpec.outputBinaryRightCaudate',SEGMENTATION_DataSink[subjectid], 'BRAINSCut.@outputBinaryRightCaudate')
                    baw200.connect(myLocalSegWF[subjectid], 'OutputSpec.outputBinaryLeftHippocampus',SEGMENTATION_DataSink[subjectid], 'BRAINSCut.@outputBinaryLeftHippocampus')
                    baw200.connect(myLocalSegWF[subjectid], 'OutputSpec.outputBinaryRightHippocampus',SEGMENTATION_DataSink[subjectid], 'BRAINSCut.@outputBinaryRightHippocampus')
                    baw200.connect(myLocalSegWF[subjectid], 'OutputSpec.outputBinaryLeftPutamen',SEGMENTATION_DataSink[subjectid], 'BRAINSCut.@outputBinaryLeftPutamen')
                    baw200.connect(myLocalSegWF[subjectid], 'OutputSpec.outputBinaryRightPutamen',SEGMENTATION_DataSink[subjectid], 'BRAINSCut.@outputBinaryRightPutamen')
                    baw200.connect(myLocalSegWF[subjectid], 'OutputSpec.outputBinaryLeftThalamus',SEGMENTATION_DataSink[subjectid], 'BRAINSCut.@outputBinaryLeftThalamus')
                    baw200.connect(myLocalSegWF[subjectid], 'OutputSpec.outputBinaryRightThalamus',SEGMENTATION_DataSink[subjectid], 'BRAINSCut.@outputBinaryRightThalamus')
                    baw200.connect(myLocalSegWF[subjectid], 'OutputSpec.outputLabelImageName', SEGMENTATION_DataSink[subjectid],'BRAINSCut.@outputLabelImageName')
                    baw200.connect(myLocalSegWF[subjectid], 'OutputSpec.outputCSVFileName', SEGMENTATION_DataSink[subjectid],'BRAINSCut.@outputCSVFileName')

                if 'FREESURFER' in WORKFLOW_COMPONENTS:
                    RunAllFSComponents=True ## A hack to avoid 26 hour run of freesurfer
                    from WorkupT1T2FreeSurfer import CreateFreeSurferWorkflow
                    myLocalFSWF[subjectid]= CreateFreeSurferWorkflow(projectid, subjectid, sessionid,"Level1_FSTest",CLUSTER_QUEUE,RunAllFSComponents)
                    baw200.connect(uidSource,'uid',myLocalFSWF[subjectid],'InputSpec.subject_id')
                    baw200.connect(PHASE_2_oneSubjWorkflow[sessionid],'OutputSpec.t1_average',myLocalFSWF[subjectid],'InputSpec.T1_files')
                    baw200.connect(PHASE_2_oneSubjWorkflow[sessionid],'OutputSpec.t2_average',myLocalFSWF[subjectid],'InputSpec.T2_files')
                    baw200.connect(PHASE_2_oneSubjWorkflow[sessionid],'OutputSpec.outputLabels',myLocalFSWF[subjectid],'InputSpec.label_file')
                    #baw200.connect(PHASE_2_oneSubjWorkflow[sessionid],'OutputSpec.outputLabels',myLocalFSWF[subjectid],'InputSpec.mask_file') #Yes, the same file as label_file!

                    ### Now define where the final organized outputs should go.
                    if RunAllFSComponents == True:
                        FS_DS[subjectid]=pe.Node(nio.DataSink(),name="FREESURFER_DS_"+str(subjectid)+"_"+str(sessionid))
                        FS_DS[subjectid].inputs.base_directory=ExperimentBaseDirectoryResults
                        FS_DS[subjectid].inputs.regexp_substitutions = [
                            ('/_uid_(?P<myuid>[^/]*)',r'/\g<myuid>')
                            ]
                        baw200.connect(myLocalFSWF[subjectid], 'OutputSpec.FreesurferOutputDirectory', FS_DS[subjectid],'FREESURFER_SUBJ.@FreesurferOutputDirectory')
                    ### Now define where the final organized outputs should go.
                    FSPREP_DataSink[subjectid]=pe.Node(nio.DataSink(),name="FREESURFER_PREP_"+str(subjectid)+"_"+str(sessionid))
                    FSPREP_DataSink[subjectid].inputs.base_directory=ExperimentBaseDirectoryResults
                    FREESURFER_PREP_PATTERNS = GenerateOutputPattern(projectid, subjectid, sessionid,'FREESURFER_PREP')
                    FSPREP_DataSink[subjectid].inputs.regexp_substitutions = FREESURFER_PREP_PATTERNS
                    print "========================="
                    print "========================="
                    print "========================="
                    print FREESURFER_PREP_PATTERNS
                    print "========================="
                    print "========================="
                    print "========================="
                    baw200.connect(myLocalFSWF[subjectid], 'OutputSpec.cnr_optimal_image', FSPREP_DataSink[subjectid],'FREESURFER_PREP.@cnr_optimal_image')

                else:
                    print "Skipping freesurfer"

    return baw200
