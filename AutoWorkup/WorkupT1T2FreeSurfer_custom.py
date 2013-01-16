#!/usr/bin/env python

from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine
# from nipype.interfaces.freesurfer import ReconAll
import fswrap

from nipype.interfaces.freesurfer.model import MS_LDA
import os

"""
    from WorkupT1T2FreeSurfer import CreateFreeSurferWorkflow
    myLocalFSWF= CreateFreeSurferWorkflow("HansFSTest")
    baw200.connect(uidSource,'uid',myLocalFSWF,'inputspec.subject_id')
    baw200.connect(SplitAvgBABC,'avgBABCT1',myLocalFSWF,'inputspec.T1_files')
"""


def MakeFreeSurferOutputDirectory(subjects_dir, subject_id):
    return subjects_dir + '/' + subject_id


def GenerateWFName(projectid, subjectid, sessionid, WFName):
    return WFName + '_' + str(subjectid) + "_" + str(sessionid) + "_" + str(projectid)


def CreateFreeSurferWorkflow_custom(projectid, subjectid, sessionid, WFname, CLUSTER_QUEUE, CLUSTER_QUEUE_LONG, RunAllFSComponents=True, RunMultiMode=True, constructed_FS_SUBJECTS_DIR='/never_use_this'):
    freesurferWF = pe.Workflow(name=GenerateWFName(projectid, subjectid, sessionid, WFname))

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['FreeSurfer_ID', 'T1_files', 'T2_files', 'subjects_dir',
                                                             'wm_prob', 'label_file', 'mask_file']), name='inputspec')
    outputsSpec = pe.Node(interface=IdentityInterface(fields=['subject_id', 'subjects_dir',
                                                              'FreeSurferOutputDirectory', 'cnr_optimal_image']), name='outputspec')

    ### HACK: the nipype interface requires that this environmental variable is set before running
    print "HACK SETTING SUBJECTS_DIR {0}".format(constructed_FS_SUBJECTS_DIR)
    os.environ['SUBJECTS_DIR'] = constructed_FS_SUBJECTS_DIR
    inputsSpec.inputs.subjects_dir = constructed_FS_SUBJECTS_DIR  # HACK

    if RunMultiMode:
        mergeT1T2 = pe.Node(interface=Merge(2), name="Merge_T1T2")
        freesurferWF.connect(inputsSpec, 'T1_files', mergeT1T2, 'in1')
        freesurferWF.connect(inputsSpec, 'T2_files', mergeT1T2, 'in2')

        # Some constants based on assumpts about the label_file from BRAINSABC
        white_label = 1
        grey_label = 2

        msLDA_GenerateWeights = pe.Node(interface=MS_LDA(), name="MS_LDA")
        MSLDA_sge_options_dictionary = {'qsub_args': '-S /bin/bash -pe smp1 1 -l h_vmem=12G,mem_free=2G -o /dev/null -e /dev/null ' + CLUSTER_QUEUE, 'overwrite': True}
        msLDA_GenerateWeights.plugin_args = MSLDA_sge_options_dictionary
        msLDA_GenerateWeights.inputs.lda_labels = [white_label, grey_label]
        msLDA_GenerateWeights.inputs.weight_file = 'weights.txt'
        msLDA_GenerateWeights.inputs.use_weights = False
        msLDA_GenerateWeights.inputs.vol_synth_file = 'synth_out.nii.gz'
        # msLDA_GenerateWeights.inputs.vol_synth_file = 'synth_out.nii.gz'
        # msLDA_GenerateWeights.inputs.shift = 0 # value to shift by

        freesurferWF.connect(mergeT1T2, 'out', msLDA_GenerateWeights, 'images')
        # freesurferWF.connect(inputsSpec,'subjects_dir',  msLDA_GenerateWeights,'subjects_dir')
        freesurferWF.connect(inputsSpec, 'label_file', msLDA_GenerateWeights, 'label_file')
        # freesurferWF.connect(inputsSpec,'mask_file',  msLDA_GenerateWeights,'mask_file') ## Mask file MUST be unsigned char
        freesurferWF.connect(msLDA_GenerateWeights, 'vol_synth_file', outputsSpec, 'cnr_optimal_image')

    freesurfer_sge_options_dictionary = {'qsub_args': '-S /bin/bash -pe smp1 1 -l h_vmem=18G,mem_free=8G -o /dev/null -e /dev/null ' + CLUSTER_QUEUE, 'overwrite': True}
    if RunAllFSComponents == True:
        print("""Run FreeSurfer ReconAll at""")
        fs_reconall = pe.Node(interface=fswrap.FSScript(), name="FS52_custom")
        fs_reconall.plugin_args = freesurfer_sge_options_dictionary
        fs_reconall.inputs.subcommand = 'autorecon'
        # fs_reconall.inputs.directive = 'all'
        # fs_reconall.inputs.fs_env_script = '' # NOTE: NOT NEEDED HERE 'FreeSurferEnv.sh'
        # fs_reconall.inputs.fs_home = ''       # NOTE: NOT NEEDED HERE
        freesurferWF.connect(inputsSpec, 'FreeSurfer_ID', fs_reconall, 'subject_id')
        if RunMultiMode:
            ## Use the output of the synthesized T1 with maximized contrast
            freesurferWF.connect(msLDA_GenerateWeights, 'vol_synth_file', fs_reconall, 'T1_files')
        else:
            ## Use the output of the T1 only image
            freesurferWF.connect(inputsSpec, 'T1_files', fs_reconall, 'T1_files')

        computeFinalDirectory = pe.Node(Function(function=MakeFreeSurferOutputDirectory, input_names=['subjects_dir', 'subject_id'], output_names=['FreeSurferOutputDirectory']), run_without_submitting=True, name="99_computeFreeSurferOutputDirectory")
        freesurferWF.connect(inputsSpec, 'subjects_dir', computeFinalDirectory, 'subjects_dir')
        freesurferWF.connect(inputsSpec, 'FreeSurfer_ID', computeFinalDirectory, 'subject_id')

        freesurferWF.connect(inputsSpec, 'label_file', fs_reconall, 'brainmask')
        freesurferWF.connect(inputsSpec, 'subjects_dir', fs_reconall, 'subjects_dir')

        freesurferWF.connect(inputsSpec, 'FreeSurfer_ID', outputsSpec, 'subject_id')
        freesurferWF.connect(inputsSpec, 'subjects_dir', outputsSpec, 'subjects_dir')
        freesurferWF.connect(computeFinalDirectory, 'FreeSurferOutputDirectory', outputsSpec, 'FreeSurferOutputDirectory')
    return freesurferWF

def CreateFreeSurferLongitudinalWorkflow(projectid, subjectid, session_ids, WFname, CLUSTER_QUEUE, CLUSTER_QUEUE_LONG, RunAllFSComponents=True, RunMultiMode=True, constructed_FS_SUBJECTS_DIR='/never_use_this', subcommand='template'):
    """ Construct the longitudinal workflow
    Step 1: Construct the within-subject cross-sectional template (using all subject's sessions)
    Step 2: Construct the longitudinal subject results (for each session individually)
    """
    freesurfer_sge_options_dictionary = {'qsub_args': '-S /bin/bash -pe smp1 1 -l h_vmem=18G,mem_free=8G -o /dev/null -e /dev/null ' + CLUSTER_QUEUE, 'overwrite': True}
    freesurferWF = pe.Workflow(name=GenerateWFName(projectid, subjectid, '', WFname))
    inputsSpec = pe.Node(interface=IdentityInterface(fields=['FreeSurfer_ID',
                                                             'subjects_dir' ]),
                         name='inputspec')
    ### HACK: the nipype interface requires that this environmental variable is set before running
    print "HACK SETTING SUBJECTS_DIR {0}".format(constructed_FS_SUBJECTS_DIR)
    os.environ['SUBJECTS_DIR'] = constructed_FS_SUBJECTS_DIR
    inputsSpec.inputs.subjects_dir = constructed_FS_SUBJECTS_DIR  # HACK
    print("""Run FreeSurfer Within Subject Template at""")
    fs_template = pe.Node(interface=fswrap.FSScript(), name="FS55_custom_cross_sectional")
    fs_template.plugin_args = freesurfer_sge_options_dictionary
    fs_template.inputs.session_ids = session_ids
    fs_template.inputs.subcommand = 'template'
    freesurferWF.connect(inputsSpec, 'FreeSurfer_ID', fs_template, 'subject_id')
    fs_longitudinal = pe.Node(interface=fswrap.FSScript(),
                              name="FS55_custom_longitudinal")
    fs_longitudinal.plugin_args = freesurfer_sge_options_dictionary
    fs_longitudinal.inputs.subcommand = 'longitudinal'
    fs_longitudinal.iterables = ("session_ids", session_ids)
    freesurferWF.connect(fs_template, 'template_id', fs_longitudinal, 'template_id')
    outputsSpec = pe.Node(interface=IdentityInterface(fields=['FreeSurferTemplateDir',
                                                              'FreeSurferLongitudinalDir']),
                          name='outputspec')
    freesurferWF.connect(fs_template, 'outDir', outputsSpec, 'FreeSurferTemplateDir')
    freesurferWF.connect(fs_longitudinal, 'outDir', outputsSpec, 'FreeSurferLongitudinalDir')
