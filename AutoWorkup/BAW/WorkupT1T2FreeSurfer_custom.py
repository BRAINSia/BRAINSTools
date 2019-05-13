"""
WorkupT1T2FreeSurfer_custom.py
=================================
Description:

Author:

Usage:

"""


import os
from builtins import str

import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.freesurfer.model import MS_LDA
from nipype.interfaces.utility import Merge, IdentityInterface

from . import fswrap


def generate_wf_name(projectid, subjectid, sessionid, WFName):
    """
    This function...

    :param projectid:
    :param subjectid:
    :param sessionid:
    :param WFName:
    :return:
    """
    return WFName + "_" + str(subjectid) + "_" + str(sessionid) + "_" + str(projectid)


def create_free_surfer_workflow_custom(
    projectid,
    subjectid,
    sessionid,
    WFname,
    CLUSTER_QUEUE,
    CLUSTER_QUEUE_LONG,
    RunAllFSComponents=True,
    RunMultiMode=True,
    constructed_FS_SUBJECTS_DIR="/never_use_this",
):
    """
    This function...

    :param projectid:
    :param subjectid:
    :param sessionid:
    :param WFname:
    :param CLUSTER_QUEUE:
    :param CLUSTER_QUEUE_LONG:
    :param RunAllFSComponents: True
    :param RunMultiMode: True
    :param constructed_FS_SUBJECTS_DIR: '/never_use_this'
    :return:
    """
    freesurferWF = pe.Workflow(
        name=generate_wf_name(projectid, subjectid, sessionid, WFname)
    )

    inputsSpec = pe.Node(
        interface=IdentityInterface(
            fields=[
                "subj_session_id",
                "T1_files",
                "T2_files",
                "subjects_dir",
                "wm_prob",
                "label_file",
                "mask_file",
            ]
        ),
        name="inputspec",
    )
    outputsSpec = pe.Node(
        interface=IdentityInterface(
            fields=["full_path_FS_output", "processed_output_name", "cnr_optimal_image"]
        ),
        name="outputspec",
    )

    ### HACK: the nipype interface requires that this environmental variable is set before running
    print(("HACK SETTING SUBJECTS_DIR {0}".format(constructed_FS_SUBJECTS_DIR)))
    os.environ["SUBJECTS_DIR"] = constructed_FS_SUBJECTS_DIR
    inputsSpec.inputs.subjects_dir = constructed_FS_SUBJECTS_DIR  # HACK

    if RunMultiMode:
        mergeT1T2 = pe.Node(interface=Merge(2), name="Merge_T1T2")
        freesurferWF.connect(inputsSpec, "T1_files", mergeT1T2, "in1")
        freesurferWF.connect(inputsSpec, "T2_files", mergeT1T2, "in2")

        # Some constants based on assumpts about the label_file from BRAINSABC
        white_label = 1
        grey_label = 2

        msLDA_GenerateWeights = pe.Node(interface=MS_LDA(), name="MS_LDA")
        MSLDA_sge_options_dictionary = {
            "qsub_args": modify_qsub_args(CLUSTER_QUEUE, 2, 1, 1),
            "overwrite": True,
        }
        msLDA_GenerateWeights.plugin_args = MSLDA_sge_options_dictionary
        msLDA_GenerateWeights.inputs.lda_labels = [white_label, grey_label]
        msLDA_GenerateWeights.inputs.weight_file = "weights.txt"
        msLDA_GenerateWeights.inputs.use_weights = False
        msLDA_GenerateWeights.inputs.vol_synth_file = "synth_out.nii.gz"
        # msLDA_GenerateWeights.inputs.vol_synth_file = 'synth_out.nii.gz'
        # msLDA_GenerateWeights.inputs.shift = 0 # value to shift by

        freesurferWF.connect(mergeT1T2, "out", msLDA_GenerateWeights, "images")
        # freesurferWF.connect(inputsSpec,'subjects_dir',  msLDA_GenerateWeights,'subjects_dir')
        freesurferWF.connect(
            inputsSpec, "label_file", msLDA_GenerateWeights, "label_file"
        )
        # freesurferWF.connect(inputsSpec,'mask_file',  msLDA_GenerateWeights,'mask_file') ## Mask file MUST be unsigned char
        freesurferWF.connect(
            msLDA_GenerateWeights, "vol_synth_file", outputsSpec, "cnr_optimal_image"
        )

    if RunAllFSComponents == True:
        print("""Run FreeSurfer ReconAll at""")
        fs_reconall = pe.Node(
            interface=fswrap.FSScript(), name="FS52_cross_" + str(sessionid)
        )
        freesurfer_sge_options_dictionary = {
            "qsub_args": modify_qsub_args(CLUSTER_QUEUE, 8, 4, 4),
            "overwrite": True,
        }
        fs_reconall.plugin_args = freesurfer_sge_options_dictionary
        fs_reconall.inputs.subcommand = "autorecon"
        # fs_reconall.inputs.directive = 'all'
        # fs_reconall.inputs.fs_env_script = '' # NOTE: NOT NEEDED HERE 'FreeSurferEnv.sh'
        # fs_reconall.inputs.fs_home = ''       # NOTE: NOT NEEDED HERE
        freesurferWF.connect(
            inputsSpec, "subj_session_id", fs_reconall, "subj_session_id"
        )
        if RunMultiMode:
            ## Use the output of the synthesized T1 with maximized contrast
            ## HACK:  REMOVE FOR NOW - NEEDS FURTHER TESTING
            ## freesurferWF.connect(msLDA_GenerateWeights, 'vol_synth_file', fs_reconall, 'T1_files')
            freesurferWF.connect(inputsSpec, "T1_files", fs_reconall, "T1_files")
            ## END HACK
        else:
            ## Use the output of the T1 only image
            freesurferWF.connect(inputsSpec, "T1_files", fs_reconall, "T1_files")

        freesurferWF.connect(inputsSpec, "label_file", fs_reconall, "brainmask")
        freesurferWF.connect(inputsSpec, "subjects_dir", fs_reconall, "subjects_dir")
        freesurferWF.connect(fs_reconall, "outDir", outputsSpec, "full_path_FS_output")
        freesurferWF.connect(
            fs_reconall, "processed_output_name", outputsSpec, "processed_output_name"
        )
    return freesurferWF


def create_free_surfer_subject_template(
    projectid,
    subjectid,
    WFname,
    CLUSTER_QUEUE,
    CLUSTER_QUEUE_LONG,
    RunAllFSComponents=True,
    RunMultiMode=True,
    constructed_FS_SUBJECTS_DIR="/never_use_this",
    subcommand="template",
):
    """ Construct the longitudinal workflow
    Step 1: Construct the within-subject cross-sectional template (using all subject's sessions)

    :param projectid:
    :param subjectid:
    :param WFname:
    :param CLUSTER_QUEUE:
    :param CLUSTER_QUEUE_LONG:
    :param RunAllFSComponents: True
    :param RunMultiMode: True
    :param constructed_FS_SUBJECTS_DIR: '/never_use_this'
    :param subcommand: 'template'
    :return:
    """
    subjectTemplate_freesurferWF = pe.Workflow(
        name=generate_wf_name(projectid, subjectid, "", WFname)
    )
    inputsSpec = pe.Node(
        interface=IdentityInterface(
            fields=["base_template_id", "subjects_dir", "list_all_subj_session_ids"]
        ),
        name="inputspec",
    )
    ### HACK: the nipype interface requires that this environmental variable is set before running
    print(("HACK SETTING SUBJECTS_DIR {0}".format(constructed_FS_SUBJECTS_DIR)))
    os.environ["SUBJECTS_DIR"] = constructed_FS_SUBJECTS_DIR
    inputsSpec.inputs.subjects_dir = constructed_FS_SUBJECTS_DIR  # HACK
    print("""Run FreeSurfer Within Subject Template at""")
    fs_template = pe.Node(
        interface=fswrap.FSScript(), name="FS52_base_" + str(subjectid)
    )
    freesurfer_sge_options_dictionary = {
        "qsub_args": modify_qsub_args(CLUSTER_QUEUE, 8, 4, 4),
        "overwrite": True,
    }
    fs_template.plugin_args = freesurfer_sge_options_dictionary
    fs_template.inputs.subcommand = "template"
    subjectTemplate_freesurferWF.connect(
        inputsSpec, "subjects_dir", fs_template, "subjects_dir"
    )
    subjectTemplate_freesurferWF.connect(
        inputsSpec, "base_template_id", fs_template, "base_template_id"
    )
    subjectTemplate_freesurferWF.connect(
        inputsSpec,
        "list_all_subj_session_ids",
        fs_template,
        "list_all_subj_session_ids",
    )

    outputsSpec = pe.Node(
        interface=IdentityInterface(
            fields=["full_path_FS_output", "processed_output_name"]
        ),
        name="outputspec",
    )
    subjectTemplate_freesurferWF.connect(
        fs_template, "outDir", outputsSpec, "full_path_FS_output"
    )
    subjectTemplate_freesurferWF.connect(
        fs_template, "processed_output_name", outputsSpec, "processed_output_name"
    )

    return subjectTemplate_freesurferWF


def create_free_surfer_longitudinal_worklfow(
    projectid,
    subjectid,
    sessionid,
    WFname,
    CLUSTER_QUEUE,
    CLUSTER_QUEUE_LONG,
    RunAllFSComponents=True,
    RunMultiMode=True,
    constructed_FS_SUBJECTS_DIR="/never_use_this",
    subcommand="longitudinal",
):
    """ Construct the longitudinal workflow
    Step 2: Construct the longitudinal subject results (for each session individually)

    :param projectid:
    :param subjectid:
    :param WFname:
    :param CLUSTER_QUEUE:
    :param CLUSTER_QUEUE_LONG:
    :param RunAllFSComponents: True
    :param RunMultiMode: True
    :param constructed_FS_SUBJECTS_DIR: '/never_use_this'
    :param subcommand: 'template'
    :return:
    """
    long_freesurferWF = pe.Workflow(
        name=generate_wf_name(projectid, subjectid, sessionid, WFname)
    )
    inputsSpec = pe.Node(
        interface=IdentityInterface(
            fields=["base_template_id", "subj_session_id", "subjects_dir"]
        ),
        name="inputspec",
    )
    ### HACK: the nipype interface requires that this environmental variable is set before running
    print(("HACK SETTING SUBJECTS_DIR {0}".format(constructed_FS_SUBJECTS_DIR)))
    os.environ["SUBJECTS_DIR"] = constructed_FS_SUBJECTS_DIR
    inputsSpec.inputs.subjects_dir = constructed_FS_SUBJECTS_DIR  # HACK

    fs_longitudinal = pe.Node(
        interface=fswrap.FSScript(), name="FS52_long_" + str(sessionid)
    )
    freesurfer_sge_options_dictionary = {
        "qsub_args": modify_qsub_args(CLUSTER_QUEUE, 8, 4, 4),
        "overwrite": True,
    }
    fs_longitudinal.plugin_args = freesurfer_sge_options_dictionary
    fs_longitudinal.inputs.subcommand = "longitudinal"
    long_freesurferWF.connect(
        inputsSpec, "subjects_dir", fs_longitudinal, "subjects_dir"
    )
    long_freesurferWF.connect(
        inputsSpec, "subj_session_id", fs_longitudinal, "subj_session_id"
    )
    long_freesurferWF.connect(
        inputsSpec, "base_template_id", fs_longitudinal, "base_template_id"
    )
    outputsSpec = pe.Node(
        interface=IdentityInterface(
            fields=["full_path_FS_output", "processed_output_name"]
        ),
        name="outputspec",
    )
    long_freesurferWF.connect(
        fs_longitudinal, "outDir", outputsSpec, "full_path_FS_output"
    )
    long_freesurferWF.connect(
        fs_longitudinal, "processed_output_name", outputsSpec, "processed_output_name"
    )

    return long_freesurferWF
