from workflow import create_logb_workflow
import sqlite3


subjects_dir = "/Shared/sinapse/CACHE/20161010_AtrophySimulation_Baseline"
connection = sqlite3.connect("/Shared/johnsonhj/HDNI/20151001_AtrophySimulation/results.db")
cursor = connection.cursor()


import os
from nipype import Workflow, DataSink, Node
base_dir = "/Shared/sinapse/CACHE/20161010_AtrophySimulation_Baseline_CACHE"
for row in cursor.execute("SELECT t1_image_file, t2_image_file, session_id from input"):
    session_id = str(row[2])
    t1_file = str(row[0])
    t2_file = str(row[1])
    tissue_classify_directory = os.path.dirname(t1_file)
    csf_file = os.path.join(tissue_classify_directory, 'POSTERIOR_CSF.nii.gz')
    posterior_files = {"CSF": csf_file}
    abc_file = os.path.join(tissue_classify_directory, "complete_brainlabels_seg.nii.gz")
    subject_directory = os.path.dirname(tissue_classify_directory)
    hncma_atlas = os.path.join(subject_directory, "WarpedAtlas2Subject", "hncma_atlas.nii.gz")
    joint_fusion_file = os.path.join(subject_directory, "JointFusion", "JointFusion_HDAtlas20_2015_dustCleaned_label.nii.gz")
    logb_wf = create_logb_workflow(name="{0}_LOGISMOSB_Workflow".format(session_id))
    wf = Workflow("AtrophySim_Baseline_{0}".format(session_id))
    datasink = Node(DataSink(), name="DataSink")
    datasink.inputs.base_directory = os.path.join(subjects_dir, session_id)
    for hemisphere in ("lh", "rh"):
        for matter in ("gm", "wm"):
            wf.connect(logb_wf, "outputspec.{0}_{1}surface_file".format(hemisphere, matter),
                       datasink, "LOGISMOSB.@{0}_{1}".format(hemisphere, matter))
    logb_wf.inputs.inputspec.t1_file = t1_file
    logb_wf.inputs.inputspec.t2_file = t2_file
    logb_wf.inputs.inputspec.posterior_files = posterior_files
    logb_wf.inputs.inputspec.hncma_atlas = hncma_atlas
    logb_wf.inputs.inputspec.joint_fusion_file = joint_fusion_file
    logb_wf.inputs.inputspec.brainlabels_file = abc_file
    wf.base_dir = base_dir
    # wf.run(plugin="SGE", plugin_args={"qsub_args": "-q HJ,all.q,COE,UI"})
    # wf.run(plugin="MultiProc", plugin_args={"n_procs": 24})
    wf.run()

