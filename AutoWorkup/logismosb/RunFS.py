from workflow import create_fs_logb_workflow_for_both_hemispheres
import sqlite3
import os
from nipype import Workflow, DataSink, Node
from nipype.interfaces.freesurfer import ReconAll
from nipype.interfaces.io import FreeSurferSource

connection = sqlite3.connect("/Shared/johnsonhj/HDNI/20151001_AtrophySimulation/results.db")
cursor = connection.cursor()

num_threads = 12

base_dir = "/Shared/sinapse/CACHE/20161010_AtrophySimulation_Baseline_CACHE"
for row in cursor.execute("SELECT t1_image_file, t2_image_file, session_id FROM input"):
    session_id = str(row[2])
    t1_file = str(row[0])
    t2_file = str(row[1])

    wf = Workflow(name="FreeSurfer_{0}".format(session_id))

    subject_directory = os.path.dirname(os.path.dirname(t1_file))

    recon_all = Node(ReconAll(), "ReconAll")
    recon_all.inputs.T1_files = [t1_file]
    recon_all.inputs.T2_file = t2_file
    recon_all.inputs.openmp = num_threads
    recon_all.inputs.subject_id = "FreeSurfer"
    recon_all.inputs.flags = "-no-isrunning"
    recon_all.inputs.subjects_dir = os.path.join("/Shared/sinapse/CACHE/20161010_AtrophySimulation_Baseline",
                                                 session_id)
    recon_all.plugin_args = plugin_args={"qsub_args": "-q HJ,UI,all.q,COE -pe smp {0}".format(num_threads),
                                         "overwrite": True}

    hncma_atlas = os.path.join(subject_directory, "WarpedAtlas2Subject", "hncma_atlas.nii.gz")

    logb = create_fs_logb_workflow_for_both_hemispheres(name="FSLOGB_{0}".format(session_id),
                                                        plugin_args={
                                                            "qsub_args": "-q HJ,UI,all.q,COE -pe smp {0}".format(
                                                                num_threads),
                                                            "overwrite": True})
    logb.inputs.inputspec.rawavg = t1_file
    logb.inputs.inputspec.t2_raw = t2_file
    logb.inputs.inputspec.aseg_presurf = os.path.join(subject_directory, "FreeSurfer", "mri",
                                                      "aseg.presurf.mgz")
    if not os.path.isfile(logb.inputs.inputspec.aseg_presurf):
        print("could not find aseg")
        import sys
        sys.exit()
    logb.inputs.inputspec.hncma_atlas = hncma_atlas

    datasink = Node(DataSink(), name="DataSink")
    datasink.inputs.base_directory = recon_all.inputs.subjects_dir
    for hemisphere in ("lh", "rh"):
        fssource = Node(FreeSurferSource(), "{0}FSSource".format(hemisphere))
        fssource.inputs.hemi = hemisphere
        wf.connect([(recon_all, fssource, [("subject_id", "subject_id"),
                                           ("subjects_dir", "subjects_dir")]),
                    (fssource, logb, [("white", "inputspec.{0}_white".format(hemisphere))])])

        for matter in ("gm", "wm"):
            wf.connect(logb, "outputspec.{0}_{1}_surf_file".format(hemisphere, matter),
                       datasink, "LOGISMOSB.FreeSurfer.@{0}_{1}".format(hemisphere, matter))

    wf.base_dir = base_dir
    wf.config['execution']['job_finished_timeout'] = 120
    #wf.run(plugin="SGEGraph", plugin_args={"qsub_args": "-q HJ,all.q,COE,UI"})
    wf.run()
