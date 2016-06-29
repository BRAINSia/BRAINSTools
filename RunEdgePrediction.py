from AutoWorkup.logismosb.maclearn.workflows import create_logismosb_machine_learning_workflow
from AutoWorkup.logismosb.maclearn.nipype_interfaces import create_identity_transform
import sys
import os
import sqlite3
import SimpleITK as sitk
import glob
from nipype import Workflow, DataSink, Node
connection = sqlite3.connect("/Shared/johnsonhj/HDNI/20151001_AtrophySimulation/results.db")
cursor = connection.cursor()

base_dir = "/Shared/sinapse/CACHE/20160819_MachineLearning_baseline_CACHE"
results_dir = os.path.join(os.path.dirname(base_dir), os.path.basename(base_dir).replace("CACHE", "Results"))
if not os.path.isdir(base_dir):
    os.mkdir(base_dir)


def create_identity_transform_file(out_file):
    transform = create_identity_transform()
    sitk.WriteTransform(transform, out_file)
    return out_file

identity_transform_file = create_identity_transform_file(os.path.join(base_dir, "identity_transform.h5"))
wm_classifier_file = "/Shared/sinapse/CACHE/20160811_Davids_MachineLearning/Classifier/white_matter_classifier.pkl"
gm_classifier_file = "/Shared/sinapse/CACHE/20160811_Davids_MachineLearning/Classifier/gray_matter_classifier.pkl"


for row in cursor.execute("SELECT t1_image_file, t2_image_file, session_id from input"):
    session_id = str(row[2])
    t1_file = str(row[0])
    t2_file = str(row[1])
    tissue_classify_directory = os.path.dirname(t1_file)
    posterior_files = dict()
    for name in ["CSF", "VB", "CRBLWM", "CRBLGM", "SURFGM", "WM"]:
        posterior_files[name] = os.path.join(tissue_classify_directory, 'POSTERIOR_{0}.nii.gz'.format(name))
        if not os.path.exists(posterior_files[name]):
            print("FILE NOT FOUND: {0}".format(posterior_files[name]))
            sys.exit()

    abc_file = os.path.join(tissue_classify_directory, "complete_brainlabels_seg.nii.gz")

    subject_directory = os.path.dirname(tissue_classify_directory)
    hncma_atlas = os.path.join(subject_directory, "WarpedAtlas2Subject", "hncma_atlas.nii.gz")
    direction_files = dict()
    for name in ["rho", "phi", "theta"]:
        direction_files[name] = os.path.join(subject_directory, "WarpedAtlas2Subject", "{0}.nii.gz".format(name))

    lh_white_surface_file = os.path.join(subject_directory, "FreeSurfer", "surf", "lh.white")
    rh_white_surface_file = os.path.join(subject_directory, "FreeSurfer", "surf", "rh.white")

    logb_wf = create_logismosb_machine_learning_workflow()
    wf = Workflow("MachineLearning_Baseline_{0}".format(session_id))
    datasink = Node(DataSink(), name="DataSink")
    datasink.inputs.base_directory = os.path.join(results_dir, session_id)
    for hemisphere in ("lh", "rh"):
        for matter in ("gm", "wm"):
            wf.connect(logb_wf, "output_spec.{0}_{1}surface_file".format(hemisphere, matter),
                       datasink, "EdgePrediction.@{0}_{1}".format(hemisphere, matter))

    logb_wf.inputs.input_spec.t1_file = t1_file
    logb_wf.inputs.input_spec.orig_t1 = t1_file
    logb_wf.inputs.input_spec.t2_file = t2_file
    logb_wf.inputs.input_spec.posteriors = posterior_files
    logb_wf.inputs.input_spec.hncma_file = hncma_atlas
    logb_wf.inputs.input_spec.abc_file = abc_file
    # logb_wf.inputs.input_spec.acpc_transform = identity_transform_file
    logb_wf.inputs.input_spec.rho = direction_files["rho"]
    logb_wf.inputs.input_spec.theta = direction_files["theta"]
    logb_wf.inputs.input_spec.phi = direction_files["phi"]
    logb_wf.inputs.input_spec.lh_white_surface_file = lh_white_surface_file
    logb_wf.inputs.input_spec.rh_white_surface_file = rh_white_surface_file
    logb_wf.inputs.input_spec.wm_classifier_file = wm_classifier_file
    logb_wf.inputs.input_spec.gm_classifier_file = gm_classifier_file
    wf.base_dir = base_dir
    # wf.run(plugin="SGE", plugin_args={"qsub_args": "-q HJ,all.q,COE,UI"})
    # wf.run(plugin="MultiProc", plugin_args={"n_procs": 24})
    wf.run()
