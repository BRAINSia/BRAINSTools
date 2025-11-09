#!/usr/bin/env python
"""
WorkupAtlasDustCleanup.py
===========================
Description:

Author:

Usage:

"""


import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.utility import Function, IdentityInterface


def run_automatic_cleanup_script(
    inFN1,
    inAtlas,
    outAtlas,
    maxIslandCount,
    useFullyConnected,
    forceLabelChange,
    noDilation,
    inFN2=None,
    includeList=None,
    excludeList=None,
):
    """
    This function...

    :param inFN1:
    :param inAtlas:
    :param outAtlast:
    :param maxIslandCount:
    :param useFullyConnected:
    :param forceLabelChange:
    :param noDilation:
    :param inFN2:
    :param includeList:
    :param excludeList:
    :return:
    """
    arguments = {
        "--inputT1Path": inFN1,
        "--inputT2Path": inFN2,
        "--inputAtlasPath": inAtlas,
        "--outputAtlasPath": outAtlas,
        "--maximumIslandVoxelCount": maxIslandCount,
        "--useFullyConnectedInConnectedComponentFilter": useFullyConnected,
        "--forceSuspiciousLabelChange": forceLabelChange,
        "--noDilation": noDilation,
        "--includeLabelsList": includeList,
        "--excludeLabelsList": excludeList,
    }

    from atlasSmallIslandCleanup import DustCleanup

    localDustCleanupObject = DustCleanup(arguments=arguments)
    localDustCleanupObject.main()

    import os

    return os.path.abspath(outAtlas)


def create_dust_cleanup_workflow(workflowFileName, onlyT1, master_config):
    """
    This function...

    :param workflowFileName:
    :param onlyT1:
    :param master_config:
    :return:
    """
    # if onlyT1:
    #    n_modality = 1
    # else:
    #    n_modality = 2
    # CLUSTER_QUEUE = master_config['queue']
    # CLUSTER_QUEUE_LONG = master_config['long_q']

    dustCleanupWF = pe.Workflow(name=workflowFileName)

    inputsSpec = pe.Node(
        interface=IdentityInterface(
            fields=[
                "subj_t1_image",  # Input T1 image
                "subj_t2_image",  # Input T2 image
                "subj_label_atlas",  # Input label atlas image
            ]
        ),
        run_without_submitting=True,
        name="inputspec",
    )

    outputsSpec = pe.Node(
        interface=IdentityInterface(
            fields=["JointFusion_HDAtlas20_2015_dustCleaned_label"]
        ),
        run_without_submitting=True,
        name="outputspec",
    )

    """
    Multimodal atlas dust cleanup if T2 exists to clean 'suspicious' dust. This stage builds islands
    using four-neighbor connectivity (useFullyConnected = False), has a max island count of 6 (instead
    of 5 as in the next stage), forces labels to change for these islands, does NOT dilate the label
    mask in order to clean all dust particles even clusters, and only suspicious (999) is cleaned.
    """
    sessionRunDustCleanupOnSuspicious = pe.Node(
        Function(
            function=run_automatic_cleanup_script,
            input_names=[
                "inFN1",
                "inFN2",
                "inAtlas",
                "outAtlas",
                "maxIslandCount",
                "useFullyConnected",
                "forceLabelChange",
                "noDilation",
                "includeList",
                "excludeList",
            ],
            output_names=["cleanedLabelImage"],
        ),
        run_without_submitting=True,
        name="sessionRunDustCleanupOnSuspicious",
    )
    dustCleanupWF.connect(
        inputsSpec, "subj_t1_image", sessionRunDustCleanupOnSuspicious, "inFN1"
    )
    if not onlyT1:
        dustCleanupWF.connect(
            inputsSpec, "subj_t2_image", sessionRunDustCleanupOnSuspicious, "inFN2"
        )
    else:
        pass

    dustCleanupWF.connect(
        inputsSpec, "subj_label_atlas", sessionRunDustCleanupOnSuspicious, "inAtlas"
    )
    sessionRunDustCleanupOnSuspicious.inputs.outAtlas = (
        "JointFusion_HDAtlas20_2015_dustCleaned_label_suspicious_maxIsland6.nii.gz"
    )
    sessionRunDustCleanupOnSuspicious.inputs.maxIslandCount = 6
    sessionRunDustCleanupOnSuspicious.inputs.useFullyConnected = False
    sessionRunDustCleanupOnSuspicious.inputs.forceLabelChange = True
    sessionRunDustCleanupOnSuspicious.inputs.noDilation = True
    sessionRunDustCleanupOnSuspicious.inputs.includeList = "999"
    sessionRunDustCleanupOnSuspicious.inputs.excludeList = None

    """
    Multimodal atlas dust cleanup if T2 exists to clean most labels after suspicious has been cleaned.
    Labels excluded from this cleaning stage may have viable isolated islands that should not be changed.
    This stage builds islands using eight-neighbor connectivity (useFullyConnected = True), has a max
    island count of 5 (instead of 6 as in the previous stage), forces labels to change for these islands,
    does dilate the label mask to avoid cleaning clustered dust, and several labels are excluded.
    """
    sessionRunDustCleanup = pe.Node(
        Function(
            function=run_automatic_cleanup_script,
            input_names=[
                "inFN1",
                "inFN2",
                "inAtlas",
                "outAtlas",
                "maxIslandCount",
                "useFullyConnected",
                "forceLabelChange",
                "noDilation",
                "includeList",
                "excludeList",
            ],
            output_names=["cleanedLabelImage"],
        ),
        run_without_submitting=True,
        name="sessionRunDustCleanup",
    )
    dustCleanupWF.connect(inputsSpec, "subj_t1_image", sessionRunDustCleanup, "inFN1")
    if not onlyT1:
        dustCleanupWF.connect(
            inputsSpec, "subj_t2_image", sessionRunDustCleanup, "inFN2"
        )
    else:
        pass

    dustCleanupWF.connect(
        sessionRunDustCleanupOnSuspicious,
        "cleanedLabelImage",
        sessionRunDustCleanup,
        "inAtlas",
    )
    sessionRunDustCleanup.inputs.outAtlas = (
        "JointFusion_HDAtlas20_2015_dustCleaned_label.nii.gz"
    )
    sessionRunDustCleanup.inputs.maxIslandCount = 5
    sessionRunDustCleanup.inputs.useFullyConnected = True
    sessionRunDustCleanup.inputs.forceLabelChange = True
    sessionRunDustCleanup.inputs.noDilation = False
    sessionRunDustCleanup.inputs.includeList = None
    sessionRunDustCleanup.inputs.excludeList = (
        "4,5,14,15,21,24,31,43,44,63,72,85,98,128,219,15000"
    )

    dustCleanupWF.connect(
        sessionRunDustCleanup,
        "cleanedLabelImage",
        outputsSpec,
        "JointFusion_HDAtlas20_2015_dustCleaned_label",
    )

    return dustCleanupWF
