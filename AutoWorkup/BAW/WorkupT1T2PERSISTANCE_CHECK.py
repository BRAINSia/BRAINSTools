#!/usr/bin/env python
"""
WorkupT1T2PERSISTANCE_CHECK.py
================================
Description:

Author:

Usage:

"""


import nipype.interfaces.io as nio  # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces.base import (
    CommandLine,
    CommandLineInputSpec,
    TraitedSpec,
    File,
    Directory,
)
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.semtools import *
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface

"""
    from WorkupT1T2PERSISTANCE_CHECK import create_persistance_check_workflow
    myLocalPERSISTANCE_CHECKWF= create_persistance_check_workflow("999999_PersistanceCheckingWorkflow")
    PERSISTANCE_CHECKWF.connect(SplitAvgBABC,'avgBABCT1',myLocalPERSISTANCE_CHECKWF,'fixedVolume')
    PERSISTANCE_CHECKWF.connect(BABC,'outputLabels',myLocalPERSISTANCE_CHECKWF,'fixedBinaryVolume')
    PERSISTANCE_CHECKWF.connect(BAtlas,'template_t1_denoised_gaussian',myLocalPERSISTANCE_CHECKWF,'movingVolume')
    PERSISTANCE_CHECKWF.connect(BAtlas,'template_brain',myLocalPERSISTANCE_CHECKWF,'movingBinaryVolume')
    PERSISTANCE_CHECKWF.connect(BLI,'outputTransformFilename',myLocalPERSISTANCE_CHECKWF,'initialTransform')
"""


def create_persistance_check_workflow(WFname):
    """ The purpose of this workflow is to debug the automatic deletion of files from the output directory.

    :param WFname:
    :return:
    """
    PERSISTANCE_CHECKWF = pe.Workflow(name=WFname)

    inputsSpec = pe.Node(
        interface=IdentityInterface(
            fields=[
                "fixedVolume",
                "fixedBinaryVolume",
                "movingVolume",
                "movingBinaryVolume",
                "initialTransform",
            ]
        ),
        name="inputspec",
    )
    PERSISTANCE_CHECKWF.connect(inputsSpec, "subject_id", fs_reconall, "subject_id")
    PERSISTANCE_CHECKWF.connect(inputsSpec, "T1_files", fs_reconall, "T1_files")

    print("DOING FILE PERSISTANCE CHECK")
    PERSISTANCE_CHECK = pe.Node(
        interface=BRAINSFit(), name="99999_PERSISTANCE_CHECK_PERSISTANCE_CHECK"
    )
    PERSISTANCE_CHECK.inputs.costMetric = "MMI"
    PERSISTANCE_CHECK.inputs.debugLevel = 10
    PERSISTANCE_CHECK.inputs.maskProcessingMode = "ROI"
    PERSISTANCE_CHECK.inputs.numberOfSamples = 1000
    PERSISTANCE_CHECK.inputs.numberOfIterations = [1500]
    PERSISTANCE_CHECK.inputs.numberOfHistogramBins = 50
    PERSISTANCE_CHECK.inputs.maximumStepLength = 0.2
    PERSISTANCE_CHECK.inputs.minimumStepLength = [0.005]
    PERSISTANCE_CHECK.inputs.transformType = ["Affine"]
    PERSISTANCE_CHECK.inputs.maxBSplineDisplacement = 7
    PERSISTANCE_CHECK.inputs.maskInferiorCutOffFromCenter = 65
    PERSISTANCE_CHECK.inputs.splineGridSize = [28, 20, 24]
    PERSISTANCE_CHECK.inputs.outputVolume = "Trial_Initializer_Output.nii.gz"
    PERSISTANCE_CHECK.inputs.outputTransform = "Trial_Initializer_Output.h5"
    PERSISTANCE_CHECK.inputs.writeOutputTransformInFloat = True

    PERSISTANCE_CHECKWF.connect(
        inputsSpec, "fixedVolume", PERSISTANCE_CHECK, "fixedVolume"
    )
    PERSISTANCE_CHECKWF.connect(
        inputsSpec, "fixedBinaryVolume", PERSISTANCE_CHECK, "fixedBinaryVolume"
    )
    PERSISTANCE_CHECKWF.connect(
        inputsSpec, "movingVolume", PERSISTANCE_CHECK, "movingVolume"
    )
    PERSISTANCE_CHECKWF.connect(
        inputsSpec, "movingBinaryVolume", PERSISTANCE_CHECK, "movingBinaryVolume"
    )
    PERSISTANCE_CHECKWF.connect(
        inputsSpec, "initialTransform", PERSISTANCE_CHECK, "initialTransform"
    )

    outputsSpec = pe.Node(
        interface=IdentityInterface(fields=["outputVolume", "outputTransform"]),
        name="outputspec",
    )

    return PERSISTANCE_CHECKWF
