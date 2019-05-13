## \author Ali Ghayoor
##
## This workflow gets the output DWI scan of CompressSensing and runs ukfTracts
##
"""
TractogrpahyWorkflow.py
============================
Description:
    The purpose of this is to...

Author:
    Ali Ghayoor
    
Usage:

"""
import os
from functools import reduce

import SimpleITK as sitk
import nipype
import nipype.interfaces.io as nio  # Data i/oS
import nipype.pipeline.engine as pe  # pypeline engine
from nipype.interfaces import ants
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


def create_tractography_workflow(WFname):
    """
    This Function takes in...

    :param WFname:
    :return:
    """
    ###### UTILITY FUNCTIONS #######
    def computer_number_of_seeds_per_vexel(inputVolume):
        """
        This Function takes in...

        :param inputVolume:
        :return:
        """
        import operator
        import SimpleITK as sitk

        inVol = sitk.ReadImage(inputVolume)
        voxelVolume = reduce(operator.mul, inVol.GetSpacing())
        # 10 seeds per voxel is used when voxel voluem is 8 mm^3.
        seedsPerVoxel = round(voxelVolume * 10 / 8)
        return int(seedsPerVoxel)

    #################################
    TractWF = pe.Workflow(name=WFname)

    inputsSpec = pe.Node(
        interface=IdentityInterface(
            fields=["DWI_Corrected_Aligned_CS", "DWIBrainMask"]
        ),
        name="inputsSpec",
    )

    outputsSpec = pe.Node(
        interface=IdentityInterface(fields=["ukfTracks"]), name="outputsSpec"
    )

    ########
    # Before running the UKF, we need to define the number of seeds per voxel
    # based on the voxel volume of the input DWI scan.
    ########

    # Step1: extract B0 from DWI volume
    EXTRACT_B0 = pe.Node(interface=extractNrrdVectorIndex(), name="EXTRACT_B0")
    EXTRACT_B0.inputs.vectorIndex = 0
    EXTRACT_B0.inputs.outputVolume = "B0_Image.nrrd"
    TractWF.connect(inputsSpec, "DWI_Corrected_Aligned_CS", EXTRACT_B0, "inputVolume")

    # Step2: Compute number of seeds per voxel
    computeNumberOfSeedsPerVoxelNode = pe.Node(
        interface=Function(
            function=computer_number_of_seeds_per_vexel,
            input_names=["inputVolume"],
            output_names=["seedsPerVoxel"],
        ),
        name="ComputeNumberOfSeedsPerVoxel",
    )
    TractWF.connect(
        EXTRACT_B0, "outputVolume", computeNumberOfSeedsPerVoxelNode, "inputVolume"
    )

    # Step3: UKF Processing
    UKFNode = pe.Node(interface=UKFTractography(), name="UKFRunRecordStates")
    UKFNode.inputs.tracts = "ukfTracts.vtp"
    UKFNode.inputs.numTensor = "2"
    UKFNode.inputs.freeWater = True  ## default False
    UKFNode.inputs.minFA = 0.06
    UKFNode.inputs.minGA = 0.06
    UKFNode.inputs.seedFALimit = 0.06
    UKFNode.inputs.Ql = 70
    UKFNode.inputs.recordLength = 2
    UKFNode.inputs.recordTensors = True
    UKFNode.inputs.recordFreeWater = True
    UKFNode.inputs.recordFA = True
    UKFNode.inputs.recordTrace = True

    TractWF.connect(inputsSpec, "DWI_Corrected_Aligned_CS", UKFNode, "dwiFile")
    TractWF.connect(inputsSpec, "DWIBrainMask", UKFNode, "maskFile")
    TractWF.connect(
        computeNumberOfSeedsPerVoxelNode, "seedsPerVoxel", UKFNode, "seedsPerVoxel"
    )
    TractWF.connect(UKFNode, "tracts", outputsSpec, "ukfTracks")

    return TractWF
