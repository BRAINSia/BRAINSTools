## \author Ali Ghayoor
##
## This workflow gets the output DWI scan of CompressSensing and runs ukfTracts
##

import os
import SimpleITK as sitk
import nipype
from nipype.interfaces import ants
from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/oS
import nipype.pipeline.engine as pe  # pypeline engine
from SEMTools import *

def CreateTractographyWorkflow(WFname):
    TractWF = pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['DWI_Corrected_Aligned_CS', 'DWIBrainMask']),
                         name='inputsSpec')

    outputsSpec = pe.Node(interface=IdentityInterface(fields=['ukfTracks']),
                          name='outputsSpec')

    # Step1: UKF Processing
    UKFNode = pe.Node(interface=UKFTractography(), name= "UKFRunRecordStates")
    UKFNode.inputs.tracts = "ukfTracts.vtk"
    UKFNode.inputs.seedsPerVoxel = 10
    UKFNode.inputs.numTensor = '2'
    UKFNode.inputs.freeWater = True ## default False
    UKFNode.inputs.minFA = 0.06
    UKFNode.inputs.minGA = 0.06
    UKFNode.inputs.seedFALimit = 0.06
    UKFNode.inputs.Ql = 70
    UKFNode.inputs.recordLength = 2

    TractWF.connect(inputsSpec, 'DWI_Corrected_Aligned_CS', UKFNode, 'dwiFile')
    TractWF.connect(inputsSpec, 'DWIBrainMask', UKFNode, 'maskFile')
    TractWF.connect(UKFNode,'tracts',outputsSpec,'ukfTracks')

    return TractWF
