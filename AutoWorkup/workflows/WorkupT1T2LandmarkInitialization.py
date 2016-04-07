#!/usr/bin/env python

from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine

import os

from nipype.interfaces.semtools.segmentation.specialized import BRAINSConstellationDetector
from nipype.interfaces.semtools.utilities.brains import BRAINSLandmarkInitializer
from nipype.interfaces.semtools.registration.brainsresample import BRAINSResample
from nipype.interfaces.semtools.segmentation.specialized import BRAINSROIAuto

from utilities.distributed import modify_qsub_args

"""
    from WorkupT1T2LandmarkInitialization import CreateLandmarkInitializeWorkflow
    myLocalLMIWF= CreateLandmarkInitializeWorkflow("LandmarkInitialize")
    landmarkInitializeWF.connect( [ (uidSource, myLocalLMIWF, [(('uid', getFirstT1, subjectDatabaseFile ), 'inputsSpec.inputVolume')] ), ])

    landmarkInitializeWF.connect( BAtlas, 'template_landmarks_50Lmks_fcsv', myLocalLMIWF,'inputsSpec.atlasLandmarkFilename')
    landmarkInitializeWF.connect( BAtlas, 'template_weights_50Lmks_wts', myLocalLMIWF,'inputsSpec.atlasWeightFilename')
    landmarkInitializeWF.connect( BAtlas, 'LLSModel_50Lmks_h5', myLocalLMIWF, 'inputspec.LLSModel')
    landmarkInitializeWF.connect( BAtlas, 'T1_50Lmks_mdl', myLocalLMIWF, 'inputspec.inputTemplateModel')

    landmarkInitializeWF.connect(BAtlas,'template_t1_denoised_gaussian',myLocalLMIWF,'inputsSpec.atlasVolume')

"""


def CreateLandmarkInitializeWorkflow(WFname, master_config, InterpolationMode, PostACPCAlignToAtlas, DoReverseInit, useEMSP=False, Debug=False):
    CLUSTER_QUEUE=master_config['queue']
    CLUSTER_QUEUE_LONG=master_config['long_q']
    landmarkInitializeWF = pe.Workflow(name=WFname)

    #############
    inputsSpec = pe.Node(interface=IdentityInterface(fields=['inputVolume',
                                                             'atlasLandmarkFilename',
                                                             'atlasWeightFilename',
                                                             'LLSModel',
                                                             'inputTemplateModel',
                                                             'atlasVolume',
                                                             'EMSP']),
                         run_without_submitting=True,
                         name='inputspec')

    #############
    outputsSpec = pe.Node(interface=IdentityInterface(fields=['outputLandmarksInACPCAlignedSpace',
                                                              'outputResampledVolume', 'outputResampledCroppedVolume',
                                                              'outputLandmarksInInputSpace',
                                                              'writeBranded2DImage',
                                                              'outputTransform', 'outputMRML', 'atlasToSubjectTransform'
                                                              ]),
                          run_without_submitting=True,
                          name='outputspec')

    ########################################################/
    # Run ACPC Detect on first T1 Image - Base Image
    ########################################################
    BCD = pe.Node(interface=BRAINSConstellationDetector(), name="BCD")
    many_cpu_BCD_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE_LONG,4,2,4), 'overwrite': True}
    BCD.plugin_args = many_cpu_BCD_options_dictionary
    ##  Use program default BCD.inputs.inputTemplateModel = T1ACPCModelFile
    # BCD.inputs.outputVolume =   "BCD_OUT" + "_ACPC_InPlace.nii.gz"                #$# T1AcpcImageList
    BCD.inputs.outputTransform = "BCD" + "_Original2ACPC_transform.h5"
    BCD.inputs.outputResampledVolume = "BCD" + "_ACPC.nii.gz"
    BCD.inputs.outputLandmarksInInputSpace = "BCD" + "_Original.fcsv"
    BCD.inputs.outputLandmarksInACPCAlignedSpace = "BCD" + "_ACPC_Landmarks.fcsv"
    BCD.inputs.writeBranded2DImage = "BCD"+"_Branded2DQCimage.png"
    # BCD.inputs.outputMRML = "BCD" + "_Scene.mrml"
    BCD.inputs.interpolationMode = InterpolationMode
    BCD.inputs.houghEyeDetectorMode = 1  # Look for dark eyes like on a T1 image, 0=Look for bright eyes like in a T2 image
    BCD.inputs.acLowerBound = 80.0  # Chop the data set 80mm below the AC PC point.

    # Entries below are of the form:
    landmarkInitializeWF.connect(inputsSpec, 'inputVolume', BCD, 'inputVolume')
    landmarkInitializeWF.connect(inputsSpec, 'atlasWeightFilename',  BCD, 'atlasLandmarkWeights')
    landmarkInitializeWF.connect(inputsSpec, 'atlasLandmarkFilename',BCD, 'atlasLandmarks')

    landmarkInitializeWF.connect(inputsSpec, 'LLSModel',             BCD, 'LLSModel')
    landmarkInitializeWF.connect(inputsSpec, 'inputTemplateModel',   BCD, 'inputTemplateModel')

    # If EMSP, pre-selected landmarks are given, force to use.
    if useEMSP:
        print("*** Use pre-selected landmark file for Landmark Detection")
        landmarkInitializeWF.connect(inputsSpec, 'EMSP', BCD, 'inputLandmarksEMSP')


    # If the atlas volume is from this subject (i.e. after template building for the longitudinal phase) then set this to True
    # Otherwise, it is probably best to let the ACPC alignment be fully defined by the landmark points themselves.
    if PostACPCAlignToAtlas:
        landmarkInitializeWF.connect(inputsSpec, 'atlasVolume',          BCD, 'atlasVolume')

    ########################################################
    # Run BLI atlas_to_subject
    ########################################################
    BLI = pe.Node(interface=BRAINSLandmarkInitializer(), name="BLI")
    BLI.inputs.outputTransformFilename = "landmarkInitializer_atlas_to_subject_transform.h5"

    landmarkInitializeWF.connect(inputsSpec, 'atlasWeightFilename', BLI, 'inputWeightFilename')
    landmarkInitializeWF.connect(inputsSpec, 'atlasLandmarkFilename', BLI, 'inputMovingLandmarkFilename')
    landmarkInitializeWF.connect(BCD, 'outputLandmarksInACPCAlignedSpace', BLI, 'inputFixedLandmarkFilename')

    ## This is for debugging purposes, and it is not intended for general use.
    if DoReverseInit == True:
        ########################################################
        # Run BLI subject_to_atlas
        ########################################################
        BLI2Atlas = pe.Node(interface=BRAINSLandmarkInitializer(), name="BLI2Atlas")
        BLI2Atlas.inputs.outputTransformFilename = "landmarkInitializer_subject_to_atlas_transform.h5"

        landmarkInitializeWF.connect(inputsSpec, 'atlasWeightFilename', BLI2Atlas, 'inputWeightFilename')
        landmarkInitializeWF.connect(inputsSpec, 'atlasLandmarkFilename', BLI2Atlas, 'inputFixedLandmarkFilename')
        landmarkInitializeWF.connect(BCD, 'outputLandmarksInInputSpace', BLI2Atlas, 'inputMovingLandmarkFilename')

        Resample2Atlas = pe.Node(interface=BRAINSResample(), name="Resample2Atlas")
        Resample2Atlas.inputs.interpolationMode = "Linear"
        Resample2Atlas.inputs.outputVolume = "subject2atlas.nii.gz"

        landmarkInitializeWF.connect(inputsSpec, 'inputVolume', Resample2Atlas, 'inputVolume')
        landmarkInitializeWF.connect(BLI2Atlas, 'outputTransformFilename', Resample2Atlas, 'warpTransform')

    if (DoReverseInit == True) and (Debug == True):
        ResampleFromAtlas = pe.Node(interface=BRAINSResample(), name="ResampleFromAtlas")
        ResampleFromAtlas.inputs.interpolationMode = "Linear"
        ResampleFromAtlas.inputs.outputVolume = "atlas2subject.nii.gz"

        landmarkInitializeWF.connect(inputsSpec, 'atlasVolume', ResampleFromAtlas, 'inputVolume')
        landmarkInitializeWF.connect(BLI, 'outputTransformFilename', ResampleFromAtlas, 'warpTransform')
        landmarkInitializeWF.connect(BCD, 'outputResampledVolume', ResampleFromAtlas, 'referenceVolume')

    BROIAUTO = pe.Node(interface=BRAINSROIAuto(), name="BROIAuto_cropped")
    many_cpu_BROIAUTO_options_dictionary = {'qsub_args': modify_qsub_args(CLUSTER_QUEUE_LONG,4,2,4), 'overwrite': True}
    BROIAUTO.plugin_args = many_cpu_BROIAUTO_options_dictionary
    BROIAUTO.inputs.outputVolume = "Cropped_BCD_ACPC_Aligned.nii.gz"
    BROIAUTO.inputs.ROIAutoDilateSize = 10
    BROIAUTO.inputs.cropOutput = True
    landmarkInitializeWF.connect(BCD, 'outputResampledVolume', BROIAUTO, 'inputVolume')

    landmarkInitializeWF.connect(BROIAUTO, 'outputVolume', outputsSpec, 'outputResampledCroppedVolume')
    landmarkInitializeWF.connect(BCD, 'outputLandmarksInACPCAlignedSpace', outputsSpec, 'outputLandmarksInACPCAlignedSpace')
    landmarkInitializeWF.connect(BCD, 'outputResampledVolume', outputsSpec, 'outputResampledVolume')
    landmarkInitializeWF.connect(BCD, 'outputLandmarksInInputSpace', outputsSpec, 'outputLandmarksInInputSpace')
    landmarkInitializeWF.connect(BCD, 'outputTransform', outputsSpec, 'outputTransform')
    landmarkInitializeWF.connect(BCD, 'outputMRML', outputsSpec, 'outputMRML')
    landmarkInitializeWF.connect(BCD, 'writeBranded2DImage', outputsSpec, 'writeBranded2DImage')
    landmarkInitializeWF.connect(BLI, 'outputTransformFilename', outputsSpec, 'atlasToSubjectTransform')

    return landmarkInitializeWF
