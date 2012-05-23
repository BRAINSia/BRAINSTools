#!/usr/bin/env python

from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine

import os

from BRAINSTools import *

"""
    from WorkupT1T2LandmarkInitialization import CreateLandmarkInitializeWorkflow
    myLocalLMIWF= CreateLandmarkInitializeWorkflow("01_LandmarkInitialize")
    landmarkInitializeWF.connect( [ (uidSource, myLocalLMIWF, [(('uid', getFirstT1, subjectDatabaseFile ), 'inputsSpec.inputVolume')] ), ])
    landmarkInitializeWF.connect( BAtlas, 'template_landmarks_31_fcsv', myLocalLMIWF,'inputsSpec.atlasLandmarkFilename')
    landmarkInitializeWF.connect( BAtlas, 'template_landmark_weights_31_csv', myLocalLMIWF,'inputsSpec.atlasWeightFilename')
    landmarkInitializeWF.connect(BAtlas,'template_t1',myLocalLMIWF,'inputsSpec.atlasVolume')
    
"""
def CreateLandmarkInitializeWorkflow(WFname,BCD_model_path,InterpolationMode,DoReverseInit=False):
    landmarkInitializeWF= pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['inputVolume',
        'atlasLandmarkFilename','atlasWeightFilename','atlasVolume']), name='InputSpec' )

    ########################################################/
    # Run ACPC Detect on first T1 Image - Base Image
    ########################################################
    BCD = pe.Node(interface=BRAINSConstellationDetector(), name="01_BCD")
    ##  Use program default BCD.inputs.inputTemplateModel = T1ACPCModelFile
    ##BCD.inputs.outputVolume =   "BCD_OUT" + "_ACPC_InPlace.nii.gz"                #$# T1AcpcImageList
    BCD.inputs.outputTransform =  "BCD" + "_Original2ACPC_transform.mat"
    BCD.inputs.outputResampledVolume = "BCD" + "_ACPC.nii.gz"
    BCD.inputs.outputLandmarksInInputSpace = "BCD" + "_Original.fcsv"
    BCD.inputs.outputLandmarksInACPCAlignedSpace = "BCD" + "_ACPC_Landmarks.fcsv"
    #BCD.inputs.outputMRML = "BCD" + "_Scene.mrml"
    BCD.inputs.interpolationMode = InterpolationMode
    BCD.inputs.houghEyeDetectorMode = 1  # Look for dark eyes like on a T1 image, 0=Look for bright eyes like in a T2 image
    BCD.inputs.acLowerBound = 80.0 # Chop the data set 80mm below the AC PC point.
    BCD.inputs.LLSModel = os.path.join(BCD_model_path,'LLSModel-2ndVersion.hdf5')
    BCD.inputs.inputTemplateModel = os.path.join(BCD_model_path,'T1-2ndVersion.mdl')

    # Entries below are of the form:
    landmarkInitializeWF.connect( inputsSpec , 'inputVolume', BCD, 'inputVolume')

    ########################################################
    # Run BLI atlas_to_subject
    ########################################################
    BLI = pe.Node(interface=BRAINSLandmarkInitializer(), name="05_BLI")
    BLI.inputs.outputTransformFilename = "landmarkInitializer_atlas_to_subject_transform.mat"

    landmarkInitializeWF.connect(inputsSpec, 'atlasWeightFilename', BLI, 'inputWeightFilename')
    landmarkInitializeWF.connect(inputsSpec, 'atlasLandmarkFilename', BLI, 'inputMovingLandmarkFilename' )
    landmarkInitializeWF.connect(BCD,'outputLandmarksInACPCAlignedSpace', BLI,'inputFixedLandmarkFilename'),
    
    ## This is for debugging purposes, and it is not intended for general use.
    if DoReverseInit == True:
        ########################################################
        # Run BLI subject_to_atlas
        ########################################################
        BLI2Atlas = pe.Node(interface=BRAINSLandmarkInitializer(), name="05_BLI2Atlas")
        BLI2Atlas.inputs.outputTransformFilename = "landmarkInitializer_subject_to_atlas_transform.mat"

        landmarkInitializeWF.connect(inputsSpec, 'atlasWeightFilename', BLI2Atlas, 'inputWeightFilename')
        landmarkInitializeWF.connect(inputsSpec, 'atlasLandmarkFilename', BLI2Atlas, 'inputFixedLandmarkFilename' )
        landmarkInitializeWF.connect(BCD,'outputLandmarksInInputSpace',BLI2Atlas,'inputMovingLandmarkFilename')
        
        Resample2Atlas=pe.Node(interface=BRAINSResample(),name="05_Resample2Atlas")
        Resample2Atlas.inputs.interpolationMode = "Linear"
        Resample2Atlas.inputs.outputVolume = "subject2atlas.nii.gz"

        landmarkInitializeWF.connect( inputsSpec , 'inputVolume', Resample2Atlas, 'inputVolume')
        landmarkInitializeWF.connect(BLI2Atlas,'outputTransformFilename',Resample2Atlas,'warpTransform')
        landmarkInitializeWF.connect(inputsSpec,'atlasVolume',Resample2Atlas,'referenceVolume')
        
    DO_DEBUG = True
    if DO_DEBUG == True:
        ResampleFromAtlas=pe.Node(interface=BRAINSResample(),name="05_ResampleFromAtlas")
        ResampleFromAtlas.inputs.interpolationMode = "Linear"
        ResampleFromAtlas.inputs.outputVolume = "atlas2subject.nii.gz"

        landmarkInitializeWF.connect( inputsSpec , 'atlasVolume', ResampleFromAtlas, 'inputVolume')
        landmarkInitializeWF.connect(BLI,'outputTransformFilename',ResampleFromAtlas,'warpTransform')
        landmarkInitializeWF.connect(BCD,'outputResampledVolume',ResampleFromAtlas,'referenceVolume')
    
    #############
    outputsSpec = pe.Node(interface=IdentityInterface(fields=['outputLandmarksInACPCAlignedSpace','outputResampledVolume','outputLandmarksInInputSpace',
            'outputTransform','outputMRML','atlasToSubjectTransform'
            ]), name='OutputSpec' )

    landmarkInitializeWF.connect(BCD,'outputLandmarksInACPCAlignedSpace',outputsSpec,'outputLandmarksInACPCAlignedSpace')
    landmarkInitializeWF.connect(BCD,'outputResampledVolume',outputsSpec,'outputResampledVolume')
    landmarkInitializeWF.connect(BCD,'outputLandmarksInInputSpace',outputsSpec,'outputLandmarksInInputSpace')
    landmarkInitializeWF.connect(BCD,'outputTransform',outputsSpec,'outputTransform')
    landmarkInitializeWF.connect(BCD,'outputMRML',outputsSpec,'outputMRML')
    landmarkInitializeWF.connect(BLI,'outputTransformFilename',outputsSpec,'atlasToSubjectTransform')

    return landmarkInitializeWF
