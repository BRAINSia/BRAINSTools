#!/usr/bin/env python

from nipype.interfaces.base import CommandLine, CommandLineInputSpec, TraitedSpec, File, Directory
from nipype.interfaces.base import traits, isdefined, BaseInterface
from nipype.interfaces.utility import Merge, Split, Function, Rename, IdentityInterface
import nipype.interfaces.io as nio   # Data i/o
import nipype.pipeline.engine as pe  # pypeline engine

from BRAINSTools import *

"""
    from WorkupT1T2TissueClassify import CreateLandmarkInitializeWorkflow
    myLocalTCWF= CreateLandmarkInitializeWorkflow("01_LandmarkInitialize")
    landmarkInitializeWF.connect( [ (uidSource, myLocalTCWF, [(('uid', getFirstT1, subjectDatabaseFile ), 'inputsSpec.inputVolume')] ), ])
    landmarkInitializeWF.connect( BAtlas, 'template_landmarks_fcsv', myLocalTCWF,'inputsSpec.inputMovingLandmarkFilename')
    landmarkInitializeWF.connect( BAtlas, 'template_landmark_weights_csv', myLocalTCWF,'inputsSpec.inputWeightFilename')
    
"""
def CreateLandmarkInitializeWorkflow(WFname,CLUSTER_QUEUE,InterpolationMode):
    landmarkInitializeWF= pe.Workflow(name=WFname)

    inputsSpec = pe.Node(interface=IdentityInterface(fields=['inputVolume',
        'inputMovingLandmarkFilename','inputWeightFilename']), name='InputSpec' )

    ########################################################
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

    landmarkInitializeWF.connect(BCD,'outputLandmarksInACPCAlignedSpace', BLI,'inputFixedLandmarkFilename'),
    
    landmarkInitializeWF.connect(inputsSpec, 'inputMovingLandmarkFilename', BLI, 'inputMovingLandmarkFilename' )
    landmarkInitializeWF.connect(inputsSpec, 'inputWeightFilename', BLI, 'inputWeightFilename')
    
    #############
    outputsSpec = pe.Node(interface=IdentityInterface(fields=['atlasToSubjectTransform','outputLabels',
            't1_corrected','t2_corrected','outputAverageImages']), name='OutputSpec' )

    tissueClassifyWF.connect(bfc_files,'t1_corrected',outputsSpec,'t1_corrected')
    tissueClassifyWF.connect(bfc_files,'t2_corrected',outputsSpec,'t2_corrected')

    tissueClassifyWF.connect(BABCext,'outputVolumes',bfc_files, 'in_files')

    tissueClassifyWF.connect(BABCext,'atlasToSubjectTransform',outputsSpec,'atlasToSubjectTransform')
    tissueClassifyWF.connect(BABCext,'outputLabels',outputsSpec,'outputLabels')

    tissueClassifyWF.connect(BABCext,'outputAverageImages',outputsSpec,'outputAverageImages')

    return landmarkInitializeWF