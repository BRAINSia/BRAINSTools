## \author Hans J. Johnson
## This file contains the code necessary to build the python module
## nodes for SEM compliant tools
##

## NOTES ON HOW TO RUN THIS.
# build with all options turned on (even the non-default options)
# export BLD_DIR=/Users/johnsonhj/src/NEP-11
# export BRAINS_SRC=${BLD_DIR}/BRAINSTools
# cd ${BRAINS_SRC}/AutoWorkup; rm -rf ${BRAINS_SRC}/AutoWorkup/SEMTools;
# python local_generate_classes.py --python_paths=${BLD_DIR}/NIPYPE --program_paths=${BLD_DIR}/bin:${PATH} --output_path=${PWD}
# for i in $(find ${BRAINS_SRC}/AutoWorkup/SEMTools  -name "*.py"); do  autopep8 --max-line-length=300 -i ${i}; done


all_known_modules_list = [
'ResampleDTILogEuclidean',
'UnbiasedNonLocalMeans',
#'BatchMake',
'BinaryMaskEditorBasedOnLandmarks',
'BRAINSCreateLabelMapFromProbabilityMaps',
#'bmCheckMemory',
#'bmGridSend',
#'bmGridStore',
#'bmSliceExtractor',
'BRAINSABC',
'BRAINSAlignMSP',
'BRAINSClipInferior',
'BRAINSConstellationDetector',
'BRAINSConstellationModeler',
'BRAINSCut',
'BRAINSDemonWarp',
'BRAINSEyeDetector',
'BRAINSFit',
'BRAINSMultiSTAPLE',
'BRAINSInitializedControlPoints',
'BRAINSLandmarkInitializer',
'BRAINSLinearModelerEPCA',
'BRAINSLmkTransform',
'BRAINSMush',
'BRAINSPosteriorToContinuousClass',
'BRAINSResample',
'BRAINSResize',
'BRAINSROIAuto',
'BRAINSSnapShotWriter',
#'BRAINSTalairach',
#'BRAINSTalairachMask',
'BRAINSTransformConvert',
'BRAINSTransformFromFiducials',
'BRAINSTrimForegroundInDirection',
'CannyEdge',
'CannySegmentationLevelSetImageFilter',
'CleanUpOverlapLabels',
#'compareTractInclusion',
'DilateImage',
'DilateMask',
'DistanceMaps',
'dtiaverage',
'dtiestim',
'DTIPrep',
'dtiprocess',
#  ERROR:  invalid name 'DTI-Reg',
'DumpBinaryTrainingVectors',
#'dwiAtlas',
'DWICompare',
'DWIConvert',
'DWISimpleCompare',
'ErodeImage',
'ESLR',
'extractNrrdVectorIndex',
'fcsv_to_hdf5',
'fiberprocess',
'fiberstats',
'fibertrack',
'FiberViewerLight',
'FindCenterOfBrain',
'FlippedDifference',
'GenerateBrainClippedImage',
'GenerateCsfClippedFromClassifiedImage',
'GenerateLabelMapFromProbabilityMap',
'GenerateSummedGradientImage',
'GenerateTestImage',
'GradientAnisotropicDiffusionImageFilter',
'gtractAnisotropyMap',
'gtractAverageBvalues',
'gtractClipAnisotropy',
'gtractConcatDwi',
'gtractCopyImageOrientation',
'gtractCoRegAnatomy',
'gtractCoregBvalues',
'gtractCostFastMarching',
#'gtractCreateGuideFiber',
#'gtractFastMarchingTracking',
#'gtractFiberTracking',
'gtractImageConformity',
'gtractInvertBSplineTransform',
'gtractInvertDisplacementField',
'gtractInvertRigidTransform',
'gtractResampleAnisotropy',
'gtractResampleB0',
'gtractResampleCodeImage',
'gtractResampleDWIInPlace',
#'gtractResampleFibers',
'gtractTensor',
'gtractTransformToDisplacementField',
'HammerAttributeCreator',
'HistogramMatchingFilter',
#'ImageCalculator',
'ImageRegionPlotter',
'insertMidACPCpoint',
'JointHistogram',
'LandmarksCompare',
'landmarksConstellationAligner',
'landmarksConstellationWeights',
'maxcurvature',
'NeighborhoodMean',
'NeighborhoodMedian',
'NoiseGenerator',
'scalartransform',
'ShuffleVectorsModule',
'SimilarityIndex',
'SphericalCoordinateGeneration',
'STAPLEAnalysis',
'TextureFromNoiseImageFilter',
'TextureMeasureFilter',
'VBRAINSDemonWarp',
#'UKFTractography'
]

launcher = ['']
import sys
import os
import shutil
import argparse

print "Running: ",' '.join(sys.argv), "\n\n"

parser = argparse.ArgumentParser()
parser.add_argument("--python_paths", dest='python_paths', type=str, help="usually just the path the nipype")
parser.add_argument("--program_paths",dest='program_paths',type=str, help="where to find the sem programs")
parser.add_argument("--output_path",  dest='output_path',  type=str, help="the location where the SEMTools directory will be generated")

args = parser.parse_args()

# Platform specific information
#     Prepend the python search paths
PYTHON_AUX_PATHS = args.python_paths
PYTHON_AUX_PATHS = PYTHON_AUX_PATHS.split(':')
PYTHON_AUX_PATHS.extend(sys.path)
sys.path = PYTHON_AUX_PATHS
from nipype.interfaces.slicer.generate_classes import generate_all_classes

#\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
#####################################################################################
#     Prepend the shell environment search paths
PROGRAM_PATHS = args.program_paths
PROGRAM_PATHS = PROGRAM_PATHS.split(':')
PROGRAM_PATHS.extend(os.environ['PATH'].split(':'))
os.environ['PATH'] = ':'.join(PROGRAM_PATHS)

SEARCH_PATHS=args.program_paths.split(':')

OUTPUT_PATH = os.path.join(args.output_path,'SEMTools')
if os.path.exists(OUTPUT_PATH):
    shutil.rmtree(OUTPUT_PATH)
os.mkdir(OUTPUT_PATH)
os.chdir(OUTPUT_PATH)

found_modules_list = list()
missing_modules_list = list()

for  test_module in all_known_modules_list:
        for currPath in SEARCH_PATHS:
                currProgPath = os.path.join(currPath,test_module)
                if os.path.exists(currProgPath):
                    found_modules_list.append(test_module)
                    break
                else:
                    missing_modules_list.append(test_module)


generate_all_classes(modules_list=found_modules_list, launcher=[])

help_file = open(os.path.join(OUTPUT_PATH,'generated.sh'),'w')
help_file.write(' '.join(sys.argv)+"\n")
help_file.close()

print "\n\nRan: ",' '.join(sys.argv),"\n\n"

print "FOUND: ", found_modules_list
for test_module in missing_modules_list:
  print("Missing Candidate Program: {prog_name}".format(prog_name=test_module))
