## \author Hans J. Johnson
## This file contains the code necessary to build the python module
## nodes for SEM compliant tools
##

## NOTES ON HOW TO RUN THIS.
# build wiht all options turned on (even the non-default options)
# export PATH=~/src/BT-build/bin:${PATH}
# export PYTHONPATH=~/src/BT-build/NIPYPE
# cd ~/src/BRAINSTools/AutoWorkup/SEMTools/; rm -rf ~/src/BRAINSTools/AutoWorkup/SEMTools/* ; python ../local_generate_classes.py
# for i in $(find ~/src/BRAINSTools/AutoWorkup/SEMTools  -name "*.py"); do  autopep8 --max-line-length=300 -i ${i}; done


all_known_modules_list = [
    #'ResampleDTILogEuclidean',
    'DTIPrep'
    'UKFTractography',
    'dtiprocess',
    'dtiestim',
    'dtiaverage',
    'BRAINSMultiSTAPLE',
    'AssignArray',
    'AverageScalarsByResampling',
    'BRAINSABC',
    'BRAINSAlignMSP',
    'BRAINSApplySurfaceLabels',
    'BRAINSAssignSurfaceFeatures',
    'BRAINSClipInferior',
    'BRAINSConstellationDetector',
    'BRAINSConstellationModeler',
    'BRAINSCreateLabelMapFromProbabilityMaps',
    'BRAINSCut',
    'BRAINSDemonWarp',
    'BRAINSEyeDetector',
    'BRAINSFit',
    'BRAINSInitializedControlPoints',
    'BRAINSLandmarkInitializer',
    'BRAINSLinearModelerEPCA',
    'BRAINSLmkTransform',
    'BRAINSMeasureSurface',
    'BRAINSMultiModeSegment',
    'BRAINSMush',
    'BRAINSPosteriorToContinuousClass',
    'BRAINSROIAuto',
    'BRAINSResample',
    'BRAINSResize',
    'BRAINSSnapShotWriter',
    'BRAINSSurfaceFlattening',
    'BRAINSSurfaceGeneration',
    'BRAINSTransformConvert',
    'BRAINSTransformFromFiducials',
    'BRAINSTrimForegroundInDirection',
    'BinaryMaskEditorBasedOnLandmarks',
    'CannyEdge',
    'CannySegmentationLevelSetImageFilter',
    'CleanUpOverlapLabels',
    'CombineLabels',
    'CompareSurfaces',
    'DWIConvert',
    'DilateImage',
    'DilateMask',
    'DistanceMaps',
    'ESLR',
    'ErodeImage',
    'FlippedDifference',
    'GenerateBrainClippedImage',
    'GenerateCsfClippedFromClassifiedImage',
    'GenerateLabelMapFromProbabilityMap',
    'GenerateSummedGradientImage',
    'GenerateTestImage',
    'GradientAnisotropicDiffusionImageFilter',
    'HammerAttributeCreator',
    'HistogramMatchingFilter',
    'IcosahedronResampler',
    'ImageRegionPlotter',
    'JointHistogram',
    'LabelMaps',
    'MultiResolutionRegistration',
    'NeighborhoodMean',
    'NeighborhoodMedian',
    'NoiseGenerator',
    'ProbabilityLabels',
    'QuadEdgeMeshClampScalars',
    'QuadEdgeMeshHistogramMatching',
    'QuadEdgeMeshPiecewiseRescale',
    'QuadEdgeMeshSimilarity',
    'RearrangeSurfaceLabels',
    'RemoveTinyLabels',
    'ResampleQuadEdgeMesh',
    'STAPLEAnalysis',
    'ShuffleVectorsModule',
    'SimilarityIndex',
    'SurfaceColor',
    'SurfaceLabelCleanUp',
    'TextureFromNoiseImageFilter',
    'TextureMeasureFilter',
    'VBRAINSDemonWarp',
    'WarpQuadEdgeMesh',
    'compareTractInclusion',
    'extractNrrdVectorIndex',
    'fcsv_to_hdf5',
    'gtractAnisotropyMap',
    'gtractAverageBvalues',
    'gtractClipAnisotropy',
    'gtractCoRegAnatomy',
    'gtractConcatDwi',
    'gtractCopyImageOrientation',
    'gtractCoregBvalues',
    'gtractCostFastMarching',
    'gtractCreateGuideFiber',
    'gtractFastMarchingTracking',
    'gtractFiberTracking',
    'gtractImageConformity',
    'gtractInvertBSplineTransform',
    'gtractInvertDisplacementField',
    'gtractInvertRigidTransform',
    'gtractResampleAnisotropy',
    'gtractResampleB0',
    'gtractResampleCodeImage',
    'gtractResampleDWIInPlace',
    'gtractResampleFibers',
    'gtractTensor',
    'gtractTransformToDisplacementField',
    'insertMidACPCpoint',
    'landmarksConstellationAligner',
    'landmarksConstellationWeights',
    'SmoothingMeshScalars'
]

launcher = ['']
import sys
import os
import shutil
import argparse

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

SEARCH_PATHS=args.program_paths.split(':')

OUTPUT_PATH = os.path.join(args.output_path,'SEMTools')
if os.path.exists(OUTPUT_PATH):
    shutil.rmtree(OUTPUT_PATH)
os.mkdir(OUTPUT_PATH)
os.chdir(OUTPUT_PATH)

found_modules_list = list()

for  test_module in all_known_modules_list:
        for currPath in SEARCH_PATHS:
                currProgPath = os.path.join(currPath,test_module)
                if os.path.exists(currProgPath):
                    found_modules_list.append(currProgPath)
                    break

print found_modules_list

generate_all_classes(modules_list=found_modules_list, launcher=[])

help_file = open(os.path.join(OUTPUT_PATH,'generated.sh'),'w')
help_file.write(' '.join(sys.argv)+"\n")
help_file.close()

print "Ran: ",' '.join(sys.argv)
