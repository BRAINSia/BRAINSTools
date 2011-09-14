#!/bin/bash



/raid0/homes/hjohnson/src/BRAINS3-Darwin-SuperBuildTest/src/bin/BRAINSInitializedControlPoints --inputVolume T1.nii.gz --splineGridSize 3,4,6 --permuteOrder 0,1,2 --outputLandmarksFile 012.fcsv --outputVolume 012.nii.gz
/raid0/homes/hjohnson/src/BRAINS3-Darwin-SuperBuildTest/src/bin/BRAINSInitializedControlPoints --inputVolume T1.nii.gz --splineGridSize 3,4,6 --permuteOrder 1,0,2 --outputLandmarksFile 102.fcsv --outputVolume 102.nii.gz
/raid0/homes/hjohnson/src/BRAINS3-Darwin-SuperBuildTest/src/bin/BRAINSInitializedControlPoints --inputVolume T1.nii.gz --splineGridSize 3,4,6 --permuteOrder 1,2,0 --outputLandmarksFile 120.fcsv --outputVolume 120.nii.gz
~/src/Slicer3-Darwin/Slicer3-build/Slicer3 [012]*.nii.gz [012]*.fcsv

