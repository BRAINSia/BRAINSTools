The BRAINSConstellationDetector project is designed for image registration, constellation (arbitrary landmark) detection, landmark visualization and manipulation for human brain MR scans. It is a subproject/plug-in for the BRAINSTools image analysis toolkit.

This project was inspired by a collaboration to identify AC-PC points from a tool written by Dr. Babak Ardekani.  This worked greatly expanded upon that original concept to allow for arbitrary landmarks and much more generalized implementation for broad class of medical images by the ITK implementors at the University of Iowa, including Kent Williams, Hans J. Johnson, Vincent A. Magnotta, Ali Ghayoor, Eun Young Kim and Greg Harris.

HISTORICAL OUTDATED DOCUMENTATION (Do not trust anything below this line)
=================================

INSTALLATION:
------------
The script buildBRAINSConstellationDetector.sh will attempt to build everything from one script on a unix platform.

README For BRAINSConstellationDetector Image Analysis Tool
Author: Wei Lu

Introduction

BRAINS:
http://www.nitrc.org/projects/brains/
BRAINSConstellationDetector:
http://www.nitrc.org/projects/brainscdetector/

Objective
This README file will guide you on how to
1.  Download the source code
2.  Build the non-GUI version
3.  Build the Qt-based GUI version
4.  Run the software
Note:
A rtf version of this README file is also available at
https://www.nitrc.org/svn/brainscdetector/branches/wei-lu/MedicalSliceViewer/brainscdetector/docs/README.rtf
You will need a NITRC account to access this file.

Download the project source code
1.  Download and install SVN client at
http://subversion.apache.org/packages.html
2.  Download BRAINSConstellationDetector source files from NITRC SVN repository by:
svn checkout https://www.nitrc.org/svn/brainscdetector/branches/wei-lu/MedicalSliceViewer

Build the non-GUI version
Use bash script (in progress)
Run the following command:
mkdir brainscdetector-build
cd brainscdetector-build
bash ../MedicalSliceViewer/brainscdetector/buildTool.sh

Use CMake
1.  Download CMake source file at
http://www.cmake.org/cmake/resources/software.html
2.  Build and install CMake (version >= 2.8)
3.  Run the following command:
3.1  For Mac OSX 64-bit
mkdir brainscdetector-build
cd brainscdetector-build
cmake –C ../MedicalSliceViewer/InitialCache-OSX64.cmake ../MedicalSliceViewer
make
3.2  For Linux 64-bit
mkdir brainscdetector-build
cd brainscdetector-build
cmake –C ../MedicalSliceViewer/InitialCache-Linux64.cmake ../MedicalSliceViewer
make

Build the Qt-based GUI version
Use bash script (in progress)
Run the following command:
mkdir brainscdetector-build
cd brainscdetector-build
bash ../MedicalSliceViewer/brainscdetector/buildTool.sh

Use CMake
1.  Download Qt libraries (version >= 4.6.2) at
http://qt.nokia.com/downloads
2.  Install the binary package
3.  Download CMake source file (version >= 2.8) at
http://www.cmake.org/cmake/resources/software.html
4.  Build and install CMake
5.  Run the following command:

5.1  For Mac OSX 64-bit
mkdir brainscdetector-build
cd brainscdetector-build
cmake –DBUILD_MEDICAL_SLICEVIEWER_GUI=ON –C ../MedicalSliceViewer/InitialCache-OSX64.cmake ../MedicalSliceViewer
5.2  For Linux 64-bit
mkdir brainscdetector-build
cd brainscdetector-build
cmake –DBUILD_MEDICAL_SLICEVIEWER_GUI=ON –C ../MedicalSliceViewer/InitialCache-Linux64.cmake ../MedicalSliceViewer
6.  Locate qmake install path for makefiles by for e.g. ccmake

7.  Run make

Run the software
There are two important executables: the modeler and the detector. The modeler helps the user to train his/her model with arbitrary landmarks. The detector will try to estimate those landmarks positions for a given new brain MR image.

Run the modeler (outdated)
Usage:
./BRAINSConstellationModeler
-i <trainingFile> or --trainingFile <trainingFile>
Where
1.  trainingFile is the input training files list filename without eye centers.
For more details on setting optional parameters please type ./BRAINSConstellationModeler -h


A sample trained model file can be found at
brainscdetector/TestData/input/newT1InputModelWithVN4.mdl
The corresponding trained llsModel file can be found at
brainscdetector/TestData/input/llsModel.txt


After running the detector, a Slicer fcsv landmark file and a Slicer mrml scene file for those named landmarks in the original physical space will be also written to the working directory.
