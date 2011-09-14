##There is a program that reads in a T1 image, determine the AC, PC, VN4, and MPJ
##points and writes out a slicer .fcsv file with 3 landmark points in it
##There are 2 files to be read in order to confirm that the landmarks are
##being picked correctly.
set volumesLogic [$::slicer3::VolumesGUI GetLogic]
set fidLogic [$::slicer3::FiducialsGUI GetLogic]

set fileNameT1  /raid0/homes/aghayoor/StandAloneBRAINSConstellationDetector-build/BRAINSConstellationDetector-build/TestSuite/TEST_LOCATION/BRAINSConstellationDetectorTest_outputVerificationScript_aligned.nii.gz
set  landmarkName /raid0/homes/aghayoor/StandAloneBRAINSConstellationDetector-build/BRAINSConstellationDetector-build/TestSuite/TEST_LOCATION/BRAINSConstellationDetectorTest_outputVerificationScript_aligned.fcsv

set volumeNode [$volumesLogic AddArchetypeVolume $fileNameT1 T1brain 0]
set selectionNode [$::slicer3::ApplicationLogic GetSelectionNode]
$selectionNode SetReferenceActiveVolumeID [$volumeNode GetID]
$::slicer3::ApplicationLogic PropagateVolumeSelection

$fidLogic LoadFiducialList $landmarkName

