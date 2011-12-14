#!/bin/bash

##
## set the BRAINSCut SRC Dir.
##
ARCH=`uname`;
if [ "$ARCH" == "Darwin" ]; then
  BRAINSBuild="/ipldev/scratch/eunyokim/src/BRAINS20111028/build-Darwin-20111028";
else
  BRAINSBuild="/scratch/PREDICT/regina/BRAINS/buildICC/";
fi

##
## Set utility script
##
utilitySRC="$BRAINSBuild/../BRAINSStandAlone/BRAINSCut/Reliability/utilities.sh"
source $utilitySRC



##
## check input arguments
##
if [ $# != 4 ]; then
  echo "Incorrect Number of Argument:: $#"  
  echo "Usage:::::"  
  echo "::::::::::"  
  echo "$0 [ShuffledListFilename] [crossValidationTargetDirectory] [Date] [ROI List Filename] "
  echo "::::::::::"  
  exit 1;
fi
## 

pseudoRandomDataList=$1;
  ## pseudoRandomDataList Ex )
  ##-----------------------
  ## /paulsen/IPIG/predict_3T_MR/site-027/0839/84201/10_AUTO.NN3Tv20110418 ManDir
  ## /paulsen/IPIG/predict_3T_MR/site-144/0619/96294/10_AUTO.NN3Tv20110418 ManDir
  ## /paulsen/IPIG/predict_3T_MR/site-048/0647/45592/10_AUTO.NN3Tv20110418 ManDir
  ##-----------------------
crossValidationTargetDirectory=$2;
Date=$3;
roiListFilename=$4;
  ## roi ListFile Ex )
  ##-----------------------
  ## l_caudate true
  ## r_caudate true
  ##-----------------------

##------------------------------------------------------------------------- ##
for testIteration in 1 2 3 4 5 6 7 8 9 10
do
   currentTargetDirectory="${crossValidationTargetDirectory}/Test${testIteration}"

   ##
   ## create qsub file
   ##
   
   for HN in 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80
   do
     QSUBFile="${currentTargetDirectory}/runApply${testIteration}${Date}${HN}.sh"
     echo "QSUBFile name is :: $QSUBFile"
     qsubHeader $QSUBFile
     currentXMLFile="${currentTargetDirectory}/${Date}_$HN.xml"
     echo " \${BRAINSBuild}/BRAINSCut --netConfiguration  ${currentXMLFile} --trainModel ">>$QSUBFile
   done

   chmod 755 $QSUBFile

done
