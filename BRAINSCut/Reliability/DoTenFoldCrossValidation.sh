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
if [ $# != 5 ]; then
  echo "Incorrect Number of Argument:: $#"  
  echo "Usage:::::"  
  echo "::::::::::"  
  echo "$0 [ShuffledListFilename] [crossValidationTargetDirectory] [Date] [ROI List Filename] [XML Script]"
  echo "::::::::::"  
  exit 1;
fi
## 
GenerateXMLEXE=$5;

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


# This function will do the 10 fold cross validation for the given list of file

if [ ! -d $crossValidationTargetDirectory ]; then
  mkdir $crossValidationTargetDirectory
fi

##------------------------------------------------------------------------- ##
function generateListOfTrainAndApply( ) {

   tN=$1;
   dataList=$2;
   outputListFilename=$3;

   curreentLineOfDataList=0;

   while read line ManDir
   do
      curreentLineOfDataList=`expr $curreentLineOfDataList + 1`;

      if [ $tN -le 5 ]; then
        endApplyIndex=`expr $tN \* 3`;
        startApplyIndex=`expr $endApplyIndex - 2`;

        dataType="Train";

        if [ $curreentLineOfDataList -ge $startApplyIndex ]; then
          if [ $curreentLineOfDataList -le $endApplyIndex ]; then
            dataType="Apply"
          fi
        fi
      fi
      if [ $tN -ge 6 ]; then
        innercurreentLineOfDataList=`expr $tN - 5`;
        innercurreentLineOfDataList=`expr $innercurreentLineOfDataList \* 4`;

        endApplyIndex=15;
        endApplyIndex=`expr $endApplyIndex  + $innercurreentLineOfDataList`;
        startApplyIndex=`expr $endApplyIndex - 3`;

        dataType="Train";

        if [ $curreentLineOfDataList -ge $startApplyIndex ]; then
          if [ $curreentLineOfDataList -le $endApplyIndex ]; then
            dataType="Apply"
          fi
        fi
      fi

      echo "$line $dataType $ManDir" >> $outputListFilename
   done < ${dataList}
}
##------------------------------------------------------------------------- ##
for testIteration in 1 2 3 4 5 6 7 8 9 10
do

   currentTargetDirectory="${crossValidationTargetDirectory}/Test${testIteration}"
   mkdir -p ${currentTargetDirectory}

   ##
   ## list file for xml
   ##
   currentListFile="${currentTargetDirectory}/${Date}.list"
   rm -f $currentListFile;
   generateListOfTrainAndApply $testIteration $pseudoRandomDataList $currentListFile

   for HN in 5 10 15 20 30 40 50 60 70 80
   do
      currentXMLFile="${currentTargetDirectory}/${Date}_$HN.xml"
      echo "from $startApplyIndex to $endApplyIndex"

      ##
      ## create xml file
      ##
      XMLGeneratorCommand="${GenerateXMLEXE} $currentListFile $roiListFilename $currentXMLFile $HN $BRAINSBuild";
      printCommandAndRun "$XMLGeneratorCommand";
   done
   
   ##
   ## create qsub file
   ##
   
   # get the machine arch 
   QSUBFile="${currentTargetDirectory}/runBRAINSCutSet${testIteration}${Date}.sh"
   echo "QSUBFile name is :: $QSUBFile"

   qsubHeader $QSUBFile

   for HN in 5
   do
     currentXMLFile="${currentTargetDirectory}/${Date}_$HN.xml"
     echo " \${BRAINSBuild}/BRAINSCut --netConfiguration  ${currentXMLFile} --createVectors --trainModel --generateProbability">>$QSUBFile
   done
   for HN in 10 15 20 30 40 50 60 70 80
   do
     currentXMLFile="${currentTargetDirectory}/${Date}_$HN.xml"
     echo " \${BRAINSBuild}/BRAINSCut --netConfiguration  ${currentXMLFile} --trainModel ">>$QSUBFile
   done

   chmod 755 $QSUBFile

done
