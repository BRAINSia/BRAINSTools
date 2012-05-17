#!/bin/bash

##
## set the BRAINSCut SRC Dir.
##
ARCH=`uname`;
if [ "$ARCH" == "Darwin" ]; then
  BRAINSBuild="/ipldev/scratch/eunyokim/src/BRAINS20111028/build-Darwin/";
else
  BRAINSBuild="/scratch/PREDICT/regina/BRAINS/buildICC-May/";
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
  echo "$0 [ShuffledListFilename] [crossValidationTargetDirectory] [Date] [ROI List Filename] [XML Script] "
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

   for HN in 20 60
   do
      currentXMLFile="${currentTargetDirectory}/${Date}_ANN$HN.xml"
      echo "from $startApplyIndex to $endApplyIndex"

      ##
      ## create xml file
      ##
      XMLGeneratorCommand="${GenerateXMLEXE} $currentListFile $roiListFilename $currentXMLFile $HN $BRAINSBuild ANN$HN";
      printCommandAndRun "$XMLGeneratorCommand";
   done


   for RF in 10 25 100
   do
      currentXMLFile="${currentTargetDirectory}/${Date}_RF$RF.xml"
      echo "from $startApplyIndex to $endApplyIndex"

      ##
      ## create xml file
      ##
      XMLGeneratorCommand="${GenerateXMLEXE} $currentListFile $roiListFilename $currentXMLFile $HN $BRAINSBuild RF$RF";
      printCommandAndRun "$XMLGeneratorCommand";
   done
   
   ##
   ## create qsub file
   ##
   
   # get the machine arch 
   QSUBFile="${currentTargetDirectory}/runBRAINSCutSet${testIteration}${Date}.sh"
   echo "QSUBFile name is :: $QSUBFile"

   qsubHeader $QSUBFile

   for HN in 20
   do
     currentXMLFile="${currentTargetDirectory}/${Date}_ANN$HN.xml"
     echo " \${BRAINSBuild}/BRAINSCut --netConfiguration  ${currentXMLFile} --createVectors --generateProbability --trainModel --applyModel">>$QSUBFile
   done

   for HN in  60 
   do
     currentXMLFile="${currentTargetDirectory}/${Date}_ANN$HN.xml"
     echo " \${BRAINSBuild}/BRAINSCut --netConfiguration  ${currentXMLFile} --applyModel --trainModel ">>$QSUBFile
   done

   for RF in  10 25 100
   do
     currentXMLFile="${currentTargetDirectory}/${Date}_RF$RF.xml"
     echo " \${BRAINSBuild}/BRAINSCut --netConfiguration  ${currentXMLFile} --applyModel --trainModel --method RandomForest --numberOfTrees $RF --randomTreeDepth 100 --NoTrainingVectorShuffling">>$QSUBFile
   done
   chmod 755 $QSUBFile

done
