
if [ $# != 4 ]; then
  echo "Incorrect Number of Argument:: $#"  
  echo "Usage::"  
  echo "$0 [ShuffledListFilename] [crossValidationTargetDirectory] [Date] [ROI List Filename]"
  exit 1;
fi

shuffledListFile=$1;
  ## shuffledListFile Ex )
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

GenerateXMLEXE="/hjohnson/HDNI/PREDICT_TRAINING/regina_ann/TrainingModels/ANNCaudate20111107/Script/GenerateXMLConfigurationFileFromList.sh"

if [ ! -d $crossValidationTargetDirectory ]; then
  mkdir $crossValidationTargetDirectory
fi


for testIteration in 1 2 3 4 5 6 7 8 9 10
do

  ##
  ## create list file
  ##

   count=0;
   currentTargetDirectory="${crossValidationTargetDirectory}/Test${testIteration}"

   if [ ! -d ${currentTargetDirectory} ]; then
     echo " Create Directory ( ${currentTargetDirectory} ) "
     mkdir ${currentTargetDirectory}
   fi

   # create list file ----
   currentListFile="${currentTargetDirectory}/${Date}.list"
   rm -f $currentListFile;
   while read line ManDir
   do
      count=`expr $count + 1`;

      if [ $testIteration -le 5 ]; then
        endApplyIndex=`expr $testIteration \* 3`;
        startApplyIndex=`expr $endApplyIndex - 2`;

        dataType="Train";

        if [ $count -ge $startApplyIndex ]; then
          if [ $count -le $endApplyIndex ]; then
            dataType="Apply"
          fi
        fi
      fi
      if [ $testIteration -ge 6 ]; then
        innercount=`expr $testIteration - 5`;
        innercount=`expr $innercount \* 4`;

        endApplyIndex=15;
        endApplyIndex=`expr $endApplyIndex  + $innercount`;
        startApplyIndex=`expr $endApplyIndex - 3`;

        dataType="Train";

        if [ $count -ge $startApplyIndex ]; then
          if [ $count -le $endApplyIndex ]; then
            dataType="Apply"
          fi
        fi
      fi

      echo "$line $dataType $ManDir" >> $currentListFile

   done < ${shuffledListFile}
   # list file ----

   for HN in 4 8 12 16 54
   do
   currentXMLFile="${currentTargetDirectory}/${Date}_$HN.xml"

   echo "from $startApplyIndex to $endApplyIndex"

  ##
  ## create xml file
  ##
  XMLGeneratorCommand="${GenerateXMLEXE} $currentListFile $roiListFilename $currentXMLFile $HN ";
  echo $XMLGeneratorCommand;
  $XMLGeneratorCommand;

  done
  #########
  ## create qsub file
  #########
  
  # get the machine arch 



  QSUBFile="${currentTargetDirectory}/runBRAINSCutSet${testIteration}${Date}.sh"
  echo "QSUBFile name is :: $QSUBFile"


  echo "#!/bin/bash">$QSUBFile
  echo "#$ -N BCut${testIteration}">>$QSUBFile
  echo "#$ -j yes">>$QSUBFile
  echo "#$ -o ${currentTargetDirectory}//runBRAINSCut${testIteration}.sh.log">>$QSUBFile
  echo "#$ -l mf=2G ">>$QSUBFile
  echo "#$ -q 64bit">>$QSUBFile

  echo " arch=\`uname\`;" >>$QSUBFile

  echo "BRAINSSRC=\"/ipldev/scratch/eunyokim/src/BRAINS20111028/build-Darwin-20111028/lib/BRAINSCut\";">>$QSUBFile

  echo "source \${BRAINSSRC}/brains3_setup.sh">>$QSUBFile
  
  for HN in 4 
  do
    currentXMLFile="${currentTargetDirectory}/${Date}_$HN.xml"
    for structure in caudate
    do
      echo " \${BRAINSSRC}/BRAINSCut --netConfiguration  ${currentXMLFile} --applyModel --createVectors --trainModel --generateProbability">>$QSUBFile
    done
  done
  for HN in 8 12 16 54
  do
    currentXMLFile="${currentTargetDirectory}/${Date}_$HN.xml"
    for structure in caudate
    do
      echo " \${BRAINSSRC}/BRAINSCut --netConfiguration  ${currentXMLFile} --trainModel --generateProbability --applyModel">>$QSUBFile
    done
  done

  chmod 755 $QSUBFile

done
