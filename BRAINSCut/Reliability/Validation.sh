##
## Second part of 10 fold cross validation
## compute SSE to find iteration number to use


WorkingDir=$1;
if [ $# != 4 ]; then
  echo "Incorrect Number of Argument"  
  echo "Usage::"  
  echo "$0 [working directory] [GenerateValidationXML Script] [MaskListFile ][HNList]"
  exit 1;
fi




maskList=$3;
HNList=$4;
# HNList=" 4 11 22 33 44 55 66"


## TODO change this to BRAINSCut build Dir

machine=`uname`
if [ "$machine" == "Darwin" ]; then
  BRAINSSRC="/ipldev/scratch/eunyokim/src/BRAINS20111028/build-Darwin-20111109/lib/"
  ## Mask list

else
  BRAINSSRC="/scratch/PREDICT/regina/BRAINS//buildICC/lib/"
  ## Mask list
fi

## Sum of squar error script
SSEScript="$BRAINSSRC/../../BRAINSStandAlone/BRAINSCut/Reliability/SSEPlot.sh";
UtilityrSRC="$BRAINSSRC/../../BRAINSStandAlone/BRAINSCut/Reliability/utilities.sh";

source $UtilityrSRC 

## Generate Validation XML file script
XMLFileEXE=$2;
for HN in $HNList
do

  for TEST in 1 2 3 4 5 6 7 8 9 10
  do
     LISTFile=`ls $WorkingDir/Test$TEST/*.list`
     echo "LISTFILE:::::$LISTFile"
     DATE=`basename $LISTFile`;
     DATE=`echo $DATE|sed 's/\.list$//'`;
     outputXML="$WorkingDir/Test$TEST/${DATE}_${HN}_Validation.xml"

     # generate validation xml file
     if [ ! -s $outputXML ]; then
       printCommandAndRun "$XMLFileEXE $LISTFile $maskList $outputXML $HN"
     fi
   done
done

## 1. Create SSE values
for TEST in 1 2 3 4 5 6 7 8 9 10
do
        for HN in $HNList
        do
          ## ---------------- QSUB -----------------------
          QSUBFile="${WorkingDir}/Test$TEST/runValidation${TEST}_HN${HN}.sh"
          echo "QSUBFile name is :: $QSUBFile"
          qsubHeader $QSUBFile      
          currentXMLFile="${WorkingDir}/Test$TEST/${DATE}_${HN}_Validation.xml"
          echo " \${BRAINSBuild}/BRAINSCut --netConfiguration  ${currentXMLFile} --computeSSEOn  --applyModel">>$QSUBFile
          chmod 755 $QSUBFile
        done
      

      
      ## ---------------- QSUB -----------------------
done

for HN in $HNList
do
   ## Analyze output 
   SSEOutputDir="$WorkingDir/SSE$HN/"
   mkdir -p $SSEOutputDir;

   SSEOutput="$SSEOutputDir/CallSSEPlotAlongTrainedModel.sh";
   rm -f $SSEOutput;
   echo "bash $SSEScript $WorkingDir/Test TrainedModels$HN  $SSEOutputDir" >> $SSEOutput;
   echo "bash $SSEScript $WorkingDir/Test TrainedModels$HN  $SSEOutputDir" 
   echo "$SSEOutput"

   chmod 755 $SSEOutput
done
