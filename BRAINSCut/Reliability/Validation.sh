##
## Second part of 10 fold cross validation
## compute SSE to find iteration number to use

WorkingDir=$1;
if [ $# != 2 ]; then
  echo "Incorrect Number of Argument"  
  echo "Usage::"  
  echo "$0 [working directory] [GenerateValidationXML Script]"
  exit 1;
fi

DATE="20111112";


## TODO change this to BRAINSCut build Dir

machine=`uname`
if [ "$machine" == "Darwin" ]; then
  BRAINSSRC="/ipldev/scratch/eunyokim/src/BRAINS20111028/build-Darwin-20111109/lib/"
  ## Mask list
  maskList="/hjohnson/HDNI/PREDICT_TRAINING/regina_ann/TrainingModels/Script/Putamen.list"
else
  BRAINSSRC="/scratch/PREDICT/regina/BRAINS//buildICC/lib/"
  ## Mask list
  maskList="/scratch/PREDICT/regina/Scirpts/Putamen.list"
fi

## Sum of squar error script
SSEScript="$BRAINSSRC/../../BRAINSStandAlone/BRAINSCut/Reliability/SSEPlot.sh";

## Generate Validation XML file script
XMLFileEXE=$2;



for HN in 4 11 22 33 44 55 66
do
  ## Analyze output 
  SSEOutputDir="$WorkingDir/AnalysisOutput_putamen_$HN/"
  mkdir -p $SSEOutputDir;

  SSEOutput="$SSEOutputDir/CallSSEPlot.sh";
  rm -f $SSEOutput;
  for TEST in 1 2 3 4 5 6 7 8 9 10
  do
     LISTFile="$WorkingDir/Test$TEST/${DATE}.list"
     outputXML="$WorkingDir/Test$TEST/${DATE}_${HN}_Validation.xml"

     # generate validation xml file
     if [ ! -s $outputXML ]; then
       cmd="$XMLFileEXE $LISTFile $maskList $outputXML $HN"
       echo $cmd
       $cmd
     fi
   done

   echo "bash $SSEScript $WorkingDir/Test TrainedModels$HN  $SSEOutputDir" >> $SSEOutput;
   chmod 755 $SSEOutput
done

for TEST in 1 2 3 4 5 6 7 8 9 10
do
      ## ---------------- QSUB -----------------------
        QSUBFile="${WorkingDir}/Test$TEST/runValidation${TEST}.sh"
        echo "QSUBFile name is :: $QSUBFile"
      
      
        echo "#!/bin/bash">$QSUBFile
        echo "#$ -N VAL${testIteration}">>$QSUBFile
        echo "#$ -j yes"         >>$QSUBFile
        echo "#$ -o $QSUBFile.log" >>$QSUBFile
        echo "#$ -l mf=2G "      >>$QSUBFile
        echo "#$ -pe smp1 1-12"  >>$QSUBFile
        echo "PLATFORM=\$(uname)">>$QSUBFile
        echo "hostname"          >>$QSUBFile
        echo "uname -a"          >>$QSUBFile
        echo "## Set global number of threads to 4 in order to minimize CPU wasting.">>$QSUBFile
        echo "export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=\$NSLOTS;">>$QSUBFile
        echo "echo \"USING NUM THREADS \${ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS}\"" >>$QSUBFile
      
        echo " arch=\`uname\`;"   >>$QSUBFile
      
        echo "BRAINSSRC=\"$BRAINSSRC\";">>$QSUBFile
      
        for HN in 4 11 22 33 44 55 66
        do
          currentXMLFile="${WorkingDir}/Test$TEST/${DATE}_${HN}_Validation.xml"
          for structure in caudate
          do
            echo " \${BRAINSSRC}/BRAINSCut --netConfiguration  ${currentXMLFile} --computeSSEOn  --applyModel">>$QSUBFile
          done

        done
      
        chmod 755 $QSUBFile
      
      ## ---------------- QSUB -----------------------
done

