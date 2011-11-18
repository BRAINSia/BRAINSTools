
ANNModelDirPrefix=$1
subDir=$2
OutputDir=$3
if [ $# != 3 ]; then
  echo "Incorrect Number of Argument:: $#"  
  echo "Usage::----------------------------------------------------"  
  echo "$0 [ANN Model Dir Prefix] [subDir] [Output Dir]"
  echo " ex)  $0 /scratch/PREDICT/regina/ANN20111112T1T2/Test TrainedModels4  /scratch/PREDICT/regina/ANN20111112T1T2/SSE04"
  echo ":----------------------------------------------------------"  
  exit 1;
fi

machine=`uname`
if [ "$machine" == "Darwin" ]; then
  BRAINSSRC="/ipldev/scratch/eunyokim/src/BRAINS20111028/"
else
  BRAINSSRC="/scratch/PREDICT/regina/BRAINS/"
fi

RScript="$BRAINSSRC/BRAINSStandAlone/BRAINSCut/Reliability/SSEPlot.R"


mkdir -p $OutputDir
SSECollectionFilename="$OutputDir/ValidationSetSSE.txt"
rm -f $SSECollectionFilename
for file in $ANNModelDirPrefix*/${subDir}/*ValidationSetSSE.txt
do
  echo "cat this file :: $file"
  cat $file >> $SSECollectionFilename;
done

bash R --slave --args $SSECollectionFilename $OutputDir/${subDir} < $RScript 
