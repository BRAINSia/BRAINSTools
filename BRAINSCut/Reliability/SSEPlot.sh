
ANNModelDirPrefix=$1
subDir=$2
OutputDir=$3
if [ $# != 3 ]; then
  echo "Incorrect Number of Argument:: $#"  
  echo "Usage::----------------------------------------------------"  
  echo "$0 [ANN Model Dir Prefix] [subDir] [Output Dir]"
  echo ":----------------------------------------------------------"  
  exit 1;
fi

RScript="/ipldev/scratch/eunyokim/src/BRAINS20111028/BRAINSStandAlone/BRAINSCut/Reliability/SSEPlot.R"

SSECollectionFilename="$OutputDir/ValidationSetSSE.txt"
rm -f $SSECollectionFilename
for file in $ANNModelDirPrefix*/${subDir}/*ValidationSetSSE.txt
do
  cat $file >> $SSECollectionFilename;
done

bash R --slave --args $SSECollectionFilename $OutputDir/${subDir} < $RScript 
