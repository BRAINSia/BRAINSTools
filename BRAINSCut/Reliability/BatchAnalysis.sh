# author: Eunyoung Regina Kim
# 2011 Nov
# ---------------------------------------------------------------------------- #
# from the list file 
# compute 'manual' and 'ann' volume and 'similarity index' along
# the different threshold with BRAINSCut Post processing step
# ---------------------------------------------------------------------------- #
inputListFilename=$1;
inputANNDir=$2;
inputHN=$3; # since output directory of ANN consists of HN
ROIName=$4

if [ $# != 4 ]; then
  echo "USAGE::: $0 [list filename] [ann dir] [HN] [roi name]"
  exit 1
fi

#
# Similarity Index execution 
SIExe="/ipldev/scratch/eunyokim/src/BRAINS20111028/build-Darwin-20111028/lib/SimilarityIndex";

#
# input
# 1. list of autoAutoWorkUp directory paired with ANN Manual directory
# 2. ANN Trail directory

#
# list file ex)
# [AutoWorkUp Dir] [Manual Dir] 
# ex)
# /paulsen/IPIG/predict_3T_MR/site-024/0706/43258/10_AUTO.NN3Tv20110418  ANN2011May04ManualCompleted
# /paulsen/MRx/PHD_120/0057/34479/10_AUTO.NN3Tv20111003 ANN20111006_ManualTrimForANN 
# /paulsen/IPIG/predict_3T_MR/site-048/0217/52712/10_AUTO.NN3Tv20110418  ANN2011May04ManualCompleted


MendatoryFileExists() #----------------------------------------------------------- #
{
  local f=$1;
  if [[ -f $f ]]; then
    return 0;
  else
    echo "File $f does not exists! "
    exit 1;
  fi
}

FileExists() #------------------------------------------------------------------- #
{
  local f=$1;
  if [[ -f $f ]]; then
    return 0;
  else
    return 1;
  fi
}

OutputDir="$inputANNDir/AnalysisOutput_${ROIName}_${inputHN}";
if ( FileExists $OutputDir ); then
  echo "file $OutputDir found"
else
    mkdir $OutputDir;
fi

RScriptListFile="${OutputDir}/RScriptList.csv"
echo "subjectID,csvFIlename,ROIName">$RScriptListFile;

while read autoWorkUpDir manDir
do
   #
   # scan information
   subjectID=( `echo $autoWorkUpDir |tr "/" "\n"|grep '[0-9]\{5\}'`)

   #
   # Compute Similarity index -------------------------------------------------- #
   SIManualROI=(`ls $autoWorkUpDir/$manDir/${subjectID}_l_${ROIName}.nii.gz`) ; 
   MendatoryFileExists $SIManualROI;
   SIAnnROI=(`ls $inputANNDir/Test*/${subjectID}_${inputHN}/ANNContinuousPredictionl_${ROIName}${subjectID}.nii.gz`)
   MendatoryFileExists $SIAnnROI

   SIThresholdInterval="0.01"
   CSVOutputFileOfCurrentScan="${OutputDir}/l_${ROIName}_${subjectID}.csv"
   
   SICommand="$SIExe --inputManualVolume   $SIManualROI \
                     --ANNContinuousVolume $SIAnnROI    \
                     --thresholdInterval   $SIThresholdInterval";

   if [[ -s $CSVOutputFileOfCurrentScan ]] ; then # if file is not empty
     echo "-------------------------------------------------------------------------"
     echo "$CSVOutputFileOfCurrentScan has data already. skip this "
     echo "-------------------------------------------------------------------------"
   else
     echo "-------------------------------------------------------------------------"
     echo $SICommand;
     $SICommand |grep ","|tee $CSVOutputFileOfCurrentScan;
     #
     # Re-format the similarity text output ----------------------------------- #
     while read line
     do
       echo "$subjectID, $line " >> ${CSVOutputFileOfCurrentScan}TEMP;
     done < $CSVOutputFileOfCurrentScan
   
     mv ${CSVOutputFileOfCurrentScan}TEMP $CSVOutputFileOfCurrentScan
   fi;

   #------------------------------------------------------------------------- #
   # make list file for RSscript
   echo "$subjectID,$CSVOutputFileOfCurrentScan,$ROIName " >> $RScriptListFile;
   echo "-------------------------------------------------------------------------"

done < $inputListFilename

#
# Relative Overlap Graph ------------------------------------------------------ #
SIList=$RScriptListFile;
SIPlotOutputPlotFilename="${RScriptListFile}.pdf";
SIPlotRScript=" /ipldev/scratch/eunyokim/src/BRAINS20111028/BRAINSStandAlone/BRAINSCut/Reliability/RelativeOverlapPlot.R"
   
SIPlotR="bash R --slave --args $SIList $SIPlotOutputPlotFilename < $SIPlotRScript "

echo $SIPlotR
$SIPlotR
   












