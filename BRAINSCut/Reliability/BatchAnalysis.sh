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

OutputDir="$inputANNDir/AnalysisOutput_${ROIName}_${inputHN}";
rm -rf $OutputDir;
mkdir $OutputDir;

RScriptListFile="${OutputDir}/RScriptList.csv"
echo "subjectID,csvFIlename,ROIName">$RScriptListFile;

while read autoWorkUpDir manDir
do
   #
   # scan information
   subjectID=( `echo $autoWorkUpDir |tr "/" "\n"|grep '[0-9]\{5\}'`)

   #
   # Compute Similarity index -------------------------------------------------- #
   SIManualROI=(`ls $autoWorkUpDir/$manDir/${subjectID}_l_${ROIName}.nii.gz`)
   SIAnnROI=(`ls $inputANNDir/Test*/${subjectID}_${inputHN}/ANNContinuousPredictionl_${ROIName}${subjectID}.nii.gz`)
   SIThresholdInterval="0.01"
   SICSVOutput="${OutputDir}/l_${ROIName}_${subjectID}.csv"
   
   SICommand="$SIExe --inputManualVolume   $SIManualROI \
                     --ANNContinuousVolume $SIAnnROI    \
                     --thresholdInterval   $SIThresholdInterval";

   echo $SICommand;
   $SICommand |grep ","|tee $SICSVOutput;
   #
   # Re-format the similarity text output --------------------------------------- #
   while read line
   do
     echo "$subjectID, $line " >> ${SICSVOutput}TEMP;
   done < $SICSVOutput
   
   mv ${SICSVOutput}TEMP $SICSVOutput

   #---------------------------------------------------------------------------- #
   # make list file for RSscript
   echo "$subjectID,$SICSVOutput,$ROIName " >> $RScriptListFile;
done < $inputListFilename

#
# Relative Overlap Graph ------------------------------------------------------- #
SIList=$RScriptListFile;
SIPlotOutputPlotFilename="${RScriptListFile}.pdf";
SIPlotRScript=" /ipldev/scratch/eunyokim/src/BRAINS20111028/BRAINSStandAlone/BRAINSCut/Reliability/RelativeOverlapPlot.R"
   
SIPlotR="bash R --slave --args $SIList $SIPlotOutputPlotFilename < $SIPlotRScript "

echo $SIPlotR
   












