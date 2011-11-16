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
BRAINSBuild=$5;

if [ $# != 5 ]; then
  echo "USAGE::: $0 [list filename] [ann dir] [HN] [roi name] [BRAINS build]"
  exit 1
fi

#
# Similarity Index execution 
SIExe="$BRAINSBuild/SimilarityIndex";

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
mkdir -p $OutputDir;

for side in l r
do
  RScriptListFile="${OutputDir}/${side}_${ROIName}_RScriptList.csv"
  echo "subjectID,csvFIlename,ROIName">$RScriptListFile;
done

while read autoWorkUpDir manDir
do
  for side in l r
  do
    RScriptListFile="${OutputDir}/${side}_${ROIName}_RScriptList.csv"
   #
   # scan information
   subjectID=( `echo $autoWorkUpDir |tr "/" "\n"|grep '[0-9]\{5\}'`)

   #
   # Compute Similarity index -------------------------------------------------- #
   SIManualROI=(`ls $autoWorkUpDir/$manDir/${subjectID}_${side}_${ROIName}.nii.gz`) ; 
   MendatoryFileExists $SIManualROI;
   SIAnnROI=(`ls $inputANNDir/Test*/${subjectID}_${inputHN}/ANNContinuousPrediction${side}_${ROIName}${subjectID}.nii.gz`)
   MendatoryFileExists $SIAnnROI

   SIThresholdInterval="0.01"
   CSVOutputFileOfCurrentScan="${OutputDir}/${side}_${ROIName}_${subjectID}.csv"
   
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
      QSUBFile="${OutputDir}/ComputeSimilarities${subjectID}_${side}_${ROIName}.sh"
      echo "QSUBFile name is :: $QSUBFile"


      echo "#!/bin/bash">$QSUBFile
      echo "#$ -N ${side}_${ROIName}_${subjectID}">>$QSUBFile
      echo "#$ -j yes"         >>$QSUBFile
      echo "#$ -o $QSUBFile.log ">>$QSUBFile
      echo "#$ -l mf=2G "      >>$QSUBFile
      echo "#$ -pe smp1 1-2"  >>$QSUBFile
      echo "PLATFORM=\$(uname)">>$QSUBFile
      echo "hostname"          >>$QSUBFile
      echo "uname -a"          >>$QSUBFile
      echo "## Set global number of threads to 4 in order to minimize CPU wasting.">>$QSUBFile
      echo "export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=\$NSLOTS;"                 >>$QSUBFile
      echo "echo \"USING NUM THREADS \${ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS}\""   >>$QSUBFile

      echo " arch=\`uname\`;"                                     >>$QSUBFile
      echo "$SICommand |grep ","|tee $CSVOutputFileOfCurrentScan" >>$QSUBFile;
      echo "#"                                                    >>$QSUBFile
      echo "# Re-format the similarity text output ----------------------------------- #" >>$QSUBFile;
      echo "while read line                                      ">>$QSUBFile
      echo "do                                                   ">>$QSUBFile
      echo "  echo \"$subjectID, \$line \" >> ${CSVOutputFileOfCurrentScan}TEMP;">>$QSUBFile
      echo "done < $CSVOutputFileOfCurrentScan                   ">>$QSUBFile
   
      echo "mv ${CSVOutputFileOfCurrentScan}TEMP $CSVOutputFileOfCurrentScan">>$QSUBFile
      chmod 755 $QSUBFile;
   fi;



   #------------------------------------------------------------------------- #
   # make list file for RSscript
   echo "$subjectID,$CSVOutputFileOfCurrentScan,$ROIName " >> $RScriptListFile;
   echo "-------------------------------------------------------------------------"

 done
done < $inputListFilename

RScript="${OutputDir}/ComputeStats.sh"
rm -f $RScript

for side in l r
do
  RScriptListFile="${OutputDir}/${side}_${ROIName}_RScriptList.csv"
  #
  # Relative Overlap Graph ------------------------------------------------------ #
  SIList=$RScriptListFile;
  SIPlotOutputPlotFilename="${RScriptListFile}.pdf";
  SIPlotRScript="$BRAINSBuild/../../BRAINSStandAlone/BRAINSCut/Reliability/RelativeOverlapPlot.R"
   
  SIPlotR="R --slave --args $SIList $SIPlotOutputPlotFilename < $SIPlotRScript "
  echo "$SIPlotR" >> $RScript ;
done














