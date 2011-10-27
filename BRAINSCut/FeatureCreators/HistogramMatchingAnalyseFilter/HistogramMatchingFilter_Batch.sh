#!/bin/sh
# 
# 2011 JAN Eun Young(Regina) Kim
#
# Read list of file and produce series of histogram for analaze
#
#

# :: 
# :: Execution Locations:: 
# :: 

bin_path="/scratch/eunyokim/src/BRAINS-Build/src/BRAINSTools/BRAINSCut-build/FeatureCreators/bin/";
src_path="/scratch/eunyokim/src/BRAINS3/BRAINSTools/BRAINSCut/FeatureCreators/HistogramMatchingFilter/"

R_cmd="HistogramMatchingFilter.R"
C_cmd="HistogramMatchingFilter"

# :: 
# :: Read in List
# :: 
# :: 
# ::
# ::
# ::
# ::

if [ $# != 4 ]; then
  echo "Wrong Usage!!! [USAGE::::]"
  echo "$0 [outputDir image_postFix referenceImage listFilename]"
  exit
fi
outputDir=$1;
image_postfix=$2;
referenceImage=$3;
listFilename=$4;

while read line
do
  #
  # Get names
  #


  scanID=(`echo $line |tr "/" "\n"|grep '[0-9]\{5\}'`);
  echo "$scanID";

  inputImageFilename="$line/${scanID}${image_postfix}";
  referenceFilename="$referenceImage";
  histogramOutputFilename="$outputDir/${scanID}_AdjustedHisgoram.nii.gz";
  outputHistogramFilename="$outputDir/${scanID}_AdjustedHisgoram.nii.gz.dat";
  histogramFigureFilename="$outputDir/${scanID}_Histogram.png";

  

  # 
  # Generate Adjusted Histogram Image
  #

  generateHistCMD="${bin_path}/${C_cmd} \
                    --inputImageFileanme ${inputImageFilename} \
                    --referenceFilename ${referenceFilename} 
                    --outputFilename ${histogramOutputFilename} 
                    --outputHistogramFilename ${outputHistogramFilename} ";
  echo "#----------------------------------------------------------------------#"
  echo ${generateHistCMD}
  echo "#----------------------------------------------------------------------#"
  #${generateHistCMD}

done < $listFilename
while read line
do
  #
  # Get names
  #


  scanID=(`echo $line |tr "/" "\n"|grep '[0-9]\{5\}'`);
  echo "$scanID";

  inputImageFilename="$line/${scanID}${image_postfix}";
  referenceFilename="$referenceImage";
  histogramOutputFilename="$outputDir/${scanID}_AdjustedHisgoram.nii.gz";
  outputHistogramFilename="$outputDir/${scanID}_AdjustedHisgoram.nii.gz.dat";
  histogramFigureFilename="$outputDir/${scanID}_Histogram.png";

  
  # 
  # Generate Histogram Figure
  #

  drawHistogram="R CMD --vanilla --slave  --args \
                 ${outputHistogramFilename} \
                 ${histogramFigureFilename} \
                 < ${src_path}/${R_cmd} "

  echo "#----------------------------------------------------------------------#"
  #echo $drawHistogram >> ${outputHistogramFilename}_BATCH.sh
  echo $drawHistogram 
  echo "#----------------------------------------------------------------------#"
  $drawHistogram 

done < $listFilename
