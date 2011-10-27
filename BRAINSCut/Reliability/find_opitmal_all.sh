program="/scratch/eunyokim/src/BRAINS-Build/src/BRAINSTools/BRAINSCut/Reliability/FindOptimalThreshold.tcl"
brains2_executable="brains2"

if [ $# != 3 ]; then
  echo "USAGE::: $0 [ANNOutputDir] [Structure] [Data List]"
  exit 1
fi
ANNOutputDir=$1;
ROI=$2;
DATALIST=$3;

ManualSubDirectory="ANN2011Apr25ManualCompleted";
#ManualSubDirectory="AmygdalaMayManualCompleted";

for SubjectANNOutputDir in $ANNOutputDir/*
do
   if [ -d "${SubjectANNOutputDir}" ]; then

      subject=( `echo $SubjectANNOutputDir |tr "/" "\n"|grep '[0-9]\{5\}'`)
      #
      # read Data List and get the data location
      #
      dataDir=( `cat $DATALIST |grep $subject` );
      dataDir="${dataDir}/${ManualSubDirectory}/"


      echo " -------------------------------------------------------------------------"
      echo "$SubjectANNOutputDir with subject of $subject"
      echo " -------------------------------------------------------------------------"
      for s in "l_${ROI}" "r_${ROI}"
      do
        ANNFilename=( `ls ${SubjectANNOutputDir}/*${s}${subject}.nii.gz` );
        ManualFilename=( `ls ${dataDir}/*${s}*.nii.gz` );
        command="\
        ${brains2_executable}  -b \
                       ${program} ${ANNFilename}  $ManualFilename \
                       |grep :|tee $ANNOutputDir/${s}RO_Vs_Threshold${subject}.txt";
        echo "$command";
        ${brains2_executable}  -b \
                               ${program} ${ANNFilename}  $ManualFilename \
                               |grep :|tee $ANNOutputDir/${s}RO_Vs_Threshold${subject}.txt
      done
    fi 
done
