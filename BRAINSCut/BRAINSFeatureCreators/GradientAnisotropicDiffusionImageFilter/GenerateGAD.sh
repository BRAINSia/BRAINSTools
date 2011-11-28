
if [ $# != 2 ]; then
  echo "Wrong Usage!!! [USAGE::::]"
  echo "$0 [List File Name] [BRAINS build]"
  exit
fi

ListFile=$1;
source_dir=$2;
GADEXE="$source_dir/GradientAnisotropicDiffusionImageFilter"

## specify source directory for execution
while read line dummy
do
  for mod in T1 T2
  do
    subject=(`echo $line |tr "/" "\n"|grep '[0-9]\{5\}'`)
    command="$GADEXE \
        --inputVolume $line/${subject}_AVG_${mod}.nii.gz \
        --outputVolume $line/${subject}_AVG_${mod}_GAD.nii.gz \
        --timeStep 0.05 \
        --numberOfIterations 5 \
        --conductance 1"
    echo "======================================================================="
    echo $command
    $command
    echo "======================================================================="
  done
      
done < ${ListFile}
