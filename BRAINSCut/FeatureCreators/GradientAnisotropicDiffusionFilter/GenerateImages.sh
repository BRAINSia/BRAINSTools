
if [ $# != 1 ]; then
  echo "Wrong Usage!!! [USAGE::::]"
  echo "$0 [List File Name]"
  exit
fi

inputVolume="_AVG_T2.nii.gz"
outputVolume="_GradAniosoDiff_AVG_T2.nii.gz"
ListFile=$1;
source_dir="/scratch/eunyokim/src/BRAINS-Build/src/bin"
GenerateEXE="$source_dir/../BRAINSTools/BRAINSCut-build/FeatureCreators/bin/GradientAnisotropicDiffusionImageFilter"

## specify source directory for execution
echo "source $source_dir/brains3_setup.sh"
source $source_dir/brains3_setup.sh

while read line
do
  subject=(`echo $line |tr "/" "\n"|grep '[0-9]\{5\}'`)
  command="$GenerateEXE \
        --inputVolume $line/${subject}${inputVolume} \
        --outputVolume $line/${subject}${outputVolume} \
        --timeStep 0.06
        --numberOfIterations 5 \
        --conductance 1"
  echo "======================================================================="
  echo $command
  $command
  echo "======================================================================="
      
done < ${ListFile}
