referenceDir="/Users/eunyokim//src/build/ReferenceAtlas-build/Atlas/Atlas_20120104/"
subjectDir="/Users/eunyokim//TestImages/_uid_PHD_DTI_THP_THP0002_THP0002_UW1/"
outputDir="/Users/eunyokim//TEST/"
probabilityMapDir="/Users/eunyokim//TestImages/ProbabilityMaps/"



BCutCMDDir="/Users/eunyokim/src/BRAINSStandAlone/AutoWorkup/BRAINSTools/"

if [ $# == 4 ]; then
  referenceDir=$1;
  subjectDir=$2;
  outputDir=$3;
  probabilityMapDir=$4;
fi

mkdir -p $outputDir

python $BCutCMDDir/BRAINSCutCMD.py \
--inputSubjectT1Filename $subjectDir/11_BABC/t1_average_BRAINSABC.nii.gz \
--inputSubjectT2Filename $subjectDir/11_BABC/t2_average_BRAINSABC.nii.gz \
--inputSubjectSGFilename /hjohnson/HDNI/PREDICT_TRAINING/regina_ann/TrainingModels/FeatureImages/GadSG/96294_average_GAD_SG.nii.gz \
--inputSubjectBrainMaskFilename $subjectDir/11_BABC/brain_label_seg.nii.gz \
--deformationFromTemplateToSubject $subjectDir/11_BABC/atlas_to_subject.mat \
--inputTemplateT1 $referenceDir/template_t1.nii.gz \
--inputTemplateBrainMask $referenceDir/template_brain.nii.gz \
--inputTemplateRhoFilename $referenceDir/spatialImages/rho.nii.gz \
--inputTemplatePhiFilename $referenceDir/spatialImages/phi.nii.gz \
--inputTemplateThetaFilename $referenceDir/spatialImages/theta.nii.gz \
--trainingVectorFilename  dummyVector.txt \
--modelFilename /Users/eunyokim//TestImages/TrainedModels/RandomForestAllSubcorticalsBalancedModel.txtD0060NT0060  \
--probabilityMapsLeftAccumben     $probabilityMapDir/l_accumben_ProbabilityMap.nii.gz \
--probabilityMapsRightAccumben    $probabilityMapDir/r_accumben_ProbabilityMap.nii.gz \
--probabilityMapsLeftCaudate      $probabilityMapDir/l_caudate_ProbabilityMap.nii.gz \
--probabilityMapsRightCaudate     $probabilityMapDir/r_caudate_ProbabilityMap.nii.gz \
--probabilityMapsLeftPutamen      $probabilityMapDir/l_putamen_ProbabilityMap.nii.gz \
--probabilityMapsRightPutamen     $probabilityMapDir/r_putamen_ProbabilityMap.nii.gz \
--probabilityMapsLeftGlobus       $probabilityMapDir/l_globus_ProbabilityMap.nii.gz \
--probabilityMapsRightGlobus      $probabilityMapDir/r_globus_ProbabilityMap.nii.gz \
--probabilityMapsLeftThalamus     $probabilityMapDir/l_thalamus_ProbabilityMap.nii.gz \
--probabilityMapsRightThalamus    $probabilityMapDir/r_thalamus_ProbabilityMap.nii.gz \
--probabilityMapsLeftHippocampus  $probabilityMapDir/l_hippocampus_ProbabilityMap.nii.gz \
--probabilityMapsRightHippocampus $probabilityMapDir/r_hippocampus_ProbabilityMap.nii.gz \
--deformationFromSubjectToTemplate dummy \
--outputBinaryLeftAccumben $outputDir/left_accumben.nii.gz \
--outputBinaryRightAccumben $outputDir/right_accumben.nii.gz \
--outputBinaryLeftCaudate $outputDir/left_caudate.nii.gz \
--outputBinaryRightCaudate $outputDir/right_caudate.nii.gz \
--outputBinaryLeftPutamen $outputDir/left_putamen.nii.gz \
--outputBinaryRightPutamen $outputDir/right_putamen.nii.gz \
--outputBinaryLeftGlobus $outputDir/left_globus.nii.gz \
--outputBinaryRightGlobus $outputDir/right_globus.nii.gz \
--outputBinaryLeftThalamus $outputDir/left_thalamus.nii.gz \
--outputBinaryRightThalamus $outputDir/right_thalamus.nii.gz \
--outputBinaryLeftHippocampus $outputDir/left_hippocampus.nii.gz \
--outputBinaryRightHippocampus $outputDir/right_hippocampus.nii.gz \
--xmlFilename $outputDir/output.xml

