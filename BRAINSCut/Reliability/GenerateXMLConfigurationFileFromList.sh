## Input : 
##     - List of Data Directory ( currentSubjectDirectory training/apply )
##     - List of Mask ( name )
##     - Output xml file name for xml file
##     - Training Model Directory ( for Registration/ProbabilityMap/Model output )

if [ $# != 4 ]; then
  echo "Incorrect Number of Argument"  
  echo "Usage::"  
  echo "$0 [ListFilename] [ROIListFilename] [OutputXMLFile] [HN]"
  exit 1;
fi
LISTofDataDirectory=$1;
LISTofMask=$2;
OutputXMLFile=$3
ANNHN=$4;

ModelDate="20110919"

ModelDirectory=`dirname $OutputXMLFile`

InputVectorDirectory="$ModelDirectory/InputVectors/"
TrainingVectorDirectory="$ModelDirectory/TrainedModels${ANNHN}/"
ProbabilityMapDirectory="$ModelDirectory/ProbabilityMaps/"

mkdir $InputVectorDirectory
mkdir $TrainingVectorDirectory
mkdir $ProbabilityMapDirectory

SpatialLocationDirectory="/hjohnson/HDNI/PREDICT_TRAINING/regina_ann/TrainingModels/3T_DataSet_2011May/ProbabilityMaps/T1_BSpline_ROIWithMay04/"


InputVectorFilename="${ModelDate}InputVector.txt"
TrainingModelFilename="${ModelDate}ANNModel.txt"
DeformTargetImage="T1";
TemplateDirectory="/ipldev//scratch/eunyokim/src/BRAINS-COMPILE/src/bin/Atlas/Atlas_20110701/"

InputVectorGradientSize="1";

ANNIteration="100"
DeformationDirectory="ANNDeformation${ModelDate}"
RegistrationID="BSpline_ROI"







##
## write header 
##
echo "<AutoSegProcessDescription> "                                                                   > $OutputXMLFile  
# 
# template
echo "  <DataSet Name=\"template\" Type=\"Atlas\" >"                                                  >> ${OutputXMLFile} 
echo "      <Image Type=\"T1\" Filename=\"${TemplateDirectory}/template_t1.nii.gz\" />"               >> ${OutputXMLFile} 
echo "      <Image Type=\"T2\" Filename=\"na\" />"                                                    >> ${OutputXMLFile} 
echo "      <Image Type=\"T1GAD\" Filename=\"na\" />"                                                    >> ${OutputXMLFile} 
echo "      <Image Type=\"T2GAD\" Filename=\"na\" />"                                                    >> ${OutputXMLFile} 
echo "      <Image Type=\"SG\" Filename=\"na\" />"                                                    >> ${OutputXMLFile} 
echo "      <SpatialLocation Type=\"rho\" Filename=\"${SpatialLocationDirectory}/rho.nii.gz\" />"     >> ${OutputXMLFile} 
echo "      <SpatialLocation Type=\"phi\" Filename=\"${SpatialLocationDirectory}/phi.nii.gz\" />"     >> ${OutputXMLFile} 
echo "      <SpatialLocation Type=\"theta\" Filename=\"${SpatialLocationDirectory}/theta.nii.gz\" />" >> ${OutputXMLFile} 
echo "  </DataSet>"                                                                                   >> ${OutputXMLFile} 

echo ""              >> ${OutputXMLFile} 



#
# Registration 
#
echo "  <RegistrationConfiguration "                                             >> ${OutputXMLFile}
echo "                     ImageTypeToUse  = \"${DeformTargetImage}\" "          >> ${OutputXMLFile}
echo "                     ID              = \"${RegistrationID}\" "             >> ${OutputXMLFile}
echo "                     BRAINSROIAutoDilateSize= \"1\" "                      >> ${OutputXMLFile}
echo "   />"                                                                     >> ${OutputXMLFile}

#
# Add ANNParams

echo "  <ANNParams        Iterations             = \"${ANNIteration}\""   >> ${OutputXMLFile} 
echo "                    MaximumVectorsPerEpoch = \"700000\""            >> ${OutputXMLFile} 
echo "                    EpochIterations        = \"100\""               >> ${OutputXMLFile} 
echo "                    ErrorInterval          = \"1\""                 >> ${OutputXMLFile} 
echo "                    DesiredError           = \"0.000001\""          >> ${OutputXMLFile} 
echo "                    NumberOfHiddenNodes    = \"${ANNHN}\""          >> ${OutputXMLFile} 
echo "                    ActivationSlope        = \"1.0\""               >> ${OutputXMLFile} 
echo "                    ActivationMinMax       = \"1.0\""               >> ${OutputXMLFile} 
echo "   />"                                                              >> ${OutputXMLFile} 

echo ""  >> ${OutputXMLFile} 

#
# Add Neural Params

echo "   <NeuralNetParams MaskSmoothingValue     = \"0.0\" "                                                  >> ${OutputXMLFile} 
echo "                    GradientProfileSize    = \"${InputVectorGradientSize}\""                            >> ${OutputXMLFile} 
echo "                    TrainingVectorFilename = \"${InputVectorDirectory}/${InputVectorFilename}\" "       >> ${OutputXMLFile} 
echo "                    TrainingModelFilename  = \"${TrainingVectorDirectory}/${TrainingModelFilename}\" "  >> ${OutputXMLFile} 
echo "                    TestVectorFilename     = \"na\" "                                                   >> ${OutputXMLFile} 
echo "                    Normalization          = \"true\" "                                                 >> ${OutputXMLFile}
echo "   />"                                                                                                  >> ${OutputXMLFile} 

echo ""  >> ${OutputXMLFile} 


echo "<ApplyModel         CutOutThresh           = \"0.05\" " >> ${OutputXMLFile} 
echo "                    MaskThresh             = \"0.4\" "  >> ${OutputXMLFile} 
echo "   />"                                                  >> ${OutputXMLFile} 
echo ""  >> ${OutputXMLFile} 




#
# Add Probability Map ( Region of Interest)
while read maskline generateVector
do
    echo "  <ProbabilityMap StructureID    = \"${maskline}\" "                                                  >> ${OutputXMLFile} 
    echo "                  Gaussian       = \"0.5\""                                                           >> ${OutputXMLFile} 
    echo "                  GenerateVector = \"${generateVector}\""                                             >> ${OutputXMLFile} 
    echo "                  Filename       = \"${ProbabilityMapDirectory}/${maskline}_ProbabilityMap.nii.gz\""  >> ${OutputXMLFile} 
    echo "   />"                                                                                                >> ${OutputXMLFile} 
    echo ""  >> ${OutputXMLFile} 
done < $LISTofMask

##
## write training/applying data set specification loop
##
# Example of line::
#   - /paulsen/IPIG/predict_3T_MR/site-073/2334/43971/10_AUTO.v020 Apply ManualMaskDir 
#   - /paulsen/IPIG/predict_3T_MR/site-073/2334/43971/10_AUTO.v020 Train ManualMaskDir

while read currentSubjectDirectory currentDataType ANNManualDirectory
do
  currentSubject=( `echo $currentSubjectDirectory|tr "/" "\n"|grep '[0-9]\{5\}'`)  

  # 
  # Add Image
  # TODO::deal with type of apply

  ANNOutputDirectory="${ModelDirectory}/${currentSubject}_${ANNHN}/"

  if [ "${currentDataType}" = Apply ]; then
    echo "  <DataSet Name=\"${currentSubject}\" Type=\"${currentDataType}\""     >> ${OutputXMLFile} 
    echo "           OutputDir=\"${ANNOutputDirectory}\" >"  >> ${OutputXMLFile}
  else
    echo "  <DataSet Name=\"${currentSubject}\" Type=\"${currentDataType}\" >"   >> ${OutputXMLFile} 
  fi

  echo   "      <Image Type=\"T1\" Filename=\"${currentSubjectDirectory}/${currentSubject}_AVG_T1.nii.gz\" />"  >> ${OutputXMLFile} 
  echo   "      <Image Type=\"T2\" Filename=\"${currentSubjectDirectory}/${currentSubject}_AVG_T2.nii.gz\" />"  >> ${OutputXMLFile} 
  echo   "      <Image Type=\"T1GAD\" Filename=\"${currentSubjectDirectory}/${currentSubject}_AVG_T1_GAD.nii.gz\" />"  >> ${OutputXMLFile} 
  echo   "      <Image Type=\"T2GAD\" Filename=\"${currentSubjectDirectory}/${currentSubject}_AVG_T2_GAD.nii.gz\" />"  >> ${OutputXMLFile} 
  echo   "      <Image Type=\"SG\" Filename=\"${currentSubjectDirectory}/ANN/${currentSubject}_Summed_Gradient.nii.gz\" />"      >> ${OutputXMLFile} 

  # 
  # including mask specification loop
  if [ "${currentDataType}" = "Train" ]; then
    while read maskline dummy
    do
      echo "      <Mask Type=\"${maskline}\" Filename=\"${currentSubjectDirectory}/${ANNManualDirectory}/${currentSubject}_${maskline}.nii.gz\" />"  >> ${OutputXMLFile} 
    done < $LISTofMask
  fi

  # 
  # generate directory
  if [ ! -d ${currentSubjectDirectory}/$DeformationDirectory/ ]; then
    mkdir ${currentSubjectDirectory}/$DeformationDirectory/
  fi
  #
  # Add Regisration parameters
  echo "      <Registration SubjToAtlasRegistrationFilename=\"${currentSubjectDirectory}/$DeformationDirectory/SubToAtlas_${currentSubject}.mat\" "  >> ${OutputXMLFile} 
  echo "                    AtlasToSubjRegistrationFilename=\"${currentSubjectDirectory}/$DeformationDirectory/AtlasToSub_${currentSubject}.mat\" "  >> ${OutputXMLFile} 
  echo "                    SubjectBinaryFilename=\"${currentSubjectDirectory}/${currentSubject}_BRAINSABC_Brain.nii.gz\" "                     >> ${OutputXMLFile}
  echo "                    AtlasBinaryFilename=\"${TemplateDirectory}/template_brain.nii.gz\" "                                         >> ${OutputXMLFile} 
  echo "                    ID=\"${RegistrationID}\" /> "                                                                                >> ${OutputXMLFile}


  echo "  </DataSet>"   >> ${OutputXMLFile} 
  echo ""               >> ${OutputXMLFile} 

done < $LISTofDataDirectory



echo "</AutoSegProcessDescription>"  >> ${OutputXMLFile} 

