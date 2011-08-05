## Input : 
##     - List of Data Directory ( directory training/apply )
##       ex)
##       /paulsen/IPIG/predict_3T_MR/site-073/2334/43971/10_AUTO.NN3Tv20101112 Apply
##       /paulsen/IPIG/predict_3T_MR/site-144/0619/96294/10_AUTO.NN3Tv20101112 Train
##       /paulsen/IPIG/predict_3T_MR/site-024/0132/38235/10_AUTO.NN3Tv20101112 Train
##       /paulsen/IPIG/predict_3T_MR/site-024/0138/84460/10_AUTO.NN3Tv20101112 Train
##       /paulsen/IPIG/predict_3T_MR/site-024/0241/37022/10_AUTO.NN3Tv20101112 Train
##       /paulsen/IPIG/predict_3T_MR/site-024/0247/34082/10_AUTO.NN3Tv20101112 Train
##       /paulsen/IPIG/predict_3T_MR/site-024/0276/91147/10_AUTO.NN3Tv20101112 Train
##       /paulsen/IPIG/predict_3T_MR/site-024/0326/10335/10_AUTO.NN3Tv20101112 Train
##       /paulsen/IPIG/predict_3T_MR/site-024/0329/77478/10_AUTO.NN3Tv20101112 Train
##       /paulsen/IPIG/predict_3T_MR/site-024/0486/87259/10_AUTO.NN3Tv20101112 Train
##       /paulsen/IPIG/predict_3T_MR/site-024/0489/40535/10_AUTO.NN3Tv20101112 Train
##
##     - List of Mask File( name )
##       ex) 
##       l_caudate
##       r_caudate
##     - Output xml file name for xml file
##     - Training Model Directory ( for Registration/ProbabilityMap/Model output )
## 

Deforming_Img="T1"
Deformation_Type="BSpline_ROI";
Atlas_Dir="/ipldev//scratch/eunyokim/src/BRAINS3/BRAINSTools/BRAINSABC/Atlas_20101105/"
HN="60";
GradSize="1";
It="1"


if [ $# != 4 ]; then
  echo "Incorrect Number of Argument"  
  echo "Usage::"  >> ${OutputXMLFile} 
  echo "$0 [File for List of Data Directory ( directory training/apply ) ]"  
  echo "                    [File for List of Mask] [Output XML file name]" } 
  echo "                    [Output Model Directory] " 
  exit 1;
else
  LISTofDataDirectory=$1
  if [ ! -f $LISTofDataDirectory ]; then
    echo "$LISTofDataDirectory : does not exists"  
    exit 1
  fi
  LISTofMask=$2
  if [ ! -f $LISTofMask ]; then
    echo "$LISTofMask : does not exists"  
    exit 1
  fi
  OutputXMLFile=$3;
  OutputModelDirectory=$4;
fi
InputVector_Dir="${OutputModelDirectory}/InputVectors/${Deforming_Img}_${Deformation_Type}/"
TrainingModel_Dir="${OutputModelDirectory}/Trainings/${Deforming_Img}_${Deformation_Type}/"
dummy_filename=(`echo ${OutputXMLFile} |tr "/" "\n"|tail -1|tr "\." "\n"|head -1`) 
##
## write header 
##
echo "<AutoSegProcessDescription> " > $OutputXMLFile  

##
## write training/applying data set specification loop
##
while read line
do
  # Example of line::
  #   - /paulsen/IPIG/predict_3T_MR/site-073/2334/43971/10_AUTO.v020 Apply
  #   - /paulsen/IPIG/predict_3T_MR/site-073/2334/43971/10_AUTO.v020 Train
  subject=( `echo $line |tr "/" "\n"|grep '[0-9]\{5\}'`)  
  directory=(`echo $line|sed "s/ .*//"`)  
  data_type=(`echo $line|sed "s/.* //"`) # space separated Apply/Train 

  # 
  # Add Image
  # TODO::deal with type of apply
  if [ "${data_type}" = Apply ]; then
    echo "  <DataSet Name=\"${subject}\" Type=\"${data_type}\""  >> ${OutputXMLFile} 
    echo "           OutputDir=\"${TrainingModel_Dir}/${It}_Apply/${subject}\" >"  >> ${OutputXMLFile} 
  else
    echo "  <DataSet Name=\"${subject}\" Type=\"${data_type}\">"  >> ${OutputXMLFile} 
  fi
  echo "      <Image Type=\"T1\" Filename=\"${directory}/${subject}_AVG_T1.nii.gz\" />"  >> ${OutputXMLFile} 
  echo "      <Image Type=\"T2\" Filename=\"${directory}/${subject}_AVG_T2.nii.gz\" />"  >> ${OutputXMLFile} 
  echo "      <Image Type=\"GMI_GM\" Filename=\"${directory}/${subject}_GMI_fAttributeV2.nii.gz\" />"  >> ${OutputXMLFile} 

  # 
  # including mask specification loop
  while read maskline
  do
    echo "      <Mask Type=\"${maskline}\" Filename=\"${directory}/ANN_train/${subject}_${maskline}.nii.gz\" />"  >> ${OutputXMLFile} 
  done < $LISTofMask
  #
  # Add Regisration parameters
  Deform_Dir="${OutputModelDirectory}/Deformations/${Deforming_Img}_${Deformation_Type}/"
  echo "      <Registration SubjToAtlasRegistrationFilename=\"${Deform_Dir}/SubToAtlas/${subject}.mat\" "  >> ${OutputXMLFile} 
  echo "                   AtlasToSubjRegistrationFilename=\"${Deform_Dir}/AtlasToSub/${subject}.mat\" "  >> ${OutputXMLFile} 
  echo "                   ID=\"${Deforming_Img}_${Deformation_Type}\" /> "  >> ${OutputXMLFile} 
  echo "  </DataSet>"  >> ${OutputXMLFile} 
  echo ""  >> ${OutputXMLFile} 

done < $LISTofDataDirectory
##
## write footer( template, annparams, etc... )
##

# 
# template
echo "  <DataSet Name=\"template\" Type=\"Atlas\" >"  >> ${OutputXMLFile} 
echo "      <Image Type=\"T1\" Filename=\"${Atlas_Dir}/template_t1.nii.gz\" />"  >> ${OutputXMLFile} 
echo "      <Image Type=\"T2\" Filename=\"na\" />"  >> ${OutputXMLFile} 
echo "      <Image Type=\"GMI_GM\" Filename=\"na\" />"  >> ${OutputXMLFile} 
echo "  </DataSet>"  >> ${OutputXMLFile} 
echo ""  >> ${OutputXMLFile} 

#
# Registration 
#
echo "  <RegistrationParams " >> ${OutputXMLFile}
echo "                     Type      = \"${Deforming_Img}_${Deformation_Type}\" ">> ${OutputXMLFile}
echo "                     Command   = \"GenerateBSplineTransform.tcl\" ">> ${OutputXMLFile}
echo "                     ImageType = \"T1\" ">> ${OutputXMLFile}
echo "                     ID        = \"${Deforming_Img}_${Deformation_Type}\" ">> ${OutputXMLFile}
echo "   />">> ${OutputXMLFile}




#
# Add Probability Map ( Region of Interest)
ProbMap_Dir="${OutputModelDirectory}/ProbabilityMaps/${Deforming_Img}_${Deformation_Type}/"
while read maskline
do
    echo "  <ProbabilityMap StructureID    = \"${maskline}\" "  >> ${OutputXMLFile} 
    echo "                  Gaussian       = \"1.0\""  >> ${OutputXMLFile} 
    echo "                  GenerateVector = \"true\""  >> ${OutputXMLFile} 
    echo "                  Filename       = \"${ProbMap_Dir}/${maskline}_ProbabilityMap.nii.gz\""  >> ${OutputXMLFile} 
    echo "                  rho            =\"${ProbMap_Dir}/rho.nii.gz\""  >> ${OutputXMLFile} 
    echo "                  phi            =\"${ProbMap_Dir}/phi.nii.gz\""  >> ${OutputXMLFile} 
    echo "                  theta            =\"${ProbMap_Dir}/theta.nii.gz\""  >> ${OutputXMLFile} 
    echo "   />"  >> ${OutputXMLFile} 
    echo ""  >> ${OutputXMLFile} 
done < $LISTofMask

#
# Add ANNParams

echo "  <ANNParams        Iterations             = \"${It}\""  >> ${OutputXMLFile} 
echo "                    MaximumVectorsPerEpoch = \"700000\""  >> ${OutputXMLFile} 
echo "                    EpochIterations        = \"100\""  >> ${OutputXMLFile} 
echo "                    ErrorInterval          = \"1\""  >> ${OutputXMLFile} 
echo "                    DesiredError           = \"0.000001\""  >> ${OutputXMLFile} 
echo "                    NumberOfHiddenNodes    = \"${HN}\""  >> ${OutputXMLFile} 
echo "                    ActivationSlope        = \"1.0\""  >> ${OutputXMLFile} 
echo "                    ActivationMinMax       = \"1.0\""  >> ${OutputXMLFile} 
echo "   />"  >> ${OutputXMLFile} 
echo ""  >> ${OutputXMLFile} 

#
# Add Neural Params

echo "   <NeuralNetParams MaskSmoothingValue     = \"0.0\" "  >> ${OutputXMLFile} 
echo "                    GradientProfileSize    = \"${GradSize}\""  >> ${OutputXMLFile} 
echo "                    TrainingVectorFilename = \"${InputVector_Dir}/${dummy_filename}_vector.txt\" "  >> ${OutputXMLFile} 
echo "                    TrainingModelFilename  = \"${TrainingModel_Dir}/${dummy_filename}_model.txt\" "  >> ${OutputXMLFile} 
echo "                    TestVectorFilename     = \"na\" "  >> ${OutputXMLFile} 
echo "   />"  >> ${OutputXMLFile} 
echo ""  >> ${OutputXMLFile} 


echo "<ApplyModel         CutOutThresh           = \"0.05\" "   >> ${OutputXMLFile} 
echo "                    CutOutGaussian         = \"0\""  >> ${OutputXMLFile} 
echo "                    MaskThresh             = \"0.3\" "  >> ${OutputXMLFile} 
echo "   />"  >> ${OutputXMLFile} 
echo ""  >> ${OutputXMLFile} 

echo "</AutoSegProcessDescription>"  >> ${OutputXMLFile} 

