# BRAINSCut
BRAINSCut is a machine-learning based segmentation tool.
Built-in models provide human brain subcortical structure segmentations.
The subcortical structures of interests include accumben necleus, caudate,
putamen, globus pallidum, thalamus, and hippocampus in left and right hemisphere.

## Stages
BRAINSCut provides four stages. Each stage is either for training or testing phase.

1. Generate Probability Maps (command line option: --generateProbability)

2. Create input vectors (command line option: --createVectors)

3. Train model (command line option: --trainModel)

4. Apply model and get segmentation (command line option --applyModel)

**Multi-modal processing, t1 and t2 weighted input, assumes that t1 and t2 are pre-aligned (co-registered) each other**
## Usage Example

### Case 1: Getting new caudate segmentation using BRAINSTools provided Model Files
Eventually, we would like to run one of following command lines:
```
# For T1 only inputs:
BRAINSCut --applyModel \
  --netConfiguration  myNetConfigurationForCaudate.xml \
  --modelFilename  T1OnlyModels/trainModelFile.txtD0060NT0060_caudate_LinearWithMask.gz \
  --method RandomForest \
  --numberOfTrees 60  --randomTreeDepth 60
# For multi-modal T1 and T2 inputs
BRAINSCut --applyModel \
  --netConfiguration  myNetConfigurationForCaudate.xml \
  --modelFilename  modelFiles/trainModelFile.txtD0060NT0060_caudate_LinearWithMask.gz \
  --method RandomForest \
  --numberOfTrees 60  --randomTreeDepth 60
```
As you notice from above command line, BRAINSTools provide two sets of models
for six subcortical structures.
One utlizes T1 image and the other utilizes T1 and T2 images.
If both T1 and T2 images are available, we recommend use both and follow the step
'multi-modal T1 and T2 inputs'. Otherwise, please follow step 't2 only inputs'.

#### Prepare Input Files
  * T1 image
  * T2 image (in the case of multimodality model)
  * SG (Sum of gradient image: in the case of multimodality model)
  * Atlas to Subject Registration (Warping) Filename (We strongly recommend ANTs)
#### Create XML
  * [Create XML] Create a working directory
  ```
  $ mkdir mySubjectTrial
  $ cd mySubjectTrial
  ```
  * [Create XML] Copy a provided xml template and modify properly.
  ```
  # For T1 only inputs:
  $ cp ./BRAINSCut/Example/caudateT1OnlyApply.xml ./myCaudateTrialOnSubject.xml
  # For multi-modal T1 and T2 inputs
  $ cp ./BRAINSCut/Example/caudateT1T2Apply.xml ./myCaudateTrialOnSubject.xml
  ```
  * [Create XML] ReplaceWithBRAINSToolsBuildDir: change this to point your BRAINSTools directory.
  Make sure all the files provided by BRAINSTools exist:
    * Atlas/Atlas_20131115/template_t1_denoised_gaussian.nii.gz
    * Atlas/Atlas_20131115/spatialImages/rho.nii.gz
    * Atlas/Atlas_20131115/spatialImages/phi.nii.gz
    * Atlas/Atlas_20131115/spatialImages/theta.nii.gz
    * Atlas/Atlas_20131115/probabilityMaps/l_caudate_ProbabilityMap.nii.gz
    * Atlas/Atlas_20131115/probabilityMaps/r_caudate_ProbabilityMap.nii.gz
  * [Create XML] Fill in \<DataSet...\> Section.
  This section describes your subject of interest. (Replace all the square bracket parts [] ).

  ```
  # For T1 only inputs:
    <DataSet Name="subject" Type="Apply"      OutputDir="./" >
    <Image Type="T1" Filename="[Your_T1_brain_MRI_directory]/[your_T1_image_filename]" />
    <Registration SubjToAtlasRegistrationFilename=""
       AtlasToSubjRegistrationFilename="[Your_registration_directory]/[your_atlas_to_subject.h5]"
       ID="BSpline_ROI" />
    </DataSet>
  # For multi-modal T1 and T2 inputs
    <DataSet Name="subject" Type="Apply"      OutputDir="./" >
    <Image Type="T1" Filename="[Your_T1_brain_MRI_directory]/[your_t1_image_filename]" />
    <Image Type="T2" Filename="[Your_T1_brain_MRI_directory]/[your_t2_image_filename]" />
    <Image Type="SG" Filename="[your_SG_image_filename]" />
    <Registration SubjToAtlasRegistrationFilename=""
      AtlasToSubjRegistrationFilename="[Your_Registration_directory]/[your_atlas_to_subject.h5]"
       ID="BSpline_ROI" />
    </DataSet>
  ```
  
  * [Create XML] *Multimodal Only* Create SG, Sum of Gradient, Image from T1 and T2 using GenerateSummedGradientImage
  ```
  $BRAINSToolsBinDir/GenerateSummedGradientImage \
    --inputVolume1 [your_T1_image_filename] \
    --inputVolume2 [your_T2_image_filename] \
    --outputFileName [your_SG_image_filename]
  ```
#### Run the command
```
# For T1 only inputs:
BRAINSCut --applyModel \
  --netConfiguration  myNetConfigurationForCaudate.xml \
  --modelFilename  T1OnlyModels/trainModelFile.txtD0060NT0060_caudate_LinearWithMask.gz \
  --method RandomForest \
  --numberOfTrees 60  --randomTreeDepth 60
# For multi-modal T1 and T2 inputs
BRAINSCut --applyModel \
  --netConfiguration  myNetConfigurationForCaudate.xml \
  --modelFilename  modelFiles/trainModelFile.txtD0060NT0060_caudate_LinearWithMask.gz \
  --method RandomForest \
  --numberOfTrees 60  --randomTreeDepth 60
```
#### Output
  * Following three outputs are segmentation:
    * subjectANNLabel_l_caudate.nii.gz : Binary file for left caudate
    * subjectANNLabel_r_caudate.nii.gz : Binary file for right caudate
    * subject_ANNLabel_seg.nii.gz : Both left and right label map
  * Rest of following was mainly for debugging:
    * ANNContinuousPredictionl_caudatesubject.nii.gz
    * subjectANNLabel_l_caudate.nii.gzdef.nii.gz
    * ANNContinuousPredictionr_caudatesubject.nii.gz
    * subjectANNLabel_r_caudate.nii.gzdef.nii.gz
    * subject_ANNLabel_seg.nii.gz_AmbiguousMap.nii.gz

### Case 2: Getting new subcortical segmentation using BRAINSTools provided Model Files *Other than caudate*
The process of getting other subcortical structures are pretty same to the Case 1. 
Make sure to use 'Example/subcorticalT1OnlyApply.xml' or  'Example/subcorticalT1T2Apply.xml'. 
The main differences are in the *Normalization*. 
```
## for accumben, putamen, globus, thalamus, and hippocampus
Normalization          = "IQR"
```
And then run the similar commend line with *proper* model file:
```
/Shared/sinapse/scratch/eunyokim/src/NamicExternal/build_Mac_201501/bin/BRAINSCut --applyModel \
  --netConfiguration  [your_xml_filename] \
  --modelFilename  [BRAINSTools_model_file_name_according_to_the_below_tables] \
  --method RandomForest \
  --numberOfTrees 60  --randomTreeDepth 60
```

### Reference for T1 only application:
ROI             | Model File | Normalization
--------------- | ---------- | ------------
l/r caudate     | T1OnlyModels/trainModelFile.txtD0060NT0060_caudate_LinearWithMask.gz | Linear
l/r accumben    | T1OnlyModels/trainModelFile.txtD0060NT0060_accumben.gz | IQR
l/r globus      | T1OnlyModels/trainModelFile.txtD0060NT0060_globus.gz | IQR
l/r putamen     | T1OnlyModels/trainModelFile.txtD0060NT0060_putamen.gz | IQR
l/r thalamus    | T1OnlyModels/trainModelFile.txtD0060NT0060_thalamus.gz | IQR
l/r hippocampus | T1OnlyModels/trainModelFile.txtD0060NT0060_hippocampus.gz | IQR

### Reference for T1 and T2 Multimodal application:
ROI             | Model File | Normalization
--------------- | ---------- | ------------
l/r caudate     | modelFiles/trainModelFile.txtD0060NT0060_caudate_LinearWithMask.gz | Linear
l/r accumben    | modelFiles/trainModelFile.txtD0060NT0060_accumben.gz | IQR
l/r globus      | modelFiles/trainModelFile.txtD0060NT0060_globus.gz | IQR
l/r putamen     | modelFiles/trainModelFile.txtD0060NT0060_putamen.gz | IQR
l/r thalamus    | modelFiles/trainModelFile.txtD0060NT0060_thalamus.gz | IQR
l/r hippocampus | modelFiles/trainModelFile.txtD0060NT0060_hippocampus.gz | IQR


## BRAINSCut Description with respect to the BRAINSTools AutoWorkUp Pipeline.

The resulting data set of bias-corrected average T1 and/or T2 images from BRAINSABC are subsequently segmented for subcortical structures using an automated segmentation framework, BRAINSCut. BRAINSCut is an extension of previous works that employs robust random forest machine learning that has been validated on multi-site MR data (Kim, Magnotta, Liu, & Johnson, 2014). The subcortical structures of interest include accumben, caudate, putamen, hippocampus, and thalamus. The result of this procedure were again visually inspected and resulted in a success rate greater than 90%. All the development processing was blinded to clinical data, such as HD gene-expansion status, gender, and age.

*References*

• Kim, E. Y., Magnotta, V. a, Liu, D., & Johnson, H. J. (2014). Stable Atlas-based Mapped Prior (STAMP) machine-learning segmentation for multicenter large-scale MRI data. Magnetic Resonance Imaging, 32(7), 832–844. doi:10.1016/j.mri.2014.04.016

• Kim, E. Y., & Johnson, H. (2010). Multi-structure segmentation of multi-modal brain images using artificial neural networks. Analysis, 7623, 76234B–76234B–12. doi:10.1117/12.844613

• Powell, S. (2006). Automated brain segmentation using neural networks. Medical Imaging 2006: Image Processing (p. 61443Q–61443Q–11). SPIE.

• Powell, S., Magnotta, V. a, Johnson, H., Jammalamadaka, V. K., Pierson, R., & Andreasen, N. C. (2008). Registration and machine learning-based automated segmentation of subcortical and cerebellar brain structures. NeuroImage, 39(1), 238–47. doi:10.1016/j.neuroimage.2007.05.063
