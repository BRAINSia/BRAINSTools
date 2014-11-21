

proc Genus0SurfaceGeneration { T1NormalizeBfc TissueClassImage Brain Ventricle LeftCaudate RightCaudate LeftPutamen
                               RightPutamen LeftThalamus RightThalamus ResultDir PatientId ScanId } {
  global BrainsConfig
  set WorkDir "$ResultDir/NewSurface"
  if { ![file isdirectory $WorkDir] } {
    file mkdir $WorkDir
    }
  set BRAINSABCAtlas $BrainsConfig(BRAINSABCAtlasDir) / template_t1.nii.gz
  set clipT1Filename $WorkDir / $
  {
    ScanId
  }

  _T1_clip.nii.gz
  set leftHemisphere $BrainsConfig(BRAINSABCAtlasDir) / template_leftHemisphere.nii.gz
  set rightHemisphere $BrainsConfig(BRAINSABCAtlasDir) / template_rightHemisphere.nii.gz
  set cerebellumMask $BrainsConfig(BRAINSABCAtlasDir) / template_cerebellum.nii.gz
  set ventricleMask $BrainsConfig(BRAINSABCAtlasDir) / ventricles.nii.gz

#Clip T1 image to the brain
  if
  {
    [Brains::Utils::CheckOutputsNewer[list $clipT1Filename] \
     [list $T1NormalizeBfc $Brain]] == false
  }

    {
    set t1Image[Brains::itk::LoadImage $T1NormalizeBfc "Signed-16bit"]
    set maskImage[Brains::itk::LoadImage $Brain "Signed-16bit"]
    set clipT1Image[Brains::itk::MaskImage $t1Image $maskImage]

    Brains::itk::SaveImage $clipT1Image $clipT1Filename
    $t1Image Delete
    $maskImage Delete
    $clipT1Image Delete
    }
  set affineTransform       $WorkDir / $ {ScanId} _EmsAtlas_Transform.mat
  set outputTransform       $WorkDir / $ {ScanId} _EmsAtlas_Warp.mhd
  set warpedLeftHemisphere  $WorkDir / $ {ScanId} _leftHemisphere_warped.nii.gz
  set warpedRightHemisphere $WorkDir / $ {ScanId} _rightHemisphere_warped.nii.gz
  set warpedCerebellum      $WorkDir / $ {ScanId} _cerebellum_warped.nii.gz
  set warpedVentricles      $WorkDir / $ {ScanId} _ventricles_warped.nii.gz
  set deformationField      $WorkDir / $ {ScanId} _Atlas_DeformationField.nii.gz
  set warpedAtlas           $WorkDir / $ {ScanId} _Atlas.nii.gz
  set                       outputList[list $affineTransform $outputTransform]
  lappend outputList $warpedCerebellum $warpedLeftHemisphere $warpedRightHemisphere
  lappend outputList $warpedVentricles $deformationField $warpedAtlas

#Register Atlas with Subject
  if {[Brains::Utils::CheckOutputsNewer $outputList \
       [list $clipT1Filename $BRAINSABCAtlas]] == false } {
#Register Atlas with Subject - Affine
    set affineTransform $WorkDir / $ {ScanId} _EmsAtlas_Transform.mat
    set TransformType   Rigid, ScaleVersor3D
    set Samples 400000
    set TranslationScale 1000
    set Iterations 1500
    set MinimumStepSize 0.05, 0.005, 0.0005
#Brains ::External::TestRunBrainsFit $clipT1Filename $BRAINSABCAtlas $affineTransform

#DemonsWarp BRAINSABCAtlas -> subjT1
    set outputTransform $WorkDir / $ {ScanId} _EmsAtlas_Warp.mhd
    set HistogramBins 256
    set MatchPoints 11
    set PyramidLevels 5
    set LevelIterations "400,200,100,10,2"
    set FilterSize "1,1,1"
    set InputPixelType "short"
    set OutputPixelType "short"
    set RegistrationType "Diffeomorphic"
    set DeformationFieldSigma 2
#Brains ::External::runBRAINSDemonsWarp $BRAINSABCAtlas $clipT1Filename \
    #$deformationField $warpedAtlas \
    #$HistogramBins $MatchPoints $PyramidLevels $LevelIterations $FilterSize \
    #$InputPixelType $OutputPixelType $RegistrationType $DeformationFieldSigma \
    #$affineTransform

#Warped Atlas based Regions
    set transformType "DeformationField"
#Brains ::External::runApplyWarp $leftHemisphere $clipT1Filename $transformType \
    #$deformationField $warpedLeftHemisphere
#Brains ::External::runApplyWarp $rightHemisphere $clipT1Filename $transformType \
    #$deformationField $warpedRightHemisphere
#Brains ::External::runApplyWarp $cerebellumMask $clipT1Filename $transformType \
    #$deformationField $warpedCerebellum
#Brains ::External::runApplyWarp $ventricleMask $clipT1Filename $transformType \
    #$deformationField $warpedVentricles
    }

#Fill in the image to eliminate ventricles, and subcortical strutures
  set ventricleMask[Brains::itk::LoadImage $warpedVentricles "Unsigned-8bit"]
  set rightCaudateMask[Brains::itk::LoadImage $RightCaudate "Unsigned-8bit"]
  set leftCaudateMask[Brains::itk::LoadImage $LeftCaudate "Unsigned-8bit"]
  set rightPutamenMask[Brains::itk::LoadImage $RightPutamen "Unsigned-8bit"]
  set leftPutamenMask[Brains::itk::LoadImage $LeftPutamen "Unsigned-8bit"]
  set rightThalamusMask[Brains::itk::LoadImage $RightThalamus "Unsigned-8bit"]
  set leftThalamusMask[Brains::itk::LoadImage $LeftThalamus "Unsigned-8bit"]
  set classImage[Brains::itk::LoadImage $TissueClassImage "Unsigned-8bit"]

#Generate a binary representation for subsequent operations
  set binaryVentricle[Brains::itk::BinaryThresholdImage $ventricleMask 1 255]
  set binaryRightCaudate[Brains::itk::BinaryThresholdImage $rightCaudateMask 1 255]
  set binaryLeftCaudate[Brains::itk::BinaryThresholdImage $leftCaudateMask 1 255]
  set binaryRightPutamen[Brains::itk::BinaryThresholdImage $rightPutamenMask 1 255]
  set binaryLeftPutamen[Brains::itk::BinaryThresholdImage $leftPutamenMask 1 255]
  set binaryRightThalamus[Brains::itk::BinaryThresholdImage $rightThalamusMask 1 255]
  set binaryLeftThalamus[Brains::itk::BinaryThresholdImage $leftThalamusMask 1 255]

  $ventricleMask Delete
  $rightCaudateMask Delete
  $leftCaudateMask Delete
  $rightPutamenMask Delete
  $leftPutamenMask Delete
  $rightThalamusMask Delete
  $leftThalamusMask Delete
  puts "Or Images"

  ##########Combine Regions to Fill Class Image ###########

  set binaryCaudate[Brains::itk::OrImage $binaryRightCaudate $binaryLeftCaudate]
  puts "binaryCaudate"
  set binaryPutamen[Brains::itk::OrImage $binaryRightPutamen $binaryLeftPutamen]
  puts "binary Putamen"
  set binaryThalamus[Brains::itk::OrImage $binaryRightThalamus $binaryLeftThalamus]
  set binaryBasalGanglia[Brains::itk::OrImage $binaryCaudate $binaryPutamen]
  set binarySubcortical[Brains::itk::OrImage $binaryThalamus $binaryBasalGanglia]
  set binaryFillRegion[Brains::itk::OrImage $binarySubcortical $binaryVentricle]
  puts "binaryFillRegion"

  $binaryRightCaudate Delete
  $binaryLeftCaudate Delete
  $binaryRightPutamen Delete
  $binaryLeftPutamen Delete
  $binaryRightThalamus Delete
  $binaryLeftThalamus Delete
  $binaryCaudate Delete
  $binaryPutamen Delete
  $binaryThalamus Delete
  $binaryBasalGanglia Delete
  $binarySubcortical Delete
  $binaryVentricle Delete

  puts "Combined Images"

#Create the Filled Class Image
  set scaledFillRegion[Brains::itk::ConstantImageMath $binaryFillRegion 230 Multiply]
  set imageList[list $scaledFillRegion $classImage]
  set filledClassImage[Brains::itk::NaryMaximumImageFilter $imageList]
  set resultClassFilename $WorkDir / $ {ScanId} _Filled_class.nii.gz
  Brains::itk::SaveImage $filledClassImage $resultClassFilename

  $scaledFillRegion Delete
  $classImage Delete
  $filledClassImage Delete

#Create a separate VTK file for each hemisphere surface
  set leftHemisphereImage $WorkDir / $ {ScanId} _leftTissueClass.nii.gz
  set rightHemisphereImage $WorkDir / $ {ScanId} _rightTissueClass.nii.gz
  set rightBinaryImage     $WorkDir / $ {ScanId} _right_binary.nii.gz
  set leftBinaryImage      $WorkDir / $ {ScanId} _left_binary.nii.gz
  set leftSurface          $ResultDir / $ {ScanId} _left_surface.vtk
  set rightSurface         $ResultDir / $ {ScanId} _right_surface.vtk
  set                      outputList[list $leftHemisphereImage $rightHemisphereImage $rightBinaryImage $leftBinaryImage
  ]
  lappend outputList $leftSurface $outputList

#Create surfaces using new BRAINSSurfaceGeneration code
  CreateGenus0BrainSurface $Brain $warpedCerebellum $leftBinaryImage \
  $warpedLeftHemisphere $resultClassFilename $leftHemisphereImage $leftSurface
  CreateGenus0BrainSurface $Brain $warpedCerebellum $rightBinaryImage \
  $warpedRightHemisphere $resultClassFilename $rightHemisphereImage $rightSurface
  return [list $leftSurface $rightSurface]

  }

#################################################################################
#NEW PROC TO CREATE A BRAIN SURFACE USING BRAINSSurfaceGeneration
#################################################################################
proc CreateGenus0BrainSurface { BrainMaskFile WarpedCerebellumFile BinaryOutputFilename \
                                HemisphereMaskFile TissueClassFile SurfaceImageFilename OutputSurfaceFilename } {

  set cerebellumMask[Brains::itk::LoadImage $WarpedCerebellumFile "Unsigned-8bit"]
  set binaryCerebellum[Brains::itk::BinaryThresholdImage $cerebellumMask 1 255]
  set notCerebellum[Brains::itk::NotImage $binaryCerebellum]
  $cerebellumMask Delete
  $binaryCerebellum Delete

  set brainMask[Brains::itk::LoadImage $BrainMaskFile "Unsigned-8bit"]
  set hemisphereMask[Brains::itk::LoadImage $HemisphereMaskFile "Unsigned-8bit"]
  set binaryBrain[Brains::itk::BinaryThresholdImage $brainMask 1 255]
  set binaryHemisphere[Brains::itk::BinaryThresholdImage $hemisphereMask 1 255]
  $brainMask Delete
  $hemisphereMask Delete

  set hemisphereRegion[Brains::itk::AndImage $binaryHemisphere $binaryBrain]
  $binaryBrain Delete
  $binaryHemisphere Delete

  set clippedRegion[Brains::itk::AndImage $hemisphereRegion $notCerebellum]
  $notCerebellum Delete
  $hemisphereRegion Delete

  set classImage[Brains::itk::LoadImage $TissueClassFile "Unsigned-8bit"]
  set hemisphereClassImage[Brains::itk::MaskImage  $classImage $clippedRegion]
  $classImage Delete
  $clippedRegion Delete

  Brains::itk::SaveImage $hemisphereClassImage "/tmp/class_hemisphere.nii.gz"

#Filter Image - Median followed by Anisotropic Diffusion
  set medianSurfaceImage[Brains::itk::MedianImageFilter $hemisphereClassImage]
  $hemisphereClassImage Delete

  set floatSurfaceImage[Brains::itk::CastImage $medianSurfaceImage "Float-single"]
  $medianSurfaceImage Delete

  set filteredSurfaceImage[Brains::itk::GradientAnisotropicDiffusionImageFilter $floatSurfaceImage 0.0625 1.0 5]
  $floatSurfaceImage Delete

#Threshold and Keep the Largest connected region
  set binary190Image[Brains::itk::BinaryThresholdImage $filteredSurfaceImage 190 255 0 "Unsigned-8bit"]
  $filteredSurfaceImage Delete

  set componentImage[Brains::itk::ConnectedComponentImage $binary190Image 0 "Unsigned-8bit"]
  $binary190Image Delete

  set relabelImage[Brains::itk::RelabelComponentImage $componentImage 500]
  $componentImage Delete

  set binaryImage[Brains::itk::BinaryThresholdImage $relabelImage 1 1 0 "Unsigned-8bit"]

#Save the 190 Binary Image for Surface generation
  Brains::itk::SaveImage $binaryImage $SurfaceImageFilename
  $binaryImage Delete

  if {[file extension $SurfaceImageFilename] == ".gz"} {
    set genus0ImageFilename[file rootname[file rootname $SurfaceImageFilename]] _genus0.nii.gz
    }
  else
    {
    set genus0ImageFilename[file rootname $SurfaceImageFilename] _genus0.nii.gz
    }
  set command "exec [file dirname [info nameofexecutable]]/GenusZeroImageFilterOriginal"
  append command " --biggestComponent"
  append command " --computeSurface"
  append command " --connectivity 18"
  append command " $SurfaceImageFilename"
  append command " $genus0ImageFilename"
  append command " $OutputSurfaceFilename"
  append command " >&@stdout"
  if {[catch
         {
         eval $ {command}
         } PrintTypeScript] } {
    puts "Error: $PrintTypeScript"
    return 1
    }

  return [list $SurfaceImageFilename $OutputSurfaceFilename]
  }

package require BrainsGlue
set baseDir / home / vince / images / surface / 011 8695
#set baseDir /Users/vince/images/surface/0118695
set T1NormalizeBfc $baseDir / 011 8695_T1.hdr
set TissueClassImage $baseDir / 011 8695_class.hdr
set BrainMask $baseDir / 011 8695_brain_trim.mask
set Ventricle $baseDir / 011 8695_ventricles.hdr
set LeftCaudate $baseDir / l_caud_cut.mask
set RightCaudate $baseDir / r_caud_cut.mask
set LeftPutamen $baseDir / l_put_cut.mask
set RightPutamen $baseDir / r_put_cut.mask
set LeftThalamus $baseDir / l_thal_cut.mask
set RightThalamus $baseDir / r_thal_cut.mask
set ResultDir $baseDir
set PatientId TMP
set ScanId 011 8695
Brains::AutoWorkup::GenusZeroSurfaceGeneration $T1NormalizeBfc $TissueClassImage \
  $BrainMask $LeftCaudate $RightCaudate $LeftPutamen \
  $RightPutamen $LeftThalamus $RightThalamus $ResultDir \
  $PatientId $ScanId
