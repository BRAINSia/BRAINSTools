---
title: Tools
nav_order: 2
---

<!-- GENERATED FILE - DO NOT EDIT.
     Regenerate with: python3 Utilities/Maintenance/generate_tool_catalog.py -->

# Tool Catalog

BRAINSTools provides the following command-line tools, grouped by the
category each declares in its Slicer module descriptor.

## Diffusion.GTRACT

### Anisotropy Map

This program will generate a scalar map of anisotropy, given a tensor representation. Anisotropy images are used for fiber tracking, but the anisotropy scalars are not defined along the path. Instead, the tensor representation is included as point data allowing all of these metrics to be computed using only the fiber tract point data. The images can be saved in any ITK supported format, but it is suggested that you use an image format that supports the definition of the image origin. This includes NRRD, NifTI, and Meta formats. These images can also be used for scalar analysis including regional anisotropy measures or VBM style analysis.

`gtractAnisotropyMap` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractAnisotropyMap.xml)

### Average B-Values

This program will directly average together the baseline gradients (b value equals 0) within a DWI scan. This is usually used after gtractCoregBvalues.

`gtractAverageBvalues` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractAverageBvalues.xml)

### B-Spline Transform Inversion

This program will invert a B-Spline transform using a thin-plate spline approximation.

`gtractInvertBSplineTransform` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractInvertBSplineTransform.xml)

### Clip Anisotropy

This program will zero the first and/or last slice of an anisotropy image, creating a clipped anisotropy image.

`gtractClipAnisotropy` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractClipAnisotropy.xml)

### Compare Tracts

This program will halt with a status code indicating whether a test tract is nearly enough included in a standard tract in the sense that every fiber in the test tract has a low enough sum of squares distance to some fiber in the standard tract modulo spline resampling of every fiber to a fixed number of points.

`compareTractInclusion` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/compareTractInclusion.xml)

### Concat DWI Images

This program will concatenate two DTI runs together.

`gtractConcatDwi` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractConcatDwi.xml)

### Copy Image Orientation

This program will copy the orientation from the reference image into the moving image. Currently, the registration process requires that the diffusion weighted images and the anatomical images have the same image orientation (i.e. Axial, Coronal, Sagittal). It is suggested that you copy the image orientation from the diffusion weighted images and apply this to the anatomical image. This image can be subsequently removed after the registration step is complete. We anticipate that this limitation will be removed in future versions of the registration programs.

`gtractCopyImageOrientation` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractCopyImageOrientation.xml)

### Coregister B-Values

This step should be performed after converting DWI scans from DICOM to NRRD format. This program will register all gradients in a NRRD diffusion weighted 4D vector image (moving image) to a specified index in a fixed image. It also supports co-registration with a T2 weighted image or field map in the same plane as the DWI data. The fixed image for the registration should be a b0 image. A mutual information metric cost function is used for the registration because of the differences in signal intensity as a result of the diffusion gradients. The full affine allows the registration procedure to correct for eddy current distortions that may exist in the data. If the eddyCurrentCorrection is enabled, relaxationFactor (0.25) and maximumStepSize (0.1) should be adjusted.

`gtractCoregBvalues` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractCoregBvalues.xml)

### Coregister B0 to Anatomy

This program will register a Nrrd diffusion weighted 4D vector image to a fixed anatomical image. Two registration methods are supported for alignment with anatomical images: Rigid and B-Spline. The rigid registration performs a rigid body registration with the anatomical images and should be done as well to initialize the B-Spline transform. The B-SPline transform is the deformable transform, where the user can control the amount of deformation based on the number of control points as well as the maximum distance that these points can move. The B-Spline registration places a low dimensional grid in the image, which is deformed. This allows for some susceptibility related distortions to be removed from the diffusion weighted images. In general the amount of motion in the slice selection and read-out directions direction should be kept low. The distortion is in the phase encoding direction in the images. It is recommended that skull stripped (i.e. image containing only brain with skull removed) images shoud be used for image co-registration with the B-Spline transform.

`gtractCoRegAnatomy` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractCoRegAnatomy.xml)

### Cost Fast Marching

This program will use a fast marching fiber tracking algorithm to identify fiber tracts from a tensor image. This program is the first portion of the algorithm. The user must first run gtractFastMarchingTracking to generate the actual fiber tracts.  This algorithm is roughly based on the work by G. Parker et al. from IEEE Transactions On Medical Imaging, 21(5): 505-512, 2002. An additional feature of including anisotropy into the vcl_cost function calculation is included.

`gtractCostFastMarching` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractCostFastMarching.xml)

### Create Displacement Field

This program will compute forward deformation from the given Transform. The size of the DF is equal to MNI space

`gtractTransformToDisplacementField` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractTransformToDisplacementField.xml)

### Create Guide Fiber

This program will create a guide fiber by averaging fibers from a previously generated tract.

`gtractCreateGuideFiber` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractCreateGuideFiber.xml)

### Extract Nrrd Index

This program will extract a 3D image (single vector) from a vector 3D image at a given vector index.

`extractNrrdVectorIndex` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/extractNrrdVectorIndex.xml)

### Fast Marching Tracking

This program will use a fast marching fiber tracking algorithm to identify fiber tracts from a tensor image. This program is the second portion of the algorithm. The user must first run gtractCostFastMarching to generate the vcl_cost image. The second step of the algorithm implemented here is a gradient descent soplution from the defined ending region back to the seed points specified in gtractCostFastMarching. This algorithm is roughly based on the work by G. Parker et al. from IEEE Transactions On Medical Imaging, 21(5): 505-512, 2002. An additional feature of including anisotropy into the vcl_cost function calculation is included.

`gtractFastMarchingTracking` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractFastMarchingTracking.xml)

### Fiber Tracking

This program implements four fiber tracking methods (Free, Streamline, GraphSearch, Guided). The output of the fiber tracking is vtkPolyData (i.e. Polylines) that can be loaded into Slicer3 for visualization. The poly data can be saved in either old VTK format files (.vtk) or in the new VTK XML format (.xml). The polylines contain point data that defines ther Tensor at each point along the fiber tract. This can then be used to rendered as glyphs in Slicer3 and can be used to define severeal scalar measures without referencing back to the anisotropy images. (1) Free tracking is a basic streamlines algorithm. This is a direct implementation of the method original proposed by Basser et al. The tracking follows the primarty eigenvector. The tracking begins with seed points in the starting region. Only those voxels above the specified anisotropy threshold in the starting region are used as seed points. Tracking terminates either as a result of maximum fiber length, low ansiotropy, or large curvature. This is a great way to explore your data. (2) The streamlines algorithm is a direct implementation of the method originally proposed by Basser et al. The tracking follows the primary eigenvector. The tracking begins with seed points in the starting region. Only those voxels above the specified anisotropy threshold in the starting region are used as seed points. Tracking terminates either by reaching the ending region or reaching some stopping criteria. Stopping criteria are specified using the following parameters: tracking threshold, curvature threshold, and max length. Only paths terminating in the ending region are kept in this method. The TEND algorithm proposed by Lazar et al. (Human Brain Mapping 18:306-321, 2003) has been instrumented. This can be enabled using the --useTend option while performing Streamlines tracking. This utilizes the entire diffusion tensor to deflect the incoming vector instead of simply following the primary eigenvector. The TEND parameters are set using the --tendF and --tendG options. (3) Graph Search tracking is the first step in the full GTRACT algorithm developed by Cheng et al. (NeuroImage 31(3): 1075-1085, 2006) for finding the tracks in a tensor image. This method was developed to generate fibers in a Tensor representation where crossing fibers occur. The graph search algorithm follows the primary eigenvector in non-ambigous regions and utilizes branching and a graph search algorithm in ambigous regions. Ambiguous tracking regions are defined based on two criteria: Branching Al Threshold (anisotropy values below this value and above the traching threshold) and Curvature Major Eigen (angles of the primary eigenvector direction and the current tracking direction). In regions that meet this criteria, two or three tracking paths are considered. The first is the standard primary eigenvector direction. The second is the seconadary eigenvector direction. This is based on the assumption that these regions may be prolate regions. If the Random Walk option is selected then a third direction is also considered. This direction is defined by a cone pointing from the current position to the centroid of the ending region. The interior angle of the cone is specified by the user with the Branch/Guide Angle parameter. A vector contained inside of the cone is selected at random and used as the third direction. This method can also utilize the TEND option where the primary tracking direction is that specified by the TEND method instead of the primary eigenvector. The parameter '--maximumBranchPoints' allows the tracking to have this number of branches being considered at a time. If this number of branch points is exceeded at any time, then the algorithm will revert back to a streamline alogrithm until the number of branches is reduced. This allows the user to constrain the computational complexity of the algorithm. (4) The second phase of the GTRACT algorithm is Guided Tracking. This method incorporates anatomical information about the track orientation using an initial guess of the fiber track. In the originally proposed GTRACT method, this would be created from the fibers resulting from the Graph Search tracking. However, in practice this can be created using any method and could be defined manually. To create the guide fiber the program gtractCreateGuideFiber can be used. This program will load a fiber tract that has been generated and create a centerline representation of the fiber tract (i.e. a single fiber). In this method, the fiber tracking follows the primary eigenvector direction unless it deviates from the guide fiber track by a angle greater than that specified by the '--guidedCurvatureThreshold' parameter. The user must specify the guide fiber when running this program.

`gtractFiberTracking` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractFiberTracking.xml)

### Image Conformity

This program will straighten out the Direction and Origin to match the Reference Image.

`gtractImageConformity` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractImageConformity.xml)

### Invert Displacement Field

This program will invert a deformatrion field. The size of the deformation field is defined by an example image provided by the user

`gtractInvertDisplacementField` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractInvertDisplacementField.xml)

### Resample Anisotropy

This program will resample a floating point image using either the Rigid or B-Spline transform. You may want to save the aligned B0 image after each of the anisotropy map co-registration steps with the anatomical image to check the registration quality with another tool.

`gtractResampleAnisotropy` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractResampleAnisotropy.xml)

### Resample B0

This program will resample a signed short image using either a Rigid or B-Spline transform. The user must specify a template image that will be used to define the origin, orientation, spacing, and size of the resampled image.

`gtractResampleB0` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractResampleB0.xml)

### Resample Code Image

This program will resample a short integer code image using either the Rigid or Inverse-B-Spline transform.  The reference image is the DTI tensor anisotropy image space, and the input code image is in anatomical space.

`gtractResampleCodeImage` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractResampleCodeImage.xml)

### Resample DWI In Place

Resamples DWI image to structural image.

`gtractResampleDWIInPlace` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractResampleDWIInPlace.xml)

### Resample Fibers

This program will resample a fiber tract with respect to a pair of deformation fields that represent the forward and reverse deformation fields.

`gtractResampleFibers` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractResampleFibers.xml)

### Rigid Transform Inversion

This program will invert a Rigid transform.

`gtractInvertRigidTransform` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractInvertRigidTransform.xml)

### Tensor Estimation

This step will convert a b-value averaged diffusion tensor image to a 3x3 tensor voxel image. This step takes the diffusion tensor image data and generates a tensor representation of the data based on the signal intensity decay, b values applied, and the diffusion difrections. The apparent diffusion coefficient for a given orientation is computed on a pixel-by-pixel basis by fitting the image data (voxel intensities) to the Stejskal-Tanner equation. If at least 6 diffusion directions are used, then the diffusion tensor can be computed. This program uses itk::DiffusionTensor3DReconstructionImageFilter. The user can adjust background threshold, median filter, and isotropic resampling.

`gtractTensor` &middot; version 5.8.0 &middot; [documentation](http://wiki.slicer.org/slicerWiki/index.php/Modules:GTRACT) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/GTRACT/Cmdline/gtractTensor.xml)

## Diffusion.Import and Export

### Diffusion-weighted DICOM Import (DWIConvert)

Converts diffusion weighted MR images in DICOM series into NRRD format for analysis in Slicer. This program has been tested on only a limited subset of DTI DICOM formats available from Siemens, GE, and Philips scanners. Work in progress to support DICOM multi-frame data. The program parses DICOM header to extract necessary information about measurement frame, diffusion weighting directions, b-values, etc, and write out a NRRD image. For non-diffusion weighted DICOM images, it loads in an entire DICOM series and writes out a single dicom volume in a .nhdr/.raw pair.

`DWIConvert` &middot; version 5.8.0 &middot; [documentation](https://slicer.readthedocs.io/en/latest/user_guide/modules/dwiconvert.html) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/DWIConvert/DWIConvert.xml)

## Diffusion.Utilities

### DWI Cleanup (BRAINS)

Remove bad gradients/volumes from DWI NRRD file.

`BRAINSDWICleanup` &middot; version 5.8.0 &middot; [documentation](https://slicer.readthedocs.io/en/latest/user_guide/modules/brainsdwicleanup.html) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSDWICleanup/BRAINSDWICleanup.xml)

## Filtering

### Edge Map Generator (BRAINS)

Inverse of Maximum Gradient Image

`GenerateEdgeMapImage` &middot; version 5.8.0 &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSSuperResolution/GenerateEdgeMap/GenerateEdgeMapImage.xml)

### Percentile Rescaling

Computes the percentile rescaling of an image

`BRAINSIntensityNormalize` &middot; version 5.8.0 &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSIntensityNormalize/BRAINSIntensityNormalize.xml)

## Quantification

### Label Statistics (BRAINS)

Compute image statistics within each label of a label map.

`BRAINSLabelStats` &middot; version 5.8.0 &middot; [documentation](https://slicer.readthedocs.io/en/latest/user_guide/modules/brainslabelstats.html) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSLabelStats/BRAINSLabelStats.xml)

## Registration

### General Registration (BRAINS)

Register a three-dimensional volume to a reference volume (Mattes Mutual Information by default). Method described in BRAINSFit: Mutual Information Registrations of Whole-Brain 3D Images, Using the Insight Toolkit, Johnson H.J., Harris G., Williams K., The Insight Journal, 2007. https://hdl.handle.net/1926/1291

`BRAINSFit` &middot; version 5.8.0 &middot; [documentation](https://slicer.readthedocs.io/en/latest/user_guide/modules/brainsfit.html) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSFit/BRAINSFit.xml)

### Registration Metric Test (BRAINS)

Compare Mattes/MSQ metric value for two input images and a possible input BSpline transform.

`PerformMetricTest` &middot; version 5.8.0 &middot; [documentation](https://slicer.readthedocs.io/en/latest/user_guide/modules/performmetrictest.html) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSFit/PerformMetricTest.xml)

### Resample Image (BRAINS)

This program collects together three common image processing tasks that all involve resampling an image volume: Resampling to a new resolution and spacing, applying a transformation (using an ITK transform IO mechanisms) and Warping (using a vector image deformation field).

`BRAINSResample` &middot; version 5.8.0 &middot; [documentation](https://slicer.readthedocs.io/en/latest/user_guide/modules/brainsresample.html) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSResample/BRAINSResample.xml)

### Resize Image (BRAINS)

This program is useful for downsampling an image by a constant scale factor.

`BRAINSResize` &middot; version 5.8.0 &middot; [documentation](https://slicer.readthedocs.io/en/latest/user_guide/modules/brainsresize.html) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSResample/BRAINSResize.xml)

## Registration.Specialized

### Fiducial Registration (BRAINS)

Computes a rigid, similarity or affine transform from a matched list of fiducials

`BRAINSTransformFromFiducials` &middot; version 5.8.0 &middot; [documentation](http://www.slicer.org/slicerWiki/index.php/Modules:TransformFromFiducials-Documentation-3.6) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSConstellationDetector/src/BRAINSTransformFromFiducials.xml)

## Segmentation.Specialized

### BRAINS Binary Mask Editor Based On Landmarks(BRAINS)

Edit a binary mask by cutting it with planes defined by three landmarks or by a landmark and an axis direction, clearing voxels on the specified side of each plane.

`BinaryMaskEditorBasedOnLandmarks` &middot; version 5.8.0 &middot; [documentation](http://www.nitrc.org/projects/brainscdetector/) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSConstellationDetector/src/BinaryMaskEditorBasedOnLandmarks.xml)

### Brain Landmark Constellation Detector (BRAINS)

This program will find the mid-sagittal plane, a constellation of landmarks in a volume, and create an AC/PC aligned data set with the AC point at the center of the voxel lattice (labeled at the origin of the image physical space.)  Part of this work is an extention of the algorithms originally described by Dr. Babak A. Ardekani, Alvin H. Bachman, Model-based automatic detection of the anterior and posterior commissures on MRI scans, NeuroImage, Volume 46, Issue 3, 1 July 2009, Pages 677-682, ISSN 1053-8119, DOI: 10.1016/j.neuroimage.2009.02.030.  (http://www.sciencedirect.com/science/article/B6WNP-4VRP25C-4/2/8207b962a38aa83c822c6379bc43fe4c)

`BRAINSConstellationDetector` &middot; version 5.8.0 &middot; [documentation](http://www.nitrc.org/projects/brainscdetector/) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSConstellationDetector/src/BRAINSConstellationDetector.xml)

### Clean Contiguous Label Map (BRAINS)

From a range of label map values, extract the largest contiguous region of those labels

`ESLR` &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSABC/brainseg/ESLR.xml)

### Create Label Map From Probability Maps (BRAINS)

Given A list of Probability Maps, generate a LabelMap.

`BRAINSCreateLabelMapFromProbabilityMaps` &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSCreateLabelMapFromProbabilityMaps/BRAINSCreateLabelMapFromProbabilityMaps.xml)

### Create best representative label map)

given a list of label map images, create a representative/average label map.

`BRAINSMultiSTAPLE` &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSMultiSTAPLE/BRAINSMultiSTAPLE.xml)

### Foreground masking (BRAINS)

This program is used to create a mask over the most prominent foreground region in an image.  This is accomplished via a combination of otsu thresholding and a closing operation.

`BRAINSROIAuto` &middot; version 5.8.0 &middot; [documentation](https://slicer.readthedocs.io/en/latest/user_guide/modules/brainsroiauto.html) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSROIAuto/BRAINSROIAuto.xml)

### Intra-subject registration, bias Correction, and tissue classification (BRAINS)

Atlas-based tissue segmentation method. This is an algorithmic extension of work performed at UNC and Utah, with significant contributions from Marcel Prastawa and the NAMIC project.

`BRAINSABC` &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSABC/brainseg/BRAINSABC.xml)

### Mask Hole Filling (BRAINS)

Cleans up a mask image by filling any holes

`BRAINSCleanMask` &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSABC/brainseg/BRAINSCleanMask.xml)

### Pure Plugs Mask

This program gets several modality image files and returns a binary mask that defines the pure plugs.

`GeneratePurePlugMask` &middot; version 5.8.0 &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSABC/brainseg/GeneratePurePlugMask.xml)

### Tissue Classification

This program will generate an 8-bit continuous tissue classified image based on BRAINSABC posterior images.

`BRAINSPosteriorToContinuousClass` &middot; version 5.8.0 &middot; [documentation](http://www.nitrc.org/plugins/mwiki/index.php/brains:BRAINSClassify) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSPosteriorToContinuousClass/BRAINSPosteriorToContinuousClass.xml)

## Utilities

### Strip Rotation (BRAINS)

Read an Image, write out same image with identity rotation matrix plus an ITK transform file

`BRAINSStripRotation` &middot; version 5.8.0 &middot; [documentation](https://slicer.readthedocs.io/en/latest/user_guide/modules/brainsstriprotation.html) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSStripRotation/BRAINSStripRotation.xml)

## Utilities.BRAINS

### ACPC-Aligned Landmark Conversion

This program converts the original landmark files to the acpc-aligned landmark files

`landmarksConstellationAligner` &middot; version 5.8.0 &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSConstellationDetector/src/landmarksConstellationAligner.xml)

### Align Mid Saggital Brain (BRAINS)

Resample an image into ACPC alignment ACPCDetect

`BRAINSAlignMSP` &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSConstellationDetector/src/BRAINSAlignMSP.xml)

### Average Fiducials

This program gets several fcsv file each one contains several landmarks with the same name but slightly different coordinates. For EACH landmark we compute the average coordination.

`GenerateAverageLmkFile` &middot; version 5.8.0 &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSConstellationDetector/landmarkStatistics/GenerateAverageLmkFile.xml)

### Brain Deface from T1/T2 image (BRAINS)

This program: 1) will deface images from a set of images.  Inputs must be ACPC aligned, and AC, PC, LE, RE provided.

`BRAINSDeface` &middot; version 5.8.0 &middot; [documentation](https://slicer.readthedocs.io/en/latest/user_guide/modules/brainsdeface.html) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSDeface/BRAINSDeface.xml)

### Brain Extraction from T1/T2 image (BRAINS)

This program: 1) generates a weighted mixture image optimizing the mean and variance and 2) produces a mask of the brain volume

`BRAINSMush` &middot; version 5.8.0 &middot; [documentation](http://mri.radiology.uiowa.edu) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSMush/BRAINSMush.xml)

### Clip Inferior of Center of Brain (BRAINS)

This program will read the inputVolume as a short int image, write the BackgroundFillValue everywhere inferior to the lower bound, and write the resulting clipped short int image in the outputVolume.

`BRAINSClipInferior` &middot; version 5.8.0 &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSConstellationDetector/src/BRAINSClipInferior.xml)

### ConstellationDetectorGUI (BRAINS)

This program provides the user with a GUI tool to view/manipulate landmarks for an input volume.

`BRAINSConstellationDetectorGUI` &middot; version 5.8.0 &middot; [documentation](http://www.nitrc.org/projects/brainscdetector/) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSConstellationDetector/gui/BRAINSConstellationDetectorGUI.xml)

### Eye Detector (BRAINS)

Locate left and right eye centers in a T1 head image using a Hough radial-voting transform seeded by the center of head mass, and write the eye-aligned resampled volume.

`BRAINSEyeDetector` &middot; version 5.8.0 &middot; [documentation](http://www.nitrc.org/projects/brainscdetector/) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSConstellationDetector/src/BRAINSEyeDetector.xml)

### Generate Landmarks Model (BRAINS)

Train up a model for BRAINSConstellationDetector

`BRAINSConstellationModeler` &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSConstellationDetector/src/BRAINSConstellationModeler.xml)

### Generate Landmarks Weights (BRAINS)

Train up a list of Weights for the Landmarks in BRAINSConstellationDetector

`landmarksConstellationWeights` &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSConstellationDetector/src/landmarksConstellationWeights.xml)

### Initialized Control Points (BRAINS)

Outputs bspline control points as landmarks

`BRAINSInitializedControlPoints` &middot; version 5.8.0 &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSInitializedControlPoints/BRAINSInitializedControlPoints.xml)

### Label Map from Probability Images

Given a list of probability maps for labels, create a discrete label map where only the highest probability region is used for the labeling.

`GenerateLabelMapFromProbabilityMap` &middot; version 5.8.0 &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSABC/brainseg/GenerateLabelMapFromProbabilityMap.xml)

### Landmark FCSV to HDF5 Converter (BRAINS)

Convert a collection of fcsv files to a HDF5 format file

`fcsv_to_hdf5` &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSConstellationDetector/src/fcsv_to_hdf5.xml)

### Landmark Linear Modeler (BRAINS)

Training linear model using EPCA. Implementation based on my MS thesis, "A METHOD FOR AUTOMATED LANDMARK CONSTELLATION DETECTION USING EVOLUTIONARY PRINCIPAL COMPONENTS AND STATISTICAL SHAPE MODELS"

`BRAINSLinearModelerEPCA` &middot; version 5.8.0 &middot; [documentation](http://www.nitrc.org/projects/brainscdetector/) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSConstellationDetector/src/BRAINSLinearModelerEPCA.xml)

### Landmark Transform (BRAINS)

This utility program estimates the affine transform to align the fixed landmarks to the moving landmarks, and then generate the resampled moving image to the same physical space as that of the reference image.

`BRAINSLmkTransform` &middot; version 5.8.0 &middot; [documentation](http://www.nitrc.org/projects/brainscdetector/) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSConstellationDetector/src/BRAINSLmkTransform.xml)

### Landmark Transform Initializer (BRAINS)

Create transformation file (*.h5) from a pair of landmarks (*fcsv) files.

`BRAINSLandmarkInitializer` &middot; version 5.8.0 &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSLandmarkInitializer/BRAINSLandmarkInitializer.xml)

### Landmarks Transformation

This program converts the original landmark file to the target landmark file using the input transform.

`BRAINSConstellationLandmarksTransform` &middot; version 5.8.0 &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSConstellationDetector/src/BRAINSConstellationLandmarksTransform.xml)

### MidACPC Landmark Insertion

This program gets a landmark fcsv file and adds a new landmark as the midpoint between AC and PC points to the output landmark fcsv file

`insertMidACPCpoint` &middot; version 5.8.0 &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSConstellationDetector/src/insertMidACPCpoint.xml)

### Segment based on rectangular region of joint histogram (BRAINS)

This tool creates binary regions based on segmenting multiple image modalitities at once.

`BRAINSMultiModeSegment` &middot; version 5.8.0 &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSMultiModeSegment/BRAINSMultiModeSegment.xml)

### Snapshot Writer (BRAINS)

Create 2D snapshot of input images. Mask images are color-coded

`BRAINSSnapShotWriter` &middot; version 5.8.0 &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSSnapShotWriter/BRAINSSnapShotWriter.xml)

### Transform Convert (BRAINS)

Convert ITK transforms to higher order transforms

`BRAINSTransformConvert` &middot; version 5.8.0 &middot; [documentation](https://slicer.readthedocs.io/en/latest/user_guide/modules/brainstransformconvert.html) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSTransformConvert/BRAINSTransformConvert.xml)

### Trim Foreground In Direction (BRAINS)

This program will trim off the neck and also air-filling noise from the inputImage.

`BRAINSTrimForegroundInDirection` &middot; version 5.8.0 &middot; [documentation](http://www.nitrc.org/projects/art/) &middot; [source](https://github.com/BRAINSia/BRAINSTools/blob/main/BRAINSConstellationDetector/src/BRAINSTrimForegroundInDirection.xml)
