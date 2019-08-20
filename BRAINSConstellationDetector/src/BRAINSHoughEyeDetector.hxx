/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
/*
 * Author: Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#include "BRAINSHoughEyeDetector.h"

namespace itk
{

template < typename TInputImage, typename TOutputImage >
typename TOutputImage::Pointer
RigidResampleInPlayByVersor3D( const typename TInputImage::ConstPointer & image,
                               VersorRigid3DTransform< double >::Pointer  versorRigid3DTfm )
{
  /** The output image will have exact the same index contents
   but with modified image info so that the index-to-physical mapping
   makes the image in the physical space aligned */
  using ResampleIPFilterType = itk::ResampleInPlaceImageFilter< TInputImage, TOutputImage >;

  typename ResampleIPFilterType::Pointer resampleIPFilter = ResampleIPFilterType::New();
  resampleIPFilter->SetInputImage( image );
  resampleIPFilter->SetRigidTransform( versorRigid3DTfm.GetPointer() );
  resampleIPFilter->Update();
  typename TOutputImage::Pointer resampledOutput = resampleIPFilter->GetOutput();
  return resampledOutput;
}

/** Isolate resample function from landamrk points **/
template < typename TInputImage, typename TOutputImage >
VersorRigid3DTransform< double >::Pointer
ResampleFromEyePoints( const typename TInputImage::PointType &    LE_Point,
                       const typename TInputImage::PointType &    RE_Point,
                       const typename TInputImage::ConstPointer & image )
{
  using InputIndexType = typename TInputImage::IndexType;
  using InputSizeType = typename TInputImage::SizeType;
  using InputRegionType = typename TInputImage::RegionType;
  using InputPointType = typename TInputImage::PointType;

  /* Transform and filter type alias */
  using VersorTransformType = VersorRigid3DTransform< double >;
  using VersorVectorType = typename VersorTransformType::OutputVectorType;

  static constexpr unsigned int Dimension = TInputImage::ImageDimension;
  /*
   * Align the image with eye centers
   */
  const InputRegionType ImageRegion = image->GetLargestPossibleRegion();
  const InputSizeType   size = ImageRegion.GetSize();

  // Compute translation
  InputPointType physicalStartLocation;
  InputPointType physicalStopLocation;
  {
    InputIndexType startIndex = ImageRegion.GetIndex();
    InputIndexType stopIndex;
    stopIndex[0] = startIndex[0] + size[0] - 1;
    stopIndex[1] = startIndex[1] + size[1] - 1;
    stopIndex[2] = startIndex[2] + size[2] - 1;
    image->TransformIndexToPhysicalPoint( startIndex, physicalStartLocation );
    image->TransformIndexToPhysicalPoint( stopIndex, physicalStopLocation );
  }

  // Space coordinate origin -> center of eye centers of input image
  VersorVectorType translation1;
  /** Center the image at a place close to the AC point.
   Here we use an approximation which is vertically at an IPD distance from center
   of eye centers, i.e., the center is at ( 0, -IPD, 0 ) in LPS coordinate system */
  VersorVectorType translation2;
  for ( unsigned int i = 0; i < Dimension; ++i )
  {
    translation1[i] = 0.5 * ( LE_Point[i] + RE_Point[i] );
    translation2[i] = 0;
  }
  const double IPD_dist = LE_Point.EuclideanDistanceTo( RE_Point );
  // https://en.wikipedia.org/wiki/Pupillary_distance  Minimum inter pupulary distance measured is 51mm for women

  // 1988 Anthropometric Survey MIN: 51mm,  Max 77mm, so add bit of margin on this stddev=3.6

  constexpr double mindistance_IPD = 51.0;
  constexpr double maxdistance_IPD = 77.0;
  constexpr double two_stddev_distance_IPD = 2.0 * 3.6;

  const double IPD = ( LE_Point.EuclideanDistanceTo( RE_Point ) );
  if ( IPD < ( mindistance_IPD - two_stddev_distance_IPD ) )
  {

    std::cerr << "ERROR:  'Left Eye' physical location must be at least 40mm to the left of the 'Right Eye': " << IPD
              << std::endl;

    std::cerr << "Right Eye: " << RE_Point << std::endl; //-27
    std::cerr << "Left Eye: " << LE_Point << std::endl;  //+31


    std::cerr << "     :   according to https://en.wikipedia.org/wiki/Pupillary_distance" << std::endl;
    exit( -1 );
  }

  if ( IPD > ( maxdistance_IPD + two_stddev_distance_IPD ) )
  {
    std::cerr << "ERROR:  'Left Eye' physical location must be at less than 86mm to the left of the 'Right Eye': "
              << IPD << std::endl;
    std::cerr << "Right Eye: " << RE_Point << std::endl; //-27
    std::cerr << "Left Eye: " << LE_Point << std::endl;  //+31

    std::cerr << "     :   according to https://en.wikipedia.org/wiki/Pupillary_distance" << std::endl;
    exit( -1 );
  }
  translation2[1] = IPD_dist;

  // Compute rotation in radian
  // about +S-axis
  VersorVectorType rotation1; // about +S-axis
  rotation1[0] = 0;
  rotation1[1] = 0;
  rotation1[2] = 1;

  // about +P-axis
  VersorVectorType rotation2; // about +P-axis
  rotation2[0] = 0;
  rotation2[1] = 1;
  rotation2[2] = 0;

  typename TOutputImage::PointType RotAngle;
  // Note: algorithm doesn't treat rotations about +L-axis
  RotAngle[0] = 0.0;
  RotAngle[1] = 0.0;
  RotAngle[2] = 0.0;


  // about +P-axis
  RotAngle[1] = std::atan( ( LE_Point[2] - RE_Point[2] ) / IPD_dist );
  // about +S-axis
  RotAngle[2] = -std::atan( ( LE_Point[1] - RE_Point[1] ) / IPD_dist );


  // Set rigid tranformation
  VersorRigid3DTransform< double >::Pointer versorRigid3DTfm = VersorRigid3DTransform< double >::New();
  versorRigid3DTfm->Translate( translation2 );
  versorRigid3DTfm->SetRotation( rotation1, -RotAngle[2] );
  versorRigid3DTfm->SetRotation( rotation2, -RotAngle[1] );
  versorRigid3DTfm->Translate( translation1 );

  return versorRigid3DTfm;
}

template < typename TInputImage, typename TOutputImage >
BRAINSHoughEyeDetector< TInputImage, TOutputImage >::BRAINSHoughEyeDetector()
{
  /** Input Parameters */
  // Note the default parameter set is designed for the IXI database
  this->m_NumberOfSpheres = 2;
  this->m_MinimumRadius = 11.;
  this->m_MaximumRadius = 13.;
  this->m_SigmaGradient = 1.;
  this->m_Variance = 1.;
  this->m_SphereRadiusRatio = 1.;
  this->m_VotingRadiusRatio = .5;
  this->m_Threshold = 10.;
  this->m_OutputThreshold = .8;
  this->m_GradientThreshold = 0.;
  this->m_NbOfThreads = 64;
  this->m_SamplingRatio = .2;
  this->m_HoughEyeDetectorMode = 1; // for T1-weighted image
  this->m_orig_lmk_CenterOfHeadMass.Fill( -9999.87654321 );

  this->m_R1 = 30;
  this->m_R2 = 120;
  this->m_Theta = 1.04719755; // 120 degrees 0.785398163;   //
                              // spread angle = 90 deg
  this->m_ResultsDir = ".";
  this->m_WritedebuggingImagesLevel = 0;

  /** Output parameters */
  this->m_AccumulatorImage = TInputImage::New();
  this->m_RoIImage = TInputImage::New();
  this->m_orig_lmk_LE.Fill( 123 );
  this->m_orig_lmk_RE.Fill( 123 );
  this->m_Failure = false;

  this->m_MaxInputPixelValue = 0;
  this->m_MinInputPixelValue = 0;

  this->m_MaxInputPixelValue = -1234;
  this->m_MinInputPixelValue = 1234;

  /** Internal parameters */
  this->m_orig2eyeFixedTransform = VersorTransformType::New();
}

template < typename TInputImage, typename TOutputImage >
void
BRAINSHoughEyeDetector< TInputImage, TOutputImage >::GenerateData()
{
  const InputImageConstPointer image( this->GetInput() );
  const InputRegionType        region = image->GetLargestPossibleRegion();
  {
    /*
     * Set RoI of the input image
     */
    this->m_RoIImage->CopyInformation( image );
    this->m_RoIImage->SetRegions( image->GetLargestPossibleRegion() );
    this->m_RoIImage->Allocate();
    this->m_RoIImage->FillBuffer( 0 );

    {
      // A unit vector pointing to the anterior direction
      typename InputPointType::VectorType unitVectorAnterior;
      unitVectorAnterior[0] = 0.;
      unitVectorAnterior[1] = -1.;
      unitVectorAnterior[2] = 0.;

      if ( ( this->m_R1 > 0 ) && ( this->m_R2 > this->m_R1 ) && ( this->m_Theta > 0 ) )
      {
        InputImageConstIterator It0( image, region );
        It0.GoToBegin();
        OutputImageIteratorWithIndex It1( this->m_RoIImage, region );
        It1.GoToBegin();

        while ( !( It0.IsAtEnd() || It1.IsAtEnd() ) )
        {
          InputPointType currPt;
          image->TransformIndexToPhysicalPoint( It1.GetIndex(), currPt );

          // Center of head mass to current vector
          const typename InputPointType::VectorType CMtoCurrVec = currPt - this->m_orig_lmk_CenterOfHeadMass;

          // posterior/anterior component of the vector
          const float CMtoCurrPA = CMtoCurrVec * unitVectorAnterior;

          // norm of CMtoCurrVec
          const float CMtoCurrNorm = CMtoCurrVec.GetNorm();

          // angle between current vector and unitVectorAnterior
          const float currTheta = acos( CMtoCurrPA / CMtoCurrNorm );

          if ( ( CMtoCurrNorm > this->m_R1 ) && ( CMtoCurrNorm < this->m_R2 ) && ( currTheta < this->m_Theta ) )
          {
            It1.Set( It0.Get() );
          }

          // save max/min pixel value info
          if ( It0.Get() < this->m_MinInputPixelValue )
          {
            this->m_MinInputPixelValue = It0.Get();
          }
          else if ( It0.Get() > this->m_MaxInputPixelValue )
          {
            this->m_MaxInputPixelValue = It0.Get();
          }
          ++It0;
          ++It1;
        }
      }
    }

    /*
     * Sphere detection with Hough transform radial voting filter
     */
    HoughFilterPointer houghFilter = HoughFilterType::New();
    if ( ( this->m_R1 > 0 ) && ( this->m_R2 > this->m_R1 ) && ( this->m_Theta > 0 ) )
    {
      houghFilter->SetInput( this->m_RoIImage );
    }
    else
    {
      std::cout << "Warning: RoI parameters are not valid"
                << "Set to ( default ) entire image" << std::endl;
      houghFilter->SetInput( image );
    }
    houghFilter->SetNumberOfSpheres( this->m_NumberOfSpheres );
    houghFilter->SetMinimumRadius( this->m_MinimumRadius );
    houghFilter->SetMaximumRadius( this->m_MaximumRadius );
    houghFilter->SetSigmaGradient( this->m_SigmaGradient );
    houghFilter->SetVariance( this->m_Variance );
    houghFilter->SetSphereRadiusRatio( this->m_SphereRadiusRatio );
    houghFilter->SetVotingRadiusRatio( this->m_VotingRadiusRatio );
    houghFilter->SetThreshold( this->m_Threshold );
    houghFilter->SetOutputThreshold( this->m_OutputThreshold );
    houghFilter->SetGradientThreshold( this->m_GradientThreshold );
    houghFilter->SetNbOfThreads( this->m_NbOfThreads );
    houghFilter->SetSamplingRatio( this->m_SamplingRatio );
    houghFilter->SetHoughEyeDetectorMode( this->m_HoughEyeDetectorMode );
    try
    {
      houghFilter->Update();
    }
    catch ( itk::ExceptionObject & excep )
    {
      std::cerr << "Failed houghFilter " << std::endl;
      std::cerr << excep << std::endl;
    }
    catch ( ... )
    {
      std::cout << "Failed on houghFilter exception occured" << std::endl;
    }
    this->m_AccumulatorImage = houghFilter->GetOutput();

    /*
     * Write debug image
     */
    if ( this->m_WritedebuggingImagesLevel > 1 )
    {
      {
        // Write debug ROI image
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetFileName( this->m_ResultsDir + "/HoughEyeROI.nii.gz" );
        writer->SetInput( this->m_RoIImage );
#if ITK_VERSION_MAJOR >= 5
        writer->SetUseCompression( true );
#endif
        try
        {
          writer->Update();
        }
        catch ( itk::ExceptionObject & excep )
        {
          std::cerr << "Cannot write the ROI image!" << std::endl;
          std::cerr << excep << std::endl;
        }
      }

      {
        // Write debug accumulator image
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetFileName( this->m_ResultsDir + "/HoughEyeAccumulator.nii.gz" );
        writer->SetInput( this->m_AccumulatorImage );
#if ITK_VERSION_MAJOR >= 5
        writer->SetUseCompression( true );
#endif
        try
        {
          writer->Update();
        }
        catch ( itk::ExceptionObject & excep )
        {
          std::cerr << "Cannot write the ROI image!" << std::endl;
          std::cerr << excep << std::endl;
        }
      }
    }

    /*
     * Computing eye centers by finding the points with highest
     * accumulated PDF
     */
    // Get some basic information of the image

    const SpheresListType spheres = houghFilter->GetSpheres();
    if ( spheres.size() < 2 )
    {
      std::cerr << "Error: The number of detected spheres is less than 2!" << std::endl;
      std::cerr << "The program will continue to run for generating some debug output for GUI corrector." << std::endl;
      std::cerr << "Make sure you set the debug level > 4." << std::endl;
      this->m_Failure = true;
      return;
    }

    InputIndexType indexEye1;
    SphereIterator itSpheres = spheres.begin();
    for ( unsigned int i = 0; i < Dimension; ++i )
    {
      indexEye1[i] = static_cast< unsigned long int >( ( *itSpheres )->GetObjectToParentTransform()->GetOffset()[i] );
    }

    ++itSpheres;

    InputIndexType indexEye2;
    for ( unsigned int i = 0; i < Dimension; ++i )
    {
      indexEye2[i] = static_cast< unsigned long int >( ( *itSpheres )->GetObjectToParentTransform()->GetOffset()[i] );
    }

    InputPointType physicalEye1;
    InputPointType physicalEye2;
    {
      image->TransformIndexToPhysicalPoint( indexEye1, physicalEye1 );
      image->TransformIndexToPhysicalPoint( indexEye2, physicalEye2 );
    }

    // We will determine the left and right of the eyes
    // Assuming that the degradation of the image does not
    // change their relative positions.
    if ( physicalEye1[0] > physicalEye2[0] ) // eye1 is on the left
    {
      for ( unsigned int i = 0; i < Dimension; ++i )
      {
        this->m_orig_lmk_LE[i] = physicalEye1[i];
        this->m_orig_lmk_RE[i] = physicalEye2[i];
      }
    }
    else
    {
      for ( unsigned int i = 0; i < Dimension; ++i )
      {
        this->m_orig_lmk_LE[i] = physicalEye2[i]; // eye2 is on the left
        this->m_orig_lmk_RE[i] = physicalEye1[i];
      }
    }

    /*
     * Metrics Collection
     */

    this->m_orig2eyeFixedTransform =
      ResampleFromEyePoints< TInputImage, TOutputImage >( this->m_orig_lmk_LE, this->m_orig_lmk_RE, image );
  }
}

template < typename TInputImage, typename TOutputImage >
void
BRAINSHoughEyeDetector< TInputImage, TOutputImage >::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << "HoughEyeDetectorMode: " << this->m_HoughEyeDetectorMode << std::endl;
  os << "Number Of Spheres: " << this->m_NumberOfSpheres << std::endl;
  os << "Minimum Radius:  " << this->m_MinimumRadius << std::endl;
  os << "Maximum Radius: " << this->m_MaximumRadius << std::endl;
  os << "Derivative Scale : " << this->m_SigmaGradient << std::endl;
  os << "Accumulator Blur Variance: " << this->m_Variance << std::endl;
  os << "Sphere Radius Ratio: " << this->m_SphereRadiusRatio << std::endl;
  os << "Voting Radius Ratio: " << this->m_VotingRadiusRatio << std::endl;
  os << "Threshold: " << this->m_Threshold << std::endl;
  os << "Output Threshold : " << this->m_OutputThreshold << std::endl;
  os << "Gradient Threshold: " << this->m_GradientThreshold << std::endl;
  os << "NbOfThreads: " << this->m_NbOfThreads << std::endl;
  os << "Sampling Ratio: " << this->m_SamplingRatio << std::endl;
  os << "Interior Radius of RoI: " << this->m_R1 << std::endl;
  os << "Exterior Radius of RoI: " << this->m_R2 << std::endl;
  os << "Spread Angle of RoI: " << this->m_Theta << std::endl;
}
} // namespace itk
