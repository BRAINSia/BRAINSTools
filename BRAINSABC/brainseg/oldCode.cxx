  {
  // TODO:  Investigate if a small gaussian filter is needed here after
  // clipping.
#if 0
  if( 0 )
    {
    typename itk::DiscreteGaussianImageFilter<TProbabilityImage, TProbabilityImage>::Pointer gaussianFilter =
      itk::DiscreteGaussianImageFilter<TProbabilityImage, TProbabilityImage>::New();
    gaussianFilter->SetInput(WarpedPriorsList[i]);
    typename itk::DiscreteGaussianImageFilter<TProbabilityImage, TProbabilityImage>::ArrayType MySigmas;
    for( size_t s = 0; s < MySigmas.Size(); s++ )
      {      // Make it 115% of the voxel spacings
             // #MySigmas[s]=1.10*WarpedPriorsList[k]->GetSpacing()[s];
      MySigmas[s] = 1.5 * WarpedPriorsList[i]->GetSpacing()[s];
      }
    gaussianFilter->SetVariance( MySigmas );
    const FloatingPrecision MaximumError = 0.01;     // Copied from example.
    gaussianFilter->SetMaximumError( MaximumError );
    gaussianFilter->Update();
    smoothProbImage = gaussianFilter->GetOutput();
    }
  if( this->m_DebugLevel > 5 )
    {
    std::stringstream CurrentEMIteration_stream("");
    CurrentEMIteration_stream << CurrentEMIteration;
    // Write the subject candidate regions

    std::ostringstream oss;
    oss << this->m_OutputDebugDir << "CANDIDIDATE_PROB_" << this->m_PriorNames[i] << "_LEVEL_"
        << CurrentEMIteration_stream.str() << ".nii.gz" << std::ends;
    std::string fn = oss.str();
    muLogMacro( << "Writing Subject Candidate Region." << fn << std::endl );
    muLogMacro( << std::endl );

    typedef itk::ImageFileWriter<TProbabilityImage> ByteWriterType;
    typename ByteWriterType::Pointer writer = ByteWriterType::New();
    writer->SetInput(smoothProbImage);
    writer->SetFileName(fn.c_str() );
    writer->UseCompressionOn();
    writer->Update();
    }
#endif
  }

/**
*/
template <class TInputImage, class TProbabilityImage>
typename ByteImageType::Pointer
EMSegmentationFilter<TInputImage, TProbabilityImage>
::ExtractNonAirRegion(const unsigned int /*CurrentIterationID*/,
                      const typename TProbabilityImage::Pointer & BackgroundPrior )
{
  ByteImageType::Pointer NonBackground;

    {
    typedef itk::ResampleImageFilter<ByteImageType, ByteImageType>
      ByteResamplerType;
    muLogMacro(<< " Computing non-air region ." << std::endl);

    typedef itk::BinaryThresholdImageFilter<TProbabilityImage, ByteImageType> ProbThresholdType;
    typename ProbThresholdType::Pointer probNonAirRegion = ProbThresholdType::New();
    probNonAirRegion->SetInput( BackgroundPrior );
    probNonAirRegion->SetInsideValue(1);
    probNonAirRegion->SetOutsideValue(0);
    probNonAirRegion->SetLowerThreshold(-1.0e10);
    probNonAirRegion->SetUpperThreshold(0.99); // Set the uppper threshold to
                                               // 99% of image
    probNonAirRegion->Update();
    NonBackground = probNonAirRegion->GetOutput();
    }

  return NonBackground;
}

typename ByteImageType::Pointer ComputeForegroundProbMask(
  const std::vector<typename TProbabilityImage::Pointer> & probList, const BoolVectorType & IsForegroundPriorVector);

template <class TInputImage, class TProbabilityImage>
typename ByteImageType::Pointer
EMSegmentationFilter<TInputImage, TProbabilityImage>
::ComputeForegroundProbMask(const std::vector<typename TProbabilityImage::Pointer> & probList,
                            const BoolVectorType & IsForegroundPriorVector )
{
  muLogMacro(<< "ComputeForegroundProbMask" << std::endl );

  const unsigned int numPriors = probList.size();

  typename ByteImageType::Pointer currForegroundMask = ByteImageType::New();
  currForegroundMask->CopyInformation(probList[0]);
  currForegroundMask->SetRegions(probList[0]->GetLargestPossibleRegion() );
  currForegroundMask->Allocate();

  const ProbabilityImageSizeType size = probList[0]->GetLargestPossibleRegion().GetSize();
  const typename ByteImageType::PixelType insideMaskValue = 100;

    {
#pragma omp parallel for
    for( LOOPITERTYPE kk = 0; kk < (LOOPITERTYPE)size[2]; kk++ )
      {
      for( LOOPITERTYPE jj = 0; jj < (LOOPITERTYPE)size[1]; jj++ )
        {
        for( LOOPITERTYPE ii = 0; ii < (LOOPITERTYPE)size[0]; ii++ )
          {
          const ProbabilityImageIndexType currIndex = {{ii, jj, kk}};
          FloatingPrecision tmp = 0.0;
          for( unsigned int iprior = 0; iprior < numPriors; iprior++ )
            {
            const bool fgflag = IsForegroundPriorVector[iprior];
            if( fgflag == true )
              {
              tmp += probList[iprior]->GetPixel(currIndex);
              }
            }
          if( tmp > 0.005 ) // Only include if the sum of the non-background
                            // priors are greater than 1/2 of 1 percent chance
            {
            currForegroundMask->SetPixel(currIndex, static_cast<typename ByteImageType::PixelType>(insideMaskValue) );
            }
          else
            {
            currForegroundMask->SetPixel(currIndex, 0);
            }
          }
        }
      }
    }
#if 1
    {
    // Pre-Dilate mask
    typedef itk::BinaryBallStructuringElement<unsigned char, 3> StructElementType;
    typedef
      itk::BinaryDilateImageFilter<ByteImageType, ByteImageType,
                                   StructElementType> DilateType;

    StructElementType structel;
    structel.SetRadius(1);
    structel.CreateStructuringElement();

    typename DilateType::Pointer dil = DilateType::New();
    dil->SetDilateValue(50);
    dil->SetKernel(structel);
    dil->SetInput(currForegroundMask);

    dil->Update();

    ByteImagePointer dilmask = dil->GetOutput();
    // a simple assignment is probably sufficient, test both ways?
    currForegroundMask = CopyImage<ByteImageType>(dilmask);
    }
#endif

#if 0
    { // Now clip to the filed of view
    if( m_NonAirRegion.IsNotNull() )
      {
      if( size != m_NonAirRegion->GetLargestPossibleRegion().GetSize() )
        {
        itkExceptionMacro(
          << "NonAirRegion mask size mismatch " << size << " != "
          << m_NonAirRegion->GetLargestPossibleRegion().GetSize() << " ." << std::endl );
        }

      unsigned long count = 0;
        {
#pragma omp parallel for reduction(+:count)
        for( LOOPITERTYPE kk = 0; kk < (LOOPITERTYPE)size[2]; kk++ )
          {
          for( LOOPITERTYPE jj = 0; jj < (LOOPITERTYPE)size[1]; jj++ )
            {
            for( LOOPITERTYPE ii = 0; ii < (LOOPITERTYPE)size[0]; ii++ )
              {
              const ProbabilityImageIndexType currIndex = {{ii, jj, kk}};
              if( ( m_NonAirRegion->GetPixel(currIndex) == 0) )
                {
                currForegroundMask->SetPixel(currIndex, 0);
                }
              else
                {
                count = count + 1;
                }
              }
            }
          }
        }
      muLogMacro(
        << "DEBUG:  Warped TemplateBrainMask size " << count << " of " << size[0] * size[1] * size[2] <<  std::endl);
      }
#if 1
      {
      // Dilate mask
      typedef itk::BinaryBallStructuringElement<unsigned char, 3> StructElementType;
      typedef
        itk::BinaryDilateImageFilter<ByteImageType, ByteImageType,
                                     StructElementType> DilateType;

      StructElementType structel;
      structel.SetRadius(1);
      structel.CreateStructuringElement();

      typename DilateType::Pointer dil = DilateType::New();
      dil->SetDilateValue(25);
      dil->SetKernel(structel);
      dil->SetInput(currForegroundMask);

      dil->Update();

      ByteImagePointer dilmask = dil->GetOutput();
      currForegroundMask = CopyImage<ByteImageType>(dilmask);
      }
#endif
    }
#endif

  if( this->m_DebugLevel > 5 )
    {                                       // DEBUG:  This code is for
                                            // debugging purposes only;
    unsigned int write_corrected_level = 0; // DEBUG:  This code is for
                                            // debugging purposes only;
    std::stringstream write_corrected_level_stream("");
    write_corrected_level_stream << write_corrected_level;
      { // DEBUG:  This code is for debugging purposes only;
      typedef itk::ImageFileWriter<ByteImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->UseCompressionOn();

      const std::string fn = this->m_OutputDebugDir + "/MASK_LEVEL_" + write_corrected_level_stream.str() + ".nii.gz";
      writer->SetInput(currForegroundMask);
      writer->SetFileName(fn.c_str() );
      writer->Update();
      muLogMacro( << "DEBUG:  Wrote image " << fn <<  std::endl );
      }
    }
  return currForegroundMask;
}
