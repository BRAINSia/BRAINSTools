#if 0
#ifndef __EMSplitPriorMST_h
#define __EMSplitPriorMST_h
template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter<TInputImage, TProbabilityImage>
::SplitPriorMST(unsigned int /*iprior*/)
{
  exit(-1);
#if 0
  const unsigned int numChannels = m_InputImages.size();
  const unsigned int numPriors = m_WarpedPriors.size();

  if( iprior >= numPriors )
    {
    itkExceptionMacro(<< "Invalid prior index");
    }

  if( m_PriorGaussianClusterCountVector[iprior] == 1 )
    {
    return; // No split necessary
    }

  unsigned int numClasses = 0;
  for( unsigned int i = 0; i < numPriors; i++ )
    {
    numClasses += m_PriorGaussianClusterCountVector[i];
    }

  muLogMacro(<< "Splitting distributions for prior " << iprior + 1 << " (MST)\n");

  InputImageSizeType size
    = m_InputImages[0]->GetLargestPossibleRegion().GetSize();

  ProbabilityImagePointer probImg = m_WarpedPriors[iprior];

  // Find max prob
  double maxP = 0.0;
    {
    // #pragma omp parallel for
    for( long kk = 0; kk < (long)size[2]; kk++ )
      {
      for( long jj = 0; jj < (long)size[1]; jj++ )
        {
        for( long ii = 0; ii < (long)size[0]; ii++ )
          {
          const ProbabilityImageIndexType currIndex = {{ii, jj, kk}};
          if( m_CurrentIterationForegroundMask->GetPixel(currIndex) == 0 )
            {
            continue;
            }
          if( probImg->GetPixel(currIndex) > maxP )
            {
            maxP = probImg->GetPixel(currIndex);
            }
          }
        }
      }
    }

  // Select uniqueClassAverageSamples by thresholding prior with value above tau
  const double tau = 0.85 * maxP;
  muLogMacro(<< "Sampling with threshold tau = " << tau << "\n");

  unsigned int numPossibleSamples = 0;
    {
    // #pragma omp parallel for
    for( long kk = 0; kk < (long)size[2]; kk++ )
      {
      for( long jj = 0; jj < (long)size[1]; jj++ )
        {
        for( long ii = 0; ii < (long)size[0]; ii++ )
          {
          const ProbabilityImageIndexType currIndex = {{ii, jj, kk}};
          if( m_CurrentIterationForegroundMask->GetPixel(currIndex) == 0 )
            {
            continue;
            }
          if( probImg->GetPixel(currIndex) >= tau )
            {
            numPossibleSamples++;
            }
          }
        }
      }
    }

  // Make a mapping of unique type names to a list of corresponding corrected
  // images indicies
  // The uniqueClassAverageSamples should only have dimension of the number of
  // unique class types,
  // and not the total number of images being corrected.
  std::map<std::string, std::list<unsigned int> > uniqueClassMappings;
    {
    for( unsigned int u = 0; u < m_InputVolumeTypes.size(); u++ )
      {
      uniqueClassMappings[m_InputVolumeTypes[u]].push_back(u);
      }
    }
  // Sample selection mask
  std::vector<QHullMSTClusteringProcess::VertexType> uniqueClassAverageSamples;
  std::vector<InputImageIndexType>                   uniqueClassAverageIndicies;
  size_t                                             numSamples = numPossibleSamples;
    {
    std::vector<unsigned char> selectMask(numPossibleSamples);
    for( unsigned int i = 0; i < numPossibleSamples; i++ )
      {
      selectMask[i] = 0;
      }

    if( numSamples > 50000 )
      {
      numSamples = 50000;
      }

    muLogMacro(<< "  Selecting " << numSamples << " / " << numPossibleSamples << "\n");

    itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer
                        rng( itk::Statistics::MersenneTwisterRandomVariateGenerator::New() );
    static unsigned int rngInitializer = 0;
    rng->Initialize(rngInitializer);
    rngInitializer++;

    if( numSamples < numPossibleSamples )
      {
      unsigned int c = 0;
      while( c < numSamples )
        {
        const unsigned int which = (unsigned int)
          rng->GetIntegerVariate(numPossibleSamples - 1);
        if( selectMask[which] != 0 )
          {
          continue;
          }
        selectMask[which] = 1;
        c++;
        }
      }
    else
      {
      for( unsigned int i = 0; i < numSamples; i++ )
        {
        selectMask[i] = 1;
        }
      }

    uniqueClassAverageSamples.clear();
    uniqueClassAverageSamples.reserve(numSamples);
    uniqueClassAverageIndicies.clear();
    uniqueClassAverageIndicies.reserve(numSamples);

    muLogMacro(<< "  Finding uniqueClassAverageSamples...\n");

    // Compute the muAccumulator
      {
      const size_t numUniqueClassMappings = uniqueClassMappings.size();
      unsigned int selectMaskIndex = 0;
        {
        // #pragma omp parallel for
        for( long kk = 0; kk < (long)size[2]; kk++ )
          {
          for( long jj = 0; jj < (long)size[1]; jj++ )
            {
            for( long ii = 0; ii < (long)size[0]; ii++ )
              {
              const ProbabilityImageIndexType currIndex = {{ii, jj, kk}};
              if( ( m_CurrentIterationForegroundMask->GetPixel(currIndex) == 0 ) ||
                  ( probImg->GetPixel(currIndex) < tau ) )
                {
                continue;
                }

              if( selectMask[selectMaskIndex] != 0 )
                {
                QHullMSTClusteringProcess::VertexType x(numUniqueClassMappings);
                unsigned int                          uniqueIndex = 0;
                for( std::map<std::string, std::list<unsigned int> >::const_iterator mi = uniqueClassMappings.begin();
                     mi != uniqueClassMappings.end();
                     mi++ )
                  {
                  const std::list<unsigned int> & currentListOfImages = mi->second;
                  const double                    invCurrentListSize = 1.0
                    / ( static_cast<double>( currentListOfImages.size() ) );
                  x[uniqueIndex] = 0.0;
                  for( std::list<unsigned int>::const_iterator li = currentListOfImages.begin();
                       li != currentListOfImages.end();
                       li++
                       )
                    {
                    const typename TInputImage::PixelType currPixel = this->m_CorrectedImages[*li]->GetPixel(currIndex);
                    x[uniqueIndex] += currPixel * invCurrentListSize;
                    }
                  uniqueIndex++;
                  }
                uniqueClassAverageSamples.push_back(x);
                uniqueClassAverageIndicies.push_back(currIndex);
                }
              ++selectMaskIndex;
              }
            }
          }
        }
      }
    }

  muLogMacro(<< "  Create MST...\n");

  QHullMSTClusteringProcess mstProc;
  mstProc.SetInputVertices(uniqueClassAverageSamples);
  mstProc.SortOn();

  muLogMacro(<< "  Allocate maps...\n");

  // TODO: Make this an std::vector
  unsigned int *mapOfSamplesIntoCluster = new unsigned int[numSamples];
    {
    for( double T = 2.0; T > 0; T -= 0.01 )
      {
      if( T < 1.0 )
        {
        itkExceptionMacro(<< "Failed clustering prior " << iprior + 1);
        }
      muLogMacro(<< "  MST clustering, T = " << T << "\n");
      const unsigned int numClusters = mstProc.GetClusters(mapOfSamplesIntoCluster, T);

      if( numClusters < m_PriorGaussianClusterCountVector[iprior] )
        {
        continue;
        }

      // Check cluster sizes
      bool sizeOK = true;
      for( unsigned int currentPriorSubClassCluster = 0;
           currentPriorSubClassCluster < m_PriorGaussianClusterCountVector[iprior];
           currentPriorSubClassCluster++ )
        {
        unsigned int numSamplesInSubClassCluster = 0;
        for( unsigned int i = 0; i < numSamples; i++ )
          {
          if( mapOfSamplesIntoCluster[i] == currentPriorSubClassCluster )
            {
            numSamplesInSubClassCluster++;
            }
          }
        if( numSamplesInSubClassCluster < 10 ) // TODO:  Determine what this
        // does?
          {
          sizeOK = false;
          }
        }
      if( sizeOK )
        {
        break;
        }
      }
    }
  // TODO:  Should be able to clear uniqueClassAverageSamples here.
  // =============Clustering Done

  // Find beginning class index for this prior
  unsigned int istart = 0;
  for( unsigned int j = 0; j < iprior; j++ )
    {
    istart += m_PriorGaussianClusterCountVector[j];
    }
  // Use the largest clusters to estimate the Gaussian parameters
  for( unsigned int currentPriorSubClassCluster = 0;
       currentPriorSubClassCluster < m_PriorGaussianClusterCountVector[iprior];
       currentPriorSubClassCluster++ )
    {
    const unsigned int iclass = istart + currentPriorSubClassCluster;

    unsigned int numSamplesInSubClassCluster = 0;
    for( unsigned int i = 0; i < numSamples; i++ )
      {
      if( mapOfSamplesIntoCluster[i] == currentPriorSubClassCluster )
        {
        numSamplesInSubClassCluster++;
        }
      }

    muLogMacro(<< " Estimating Gaussian parameters for prior " << iprior + 1 << ":class " << iclass + 1);
    muLogMacro(<< " with " << numSamplesInSubClassCluster << " uniqueClassAverageSamples\n");

    QHullMSTClusteringProcess::VertexType muAccumulator( m_InputVolumeTypes.size() );
    std::fill( muAccumulator.begin(), muAccumulator.end(), 0.0F );
    for( unsigned int ichan = 0; ichan < numChannels; ichan++ )
      {
      muAccumulator[ichan] = 0.0;
      for( unsigned int i = 0; i < numSamples; i++ )
        {
        if( mapOfSamplesIntoCluster[i] == currentPriorSubClassCluster )
          {
          const typename TInputImage::PixelType currPixel
            = this->m_CorrectedImages[ichan]->GetPixel(uniqueClassAverageIndicies[i]);
          muAccumulator[ichan] += currPixel;
          }
        }
      muAccumulator[ichan] /= static_cast<double>( numSamplesInSubClassCluster );
      this->m_ListOfClassStatistics[iclass].m_Means[ichan] = muAccumulator[ichan];
      }
    // Compute symmetric covariance matrix.
    MatrixType tmpSubClassSampleCovarianceBetweenChannels(numChannels, numChannels);
    for( unsigned int channel1Index = 0; channel1Index < numChannels; channel1Index++ )
      {
      const double channel1mu = this->m_ListOfClassStatistics[iclass].m_Means[channel1Index];
      for( unsigned int channel2Index = channel1Index; channel2Index < numChannels; channel2Index++ )
        {
        const double channel2mu = this->m_ListOfClassStatistics[iclass].m_Means[channel2Index];
        double       tmpSubClassClusterCovarianceBetweenChannels12 = 0.0;
        for( unsigned int i = 0; i < numSamples; i++ )
          {
          if( mapOfSamplesIntoCluster[i] == currentPriorSubClassCluster )
            {
            const typename TInputImage::PixelType channel1Pixel
              = this->m_CorrectedImages[channel1Index]->GetPixel(uniqueClassAverageIndicies[i]);
            const double diff1 = static_cast<double>( channel1Pixel ) - channel1mu;

            const typename TInputImage::PixelType channel2Pixel
              = this->m_CorrectedImages[channel2Index]->GetPixel(uniqueClassAverageIndicies[i]);
            const double diff2 = static_cast<double>( channel2Pixel ) - channel2mu;
            tmpSubClassClusterCovarianceBetweenChannels12 += ( diff1 * diff2 );
            }
          }
        tmpSubClassClusterCovarianceBetweenChannels12
          /= ( static_cast<double>( numSamplesInSubClassCluster ) - 1.0 + vnl_math::eps );
        if( channel1Index == channel2Index )
          {
          tmpSubClassClusterCovarianceBetweenChannels12 += 1e-20;
          }

        tmpSubClassSampleCovarianceBetweenChannels(channel1Index,
                                                   channel2Index) = tmpSubClassClusterCovarianceBetweenChannels12;
        tmpSubClassSampleCovarianceBetweenChannels(channel2Index,
                                                   channel1Index) = tmpSubClassClusterCovarianceBetweenChannels12;
        }
      }
    // Scale the covariance up a bit for softer initialization
    tmpSubClassSampleCovarianceBetweenChannels *= 1.2;
    this->m_ListOfClassStatistics[iclass].m_Covariance = tmpSubClassSampleCovarianceBetweenChannels;
    }
  for( unsigned int iclass = 0; iclass < numClasses; iclass++ )
    {
    const double detcov = vnl_determinant(this->m_ListOfClassStatistics[iclass].m_Covariance);
    if( detcov <= 0.0 )
      {
      itkExceptionMacro(<< "Determinant of covariance for class " << iclass
                        << " is <= 0.0 (" << detcov << "), covariance matrix:\n"
                        << this->m_ListOfClassStatistics[iclass].m_Covariance);
      }
    }

  uniqueClassAverageSamples.clear();
  uniqueClassAverageIndicies.clear(); // Free the memory.

  delete[] mapOfSamplesIntoCluster;

  // Special case for background, always set the "darkest" mean
  // to be the last and set it to zero
  if( iprior == ( numPriors - 1 ) )
    {
    unsigned int imin = istart;
    double       minm = this->m_ListOfClassStatistics[istart].m_Means.squared_magnitude();
    for( unsigned int currentPriorSubClassCluster = 1;
         currentPriorSubClassCluster < m_PriorGaussianClusterCountVector[iprior];
         currentPriorSubClassCluster++ )
      {
      unsigned int iclass = istart + currentPriorSubClassCluster;
      double       mag = this->m_ListOfClassStatistics[iclass].m_Means.squared_magnitude();
      if( mag < minm )
        {
        imin = iclass;
        minm = mag;
        }
      }

    if( imin != ( numClasses - 1 ) )
      {
      muLogMacro(
        << "  Replacing " << this->m_ListOfClassStatistics[imin].m_Means << " with zero\n");
      VectorType v = this->m_ListOfClassStatistics[numClasses - 1].m_Means;
      this->m_ListOfClassStatistics[imin].m_Means = v;
      }
    for( unsigned int ichan = 0; ichan < numChannels; ichan++ )
      {
      this->m_ListOfClassStatistics[numClasses - 1].m_Means[ichan] = 0;
      }
    }
#endif
}

#endif
#endif
