#ifndef __itkICCIterativeInverseDisplacementFieldImageFilter_txx
#define __itkICCIterativeInverseDisplacementFieldImageFilter_txx

#include "itkICCIterativeInverseDisplacementFieldImageFilter.h"
#include "itkProgressReporter.h"
#include <vcl_algorithm.h>

namespace itk
{
// ----------------------------------------------------------------------------
// Constructor
template <class TInputImage, class TOutputImage>
ICCIterativeInverseDisplacementFieldImageFilter<TInputImage,
                                                TOutputImage>::ICCIterativeInverseDisplacementFieldImageFilter()
{
  m_NumberOfIterations = 10000;
  m_StopValue = 1e-7F;
  m_Time = 0;
}

// ----------------------------------------------------------------------------
template <class TInputImage, class TOutputImage>
void ICCIterativeInverseDisplacementFieldImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  TimeType time;

  time.Start(); // time measurement

  InputImageConstPointer inputPtr = this->GetInput(0);
  OutputImagePointer     outputPtr = this->GetOutput(0);

  // some checks
  if( inputPtr.IsNull() )
    {
    itkExceptionMacro("\n Input is missing.");
    }
  if( !TInputImage::ImageDimension == TOutputImage::ImageDimension )
    {
    itkExceptionMacro("\n Image Dimensions must be the same.");
    }

  outputPtr->SetRegions(inputPtr->GetRequestedRegion() );
  outputPtr->SetOrigin(inputPtr->GetOrigin() );
  outputPtr->SetSpacing(inputPtr->GetSpacing() );
  outputPtr->SetDirection(inputPtr->GetDirection() );
  outputPtr->Allocate();

  // Pertubation should be changed such that at the 100'th percentile the step size is less less than the epsilon

  this->ComputeInverse(inputPtr, outputPtr);

  time.Stop();
  m_Time = time.GetMean();
}

template <class TInputImage, class TOutputImage>
void
ICCIterativeInverseDisplacementFieldImageFilter<TInputImage, TOutputImage>
::ComputeInverse(InputImageConstPointer& inputPtr, OutputImagePointer& outputPtr)
{
  ThreadStruct str;

  str.Filter = this;
  str.inputPtr = inputPtr;
  str.outputPtr = outputPtr;

  this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads() );
  this->GetMultiThreader()->SetSingleMethod(this->ComputeInverseThreaderCallback, &str);
  this->GetMultiThreader()->SingleMethodExecute();
}

template <class TInputImage, class TOutputImage>
ITK_THREAD_RETURN_TYPE
ICCIterativeInverseDisplacementFieldImageFilter<TInputImage, TOutputImage>
::ComputeInverseThreaderCallback(void * arg)
{
  ThreadStruct * str;
  int            total, threadId, threadCount;

  threadId = ( (MultiThreader::ThreadInfoStruct *)(arg) )->ThreadID;
  threadCount = ( (MultiThreader::ThreadInfoStruct *)(arg) )->NumberOfThreads;

  str = (ThreadStruct *)( ( (MultiThreader::ThreadInfoStruct *)(arg) )->UserData);

  // Execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
  ThreadRegionType splitRegion;
  total = str->Filter->SplitRequestedRegion(threadId, threadCount,
                                            splitRegion);

  if( threadId < total )
    {
    str->Filter->ThreadedComputeInverse(str->inputPtr, str->outputPtr, splitRegion, threadId);
    }

  return ITK_THREAD_RETURN_VALUE;
}

template <class TInputImage, class TOutputImage>
void
ICCIterativeInverseDisplacementFieldImageFilter<TInputImage, TOutputImage>
::ThreadedComputeInverse(InputImageConstPointer& inputPtr, OutputImagePointer& outputPtr,
                         const ThreadRegionType & regionToProcess,
                         int)
{
  OutputIterator     OutputIt = OutputIterator(outputPtr, regionToProcess);
  InputConstIterator InputIt = InputConstIterator(inputPtr, regionToProcess);

  OutputImagePixelType outputValue;

  InputImageSizeType size = inputPtr->GetLargestPossibleRegion().GetSize();

  const float fnx = static_cast<float>(size[0]);
  const float fny = static_cast<float>(size[1]);
  const float fnz = static_cast<float>(size[2]);

  InputImageSpacingType spacing = inputPtr->GetSpacing();
  const float           fEpsilon1 = m_StopValue * spacing[0];
  const float           fEpsilon2 = m_StopValue * spacing[1];
  const float           fEpsilon3 = m_StopValue * spacing[2];

  unsigned int OnePercentIncrementOfMaxIterations = m_NumberOfIterations / 110; // NOTE:  We want to do 10% extra effort
  const float  MinEpsilon = vcl_min(fEpsilon1, vcl_min(fEpsilon2, fEpsilon3) );
  const float  PerturbationScaleFactor = 0.5 * vcl_pow(MinEpsilon, 0.01F);

  InputIt.GoToBegin();
  OutputIt.GoToBegin();
  while( !OutputIt.IsAtEnd() )
    {
    InputImageIndexType index = InputIt.GetIndex();
    const float         fBasePnt3 = static_cast<float>(index[2]) / fnz;
    const float         fBasePnt2 = static_cast<float>(index[1]) / fny;
    const float         fBasePnt1 = static_cast<float>(index[0]) / fnx;

    unsigned int iIteration = 0;
    float        fDelta1, fDelta2, fDelta3;
    float        fDispU1, fDispU2, fDispU3;

      {
      /* Make an inital fGuess of where the next ConstPixel should come from */
      float fDestPnt1 = fBasePnt1;
      float fDestPnt2 = fBasePnt2;
      float fDestPnt3 = fBasePnt3;

      unsigned int current_check_level = iIteration;
      float        perturbation = 1.0F;

      do
        {
        do
          {
            {
            float fullDestPnt1 = fnx * fBasePnt1;
            float fullDestPnt2 = fny * fBasePnt2;
            float fullDestPnt3 = fnz * fBasePnt3;

            typename InputImageType::PixelType pixel =
              TrilinearInterpolationFast(
                fullDestPnt1,
                fullDestPnt2,
                fullDestPnt3,
                size);
            fDispU1 = pixel[0];
            fDispU2 = pixel[1];
            fDispU3 = pixel[2];
            }
            {
            /*This is the actual fGuessU1(x) = pfU1+fDestPnt1 */
            const float fGuessU1 = fDispU1 + fDestPnt1;
            fDelta1 = fBasePnt1 - fGuessU1;
            fDestPnt1 += (fDelta1) * perturbation;
            /*This is the actual fGuessU2(x) = pfU2+fDestPnt2 */
            const float fGuessU2 = fDispU2 + fDestPnt2;
            fDelta2 = fBasePnt2 - fGuessU2;
            fDestPnt2 += (fDelta2) * perturbation;
            /*This is the actual fGuessU3(x) = pfU3+fDestPnt3 */
            const float fGuessU3 = fDispU3 + fDestPnt3;
            fDelta3 = fBasePnt3 - fGuessU3;
            fDestPnt3 += (fDelta3) * perturbation;
            }
          iIteration++;
          }
        while( ( (vcl_abs(fDelta1) > fEpsilon1) || (vcl_abs(fDelta2) > fEpsilon2) ||
                 (vcl_abs(fDelta3) > fEpsilon3) ) && (iIteration < current_check_level) );

        perturbation *= PerturbationScaleFactor;  // Reduce stepsize in attemp to get out of local minimum. // 0.5
                                                  // decreases too quickly
        current_check_level = current_check_level + OnePercentIncrementOfMaxIterations;
        }
      while( ( (vcl_abs(fDelta1) > fEpsilon1) || (vcl_abs(fDelta2) > fEpsilon2) ||
               (vcl_abs(fDelta3) > fEpsilon3) ) && (iIteration < m_NumberOfIterations) );
      }

    outputValue[0] = -fDispU1;
    outputValue[1] = -fDispU2;
    outputValue[2] = -fDispU3;

    OutputIt.Set( outputValue );

    ++InputIt;
    ++OutputIt;
    }   // end while loop
}

template <class TInputImage, class TOutputImage>
typename ICCIterativeInverseDisplacementFieldImageFilter<TInputImage, TOutputImage>
::OutputImagePixelType
ICCIterativeInverseDisplacementFieldImageFilter<TInputImage, TOutputImage>
::TrilinearInterpolationFast(float& fDesiredX, float& fDesiredY, float& fDesiredZ, InputImageSizeType size)
{
#if defined(__sgi)
  const float fXtemp = floorf( fDesiredX );
  const float fYtemp = floorf( fDesiredY );
  const float fZtemp = floorf( fDesiredZ );
#else
  const float fXtemp = rintf( fDesiredX - 0.5F);
  const float fYtemp = rintf( fDesiredY - 0.5F);
  const float fZtemp = rintf( fDesiredZ - 0.5F);
#endif
    {
    const int iFlooredX = static_cast<int>(fXtemp);
    const int iFlooredY = static_cast<int>(fYtemp);
    const int iFlooredZ = static_cast<int>(fZtemp);

    const int nx = size[0];
    const int ny = size[1];
    const int nz = size[2];
    /*********************/
    /*Start Interpolation */
    /*********************/
    /* If in the Middle of the Image */
    if( ( iFlooredX > 0 ) && ( iFlooredX < nx - 1 ) && ( iFlooredY > 0 )
        && ( iFlooredY < ny - 1 ) && ( iFlooredZ > 0 ) && ( iFlooredZ < nz - 1 ) )
      {
      const float         fXDelta = fDesiredX - fXtemp;
      InputImageIndexType a000 = {{0, 0, 0}};
      a000[0] = iFlooredX;
      a000[1] = iFlooredY;
      a000[2] =  iFlooredZ;
      const PixelType     a000_eval = this->GetInput(0)->GetPixel(a000);
      InputImageIndexType a001;
      a001[0] = a000[0] + 1;
      a001[1] = a000[1];
      a001[2] =  a000[2];
      const PixelType     a001_eval = this->GetInput(0)->GetPixel(a001);
      const PixelType     b00_eval = fXDelta * ( a001_eval -  a000_eval  ) + ( a000_eval );
      InputImageIndexType a010;
      a010[0] = a000[0];
      a010[1] = a000[1] + 1;
      a010[2] =  a000[2];
      const PixelType     a010_eval = this->GetInput(0)->GetPixel(a010);
      InputImageIndexType a011;
      a011[0] = a010[0] + 1;
      a011[1] = a010[1];
      a011[2] =  a010[2];
      const PixelType     b01_eval = fXDelta * ( this->GetInput(0)->GetPixel(a011) -  a010_eval  ) + ( a010_eval );
      const float         fYDelta = fDesiredY - fYtemp;
      const PixelType     c0_eval = fYDelta * ( b01_eval - b00_eval ) + b00_eval;
      InputImageIndexType a100;
      a100[0] = a000[0];
      a100[1] = a000[1];
      a100[2] =  a000[2] + 1;
      const PixelType     a100_eval = this->GetInput(0)->GetPixel(a100);
      InputImageIndexType a101;
      a101[0] = a100[0] + 1;
      a101[1] = a100[1];
      a101[2] =  a100[2];
      const PixelType     b10_eval = fXDelta * ( this->GetInput(0)->GetPixel(a101) -  a100_eval  ) + ( a100_eval );
      InputImageIndexType a110;
      a110[0] = a010[0];
      a110[1] = a010[1];
      a110[2] =  a010[2] + 1;
      const PixelType     a110_eval = this->GetInput(0)->GetPixel(a110);
      InputImageIndexType a111;
      a111[0] = a011[0];
      a111[1] = a011[1];
      a111[2] =  a011[2] + 1;
      const PixelType b11_eval = fXDelta * ( this->GetInput(0)->GetPixel(a111) -  a110_eval  ) + ( a110_eval );
      const PixelType c1_eval = fYDelta * ( b11_eval - b10_eval ) + b10_eval;
      const float     fZDelta = fDesiredZ - fZtemp;
      return fZDelta * ( c1_eval - c0_eval ) + c0_eval;
      }

    else
      {
      /* If on a boundary */
      const float fXDelta = fDesiredX - fXtemp;
      const int   iNextX = iFlooredX + 1;
      const int   iNextY = iFlooredY + 1;
      const int   iNextZ = iFlooredZ + 1;
      /*Circular Boundaries */
      InputImageIndexType index;
      index = BoundaryIndexing(iFlooredX, iFlooredY, iFlooredZ, nx, ny, nz );
      const PixelType a000_eval = this->GetInput(0)->GetPixel(index);
      index = BoundaryIndexing(iNextX, iFlooredY, iFlooredZ, nx, ny, nz );
      const PixelType a001_eval = this->GetInput(0)->GetPixel(index);
      const PixelType b00_eval = fXDelta * ( ( a001_eval ) - ( a000_eval ) ) + ( a000_eval );
      index = BoundaryIndexing(iFlooredX, iNextY, iFlooredZ, nx, ny, nz );
      const PixelType a010_eval = this->GetInput(0)->GetPixel(index);
      index = BoundaryIndexing(iNextX, iNextY, iFlooredZ, nx, ny, nz );
      const PixelType a011_eval = this->GetInput(0)->GetPixel(index);
      const PixelType b01_eval = fXDelta * ( ( a011_eval ) - ( a010_eval ) ) + ( a010_eval );
      const float     fYDelta = fDesiredY - fYtemp;
      const PixelType c0_eval = fYDelta * ( b01_eval - b00_eval ) + b00_eval;
      index = BoundaryIndexing(iFlooredX, iFlooredY, iNextZ, nx, ny, nz );
      const PixelType a100_eval = this->GetInput(0)->GetPixel(index);
      index = BoundaryIndexing(iNextX, iFlooredY, iNextZ, nx, ny, nz );
      const PixelType a101_eval = this->GetInput(0)->GetPixel(index);
      const PixelType b10_eval = fXDelta * ( ( a101_eval ) - ( a100_eval ) ) + ( a100_eval );
      index = BoundaryIndexing(iFlooredX, iNextY, iNextZ, nx, ny, nz );
      const PixelType a110_eval = this->GetInput(0)->GetPixel(index);
      index = BoundaryIndexing(iNextX, iNextY, iNextZ, nx, ny, nz );
      const PixelType a111_eval = this->GetInput(0)->GetPixel(index);
      const PixelType b11_eval = fXDelta * ( ( a111_eval ) - ( a110_eval ) ) + ( a110_eval );
      const PixelType c1_eval = fYDelta * ( b11_eval - b10_eval ) + b10_eval;
      const float     fZDelta = fDesiredZ - fZtemp;
      return fZDelta * ( c1_eval - c0_eval ) + c0_eval;
      }
    }
}

template <class TInputImage, class TOutputImage>
typename ICCIterativeInverseDisplacementFieldImageFilter<TInputImage, TOutputImage>
::InputImageIndexType
ICCIterativeInverseDisplacementFieldImageFilter<TInputImage, TOutputImage>
::BoundaryIndexing(int iDesiredX,
                   int iDesiredY,
                   int iDesiredZ,
                   const int iMaxX,
                   const int iMaxY,
                   const int iMaxZ )
{
  iDesiredZ = ( ( iDesiredZ > -1 ) ?
                ( ( iDesiredZ < iMaxZ ) ? iDesiredZ : ( iMaxZ - 1 ) ) : ( 0 ) );
  iDesiredY = ( ( iDesiredY > -1 ) ?
                ( ( iDesiredY < iMaxY ) ? iDesiredY : ( iMaxY - 1 ) ) : ( 0 ) );
  iDesiredX =  ( (iDesiredX > -1 ) ?
                 ( ( iDesiredX < iMaxX ) ? iDesiredX : ( iMaxX - 1 ) ) : ( 0 ) );
  InputImageIndexType rval;
  rval[0] = iDesiredX;
  rval[1] = iDesiredY;
  rval[2] = iDesiredZ;
  return rval;
}

// ----------------------------------------------------------------------------
template <class TInputImage, class TOutputImage>
void ICCIterativeInverseDisplacementFieldImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Number of iterations: " << m_NumberOfIterations << std::endl;
  os << indent << "Stop value:           " << m_StopValue << " mm" << std::endl;
  os << indent << "Elapsed time:         " << m_Time << " sec" << std::endl;
  os << std::endl;
}
} // end namespace itk

#endif
