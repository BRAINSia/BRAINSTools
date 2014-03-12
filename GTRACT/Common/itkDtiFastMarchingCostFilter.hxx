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
/*=========================================================================

 Program:   GTRACT (Guided Tensor Restore Anatomical Connectivity Tractography)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2006/03/29 14:53:40 $
 Version:   $Revision: 1.9 $

   Copyright (c) University of Iowa Department of Radiology. All rights reserved.
   See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
   for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkDtiFastMarchingCostFilter_hxx
#define __itkDtiFastMarchingCostFilter_hxx

#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkNeighborhoodAlgorithm.h>
#include <itkImageToImageFilter.h>
#include <itkIOCommon.h>
#include <itkVector.h>
#include <itkListSample.h>
#include <vnl/vnl_vector.h>
#include <itkFixedArray.h>
#include "itkMetaDataObject.h"
#include <itkDiffusionTensor3D.h>

#include "itkDtiFastMarchingCostFilter.h"

#include <iostream>

#include "itkNumericTraits.h"
#include "vnl/vnl_math.h"
#include <algorithm>

namespace itk
{
/*
 *
 */
template <class TLevelSet, class TTensorImage>
DtiFastMarchingCostFilter<TLevelSet, TTensorImage>
::DtiFastMarchingCostFilter() :
  m_TrialHeap()
{
  this->ProcessObject::SetNumberOfRequiredInputs(0);

  OutputSizeType outputSize;
  outputSize.Fill( 16 );
  typename LevelSetImageType::IndexType outputIndex;
  outputIndex.Fill( 0 );

  m_OutputRegion.SetSize( outputSize );
  m_OutputRegion.SetIndex( outputIndex );

  m_OutputOrigin.Fill( 0.0 );
  m_OutputSpacing.Fill( 1.0 );
  m_OverrideOutputInformation = false;

  m_AlivePoints = NULL;
  m_TrialPoints = NULL;
  m_ProcessedPoints = NULL;

  m_LabelImage = LabelImageType::New();
  m_AnisotropyImage = AnisotropyImageType::New();
  m_OutputSpeedImage = OutputSpeedImageType::New();
  m_EigenvectorImage =  EigenvectorImageType::New();

  typedef typename LevelSetImageType::PixelType PixType;
  m_LargeValue    = static_cast<PixType>( NumericTraits<PixType>::max() / 2.0 );
  // m_LargeValue    = static_cast<PixelType>( 10000.0 );
  m_StoppingValue = static_cast<double>( m_LargeValue );
  m_AnisotropyWeight = static_cast<float>( 0.0);
  m_CollectPoints = false;
  m_NormalizationFactor = 1.0;
}

/*
 *
 */
template <class TLevelSet, class TTensorImage>
void
DtiFastMarchingCostFilter<TLevelSet, TTensorImage>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Alive points: " << m_AlivePoints.GetPointer() << std::endl;
  os << indent << "Trial points: " << m_TrialPoints.GetPointer() << std::endl;
  os << indent << "Stopping value: " << m_StoppingValue << std::endl;
  os << indent << "Large Value: "
     << static_cast<typename NumericTraits<PixelType>::PrintType>(m_LargeValue)
     << std::endl;
  os << indent << "Normalization Factor: " << m_NormalizationFactor << std::endl;
  os << indent << "Collect points: " << m_CollectPoints << std::endl;
  os << indent << "OverrideOutputInformation: ";
  os << m_OverrideOutputInformation << std::endl;
  os << indent << "OutputRegion: " << m_OutputRegion << std::endl;
  os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;
  os << indent << "OutputOrigin:  " << m_OutputOrigin << std::endl;
  os << indent << "OutputDirection:  " << m_OutputDirection << std::endl;
}

/*
 *
 */
template <class TLevelSet, class TTensorImage>
void
DtiFastMarchingCostFilter<TLevelSet, TTensorImage>
::GenerateOutputInformation()
{
  // copy output information from input image
  Superclass::GenerateOutputInformation();

  // use user-specified output information

  if( this->GetInput() == NULL || m_OverrideOutputInformation )
    {
    LevelSetPointer output = this->GetOutput();
    output->SetLargestPossibleRegion( m_OutputRegion );
    output->SetSpacing( m_OutputSpacing );
    output->SetOrigin( m_OutputOrigin );
    output->SetDirection( m_OutputDirection );
    }
}

/*
 *
 */
template <class TLevelSet, class TTensorImage>
void
DtiFastMarchingCostFilter<TLevelSet, TTensorImage>
::EnlargeOutputRequestedRegion(
  DataObject *output )
{
  // enlarge the requested region of the output
  // to the whole data set
  TLevelSet *imgData;

  imgData = dynamic_cast<TLevelSet *>( output );
  if( imgData )
    {
    imgData->SetRequestedRegionToLargestPossibleRegion();
    }
  else
    {
    // Pointer could not be cast to TLevelSet *
    itkWarningMacro( << "itk::FastMarchingImageFilter"
                     << "::EnlargeOutputRequestedRegion cannot cast "
                     << typeid( output ).name() << " to "
                     << typeid( TLevelSet * ).name() );
    }
}

/*
 *
 */
template <class TLevelSet, class TTensorImage>
void DtiFastMarchingCostFilter<TLevelSet, TTensorImage>
// ::Initialize( const TensorImageType * tensorImage, LevelSetImageType * output
// )
::GenerateData()
{
  LevelSetPointer         output      = this->GetOutput();
  TensorImageConstPointer tensorImage = this->GetInput();

  if( m_CollectPoints )
    {
    m_ProcessedPoints = NodeContainer::New();
    }

  output->SetBufferedRegion( output->GetRequestedRegion() );

  // Set Meta Data Orientation Information
  output->SetMetaDataDictionary( tensorImage->GetMetaDataDictionary() );
  output->Allocate();

  // cache some buffered region information
  m_BufferedRegion = output->GetBufferedRegion();
  m_StartIndex = m_BufferedRegion.GetIndex();
  m_LastIndex = m_StartIndex + m_BufferedRegion.GetSize();
  typename LevelSetImageType::OffsetType offset;
  offset.Fill( 1 );
  m_LastIndex -= offset;

  // allocate memory for the PointTypeImage
  m_LabelImage->CopyInformation( output );
  m_LabelImage->SetBufferedRegion(
    output->GetBufferedRegion() );
  m_LabelImage->Allocate();

  // allocate memory for OutputSpeedImage
  m_OutputSpeedImage->CopyInformation( output );
  m_OutputSpeedImage->SetBufferedRegion(
    output->GetBufferedRegion() );
  m_OutputSpeedImage->SetMetaDataDictionary( output->GetMetaDataDictionary() );
  m_OutputSpeedImage->Allocate();

  /*Read Tensor Image and set principal eigenvector image*/
  m_EigenvectorImage->SetRegions( tensorImage->GetLargestPossibleRegion() );
  m_EigenvectorImage->SetSpacing( tensorImage->GetSpacing() );
  m_EigenvectorImage->SetOrigin( tensorImage->GetOrigin() );
  m_EigenvectorImage->SetDirection( tensorImage->GetDirection() );
  m_EigenvectorImage->Allocate();

  typedef itk::ImageRegionIterator<EigenvectorImageType> EigIteratorType;
  EigIteratorType eigIt( m_EigenvectorImage, m_EigenvectorImage->GetRequestedRegion() );

  typedef itk::ImageRegionConstIterator<TensorImageType> TensorImageIteratorType;
  TensorImageIteratorType tensorIt( tensorImage, tensorImage->GetRequestedRegion() );
  for( eigIt.GoToBegin(), tensorIt.GoToBegin(); !eigIt.IsAtEnd() && !tensorIt.IsAtEnd(); ++eigIt, ++tensorIt )
    {
    EigenvectorPixelType principalEigenvector;
    TensorImageIndexType tensorIndex = tensorIt.GetIndex();

    // Define type of tensor pixel
    const unsigned int tensElements = 6;
    typedef itk::Vector<float, tensElements> VectorTensorPixelType;
    VectorTensorPixelType tensor;

    TensorImagePixelType tensorPixel =  tensorImage->GetPixel(tensorIndex);
    for( unsigned int i = 0; i < tensElements; i++ )
      {
      tensor[i] = tensorPixel[i];
      }

    if( tensor.GetNorm() != 0 )
      {
      typename TensorImagePixelType::EigenValuesArrayType   eigenValues;
      typename TensorImagePixelType::EigenVectorsMatrixType eigenVectors;
      tensorPixel.ComputeEigenAnalysis(eigenValues, eigenVectors);
      for( unsigned int i = 0; i < dimension; i++ )
        {
        principalEigenvector[i] = eigenVectors[( dimension - 1 )][i];
        }
      }

    else
      {
      for( unsigned int i = 0; i < dimension; i++ )
        {
        principalEigenvector[i] = 0.0;
        }
      }

    eigIt.Set( principalEigenvector );
    }

  // set all output value to Large Value
  typedef ImageRegionIterator<LevelSetImageType>
    OutputIterator;

  OutputIterator outIt( output, output->GetBufferedRegion() );

  PixelType outputPixel;
  outputPixel = m_LargeValue;
  for( outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt )
    {
    outIt.Set( outputPixel );
    }

  // set all speed output values to  0
  typedef ImageRegionIterator<OutputSpeedImageType> OutputSpeedIterator;

  OutputSpeedIterator speedIt( m_OutputSpeedImage, m_OutputSpeedImage->GetBufferedRegion() );

  PixelType outputSpeedPixel;
  outputSpeedPixel = 0.0;
  for( speedIt.GoToBegin(); !speedIt.IsAtEnd(); ++speedIt )
    {
    speedIt.Set( outputSpeedPixel );
    }

  typedef ImageRegionIterator<LabelImageType> LabelIterator;
  LabelIterator typeIt( m_LabelImage, m_LabelImage->GetBufferedRegion() );

  // process input alive points
  AxisNodeType node;

  if( m_AlivePoints )
    {
    typename NodeContainer::ConstIterator pointsIter;
    pointsIter = m_AlivePoints->Begin();
    typename NodeContainer::ConstIterator pointsEnd = m_AlivePoints->End();
    for( ; pointsIter != pointsEnd; ++pointsIter )
      {
      // get node from alive points container
      node = pointsIter.Value();

      // check if node index is within the output level set
      if( !m_BufferedRegion.IsInside( node.GetIndex() ) )
        {
        continue;
        }

      // check if node is valid and not outside of brain
      EigenvectorPixelType eigPixel = m_EigenvectorImage->GetPixel( node.GetIndex() );

      float checkEigValue = 0.0;
      bool  pass = false;
      for( unsigned int i = 0; i < dimension; i++ )
        {
        checkEigValue = eigPixel[i];
        if( ( checkEigValue > 0.0 ) || ( checkEigValue < 0.0 ) )
          {
          pass = true;
          }
        }

      if( !pass )
        {
        std::cout << "Warning seed point is outside brain region: " <<  node.GetIndex() << std::endl;
        continue;  // out of brain region
        }
      // set all points type to FarPoint
      for( typeIt.GoToBegin(); !typeIt.IsAtEnd(); ++typeIt )
        {
        typeIt.Set( FarPoint );
        }

      // make sure the heap is empty
      while( !m_TrialHeap.empty() )
        {
        m_TrialHeap.pop();
        }

      // make this an alive point
      m_LabelImage->SetPixel( node.GetIndex(), AlivePoint );

      outputPixel = node.GetValue();
      output->SetPixel( node.GetIndex(), outputPixel );

      // Set speed of initial seed points to F=|eigenvector(AlivePoint)|=1.0
      outputSpeedPixel = 1.0;

      /*Scale by FA*/
      AnisotropyImagePixelType aniso = 0.0;

      if( m_AnisotropyWeight > 0 )
        {
        aniso =  m_AnisotropyImage->GetPixel( node.GetIndex() );
        }
      outputSpeedPixel = outputSpeedPixel * ( 1 - m_AnisotropyWeight ) + aniso * m_AnisotropyWeight;

      // Normalize Speed
      PixelType normOutputSpeedPixel = outputSpeedPixel / m_NormalizationFactor;

      m_OutputSpeedImage->SetPixel( node.GetIndex(),  normOutputSpeedPixel );
      this->InitializeTrialPoints( /* tensorImage, */ node.GetIndex(), output );
      // this->UpdateNeighbors( tensorImage, node.GetIndex(), output );
      this->UpdateFront( /* tensorImage,*/ output );
      }
    }
}

/*
 *
 */
template <class TLevelSet, class TTensorImage>
void DtiFastMarchingCostFilter<TLevelSet, TTensorImage>
::UpdateFront( /* const TensorImageType * tensorImage, */ LevelSetImageType *output )
{
  // process points on the heap
  AxisNodeType node;
  double       currentValue;
  double       oldProgress = 0;

  this->UpdateProgress( 0.0 ); // Send first progress event

  while( !m_TrialHeap.empty() )
    {
    // get the node with the smallest value
    node = m_TrialHeap.top();
    m_TrialHeap.pop();

    // check if node index is within the output level set
    if( !m_BufferedRegion.IsInside( node.GetIndex() ) )
      {
      continue;
      }

    // check if node index is within brain region
    // const EigenvectorPixelType eigPixel = m_EigenvectorImage->GetPixel( node.GetIndex() );

    // does this node contain the current value ?
    currentValue = (double)output->GetPixel( node.GetIndex() );

    if( node.GetValue() != currentValue )
      {
      continue;
      }

    // is this node already alive ?
    if( m_LabelImage->GetPixel( node.GetIndex() ) != TrialPoint )
      {
      continue;
      }

    if( currentValue > m_StoppingValue )
      {
      break;
      }

    if( m_CollectPoints )
      {
      m_ProcessedPoints->InsertElement( m_ProcessedPoints->Size(), node );
      }

    // set this node as alive
    m_LabelImage->SetPixel( node.GetIndex(), AlivePoint );

    // update its neighbors
    this->UpdateNeighbors( /* tensorImage, */ node.GetIndex(), output );

    // Send events every certain number of points.
    const double newProgress = currentValue / m_StoppingValue;
    if( newProgress - oldProgress > 0.01 )   // update every 1%
      {
      this->UpdateProgress( newProgress );
      oldProgress = newProgress;
      if( this->GetAbortGenerateData() )
        {
        this->InvokeEvent( AbortEvent() );
        this->ResetPipeline();
        ProcessAborted e(__FILE__, __LINE__);

        e.SetDescription("Process aborted.");
        e.SetLocation(ITK_LOCATION);
        throw e;
        }
      }
    }
}

/*
 *
 */

template <class TLevelSet, class TTensorImage>
double
DtiFastMarchingCostFilter<TLevelSet, TTensorImage>
::InitializeTrialPoints(
  // const TensorImageType * tensorImage,
  const IndexType & index,
  LevelSetImageType *output )

{
  IndexType    neighIndex = index; // index of input alive point
  AxisNodeType node;

  // typedef vnl_vector_fixed<float,dimension> TVector;
  TVector distance, normal;

  distance.fill(0);
  normal.fill(0);
  OutputSpacingType spacing = this->GetOutput()->GetSpacing();
  PixelType         outputPixel,  priorOutputPixel;
  double            solution(0.0);

  double outputSpeedPixel;

  // make sure the heap is empty
  while( !m_TrialHeap.empty() )
    {
    m_TrialHeap.pop();
    }

  // Get complete neighborhood of alive point to process as trial points

  ConstNeighborhoodIteratorType::RadiusType radius;
  radius.Fill(1);
  ConstNeighborhoodIteratorType eigNeighborIt( radius, m_EigenvectorImage,
                                               m_EigenvectorImage->GetLargestPossibleRegion() );

  eigNeighborIt.SetLocation(neighIndex); // set center at alive point index

  EigenvectorImageType::IndexType           eigIndex = eigNeighborIt.GetIndex();
  ConstNeighborhoodIteratorType::OffsetType eigoffset;
  eigoffset.Fill(0);
  for( unsigned i = 0; i < eigNeighborIt.Size(); i++ )
    {
    eigIndex = eigNeighborIt.GetIndex(i); // eigIndex is new trial point

    if( !m_BufferedRegion.IsInside(eigIndex ) )   // out of image region
      {
      continue;
      }

    // bool pass=false;
    // float checkEigValue=0.0;

    // const EigenvectorPixelType eigPixel = m_EigenvectorImage->GetPixel(eigIndex);

    /*
    for(unsigned int j=0; j<dimension; j++)
    {
      checkEigValue=eigPixel[j];
      if ( ( checkEigValue > 0.0) || (checkEigValue < 0.0 ) )
        pass=true;
    }

    if  (!pass)
    {
      continue; //out of brain region
    }
    */

    if( m_LabelImage->GetPixel( eigIndex ) != AlivePoint )
      {
      // calculate distance
      eigoffset = eigNeighborIt.GetOffset(i);
      for( unsigned int k = 0; k < dimension; k++ )
        {
        distance[k] = eigoffset[k] * vnl_math_abs(spacing[k]);
        }

      // Set speed of initial trial points to F=|dot
      // product(eigenvector(TrialPoint),eigenvector(TrialPoint)|
      TVector trialEigvalue = m_EigenvectorImage->GetPixel(eigIndex).GetVnlVector();
      normal = distance;
      normal = normal.normalize();

      // outputSpeedPixel =
      // (vnl_math_abs(dot_product(trialEigvalue,trialEigvalue)) );
      outputSpeedPixel = ( vnl_math_abs( dot_product(trialEigvalue, normal) ) );

      /*Scale by FA*/
      AnisotropyImagePixelType aniso = 0.0;

      if( m_AnisotropyWeight > 0 )
        {
        aniso =  m_AnisotropyImage->GetPixel( node.GetIndex() );
        }
      outputSpeedPixel = outputSpeedPixel * ( 1 - m_AnisotropyWeight ) + aniso * m_AnisotropyWeight;

      // Normalize Speed
      PixelType normOutputSpeedPixel = outputSpeedPixel / m_NormalizationFactor;

      // Calculate vcl_cost (time): neighbor time (of selected Alive point) +
      // distance/speed (of Trial point)
      if( normOutputSpeedPixel > 0.0 )   // not zero
        {
        double neighTime = output->GetPixel( neighIndex );
        double trialTime = ( ( distance.magnitude() ) / normOutputSpeedPixel );
        solution = neighTime + trialTime;
        }
      else
        {
        solution = m_LargeValue;
        }

      solution = static_cast<PixelType>(solution);

      priorOutputPixel = output->GetPixel( eigIndex); // Previous time of trial
                                                      // point

      if( ( solution < priorOutputPixel ) & ( solution < m_LargeValue ) )
        {
        // write solution to m_OutputLevelSet
        outputPixel = solution;
        output->SetPixel( eigIndex, outputPixel );

        // write output speed of trial point to m_OutputSpeedImage
        m_OutputSpeedImage->SetPixel( eigIndex, normOutputSpeedPixel );

        // insert point into trial heap
        m_LabelImage->SetPixel( eigIndex, TrialPoint );
        node.SetValue( outputPixel );
        node.SetIndex( eigIndex );
        m_TrialHeap.push( node );
        }
      }
    }
  return solution;
}

/*
 *
 */
template <class TLevelSet, class TTensorImage>
void
DtiFastMarchingCostFilter<TLevelSet, TTensorImage>
::UpdateNeighbors(
  //  const TensorImageType * tensorImage,
  const IndexType & index,
  LevelSetImageType *output )

{
  IndexType neighIndex = index; // index of most recent alive point

  // Get complete neighborhood of alive point to process as trial points
  typedef itk::ConstNeighborhoodIterator<EigenvectorImageType> ConstNIterType;

  ConstNIterType::RadiusType radius;
  radius.Fill(1);
  ConstNIterType eigNeighborIt( radius, m_EigenvectorImage, m_EigenvectorImage->GetLargestPossibleRegion() );

  eigNeighborIt.SetLocation(neighIndex); // set center at alive point index

  EigenvectorImageType::IndexType eigIndex = eigNeighborIt.GetIndex();
  for( unsigned i = 0; i < eigNeighborIt.Size(); i++ )
    {
    eigIndex = eigNeighborIt.GetIndex(i); // eigIndex is new trial point

    if( !m_BufferedRegion.IsInside(eigIndex ) )   // out of image region
      {
      continue;
      }

    bool                 pass = false;
    float                checkEigValue = 0.0;
    EigenvectorPixelType eigPixel = m_EigenvectorImage->GetPixel(eigIndex);
    for( unsigned int j = 0; j < dimension; j++ )
      {
      checkEigValue = eigPixel[j];
      if( ( checkEigValue > 0.0 ) || ( checkEigValue < 0.0 ) )
        {
        pass = true;
        }
      }

    if( !pass )
      {
      continue;   // out of brain region
      }

    if( m_LabelImage->GetPixel( eigIndex ) != AlivePoint )
      {
      this->UpdateValue( /* tensorImage, */ eigIndex, output );
      // this->ModifiedUpdateValue( tensorImage, eigIndex, output );
      // this->SimpleModifiedUpdateValue( tensorImage, neighIndex, eigIndex,
      // output );
      }
    }
}

/*
 *
 */

template <class TLevelSet, class TTensorImage>
double
DtiFastMarchingCostFilter<TLevelSet, TTensorImage>

::UpdateValue(
  // const TensorImageType * tensorImage,
  const IndexType & index,
  LevelSetImageType *output )
{
  IndexType neighIndex = index; // index of trial point

  PixelType outputPixel;
  PixelType priorOutputPixel;

  double outputSpeedPixel(0.0);

  double neighSpeedPixel(0.0);

  AxisNodeType node;

  // typedef vnl_vector_fixed<float,dimension> TVector;
  typedef std::list<TVector> VectorListType;
  VectorListType    offsetList;
  TVector           sum; sum.fill(0);
  TVector           neighOffset;
  double            solution;
  OutputSpacingType spacing = this->GetOutput()->GetSpacing();

  // Get complete neighborhood of trial point to calculate normal

  typedef itk::ConstNeighborhoodIterator<EigenvectorImageType> ConstNIterType;
  ConstNIterType::RadiusType radius;
  radius.Fill(1);
  ConstNIterType eigNeighborIt( radius, m_EigenvectorImage, m_EigenvectorImage->GetLargestPossibleRegion() );

  eigNeighborIt.SetLocation(neighIndex);   // set center at trial point index
  ConstNIterType::OffsetType eigoffset;
  eigoffset.Fill(0);
  offsetList.clear();
  for( unsigned i = 0; i < eigNeighborIt.Size(); i++ )
    {
    EigenvectorImageType::IndexType eigIndex = eigNeighborIt.GetIndex(i);

    if( !m_BufferedRegion.IsInside( eigIndex ) )
      {
      continue;
      }

    if( m_LabelImage->GetPixel( eigIndex ) == AlivePoint )
      {
      eigoffset = eigNeighborIt.GetOffset(i);
      for( unsigned int k = 0; k < dimension; k++ )
        {
        neighOffset[k] = eigoffset[k] * vnl_math_abs(spacing[k]);
        }

      // Normal calculation based on neighborhood
      offsetList.push_back(neighOffset);
      sum += neighOffset;       // sum all offsets
      }
    }

  TVector normal = sum.normalize();

  /* Find Alive point that is closet to direction "-normal"
     by computing vector angle: vcl_cos theta= dot product (v1,normal)/|v1|  */

  TVector rOffset; rOffset.fill(0);  // alive point in direction -normal
  neighOffset.fill(0);
  neighOffset = offsetList.back();
  offsetList.pop_back();
  float vcl_cosAngle1 = ( dot_product(normal, neighOffset) ) / ( neighOffset.magnitude() );
  rOffset = neighOffset;

  while( !offsetList.empty() )
    {
    neighOffset.fill(0);
    neighOffset = offsetList.back();
    offsetList.pop_back();
    float vcl_cosAngle2 = ( dot_product(normal, neighOffset) ) / ( neighOffset.magnitude() );

    if( vcl_cosAngle2 > vcl_cosAngle1 )
      {
      vcl_cosAngle1 = vcl_cosAngle2;
      rOffset = neighOffset;
      }
    }

  // reset neighbor offset so it is alive point in direction -normal
  neighOffset = rOffset;
  // convert neighbor Offsets to neighbor Index (previously index of trial
  // point)
  for( int i = 0; i < dimension; i++ )
    {
    long int val
      = static_cast<long int>( neighOffset[i]
                               / vnl_math_abs(spacing[i]) );
    neighIndex[i] += val;
    }

  // Compute Speed: F(r)=min[F(AlivePoint), |dot product(principal
  // eigenvector(AlivePoint),normal)| ]
  TVector neighEigvalue = m_EigenvectorImage->GetPixel(neighIndex).GetVnlVector();
  double  trialSpeedPixel = ( vnl_math_abs( dot_product(neighEigvalue, normal) ) );

  /*Scale trialSpeed by FA, Fractional Anisotropy of trial point*/
  AnisotropyImagePixelType aniso = 0.0;
  if( m_AnisotropyWeight > 0 )
    {
    aniso = m_AnisotropyImage->GetPixel(index);
    }
  trialSpeedPixel = trialSpeedPixel * ( 1 - m_AnisotropyWeight ) + aniso * m_AnisotropyWeight;

  // Normalize Trial Speed
  if( trialSpeedPixel > 0.0 )
    {
    trialSpeedPixel /= m_NormalizationFactor;

    // Compute Speed: F(r)=min[F(AlivePoint), F(TrialPoint ]
    neighSpeedPixel = m_OutputSpeedImage->GetPixel( neighIndex );
    outputSpeedPixel = vnl_math_min( neighSpeedPixel, trialSpeedPixel );

    // Calculate vcl_cost (time): neighbor time (of selected Alive point) +
    // distance/speed (of Trial point)
    double neighTime = output->GetPixel( neighIndex );
    double trialTime = ( ( neighOffset.magnitude() ) / outputSpeedPixel );
    solution = neighTime + trialTime; // Total time
    }
  else
    {
    solution = m_LargeValue; // else trialSpeedPixel is zero
    }

  priorOutputPixel = output->GetPixel( index ); // Previous time of trial point

  solution = static_cast<PixelType>(solution);

  if( ( solution < priorOutputPixel ) & ( solution < m_LargeValue ) )
    {
    // write solution to m_OutputLevelSet
    outputPixel = solution;
    output->SetPixel( index, outputPixel );

    // write output speed of trial point to m_OutputSpeedImage
    m_OutputSpeedImage->SetPixel( index, static_cast<PixelType>(outputSpeedPixel) );

    // insert Trial point into trial heap
    m_LabelImage->SetPixel( index, TrialPoint );
    node.SetValue( outputPixel );
    node.SetIndex( index );
    m_TrialHeap.push( node );
    }

  return solution;
}

/*
 *
 */

template <class TLevelSet, class TTensorImage>
double
DtiFastMarchingCostFilter<TLevelSet, TTensorImage>

::ModifiedUpdateValue(
  const TensorImageType *tensorImage,
  const IndexType & index,
  LevelSetImageType *output )
{
  IndexType neighIndex = index; // index of trial point

  PixelType                       outputPixel;
  PixelType                       priorOutputPixel;
  double                          outputSpeedPixel = 0.0;
  double                          neighSpeedPixel = 0.0;
  EigenvectorImageType::IndexType eigIndex;
  double                          trialSpeedPixel = -1.0; //
                                                          // trialSpeedPixel>=0.0;

  AxisNodeType node;

  // typedef vnl_vector_fixed<float,dimension> TVector;
  TVector           aliveOffset;
  double            solution;
  OutputSpacingType spacing = this->GetOutput()->GetSpacing();

  // Get complete neighborhood of trial point to calculate normal

  typedef itk::ConstNeighborhoodIterator<EigenvectorImageType> ConstNIterType;
  ConstNIterType::RadiusType radius;
  radius.Fill(1);
  ConstNIterType eigNeighborIt( radius, m_EigenvectorImage, m_EigenvectorImage->GetLargestPossibleRegion() );

  eigNeighborIt.SetLocation(neighIndex);   // set center at trial point index
  ConstNIterType::OffsetType eigoffset;
  eigoffset.Fill(0);
  eigIndex = eigNeighborIt.GetIndex();
  for( unsigned i = 0; i < eigNeighborIt.Size(); i++ )
    {
    eigIndex = eigNeighborIt.GetIndex(i);

    if( !m_BufferedRegion.IsInside(eigIndex ) )
      {
      continue;
      }

    if( m_LabelImage->GetPixel( eigIndex ) == AlivePoint )
      {
      eigoffset = eigNeighborIt.GetOffset(i);
      TVector neighOffset, normal;
      neighOffset.fill(0); normal.fill(0);
      // New Normal calculation n=r-r', where r' is AlivePoint in neighborhood
      // of r, TrialPoint
      for( unsigned int k = 0; k < dimension; k++ )  // convert offset to vector
        {
        neighOffset[k] = eigoffset[k] * vnl_math_abs(spacing[k]);
        }
      normal = neighOffset;
      normal = normal.normalize();

      // Compute Speed: F(r)= |dot product(principal
      // eigenvector(AlivePoint),normal)|
      TVector aliveEigvalue = m_EigenvectorImage->GetPixel(eigIndex).GetVnlVector();
      double  maxSpeedPixel = vnl_math_abs( dot_product(aliveEigvalue, normal) );

      // Get max Speedvalue
      if( maxSpeedPixel > trialSpeedPixel )
        {
        trialSpeedPixel = maxSpeedPixel;
        aliveOffset = neighOffset;               // offset of AlivePoint
                                                 // selected
        }
      }   // end of "if" Alive Point
    } // end neighborhood "for" loop
  // if (trialSpeedPixel == 0.0)
  // {trialSpeedPixel /= m_NormalizationFactor;}
  /*Convert alive point Offsets to neighbor Index
  (before was Trial Point index, now it is Alive Point index)*/
  for( int i = 0; i < dimension; i++ )
    {
    long int val
      = static_cast<long int>( aliveOffset[i] / vnl_math_abs(spacing[i]) );
    neighIndex[i] += val;
    }

  /*Scale trialSpeed by FA, Fractional Anisotropy of trial point*/
  TensorImagePixelType tensorPixel =  tensorImage->GetPixel(index);
  float                ai = tensorPixel.GetFractionalAnisotropy();
  trialSpeedPixel *= ai;

  // Normalize Trial Speed
  if( trialSpeedPixel > 0.0 )
    {
    trialSpeedPixel /= m_NormalizationFactor;

    // Compute Speed: F(r)=min[F(AlivePoint), F(TrialPoint ]
    neighSpeedPixel = m_OutputSpeedImage->GetPixel( neighIndex );
    outputSpeedPixel = vnl_math_min( neighSpeedPixel, trialSpeedPixel );

    // Calculate vcl_cost (time): neighbor time (of selected Alive point) +
    // distance/speed (of Trial point)
    double neighTime = output->GetPixel( neighIndex );
    double trialTime = ( ( aliveOffset.magnitude() ) / outputSpeedPixel );
    solution = neighTime + trialTime; // Total time
    }
  else
    {
    solution = m_LargeValue;
    }

  priorOutputPixel = output->GetPixel( index ); // Previous time of trial point
  solution = static_cast<PixelType>(solution);

  if( ( solution < priorOutputPixel ) & ( solution < m_LargeValue ) )
    {
    // write solution to m_OutputLevelSet
    outputPixel = solution;
    output->SetPixel( index, outputPixel );

    // write output speed of trial point to m_OutputSpeedImage
    m_OutputSpeedImage->SetPixel( index, static_cast<PixelType>(outputSpeedPixel) );

    // insert trial point into trial heap
    m_LabelImage->SetPixel( index, TrialPoint );
    node.SetValue( outputPixel );
    node.SetIndex( index );
    m_TrialHeap.push( node );
    }

  return solution;
}

/*
 *
 */
} // namespace itk

#endif
