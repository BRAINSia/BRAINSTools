/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkScalarQuadEdgeMeshToListAdaptorTest1.cxx,v $
  Language:  C++
  Date:      $Date: 2008-03-10 19:46:31 $
  Version:   $Revision: 1.37 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkVector.h"
#include "itkListSample.h"
#include "itkHistogram.h"
#include "itkQuadEdgeMesh.h"

#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkScalarQuadEdgeMeshToListAdaptor.h"
#include "itkSampleToHistogramFilter.h"

int main( int argc, char * argv [] )
{
  if( argc < 5 )
    {
    std::cerr << "Missing command line arguments" << std::endl;
    std::cerr << "Usage :  QuadEdgeMeshHistogram  inputMeshFileName " << std::endl;
    std::cerr << "NumberOfBins Min Max" << std::endl;
    return -1;
    }

  typedef float MeshPixelType;
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> MeshType;

  typedef itk::QuadEdgeMeshVTKPolyDataReader<MeshType> InputReaderType;

  InputReaderType::Pointer inputMeshReader = InputReaderType::New();
  inputMeshReader->SetFileName( argv[1] );
  inputMeshReader->Update();

  typedef itk::ScalarQuadEdgeMeshToListAdaptor<MeshType> AdaptorType;

  AdaptorType::Pointer adaptor = AdaptorType::New();

  adaptor->SetMesh(  inputMeshReader->GetOutput() );

  try
    {
    adaptor->Compute();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  typedef AdaptorType::MeasurementType                       MeasurementType;
  typedef AdaptorType::MeasurementVectorType                 MeasurementVectorType;
  typedef itk::Statistics::ListSample<MeasurementVectorType> ListSampleType;

  typedef itk::Statistics::Histogram<MeasurementType,
                                     itk::Statistics::DenseFrequencyContainer2> HistogramType;

  typedef itk::Statistics::SampleToHistogramFilter<
      ListSampleType, HistogramType> HistogramFilterType;

  HistogramFilterType::Pointer histogramFilter = HistogramFilterType::New();

  typedef HistogramFilterType::HistogramType HistogramType;

  histogramFilter->SetInput( adaptor->GetSample() );

  histogramFilter->SetMarginalScale( 10.0 );

  typedef HistogramFilterType::HistogramSizeType HistogramSizeType;

  MeasurementVectorType value;

  unsigned int numberOfScalarComponents =
    itk::Statistics::MeasurementVectorTraits::GetLength( value );

  HistogramSizeType histogramSize( numberOfScalarComponents );
  histogramSize[0] = atoi( argv[2] );

  histogramFilter->SetHistogramSize( histogramSize );

  typedef HistogramFilterType::HistogramMeasurementVectorType HistogramMeasurementVectorType;

  HistogramMeasurementVectorType histogramBinMinimum( numberOfScalarComponents );
  histogramBinMinimum.Fill(  atof(argv[3]) - 0.5 );

  HistogramMeasurementVectorType histogramBinMaximum( numberOfScalarComponents );
  histogramBinMaximum.Fill(  atof(argv[4]) - 0.5 );

  histogramFilter->SetHistogramBinMinimum( histogramBinMinimum );
  histogramFilter->SetHistogramBinMaximum( histogramBinMaximum );

  histogramFilter->Update();

  HistogramType::ConstPointer histogram = histogramFilter->GetOutput();

  const unsigned int histogramTotalNumberOfBins = histogram->Size();

  std::cout << "Total Number of Bins " << histogramTotalNumberOfBins << std::endl;
  for( unsigned int bin = 0; bin < histogramTotalNumberOfBins; bin++ )
    {
    std::cout << "bin = " << bin << " frequency = ";
    std::cout << histogram->GetFrequency( bin, 0 ) << std::endl;
    }

  return 0;
}
