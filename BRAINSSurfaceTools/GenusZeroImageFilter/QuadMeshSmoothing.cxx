// remove the shortest edge at each iteration
// until the given criterion is fulfilled.
// NumberOfIters: 5

#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshParamMatrixCoefficients.h"

using Coord = double;
constexpr unsigned int Dimension = 3;

// Declaration of the type of Mesh
using MeshType = itk::QuadEdgeMesh<Coord, Dimension>;

#include "itkSmoothingQuadEdgeMeshFilter.h"
using SmoothingType = itk::SmoothingQuadEdgeMeshFilter<MeshType, MeshType>;

#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "QuadMeshSmoothingCLP.h"

int main( int argc, char * argv [] )
{

  PARSE_ARGS;

  // Here read a mesh from a file
  using ReaderType = itk::QuadEdgeMeshVTKPolyDataReader<MeshType>;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputSurface );

  // Note that any other coefficients could have been used
  itk::OnesMatrixCoefficients<MeshType> coeff0;

  SmoothingType::Pointer filter = SmoothingType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetNumberOfIterations( numberOfIterations );
  filter->SetRelaxationFactor( relaxationFactor );
  if( delaunayConforming )
    {
    filter->SetDelaunayConforming( true );
    }
  filter->SetCoefficientsMethod( &coeff0 );
  filter->Update();

  using WriterType = itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( outputSurface );
  writer->Update();

  return EXIT_SUCCESS;
}
