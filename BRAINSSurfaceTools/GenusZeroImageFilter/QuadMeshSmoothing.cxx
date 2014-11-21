//remove the shortest edge at each iteration 
//until the given criterion is fulfilled.
//NumberOfIters: 5

#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshParamMatrixCoefficients.h"

typedef double Coord;
const unsigned int Dimension = 3;

// Declaration of the type of Mesh
typedef itk::QuadEdgeMesh< Coord, Dimension > MeshType;

#if ITK_VERSION_MAJOR < 4
#include "itkQuadEdgeMeshSmoothing.h"
typedef itk::QuadEdgeMeshSmoothing< MeshType, MeshType > SmoothingType;
#else
#include "itkSmoothingQuadEdgeMeshFilter.h"
typedef itk::SmoothingQuadEdgeMeshFilter< MeshType, MeshType > SmoothingType;
#endif

#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "QuadMeshSmoothingCLP.h"

int main( int argc, char * argv [] )
{

  PARSE_ARGS;

  // Here read a mesh from a file
  typedef itk::QuadEdgeMeshVTKPolyDataReader< MeshType >  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputSurface );

  // Note that any other coefficients could have been used
  itk::OnesMatrixCoefficients< MeshType > coeff0;

  SmoothingType::Pointer filter = SmoothingType::New( );
  filter->SetInput( reader->GetOutput() );
  filter->SetNumberOfIterations( numberOfIterations );
  filter->SetRelaxationFactor( relaxationFactor );
  if ( delaunayConforming )
  {
    filter->SetDelaunayConforming( true );
  }
  filter->SetCoefficientsMethod( &coeff0 );
  filter->Update();

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter< MeshType >  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( outputSurface );
  writer->Update(); 

  return EXIT_SUCCESS;
}
