// remove the shortest edge at each iteration
// until the given criterion is fulfilled.

#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshDecimationCriteria.h"

typedef double Coord;
const unsigned int Dimension = 3;

// Declaration of the type of Mesh
typedef itk::QuadEdgeMesh<Coord, Dimension> MeshType;
// Declaration of the stopping criterion
// By default the cost function is to be minimized
typedef itk::NumberOfFacesCriterion<MeshType> CriterionType;

#if ITK_VERSION_MAJOR < 4
#include "itkQuadEdgeMeshSquaredEdgeLengthDecimation.h"
typedef itk::QuadEdgeMeshSquaredEdgeLengthDecimation
  <MeshType,
   MeshType,
   CriterionType> DecimationType;
#else
#include "itkSquaredEdgeLengthDecimationQuadEdgeMeshFilter.h"
typedef itk::SquaredEdgeLengthDecimationQuadEdgeMeshFilter
  <MeshType,
   MeshType,
   CriterionType> DecimationType;
#endif

#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

#include "QuadMeshDecimationCLP.h"

int main( int argc, char * argv [] )
{

  PARSE_ARGS;

  // Here read a mesh from a file
  typedef itk::QuadEdgeMeshVTKPolyDataReader<MeshType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputSurface );

  // Here we assume that the user wants a mesh with N faces in output.
  // N is given by user
  CriterionType::Pointer criterion = CriterionType::New();
  if( topologyChange )
    {
    criterion->SetTopologicalChange( true );
    }
  else
    {
    criterion->SetTopologicalChange( false );
    }
  criterion->SetNumberOfElements( numberOfElements );

  DecimationType::Pointer decimate = DecimationType::New();
  decimate->SetInput( reader->GetOutput() );
  decimate->SetCriterion( criterion );

  // Here we allow the relocation procedure
  decimate->Update();

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( decimate->GetOutput() );
  writer->SetFileName(outputSurface);
  writer->Update();

  return EXIT_SUCCESS;
}
