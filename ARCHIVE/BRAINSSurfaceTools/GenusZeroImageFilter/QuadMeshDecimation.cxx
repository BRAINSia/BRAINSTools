// remove the shortest edge at each iteration

// until the given criterion is fulfilled.

#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshDecimationCriteria.h"

using Coord = double;
constexpr unsigned int Dimension = 3;

// Declaration of the type of Mesh
using MeshType = itk::QuadEdgeMesh<Coord, Dimension>;
// Declaration of the stopping criterion
// By default the cost function is to be minimized
using CriterionType = itk::NumberOfFacesCriterion<MeshType>;

#include "itkSquaredEdgeLengthDecimationQuadEdgeMeshFilter.h"
typedef itk::SquaredEdgeLengthDecimationQuadEdgeMeshFilter<MeshType, MeshType, CriterionType> DecimationType;

#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "QuadMeshDecimationCLP.h"

int
main(int argc, char * argv[])
{

  PARSE_ARGS;

  // Here read a mesh from a file
  using ReaderType = itk::QuadEdgeMeshVTKPolyDataReader<MeshType>;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputSurface);

  // Here we assume that the user wants a mesh with N faces in output.
  // N is given by user
  CriterionType::Pointer criterion = CriterionType::New();
  if (topologyChange)
  {
    criterion->SetTopologicalChange(true);
  }
  else
  {
    criterion->SetTopologicalChange(false);
  }
  criterion->SetNumberOfElements(numberOfElements);

  DecimationType::Pointer decimate = DecimationType::New();

  decimate->SetInput(reader->GetOutput());

  decimate->SetCriterion(criterion);

  // Here we allow the relocation procedure

  decimate->Update();

  using WriterType = itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(decimate->GetOutput());
  writer->SetFileName(outputSurface);
  writer->Update();

  return EXIT_SUCCESS;
}
