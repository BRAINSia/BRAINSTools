/*=========================================================================

 Program:   BRAINS (Brain Research: Analysis of Images, Networks, and Systems)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2011/07/21 08:40:00 $
 Version:   $Revision: 1.9 $

   Copyright (c) University of Iowa Department of Radiology. All rights reserved.
   See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
   for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "CompareSurfacesCLP.h"

#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"

#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkCellData.h"

#include <vector>
#include "vtkMath.h"
#include <vtkSmartPointer.h>

int main( int argc, char *argv[] )
{
  PARSE_ARGS;

  if( (inputSurfaceFile == "") || (refSurfaceFile == "") )
    {
    std::cerr << "No input surfaces specified" << std::endl;
    return EXIT_FAILURE;
    }
  std::cout << "Compare " << inputSurfaceFile << std::endl;
  std::cout << " with " << refSurfaceFile << std::endl;
  if( vertexLocation )
    {
    std::cout << "\tVertex locations" << std::endl;
    std::cout << "\twith tolerance: " << tolerance << std::endl;
    }
  if( scalarArray )
    {
    std::cout << "\tScalar values" << std::endl;
    std::cout << "\twith tolerance: " << tolerance << std::endl;
    }
  std::cout << "------------------------------------------------------" << std::endl;

  // read the input surface
  vtkSmartPointer<vtkPolyDataReader> inputReader =
    vtkSmartPointer<vtkPolyDataReader>::New();

  inputReader->SetFileName( inputSurfaceFile.c_str() );
  inputReader->Update();

  vtkPolyData *inputSurface = inputReader->GetOutput();

  // read the reference surface
  vtkSmartPointer<vtkPolyDataReader> refReader =
    vtkSmartPointer<vtkPolyDataReader>::New();

  refReader->SetFileName( refSurfaceFile.c_str() );
  refReader->Update();

  vtkPolyData *refSurface = refReader->GetOutput();

  unsigned int numberOfPoints = inputSurface->GetNumberOfPoints();
  if( numberOfPoints != refSurface->GetNumberOfPoints() )
    {
    std::cout << "The number of points are different between two surfaces." << std::endl;
    std::cout << "Test failed." << std::endl;
    return EXIT_FAILURE;
    }

  // if compare the vertex location
  double pInput[3];
  double pRef[3];
  if( vertexLocation )
    {
    for( unsigned int i = 0; i < numberOfPoints; i++ )
      {
      inputSurface->GetPoint(i, pInput);
      refSurface->GetPoint(i, pRef);

      // calculate the distance between corresponding points
      double disP = vtkMath::Distance2BetweenPoints(pInput, pRef);

      if( disP > tolerance )
        {
        std::cout << "Distance at point " << i << " is bigger than tolerance." << std::endl;
        std::cout << "Test failed." << std::endl;
        return EXIT_FAILURE;
        }
      }
    }

  if( scalarArray )
    {
    // compare number of arrays
    int numArrays = inputSurface->GetPointData()->GetNumberOfArrays();
    if( refSurface->GetPointData()->GetNumberOfArrays() != numArrays )
      {
      std::cout << "The number of arrays are different between two surfaces." << std::endl;
      std::cout << "Test failed." << std::endl;
      return EXIT_FAILURE;
      }
    for( int array_i = 0; array_i < numArrays; array_i++ )
      {
      // compare array name
      const char *arrayName_input = inputSurface->GetPointData()->GetArrayName(array_i);
      const char *arrayName_ref = refSurface->GetPointData()->GetArrayName(array_i);
      if( strcmp(arrayName_input, arrayName_ref) != 0 )
        {
        std::cout << "The " << array_i << "th array name is different between two surfaces." << std::endl;
        std::cout << "Test failed." << std::endl;
        return EXIT_FAILURE;
        }

      // compare scalar values of each array
      vtkDataArray *inputArray = inputSurface->GetPointData()->GetArray(array_i);
      vtkDataArray *refArray = refSurface->GetPointData()->GetArray(array_i);
      // compare number of tuples
      if( inputArray->GetNumberOfTuples() != refArray->GetNumberOfTuples() )
        {
        std::cout << "The number of tuples at " << array_i << "th array are different between two surfaces."
                  << std::endl;
        std::cout << "Test failed." << std::endl;
        return EXIT_FAILURE;
        }
      // compare number of components
      if( inputArray->GetNumberOfComponents() != refArray->GetNumberOfComponents() )
        {
        std::cout << "The number of components at " << array_i << "th array are different between two surfaces."
                  << std::endl;
        std::cout << "Test failed." << std::endl;
        return EXIT_FAILURE;
        }
      // compare each scalar value in that array
      for( int i = 0; i < inputArray->GetNumberOfTuples(); i++ )
        {
        for( int c_i = 0; c_i < inputArray->GetNumberOfComponents(); c_i++ )
          {
          double scalarDist = fabs( inputArray->GetComponent(i, c_i) - refArray->GetComponent(i, c_i) );
          if( scalarDist > tolerance )
            {
            std::cout << "The scalar values of " << c_i << "th component at " << i << "th tuple in " << array_i
                      << "th array are different between two surfaces." << std::endl;
            std::cout << "Test failed." << std::endl;
            return EXIT_FAILURE;
            }
          }
        }
      }
    }

  if( cellData )
    {
    // compare the number of cellData arrays
    int numCellDatas = inputSurface->GetCellData()->GetNumberOfArrays();
    if( refSurface->GetCellData()->GetNumberOfArrays() != numCellDatas )
      {
      std::cout << "The number of cellDatas are different between two surfaces." << std::endl;
      std::cout << "Test failed." << std::endl;
      return EXIT_FAILURE;
      }
    for( int array_i = 0; array_i < numCellDatas; array_i++ )
      {
      // compare cellData name
      const char *arrayName_input = inputSurface->GetCellData()->GetArrayName(array_i);
      const char *arrayName_ref = refSurface->GetCellData()->GetArrayName(array_i);
      if( strcmp(arrayName_input, arrayName_ref) != 0 )
        {
        std::cout << "The " << array_i << "th array name is different between two surfaces." << std::endl;
        std::cout << "Test failed." << std::endl;
        return EXIT_FAILURE;
        }

      // compare scalar values of each array
      vtkDataArray *inputArray = inputSurface->GetCellData()->GetArray(array_i);
      vtkDataArray *refArray = refSurface->GetCellData()->GetArray(array_i);
      // compare number of tuples
      if( inputArray->GetNumberOfTuples() != refArray->GetNumberOfTuples() )
        {
        std::cout << "The number of tuples at " << array_i << "th array are different between two surfaces."
                  << std::endl;
        std::cout << "Test failed." << std::endl;
        return EXIT_FAILURE;
        }
      // compare number of components
      if( inputArray->GetNumberOfComponents() != refArray->GetNumberOfComponents() )
        {
        std::cout << "The number of components at " << array_i << "th array are different between two surfaces."
                  << std::endl;
        std::cout << "Test failed." << std::endl;
        return EXIT_FAILURE;
        }
      // compare each scalar value in that array
      for( int i = 0; i < inputArray->GetNumberOfTuples(); i++ )
        {
        for( int c_i = 0; c_i < inputArray->GetNumberOfComponents(); c_i++ )
          {
          double scalarDist = fabs( inputArray->GetComponent(i, c_i) - refArray->GetComponent(i, c_i) );
          if( scalarDist > tolerance )
            {
            std::cout << "The scalar values of " << c_i << "th component at " << i << "th tuple in " << array_i
                      << "th array are different between two surfaces." << std::endl;
            std::cout << "Test failed." << std::endl;
            return EXIT_FAILURE;
            }
          }
        }
      }
    }

  std::cout << "Test succeed!" << std::endl;
  return EXIT_SUCCESS;
}
