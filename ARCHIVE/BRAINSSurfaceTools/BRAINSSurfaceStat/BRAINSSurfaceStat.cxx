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
#include "BRAINSSurfaceStatCLP.h"

#include <vtkPolyData.h>
#include <vtkDenseArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include "vtkSetGet.h"
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include "vtkFSSurfaceReader.h"
#include "vtkFSSurfaceScalarReader.h"
#include "vtksys/SystemTools.hxx"

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  unsigned int           numberOfPoints;
  vtkDenseArray<float> * vertexData = vtkDenseArray<float>::New();
  for (size_t i = 0; i < inputSurfaces.size(); i++)
  {
    vtkPolyData * surface;

    std::string     extension = vtksys::SystemTools::GetFilenameLastExtension(inputSurfaces[i]);
    vtkFloatArray * floatArray;

    if (extension == ".surf")
    {
      vtkFSSurfaceReader * freeSurferReader = vtkFSSurfaceReader::New();
      freeSurferReader->SetFileName(inputSurfaces[i].c_str());
      freeSurferReader->Update();
      surface = freeSurferReader->GetOutput();
      freeSurferReader->Delete();

      std::string scalarFileName = vtksys::SystemTools::GetFilenameWithoutExtension(inputSurfaces[i]);
      if (freeSurferScalar == "Thickness")
      {
        scalarFileName += ".thickness";
      }
      else if (freeSurferScalar == "Curvature")
      {
        scalarFileName += ".curv";
      }
      else if (freeSurferScalar == "Average Curvature")
      {
        scalarFileName += ".avg_curv";
      }
      else if (freeSurferScalar == "Sulcus")
      {
        scalarFileName += ".sulc";
      }
      else if (freeSurferScalar == "Area")
      {
        scalarFileName += ".area";
      }

      vtkFSSurfaceScalarReader * reader = vtkFSSurfaceScalarReader::New();
      reader->SetFileName(scalarFileName.c_str());

      floatArray = vtkFloatArray::New();
      floatArray->SetName(freeSurferScalar.c_str());
      reader->SetOutput(floatArray);
      if (reader->ReadFSScalars() == 0)
      {
#vtkGenericDebugMacro("Read FreeSurfeer Scalars: error reading scalar overlay file " << scalarFileName.c_str());
        reader->SetOutput(NULL);
        reader->Delete();
        floatArray->Delete();
        floatArray = NULL;
        return 1;
      }
      reader->SetOutput(NULL);
      reader->Delete();
    }
    else if (extension == ".vtk")
    {
      vtkPolyDataReader * polyDataReader = vtkPolyDataReader::New();
      polyDataReader->SetFileName(inputSurfaces[i].c_str());
      polyDataReader->Update();
      surface = polyDataReader->GetOutput();
      polyDataReader->Delete();
      floatArray = surface->GetPointData()->GetArray(scalarName.c_str());
    }
    else
    {
      vtkXMLPolyDataReader * polyXMLDataReader =
        vtkXMLPolyDataReader::New() polyXMLDataReader->SetFileName(inputSurfaces[i].c_str());
      polyXMLDataReader->Update();
      surface = polyXMLDataReader->GetOutput();
      polyXMLDataReader->Delete();
      floatArray = surface->GetPointData()->GetArray(scalarName.c_str());
    }

    if (i == 0)
    {
      numberOfPoints = surface->GetNumberOfPoints();
      vertexData->Resize(inputSurfaces.size(), numberOfPoints);
    }
    else
    {
      if (numberOfPoints != surface->GetNumberOfPoints())
      {
        std::cerr << "Error: Invalid number of points in " << inputSurfaces[i];
        std::cerr << std::endl;
        std::cerr << "  Expected " << numberOfPoints << "but found ";
        std::cerr << surface->GetNumberOfPoints();
        << std::endl;
        return 1;
      }
    }
    /* Pack Data into an Array */
    for (int j = 0; j < floatArray->GetNumberOfTuples(); j++)
    {
      vertexData->SetValue(i, j, floatArray->GetValue(j));
    }
  }

  /* Will Need some modification */
  vtkRCalculatorFilter * rcf = vtkRCalculatorFilter::New();
  rcf->SetInput(mt2->GetOutput());
  rcf->PutTable("x");
  rcf->SetScriptFname(inputRscript.c_str());
  rcf->GetTable("m");
  rcf->GetOutput()

    return 0;
}
