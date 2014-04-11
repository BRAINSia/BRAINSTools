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

  Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   vtkITK
  Module:    $HeadURL$
  Date:      $Date$
  Version:   $Revision$

==========================================================================*/

#ifndef __vtkITKArchetypeImageSeriesScalarReader_h
#define __vtkITKArchetypeImageSeriesScalarReader_h

#include "vtkITKArchetypeImageSeriesReader.h"

#include "itkImageFileReader.h"

class VTK_ITK_EXPORT vtkITKArchetypeImageSeriesScalarReader : public vtkITKArchetypeImageSeriesReader
{
public:
  static vtkITKArchetypeImageSeriesScalarReader * New();

  vtkTypeMacro(vtkITKArchetypeImageSeriesScalarReader, vtkITKArchetypeImageSeriesReader);
  void PrintSelf(ostream& os, vtkIndent indent);

protected:
  vtkITKArchetypeImageSeriesScalarReader();
  ~vtkITKArchetypeImageSeriesScalarReader();

  void ExecuteData(vtkDataObject *data);

  static void ReadProgressCallback(itk::ProcessObject* obj, const itk::ProgressEvent &, void* data);

  /// private:
};

#endif
