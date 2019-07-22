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
#include "WriteIdentityTransformCLP.h"
#include "itkTransformFileWriter.h"
#include "itkAffineTransform.h"

int
main( int argc, char ** argv )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  using TransformationType = itk::AffineTransform< double, 3 >;

  TransformationType::Pointer transform = TransformationType::New();

  transform->SetIdentity();

  using WriterType = itk::TransformFileWriter;
  WriterType::Pointer writer = WriterType::New();

  writer->SetFileName( "./identityTransformation.mat" );
  writer->SetInput( transform );
        writer->SetUseCompression( true );
  writer->Update();
}
