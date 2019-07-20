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
#ifndef __EMSParameters_h
#define __EMSParameters_h

#include "itkObject.h"
#include "itkObjectFactory.h"

#include <iostream>
#include <string>

#include <vector>

/**
 * \class EMSParameters
 */
class EMSParameters : public itk::Object
{
public:
  using Self = EMSParameters;
  using Pointer = itk::SmartPointer< Self >;
  using ConstPointer = itk::SmartPointer< const Self >;

  itkNewMacro( Self );

  // Make sure all values are OK
  virtual bool
  CheckValues() const;

  virtual void
  PrintSelf( std::ostream & os, itk::Indent ) const override;

  itkGetConstMacro( Suffix, std::string );
  itkSetMacro( Suffix, std::string );

  itkGetConstMacro( AtlasDirectory, std::string );
  itkSetMacro( AtlasDirectory, std::string );

  itkGetConstMacro( AtlasOrientation, std::string );
  itkSetMacro( AtlasOrientation, std::string );

  itkGetConstMacro( DoAtlasWarp, bool );
  itkSetMacro( DoAtlasWarp, bool );

  itkGetConstMacro( OutputDirectory, std::string );
  itkSetMacro( OutputDirectory, std::string );

  itkGetConstMacro( OutputFormat, std::string );
  itkSetMacro( OutputFormat, std::string );

  void
  AddImage( std::string s, std::string orientation );

  void
  ClearImages();

  std::vector< std::string >
  GetImages()
  {
    return m_Images;
  }

  std::vector< std::string >
  GetImageOrientations()
  {
    return m_ImageOrientations;
  }

  itkGetConstMacro( FilterMethod, std::string );
  itkSetMacro( FilterMethod, std::string );

  itkGetConstMacro( FilterIterations, unsigned int );
  itkSetMacro( FilterIterations, unsigned int );

  itkGetConstMacro( FilterTimeStep, float );
  itkSetMacro( FilterTimeStep, float );

  itkGetConstMacro( MaxBiasDegree, unsigned int );
  itkSetMacro( MaxBiasDegree, unsigned int );

  itkGetConstMacro( AtlasWarpGridX, unsigned int );
  itkSetMacro( AtlasWarpGridX, unsigned int );
  itkGetConstMacro( AtlasWarpGridY, unsigned int );
  itkSetMacro( AtlasWarpGridY, unsigned int );
  itkGetConstMacro( AtlasWarpGridZ, unsigned int );
  itkSetMacro( AtlasWarpGridZ, unsigned int );

  itkGetConstMacro( Prior1, float );
  itkSetMacro( Prior1, float );
  itkGetConstMacro( Prior2, float );
  itkSetMacro( Prior2, float );
  itkGetConstMacro( Prior3, float );
  itkSetMacro( Prior3, float );
  itkGetConstMacro( Prior4, float );
  itkSetMacro( Prior4, float );

  itkGetConstMacro( AtlasLinearMapType, std::string );
  itkSetMacro( AtlasLinearMapType, std::string );

  itkGetConstMacro( ImageLinearMapType, std::string );
  itkSetMacro( ImageLinearMapType, std::string );

  virtual ~EMSParameters() = default;

protected:
  EMSParameters();

  std::string m_Suffix;

  std::string m_AtlasDirectory;
  std::string m_AtlasOrientation;

  bool m_DoAtlasWarp;

  unsigned int m_AtlasWarpGridX;
  unsigned int m_AtlasWarpGridY;
  unsigned int m_AtlasWarpGridZ;

  std::string m_OutputDirectory;
  std::string m_OutputFormat;

  std::vector< std::string > m_Images;
  std::vector< std::string > m_ImageOrientations;

  std::string  m_FilterMethod;
  unsigned int m_FilterIterations;
  float        m_FilterTimeStep;

  unsigned int m_MaxBiasDegree;

  float m_Prior1;
  float m_Prior2;
  float m_Prior3;
  float m_Prior4;

  std::string m_AtlasLinearMapType;
  std::string m_ImageLinearMapType;
};

#endif
