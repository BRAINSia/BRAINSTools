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
/*
 * Author: Hans J. Johnson, Wei Lu
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2010
 */

#ifndef _landmarksDataSet_h
#define _landmarksDataSet_h
#include "itkMacro.h"
#include "itkPoint.h"
#include <map>
#include <string>
#include <iostream>
#include <cstdlib>

/*
 * This class contains a single set of landmarks and it's associated intensity image file.
 * This is usually used to hold information like the following:
 @TEST_DATA_DIR@/A.nii.gz
 RP 128.25 -105.341 115.861
 AC 128.25 -127.356 127.553
 PC 128.25 -97.3915 127.251
 END
 * that is read in from a different class.
 */

class landmarksDataSet :  // public PointMapType
  public std::map<std::string, itk::Point<double, 3> >
{
private:
  typedef enum { eof, bad, shortLine } err_flags;
public:
  typedef itk::Point<double, 3>            PointType;
  typedef std::map<std::string, PointType> Superclass;

  landmarksDataSet()
  {
  }

  // Constructor for use by
  // landmarksConstellationTrainingDefinitionIO::ReadFile(const std::string)
  //
  landmarksDataSet(const std::string & TrainingImageFilename, const std::string & LandmarksFilename)
  {
    this->SetImageFilename(TrainingImageFilename);
    this->SetLandmarkFilename(LandmarksFilename);

    std::ifstream singleLandmarkFile( LandmarksFilename.c_str() );
    if( !singleLandmarkFile.is_open() )
      {
      itkGenericExceptionMacro(<< "Input model file cannot open! :" << LandmarksFilename);
      }
    landmarksDataSet::PointType tempPoint;
    char                        linebuf[300];

    // Read in RP, AC, PC, VN4, and LE, RE, and other "new" landmarks
    while( singleLandmarkFile.getline(linebuf, 300) )
      {
      // skip comment lines and newlines
      while( linebuf[0] == '#' || linebuf[0] == '\0' )
        {
        singleLandmarkFile.getline(linebuf, 300);
        }

      std::string            CSVs(linebuf);
      std::string::size_type posA = CSVs.find(',');
      if( posA == std::string::npos )
        {
        throw landmarksDataSet::shortLine;
        }
      std::string            PointName( CSVs.substr(0, posA - 0) );
      std::string::size_type posB = 0;
      for( unsigned int j = 0; j < 3; ++j )  // Image Dimension = 3
        {
        posB = CSVs.find(',', posA + 1);
        if( posB == std::string::npos )
          {
          throw landmarksDataSet::shortLine;
          }
        std::string datum = CSVs.substr( posA + 1, posB - ( posA + 1 ) );
        tempPoint[j] = atof( datum.c_str() );
        posA = posB;
        }
      // NOTE:  RAS is was slicer requires, but ITK is internall LPS, so we need
      // to negate the first two landmark points
      tempPoint[0] *= -1;
      tempPoint[1] *= -1;
      this->SetNamedPoint(PointName, tempPoint);
      }

    // RP, AC, PC, VN4, LE and RE are six "must-have" landmarks
    if( ( this->find("RP") == this->end() )
        || ( this->find("AC") == this->end() )
        || ( this->find("PC") == this->end() )
        || ( this->find("VN4") == this->end() )
        || ( this->find("LE") == this->end() )
        || ( this->find("RE") == this->end() ) )
      {
      itkGenericExceptionMacro(<< "Essential landmarks are missing!"
                               << "Make sure the RP (MPJ), AC, PC, 4VN (4th ventricle notch), "
                               << "LE (left eye centers), and RE (right eye centers) are all defined.");
      }

    singleLandmarkFile.close();
  }

  PointType GetNamedPoint(const std::string & NamedPoint) const
  {
    landmarksDataSet::const_iterator it = this->find(NamedPoint);

    //
    // originally, no check at all for missing points.  Next best thing,
    // throw an exception
    if( it == this->end() )
      {
      throw;
      }
    return it->second;
  }

  void SetNamedPoint(const std::string & NamedPoint, const PointType & NewPoint)
  {
    ( *this )[NamedPoint] = NewPoint;
  }

  void SetImageFilename(const std::string & NewImageFilename)
  {
    m_ImageFilename = NewImageFilename;
  }

  std::string GetImageFilename(void) const
  {
    return m_ImageFilename;
  }

  void SetLandmarkFilename(const std::string & NewImageFilename)
  {
    m_LandmarkFilename = NewImageFilename;
  }

  std::string GetLandmarkFilename(void) const
  {
    return m_LandmarkFilename;
  }

private:
  std::string m_ImageFilename;
  std::string m_LandmarkFilename;
};

inline
std::ostream & operator<<(std::ostream & os, const landmarksDataSet & ds)
{
  os << ds.GetImageFilename() << " "
     << ds.GetNamedPoint("RP") << std::endl
     << ds.GetNamedPoint("AC") << std::endl
     << ds.GetNamedPoint("PC") << std::endl
     << ds.GetNamedPoint("VN4") << std::endl
     << ds.GetNamedPoint("LE") << std::endl
     << ds.GetNamedPoint("RE") << std::endl
     << std::endl;
  return os;
}

#endif // #ifndef _landmarksDataSet_h
