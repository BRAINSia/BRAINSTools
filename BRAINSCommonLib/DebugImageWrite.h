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
#ifndef __DebugImageWrite_h
#define __DebugImageWrite_h
#include "BRAINSCommonLib.h"
#if defined( BRAINS_DEBUG_IMAGE_WRITE )
#  include "itkIO.h"
#  include "itksys/SystemTools.hxx"

#  define DEFINE_DEBUG_IMAGE_COUNTER                                                                                   \
    namespace DebugImageWrite                                                                                          \
    {                                                                                                                  \
    int fileSequenceNumber = 0;                                                                                        \
    }

namespace DebugImageWrite
{
extern int fileSequenceNumber;
inline std::string
twodigits( unsigned int x )
{
  std::string rval;

  rval = '0' + static_cast< char >( x / 10 );
  rval += '0' + static_cast< char >( x % 10 );
  return rval;
}

extern int fileSequenceNumber;

template < typename ImageType >
void
DebugOutput( int LINE, const char * FILE, typename ImageType::Pointer img, int imageIndex = -1, const char * name = 0 )
{
  std::string fname( FILE );

  fname = itksys::SystemTools::GetFilenameName( fname );
  std::stringstream filename;
  filename << "DBG_";
  filename << twodigits( fileSequenceNumber );
  filename << "_";
  filename << name;
  if ( imageIndex != -1 )
  {
    filename << twodigits( fileSequenceNumber );
  }
  filename << "_" << fname << "_" << LINE << ".nii.gz";
  std::cerr << "Writing " << filename.str() << " " << img.GetPointer() << std::endl;
  itkUtil::WriteImage< ImageType >( img, filename.str() );
  ++fileSequenceNumber;
}
} // namespace DebugImageWrite

#  define DebugOutput( imageType, image )                                                                              \
    DebugImageWrite::DebugOutput< imageType >( __LINE__, __FILE__, image, -1, #image )
#  define DebugOutputN( imageType, image, N, name )                                                                    \
    DebugImageWrite::DebugOutput< imageType >( __LINE__, __FILE__, image, N, #name )
#  define DebugOutputWName( imageType, image, name )                                                                   \
    DebugImageWrite::DebugOutput< imageType >( __LINE__, __FILE__, image, -1, #name )

#else
#  define DEFINE_DEBUG_IMAGE_COUNTER
#  define DebugOutput( imageType, image )
#  define DebugOutputN( imageType, image, N, name )
#  define DebugOutputWName( imageType, image, name )
#endif // BRAINS_DEBUG_IMAGE_WRITE
#endif // __DebugImageWrite_h
