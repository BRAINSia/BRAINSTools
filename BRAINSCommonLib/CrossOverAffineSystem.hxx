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
#ifndef __CrossOverAffineSystem_hxx
#define __CrossOverAffineSystem_hxx

#include "itkNumericTraits.h"
#include "CrossOverAffineSystem.h"

/**
 * Constructor
 */
template < typename TCoordinateType, unsigned int NDimensions >
CrossOverAffineSystem< TCoordinateType, NDimensions >::CrossOverAffineSystem()
  : m_InhaleEncodeConversion()
  , m_InhaleDecodeConversion()
  , m_ExhaleEncodeConversion()
  , m_ExhaleDecodeConversion()
{
  m_InhaleEncodeConversion = AffineTransformType::New();
  m_InhaleDecodeConversion = AffineTransformType::New();
  m_ExhaleEncodeConversion = AffineTransformType::New();
  m_ExhaleDecodeConversion = AffineTransformType::New();

  m_InhaleEncodeConversion->SetIdentity();
  m_InhaleDecodeConversion->SetIdentity();
  m_ExhaleEncodeConversion->SetIdentity();
  m_ExhaleDecodeConversion->SetIdentity();
}

/**
 * Destructor
 */
template < typename TCoordinateType, unsigned int NDimensions >
CrossOverAffineSystem< TCoordinateType, NDimensions >::~CrossOverAffineSystem()
{
  return;
}

/**
 * Utility function not provided in general vector implementations.
 */
template < typename TCoordinateType, unsigned int NDimensions >
typename CrossOverAffineSystem< TCoordinateType, NDimensions >::VectorType
Reciprocal( const typename CrossOverAffineSystem< TCoordinateType, NDimensions >::VectorType & Scale )
{
  typedef typename CrossOverAffineSystem< TCoordinateType, NDimensions >::VectorType VectorType;
  VectorType                                                                         ReciprocalScale( Scale );
  for ( unsigned int i = 0; i < NDimensions; ++i )
  {
    ReciprocalScale[i] = 1.0 / ReciprocalScale[i];
  }
  return ReciprocalScale;
}

// #define VERBOSE_PRINTING 1

/**
 * Build up the Inhale and Exhale transform enclosures by wrapping as follows:
 * (ResliceScaleReciprocal * InhaleEncode) * T * (InhaleDecode * StandardScale)
 * (ResliceScale * ExhaleEncode) * T * (ExhaleDecode * StandardScaleReciprocal)
 * So in EstablishCrossOverSystemForAir16() in TransformAdaptor.hxx,
 * the two arguments are to be given from the 'inhale' point of view.
 */
template < typename TCoordinateType, unsigned int NDimensions >
void
CrossOverAffineSystem< TCoordinateType, NDimensions >::EncloseInScaling( const VectorType & EncodeScale,
                                                                         const VectorType & DecodeScale )
{
#ifdef VERBOSE_PRINTING
  std::cout << "Encode Scale: " << EncodeScale << std::endl;
  std::cout << "Decode Scale: " << DecodeScale << std::endl;
#endif

  const bool EncodeApplyUpstream = false;
  const bool DecodeApplyUpstream = true;

  m_InhaleEncodeConversion->Scale( EncodeScale, EncodeApplyUpstream );
  m_InhaleDecodeConversion->Scale( DecodeScale, DecodeApplyUpstream );

  VectorType ReciprocalEncodeScale( Reciprocal< TCoordinateType, NDimensions >( EncodeScale ) );

  VectorType ReciprocalDecodeScale( Reciprocal< TCoordinateType, NDimensions >( DecodeScale ) );

  m_ExhaleEncodeConversion->Scale( ReciprocalEncodeScale, EncodeApplyUpstream );
  m_ExhaleDecodeConversion->Scale( ReciprocalDecodeScale, EncodeApplyUpstream );
}

/**
 * Build up the Inhale and Exhale transform enclosures by wrapping as follows:
 * (-ResliceShift * InhaleEncode) * T * (InhaleDecode * StandardShift)
 * (ResliceShift * ExhaleEncode) * T * (ExhaleDecode * -StandardShift)
 * So in EstablishCrossOverSystemForAir16() in TransformAdaptor.hxx,
 * the two arguments are to be given from the 'inhale' point of view.
 */
template < typename TCoordinateType, unsigned int NDimensions >
void
CrossOverAffineSystem< TCoordinateType, NDimensions >::EncloseInTranslation( const VectorType & EncodeShift,
                                                                             const VectorType & DecodeShift )
{
#ifdef VERBOSE_PRINTING
  std::cout << "Encode Shift: " << EncodeShift << std::endl;
  std::cout << "Decode Shift: " << DecodeShift << std::endl;
#endif

  const bool EncodeApplyUpstream = false;
  const bool DecodeApplyUpstream = true;

  m_InhaleEncodeConversion->Translate( EncodeShift, EncodeApplyUpstream );
  m_InhaleDecodeConversion->Translate( DecodeShift, DecodeApplyUpstream );

  m_ExhaleEncodeConversion->Translate( -EncodeShift, EncodeApplyUpstream );
  m_ExhaleDecodeConversion->Translate( -DecodeShift, DecodeApplyUpstream );
#ifdef VERBOSE_PRINTING
  std::cout << " -- InhaleEncoder Shift: " << m_InhaleEncodeConversion->GetOffset() << std::endl;
  std::cout << " -- InhaleDecoder Shift: " << m_InhaleDecodeConversion->GetOffset() << std::endl;
  std::cout << " -- ExhaleEncoder Shift: " << m_ExhaleEncodeConversion->GetOffset() << std::endl;
  std::cout << " -- ExhaleDecoder Shift: " << m_ExhaleDecodeConversion->GetOffset() << std::endl;
#endif
}

template < typename TCoordinateType, unsigned int NDimensions >
void
CrossOverAffineSystem< TCoordinateType, NDimensions >::EncloseInCentering( const PointType & EncodeCenter,
                                                                           const PointType & DecodeCenter )
{
#ifdef VERBOSE_PRINTING
  std::cout << "Encode Center: " << EncodeCenter << std::endl;
  std::cout << "Decode center: " << DecodeCenter << std::endl;
#endif

  m_InhaleEncodeConversion->SetCenter( EncodeCenter );
  m_InhaleDecodeConversion->SetCenter( DecodeCenter );

  m_ExhaleEncodeConversion->SetCenter( EncodeCenter );
  m_ExhaleDecodeConversion->SetCenter( DecodeCenter );
#ifdef VERBOSE_PRINTING
  std::cout << " -- InhaleEncoder Shift: " << m_InhaleEncodeConversion->GetOffset() << std::endl;
  std::cout << " -- InhaleDecoder Shift: " << m_InhaleDecodeConversion->GetOffset() << std::endl;
  std::cout << " -- ExhaleEncoder Shift: " << m_ExhaleEncodeConversion->GetOffset() << std::endl;
  std::cout << " -- ExhaleDecoder Shift: " << m_ExhaleDecodeConversion->GetOffset() << std::endl;
#endif
}

/**
 * Build up the Inhale and Exhale transform enclosures by wrapping in
 * two entire transforms in like fashion to the more useful routines above.
 * The two arguments are to be given from the 'inhale' point of view.
 */
template < typename TCoordinateType, unsigned int NDimensions >
void
CrossOverAffineSystem< TCoordinateType, NDimensions >::EncloseInAffineTransforms( AffineTransformPointer Encode,
                                                                                  AffineTransformPointer Decode )
{
  const bool EncodeApplyUpstream = false;
  const bool DecodeApplyUpstream = true;

  m_InhaleEncodeConversion->Compose( Encode, EncodeApplyUpstream );
  m_InhaleDecodeConversion->Compose( Decode, DecodeApplyUpstream );

  AffineTransformPointer EncodeInverse = AffineTransformType::New();
  AffineTransformPointer DecodeInverse = AffineTransformType::New();
  Encode->GetInverse( EncodeInverse );
  Decode->GetInverse( DecodeInverse );

  m_ExhaleEncodeConversion->Compose( EncodeInverse, EncodeApplyUpstream );
  m_ExhaleDecodeConversion->Compose( DecodeInverse, DecodeApplyUpstream );
}

#endif
