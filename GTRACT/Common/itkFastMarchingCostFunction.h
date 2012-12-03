/*=========================================================================

 Program:   GTRACT (Guided Tensor Restore Anatomical Connectivity Tractography)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2006/03/29 14:53:40 $
 Version:   $Revision: 1.9 $

   Copyright (c) University of Iowa Department of Radiology. All rights reserved.
   See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
   for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkFastMarchingCostFunction_h
#define __itkFastMarchingCostFunction_h

#include "itkSingleValuedCostFunction.h"

#include <vnl/vnl_math.h>
#include <iostream>
#include <fstream>
#include <itkImage.h>
#include <itkObject.h>
#include <itkImageToImageFilter.h>
#include <itkIOCommon.h>
#include "itkProcessObject.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "GtractTypes.h"

#include <map>
#include <string>

#include "itkIndex.h"
#include "vnl/vnl_math.h"

namespace itk
{
class FastMarchingCostFunction : public SingleValuedCostFunction
{
public:
  /** Standard class typedefs */
  typedef FastMarchingCostFunction Self;
  typedef SingleValuedCostFunction Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(FastMarchingCostFunction, SingleValuedCostFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Some convenient typedefs. */
  typedef itk::Image<float, 3>         CostImageType;
  typedef CostImageType::Pointer       CostImagePointer;
  typedef CostImageType::ConstPointer  CostImageConstPointer;
  typedef CostImageType::RegionType    CostImageRegionType;
  typedef CostImageType::SizeType      CostImageSizeType;
  typedef CostImageType::SpacingType   CostImageSpacingType;
  typedef CostImageType::PointType     CostImagePointType;
  typedef CostImageType::PixelType     CostImagePixelType;
  typedef CostImageType::IndexType     CostImageIndexType;
  typedef CostImageType::DirectionType CostImageDirectionType;

  typedef itk::LinearInterpolateImageFunction<CostImageType, float> CostIPType;      //
                                                                                     //
                                                                                     // ScalarIPType;
  typedef  CostIPType::Pointer CostIPTypePointer;

  /** ImageDimension constants */
  itkStaticConstMacro(CostImageDimension, unsigned int, 3);

  /*  A position in the optimization space. */
  typedef Superclass::ParametersType ParametersType;

  /*  vcl_cost function derivative (gradient). */
  typedef Superclass::DerivativeType DerivativeType;

  /*  vcl_cost function value. */
  typedef Superclass::MeasureType MeasureType;

  /** Set/Get input Cost Image  */
  itkSetObjectMacro( CostImage, CostImageType );
  itkGetConstObjectMacro( CostImage, CostImageType );

  // Returns dimension of image
  unsigned int GetNumberOfParameters() const;

  /** This method returns the value of the vcl_cost function for
    * the specified parameters, or position. */
  MeasureType GetValue( const ParametersType & parameters ) const;

  /** This method returns the derivative of the vcl_cost function corresponding
    * to the specified parameters.   */
  void GetDerivative( const ParametersType & parameters, DerivativeType & derivative ) const;

protected:
  FastMarchingCostFunction();
  // virtual ~FastMarchingCostFunction(){};
  ~FastMarchingCostFunction()
  {
  }

  void PrintSelf(std::ostream & os, Indent indent) const;

  // void SetMetaDataHeader();

  CostImagePointer    m_CostImage;
  CostImageRegionType m_BufferedRegion;
  CostIPTypePointer   m_CostIP;
private:

  FastMarchingCostFunction(const Self &); // purposely not implemented
  void operator=(const Self &);           // purposely not implemented
};                                        // end of class
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFastMarchingCostFunction.cxx"
#endif

#endif
