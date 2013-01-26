/*=================================================
  Program: itkBrains2LandmarkReader.h
  Date:    2009-10-22 10:26
  Version: 1.0
  Author:  Yongqiang Zhao
==================================================*/

#ifndef __itkBrains2LandmarkReader_h
#define __itkBrains2LandmarkReader_h

#include "itkConfigure.h"

#include "itkLightProcessObject.h"
#include "itkPointSet.h"
#include "itkImage.h"
#include <fstream>
#include <string>

// Notice: The reference image should have the identity direction.
// Otherwise, the result of GetPhysicalPoint() is not correct.

namespace itk
{
template <class TPixelType, unsigned Dimension>
class Brains2LandmarkReader : public LightProcessObject
{
public:
  typedef Brains2LandmarkReader    Self;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  typedef float PixelType;
  typedef float PointType;

  typedef PointSet<PointType, Dimension>  PointSetType;
  typedef PointSet<TPixelType, Dimension> InputPointSetType;

  typedef Image<PixelType, Dimension> ImageType;
  typedef typename ImageType::Pointer ImagePointer;

  itkNewMacro(Self);
  typedef Object Superclass;
  itkTypeMacro(Brains2LandmarkReader, LightProcessObject);

  itkSetStringMacro(FileName);
  itkGetStringMacro(FileName);

  PointSetType * GetPointSet()
  {
    return m_PointSet;
  }

  PointSetType * GetPhysicalPointSet()
  {
    return m_PPointSet;
  }

  void SetReferenceImage(ImageType * ig)
  {
    m_ReferenceImage = ig;
  }

  void Update();

protected:
  Brains2LandmarkReader(const Self &);
  Brains2LandmarkReader & operator=(const Self &);

  Brains2LandmarkReader();
  virtual ~Brains2LandmarkReader()
  {
  };
//   void GenerateData();
private:
  std::string m_FileName;
  typename InputPointSetType::Pointer  m_PointSet;
  typename PointSetType::Pointer  m_PPointSet;
  ImagePointer m_ReferenceImage;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBrains2LandmarkReader.hxx"
#endif

#endif // _itkBrains2LandmarkReader_h
