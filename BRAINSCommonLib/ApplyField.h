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
/**
  * \defgroup AF Apply Field
  * \ingroup Reg
  */
#ifndef __ApplyField_h
#define __ApplyField_h

#include "itkObjectFactory.h"
#include "itkObject.h"

namespace itk
{
template <typename TDisplacementField, typename TInputImage,
          typename TOutputImage>
class ApplyField : public Object
{
public:
  typedef ApplyField               Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  itkTypeMacro(MIMApplication, Object);

  itkNewMacro(Self);

  typedef TInputImage                         InputImageType;
  typedef TOutputImage                        OutputImageType;
  typedef typename OutputImageType::PixelType PixelType;
  typedef typename InputImageType::Pointer    ImagePointer;

  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

  itkSetObjectMacro(InputImage, InputImageType);
  itkGetConstObjectMacro(InputImage, InputImageType);
  itkGetConstObjectMacro(OutputImage, OutputImageType);

  /** Set/Get value to replace thresholded pixels. Pixels that lie *
    *  within Lower and Upper (inclusive) will be replaced with this
    *  value. The default is 1. */
  itkSetMacro(DefaultPixelValue,  PixelType);
  itkGetMacro(DefaultPixelValue,  PixelType);

  itkSetObjectMacro(DisplacementField, TDisplacementField);

  void Execute();

  void ReleaseDataFlagOn();

protected:
  ApplyField();
  ~ApplyField() override;
private:
  typename InputImageType::Pointer m_InputImage;
  typename OutputImageType::Pointer m_OutputImage;
  typename TDisplacementField::Pointer m_DisplacementField;
  PixelType m_DefaultPixelValue;
};
}
#ifndef ITK_MANUAL_INSTANTIATION
#include "ApplyField.hxx"
#endif
#endif
