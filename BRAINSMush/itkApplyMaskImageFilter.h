#ifndef __itkApplyMaskImageFilter_h
#define __itkApplyMaskImageFilter_h

#include <itkImageToImageFilter.h>

namespace itk
{
template <class TInputImage, class TOutputImage>
class ITK_EXPORT ApplyMaskImageFilter :
  public         ImageToImageFilter<TInputImage, TOutputImage>
{
public:    /* define methods available to everyone */

  /** Standard class typedefs */
  typedef ApplyMaskImageFilter                          Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** method for creation through the object factory */
  itkNewMacro(Self);

  /** run-time type information (and related methods) */
  itkTypeMacro(ApplyMaskImageFilter, ImageToImageFilter);

  /** ImageDimension enumeration */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** typedef to describe the output image region type */
  typedef typename TOutputImage::RegionType OutputImageRegionType;

  /** inherited typedefs */
  typedef typename Superclass::InputImageType         InputImageType;
  typedef typename Superclass::InputImagePointer      InputImagePointer;
  typedef typename Superclass::InputImageConstPointer InputImageConstPointer;
  typedef typename Superclass::OutputImageType        OutputImageType;
  typedef typename Superclass::OutputImagePointer     OutputImagePointer;

  /** pixel related typedefs */
  typedef typename InputImageType::PixelType  InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;

  itkSetMacro(InvertMask, bool);
  itkGetMacro(InvertMask, bool);

  /** Set/get the mask to be applied to the image */
  void SetMaskImage(const InputImageType *reference);

  const InputImageType * GetMaskImage(void);

protected: /* define methods available only to related classes */

  ApplyMaskImageFilter();
  ~ApplyMaskImageFilter()
  {
  }

  void PrintSelf(std::ostream & os, Indent indent) const;

  void GenerateData();

private:  /* define methods available only to this class */

  bool m_InvertMask;

  ApplyMaskImageFilter(const Self &);   // purposely not implemented
  void operator=(const Self &);         // purposely not implemented
};
}   // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkApplyMaskImageFilter.hxx"
#endif

#endif
