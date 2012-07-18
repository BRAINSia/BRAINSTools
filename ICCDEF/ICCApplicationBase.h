#ifndef _ICCApplicationBase_h
#define _ICCApplicationBase_h

#include "itkObjectFactory.h"
#include "itkObject.h"

namespace itk
{
/** \class ICCApplicationBase
 *
 * This class ties together an input parser, a preprocessor,
 * a registrator components to
 * form a deformable registration/atlas segmentation application.
 *
 */
template <typename TPreprocessor,
          typename TRegistrator>
class ICCApplicationBase : public Object
{
public:

  /** Standard class typedefs. */
  typedef ICCApplicationBase       Self;
  typedef Object                   Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MIMApplication, Object);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Preprocessor type. */
  typedef TPreprocessor                      PreprocessorType;
  typedef typename PreprocessorType::Pointer PreprocessorPointer;

  /** Registrator type. */
  typedef TRegistrator                      RegistratorType;
  typedef typename RegistratorType::Pointer RegistratorPointer;

  /**Set Debug mode*/
  itkSetMacro(OutDebug, bool);
  itkGetConstMacro(OutDebug, bool);

  RegistratorType * GetRegistratorType(void)
  {
    return m_Registrator;
  }

  /** Execute the application. */
  virtual void Execute();

protected:

  ICCApplicationBase();
  virtual ~ICCApplicationBase()
  {
  }

  /*** Initialize the preprocessor */
  virtual void InitializePreprocessor()
  {
  }

  /*** Initialize the registrator  */
  virtual void InitializeRegistrator()
  {
  }

  PreprocessorPointer m_Preprocessor;
  RegistratorPointer  m_Registrator;
  bool                m_OutDebug;
};
}   // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "ICCApplicationBase.txx"
#endif

#endif
