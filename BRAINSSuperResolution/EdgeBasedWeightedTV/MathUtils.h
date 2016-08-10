//
// Created by Johnson, Hans J on 7/15/16.
//
#ifndef SR_MATHUTILS_H_H
#define SR_MATHUTILS_H_H

#include <itkBinaryFunctorImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkPowImageFilter.h>
#include <itkExpImageFilter.h>
#include <itkSqrtImageFilter.h>
#include <itkSquareImageFilter.h>

#if 0

namespace Functor
{
/**
 * \class AddScaledVersions
 * \brief  Implement
 * \ingroup ITKImageIntensity
 */
    template< typename TInput1, typename TInput2 = TInput1, typename TOutput = TInput1 >
    class AddScaledVersions
    {
    public:
      AddScaledVersions() {}
      ~AddScaledVersionsMult() {}
      bool operator!=(const AddScaledVersionsMult &) const
      {
        return false;
      }

      bool operator==(const AddScaledVersionsMult & other) const
      {
        return !( *this != other );
      }

      inline TOutput operator()(const TInput1 & A, const TInput2 & B) const
      { return static_cast<TOutput>( m_Ascaler*A + m_Bscaler*B ); }
      private:
         double m_Ascaler;
         double m_Bscaler;
         double m_OveralScaler;
    };
}


template< typename TInputImage1, typename TInputImage2 = TInputImage1, typename TOutputImage = TInputImage1 >
class AddScaledVersions:
    public
    BinaryFunctorImageFilter< TInputImage1, TInputImage2, TOutputImage,
        Functor::Mult<
        typename TInputImage1::PixelType,
        typename TInputImage2::PixelType,
        typename TOutputImage::PixelType >   >
#endif



//If image types don't match, always return false!
template<typename TPInputImage1, typename TPInputImage2>
bool isSameImage(TPInputImage1 , TPInputImage2 ) {
  return false;
}

template<typename TPInputImage1>
bool isSameImage(TPInputImage1 inImg1, TPInputImage1 inImg2) {
  return (inImg1.GetPointer() == inImg2.GetPointer());
}

/* Macro for generating image to image operations */
#define MAKE_MATH_OPERATION_II(OPERATION)                                                  \
{                                                                                          \
  typedef itk::OPERATION##ImageFilter<typename TPInputImage1::ObjectType,                  \
                                      typename TPInputImage2::ObjectType,                  \
                                      typename TPOutputImage::ObjectType> OPERATION##Type; \
  typename OPERATION##Type::Pointer do##OPERATION = OPERATION##Type::New();                \
  do##OPERATION->SetInput1(inImg1);                                                        \
  do##OPERATION->SetInput2(inImg2);                                                        \
  if ( inImg1.IsNotNull() && isSameImage(outImg, inImg1 ) )                                \
  {                                                                                        \
  do##OPERATION->InPlaceOn();                                                              \
  }                                                                                        \
  else                                                                                     \
  {                                                                                        \
  do##OPERATION->GraftOutput(outImg);                                                      \
  }                                                                                        \
  do##OPERATION->Update();                                                                 \
  outImg=do##OPERATION->GetOutput();                                                       \
}

/*
 * This uberfunction is designed to make doing Image math operations easier
 * NOTE: Returning the outImg is for convenience of chaning operations together
 * outImg = opII('*',in1,in2,outImg);    outImg = in1 * in2;
 *
 */
template<typename TPInputImage1, typename TPInputImage2, typename TPOutputImage>
TPOutputImage opII(
    TPOutputImage outImg,
    TPInputImage1 inImg1,
    const char op,
    TPInputImage2 inImg2) {
  switch (op) {
    case '+': MAKE_MATH_OPERATION_II(Add);
      break;
    case '-': MAKE_MATH_OPERATION_II(Subtract);
      break;
    case '*': MAKE_MATH_OPERATION_II(Multiply);
      break;
    default:
      std::cout << "ERROR: Invalid operator" << std::endl;
      assert(0 == 1); //ERROR bad operator
      break;
  }
  return outImg;
}

/*
 * This uberfunction is designed to make doing Image math operations easier
 * NOTE: Returning the outImg is for convenience of chaning operations together
 * outImg = opII('*',in1,in2,outImg);    outImg = in1 * in2;
 *
 */
template<typename TPInputImage1, typename TPInputImage2, typename TPOutputImage>
TPOutputImage opII_scalar(
    TPOutputImage outImg,
    TPInputImage1 inImg1,
    const char op,
    TPInputImage2 inImg2) {
  switch (op) {
    case '+':
    case '-':
    case '*':
      opII(outImg, inImg1, op, inImg2);
      break;
    case '/': MAKE_MATH_OPERATION_II(Divide);
      break;
    default:
      std::cout << "ERROR: Invalid operator" << std::endl;
      assert(0 == 1); //ERROR bad operator
      break;
  }
  return outImg;
}
/*
 * This uberfunction is designed to make doing Image math operations easier
 * NOTE: Returning the outImg is for convenience of chaning operations together
 * outImg = opII('*',in1,in2,outImg);    outImg = in1 * in2;
 *
 */
//Specialization for CV vector multiplication
template<typename TPInputImage1, typename TPInputImage2, typename TPOutputImage>
TPOutputImage opII_CVmult(
    TPOutputImage outImg,
    TPInputImage1 inImg1,
    const char op,
    TPInputImage2 inImg2) {
  switch (op) {
    case '*': MAKE_MATH_OPERATION_II(Multiply);
      break;
    default:
      std::cout << "ERROR: Invalid operator" << std::endl;
      assert(0 == 1); //ERROR bad operator
      break;
  }
  return outImg;
}

/* Macro for generating image to image operations */
#define MAKE_MATH_OPERATION_IC(OPERATION)                                                  \
{                                                                                          \
  typedef itk::OPERATION##ImageFilter<typename TPInputImage1::ObjectType,                  \
                                      typename TPInputImage1::ObjectType,                  \
                                      typename TPOutputImage::ObjectType> OPERATION##Type; \
  typename OPERATION##Type::Pointer do##OPERATION = OPERATION##Type::New();                \
  do##OPERATION->SetInput1(inImg1);                                                        \
  do##OPERATION->SetConstant2(constant2);                                                  \
  if ( isSameImage(outImg, inImg1 ) )                                                      \
  {                                                                                        \
  do##OPERATION->InPlaceOn();                                                              \
  }                                                                                        \
  else                                                                                     \
  {                                                                                        \
  do##OPERATION->GraftOutput(outImg);                                                      \
  }                                                                                        \
  do##OPERATION->Update();                                                                 \
  outImg=do##OPERATION->GetOutput();                                                       \
}


/* Macro for generating image to image operations */
#define MAKE_UNARY_MATH_OPERATION_IC(OPERATION)                                            \
{                                                                                          \
  typedef itk::OPERATION##ImageFilter<typename TPInputImage1::ObjectType,                  \
                                      typename TPOutputImage::ObjectType> OPERATION##Type; \
  typename OPERATION##Type::Pointer do##OPERATION = OPERATION##Type::New();                \
  do##OPERATION->SetInput(inImg1);                                                         \
  if ( isSameImage(outImg, inImg1 ) )                                                      \
  {                                                                                        \
  do##OPERATION->InPlaceOn();                                                              \
  }                                                                                        \
  else                                                                                     \
  {                                                                                        \
  do##OPERATION->GraftOutput(outImg);                                                      \
  }                                                                                        \
  do##OPERATION->Update();                                                                 \
  outImg=do##OPERATION->GetOutput();                                                       \
}

/*
 * This uberfunction is designed to make doing Image math operations easier
 * NOTE: Returning the outImg is for convenience of chaning operations together
 * outImg = opII('*',in1,in2,outImg);    outImg = in1 * in2;
 *
 */
/* outImg = inImg1 'op' constant2 */
template<typename TPInputImage1, typename TPOutputImage>
TPOutputImage opIC(TPOutputImage outImg,
                   TPInputImage1 inImg1,
                   const char op,
                   const typename TPInputImage1::ObjectType::PixelType constant2) {
  switch (op) {
    case '+': MAKE_MATH_OPERATION_IC(Add); // I+C
      break;
    case '-': MAKE_MATH_OPERATION_IC(Subtract); //I-C
      break;
    case '*': MAKE_MATH_OPERATION_IC(Multiply); // I*C
      break;
    default:
      std::cout << "ERROR: Invalid operator: " << op << std::endl;
      assert(0 == 1); //ERROR bad operator
      break;
  }
  return outImg;
}

template<typename TPInputImage1, typename TPOutputImage>
TPOutputImage opIC(TPOutputImage outImg,
                   const typename TPInputImage1::ObjectType::PixelType constant2,
                   const char op,
                   TPInputImage1 inImg1) {
  return opIC(outImg, inImg1, op, constant2);
};



/*
 * This uberfunction is designed to make doing Image math operations easier
 * NOTE: Returning the outImg is for convenience of chaning operations together
 * outImg = opII('*',in1,in2,outImg);    outImg = in1 * in2;
 *
 */
// Specialization for extra math operations that only apply
// to scalar images.  For example std::sqrt(CovariantVector) is not valid.
template<typename TPInputImage1, typename TPOutputImage>
TPOutputImage opIC_scalar(
    TPOutputImage outImg,
    const std::string opS,
    TPInputImage1 inImg1) {
  const char op = opS[0];
  switch (op) {
    case 's': MAKE_UNARY_MATH_OPERATION_IC(Sqrt);
      break;
    case 'e': MAKE_UNARY_MATH_OPERATION_IC(Exp);
      break;
    default:
      std::cout << "ERROR: Invalid operator: " << op << std::endl;
      assert(0 == 1); //ERROR bad operator
      break;
  }
  return outImg;
}

/*
 * This uberfunction is designed to make doing Image math operations easier
 * NOTE: Returning the outImg is for convenience of chaning operations together
 * outImg = opII('*',in1,in2,outImg);    outImg = in1 * in2;
 *
 */
// Specialization for extra math operations that only apply
// to scalar images.  For example std::sqrt(CovariantVector) is not valid.
template<typename TPInputImage1, typename TPOutputImage>
TPOutputImage opIC_scalar(
    TPOutputImage outImg,
    TPInputImage1 inImg1,
    const char op,
    const typename TPInputImage1::ObjectType::PixelType constant2) {
  switch (op) {
    case '+':
    case '-':
    case '*':
      opIC(outImg, inImg1, op, constant2);
      break;
    case '^':
      if (constant2 == 0.5F) {
        opIC_scalar(outImg, "sqrt", inImg1);
      }
      else {
        MAKE_MATH_OPERATION_IC(Pow);
      }
      break;
    default:
      std::cout << "ERROR: Invalid operator: " << op << " " << __FILE__ << " " << __LINE__ << std::endl;
      assert(0 == 1); //ERROR bad operator
      break;
  }
  return outImg;
}

//Use Commutative property for easing migration
template<typename TPInputImage1, typename TPOutputImage>
TPOutputImage opIC_scalar(
    TPOutputImage outImg,
    const typename TPInputImage1::ObjectType::PixelType constant2,
    const char op,
    TPInputImage1 inImg1) {
  return opIC_scalar(outImg, inImg1, op, constant2);//
};


#endif //SR_MATHUTILS_H_H
