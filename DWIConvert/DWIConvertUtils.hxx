
#ifndef DWIConvertUtils_HXX
#define DWIConvertUtils_HXX

template <typename TArg>
int
CheckArg(const char * argName, const TArg & argVal, const TArg & emptyVal)
{
  if (argVal == emptyVal)
  {
    std::cerr << "Missing argument " << argName << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

template <typename TImage>
int
WriteVolume(const TImage * img, const std::string & fname)
{
  typename itk::ImageFileWriter<TImage>::Pointer imgWriter = itk::ImageFileWriter<TImage>::New();

  imgWriter->SetInput(img);
  imgWriter->SetFileName(fname.c_str());
  try
  {
    imgWriter->Update();
  }
  catch (itk::ExceptionObject & excp)
  {
    std::cerr << "Exception thrown while writing " << fname << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

template <typename TImage>
int
ReadScalarVolume(typename TImage::Pointer & img, const std::string & fname, bool allowLossyConversion)
{
  typename itk::ImageFileReader<TImage>::Pointer imgReader = itk::ImageFileReader<TImage>::New();

  imgReader->SetFileName(fname.c_str());
  try
  {
    imgReader->Update();
  }
  catch (itk::ExceptionObject & excp)
  {
    std::cerr << "Exception thrown while reading " << fname << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  img = imgReader->GetOutput();

  {
    itk::ImageIOBase *                imageIO = imgReader->GetModifiableImageIO();
    itk::ImageIOBase::IOComponentEnum ioType = itk::ImageIOBase::MapPixelType<typename TImage::PixelType>::CType;
    if (!allowLossyConversion)
    {
      if (imageIO->GetComponentType() != ioType)
      {
        std::cerr << "Error: ReadVolume: Unsupported source pixel type." << std::endl
                  << "  Input volume:  " << imageIO->GetComponentTypeAsString(imageIO->GetComponentType()) << std::endl
                  << "  Output volume: " << imageIO->GetComponentTypeAsString(ioType) << std::endl
                  << "The only supported output type is <short>. "
                  << "You may consider using allowLossyConversion option." << std::endl
                  << "However, use this option with caution! "
                  << "Conversion from images of a different type may cause data loss due to rounding or truncation."
                  << std::endl;
        return EXIT_FAILURE;
      }
    }
    if (imageIO->GetComponentType() != ioType)
    {
      using DoubleImageType = itk::Image<double, TImage::ImageDimension>;
      using DoubleImageReaderType = itk::ImageFileReader<DoubleImageType>;
      typename DoubleImageReaderType::Pointer doubleReader = DoubleImageReaderType::New();
      doubleReader->SetFileName(fname.c_str());
      try
      {
        imgReader = nullptr; // Throw away existing reader (save memory)
        img = nullptr;       // Throw away existing version of image (save memory)
        doubleReader->Update();
      }
      catch (itk::ExceptionObject & excp)
      {
        std::cerr << "Exception thrown while reading " << fname << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
      }
      using RescaleIntensityType = itk::RescaleIntensityImageFilter<DoubleImageType, TImage>;
      typename RescaleIntensityType::Pointer rescaler = RescaleIntensityType::New();
      rescaler->SetInput(doubleReader->GetOutput());
      rescaler->SetOutputMinimum(itk::NumericTraits<typename TImage::PixelType>::Zero);
      rescaler->SetOutputMaximum(itk::NumericTraits<typename TImage::PixelType>::max());
      rescaler->Update();
      img = rescaler->GetOutput();
    }
  }
  return EXIT_SUCCESS;
}

template <typename TImage>
int
ReadVectorVolume(typename TImage::Pointer & img, const std::string & fname, bool allowLossyConversion)
{
  typename itk::ImageFileReader<TImage>::Pointer imgReader = itk::ImageFileReader<TImage>::New();

  imgReader->SetFileName(fname.c_str());
  try
  {
    imgReader->Update();
  }
  catch (itk::ExceptionObject & excp)
  {
    std::cerr << "Exception thrown while reading " << fname << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  img = imgReader->GetOutput();

  {
    itk::ImageIOBase *                imageIO = imgReader->GetModifiableImageIO();
    itk::ImageIOBase::IOComponentEnum ioType = itk::ImageIOBase::MapPixelType<typename TImage::PixelType>::CType;
    if (!allowLossyConversion)
    {
      if (imageIO->GetComponentType() != ioType)
      {
        std::cerr << "Error: ReadVolume: Unsupported source pixel type." << std::endl
                  << "  Input volume:  " << imageIO->GetComponentTypeAsString(imageIO->GetComponentType()) << std::endl
                  << "  Output volume: " << imageIO->GetComponentTypeAsString(ioType) << std::endl
                  << "The only supported output type is <short>. "
                  << "You may consider using allowLossyConversion option." << std::endl
                  << "However, use this option with caution! "
                  << "Conversion from images of a different type may cause data loss due to rounding or truncation."
                  << std::endl;
        return EXIT_FAILURE;
      }
    }
    if (imageIO->GetComponentType() != ioType)
    {
      using DoubleImageType = itk::VectorImage<double, TImage::ImageDimension>;
      using DoubleImageReaderType = itk::ImageFileReader<DoubleImageType>;
      typename DoubleImageReaderType::Pointer doubleReader = DoubleImageReaderType::New();
      doubleReader->SetFileName(fname.c_str());
      try
      {
        imgReader = nullptr; // Throw away existing reader (save memory)
        img = nullptr;       // Throw away existing version of image (save memory)
        doubleReader->Update();
      }
      catch (itk::ExceptionObject & excp)
      {
        std::cerr << "Exception thrown while reading " << fname << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
      }
      typename DoubleImageType::Pointer dimg = doubleReader->GetOutput();
      double                            dmin = itk::NumericTraits<double>::max();
      double                            dmax = itk::NumericTraits<double>::min();
      const size_t                      vectorSize = dimg->GetVectorLength();
      {
        itk::ImageRegionConstIterator<DoubleImageType> it(dimg, dimg->GetLargestPossibleRegion());
        for (it.GoToBegin(); it.IsAtEnd(); ++it)
        {
          const typename DoubleImageType::PixelType & tmp = it.Get();
          for (size_t i = 0; i < vectorSize; ++i)
          {
            const double & value = tmp[i];
            dmin = std::min(dmin, value);
            dmax = std::max(dmax, value);
          }
        }
      }
      const double scale = (itk::NumericTraits<typename TImage::PixelType::ValueType>::max() - 1) / (dmax - dmin);

      itk::ImageRegionIterator<DoubleImageType> it(dimg, dimg->GetLargestPossibleRegion());
      itk::ImageRegionIterator<TImage>          outit(img, img->GetLargestPossibleRegion());
      for (it.GoToBegin(), outit.GoToBegin(); it.IsAtEnd(); ++it, ++outit)
      {
        typename DoubleImageType::PixelType tmp = it.Get();
        typename TImage::PixelType          outtmp = outit.Get();
        for (size_t i = 0; i < vectorSize; ++i)
        {
          const double & value = tmp[i];
          outtmp[i] = (value - dmin) * scale;
        }
        outit.Set(outtmp);
      }
    }
  }
  return EXIT_SUCCESS;
}


template <typename TImage>
int
RecoverMeasurementFrame(const TImage * img, DWIMetaDataDictionaryValidator::RotationMatrixType & MeasurementFrame)
{

  DWIMetaDataDictionaryValidator myDWIValidator;
  myDWIValidator.SetMetaDataDictionary(img->GetMetaDataDictionary());
  MeasurementFrame = myDWIValidator.GetMeasurementFrame();
  return EXIT_SUCCESS;
}

template <typename TImage>
int
RecoverBVectors(const TImage * img, DWIMetaDataDictionaryValidator::GradientTableType & bVecs)
{
  bVecs.clear();

  DWIMetaDataDictionaryValidator myDWIValidator;
  myDWIValidator.SetMetaDataDictionary(img->GetMetaDataDictionary());
  bVecs = myDWIValidator.GetGradientTable();
  if (bVecs.empty())
  {
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

template <typename TImage>
int
RecoverBValue(const TImage * img, double & val)
{
  DWIMetaDataDictionaryValidator myDWIValidator;
  myDWIValidator.SetMetaDataDictionary(img->GetMetaDataDictionary());
  val = myDWIValidator.GetBValue();
  return EXIT_SUCCESS;
}

template <typename TImage>
int
RecoverBValues(const TImage *                                            inputVol,
               const DWIMetaDataDictionaryValidator::GradientTableType & bVectors,
               std::vector<double> &                                     bValues)
{
  bValues.clear();

  double BValue;
  if (RecoverBValue(inputVol, BValue) != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }
  for (const auto & bVector : bVectors)
  {
    double norm = std::sqrt((bVector[0] * bVector[0]) + (bVector[1] * bVector[1]) + (bVector[2] * bVector[2]));
    if (std::abs(1.0 - norm) < 0.05) // Asssume norm is 1 if with 5% of 1 to stablize NRRD->FSL->NRRD numerical
                                     // stability.
    {
      norm = 1.0;
    }
    // bval_i = (G_norm)^2 * bval_max
    double bval = norm * BValue;
    if (std::abs(bval - itk::Math::Round<double>(bval)) < 1e-2)
    {
      bval = itk::Math::Round<double>(bval);
    }
    bValues.push_back(bval);
  }
  return EXIT_SUCCESS;
}


template <typename TValue>
bool
CloseEnough(const TValue & a, const TValue & b, double magdiv)
{
  double averageMag = (std::fabs(static_cast<double>(a)) + std::fabs(static_cast<double>(b))) / 2.0;
  double diff = std::fabs(static_cast<double>(a) - static_cast<double>(b));

  // case one -- both near zero
  if (averageMag < 0.000001)
  {
    return true;
  }
  // case 2 -- diff > average / 100000;
  if (diff > (averageMag / magdiv))
  {
    return false;
  }
  return true;
}

template <typename TVal>
bool
CloseEnough(const std::vector<TVal> & a, const std::vector<TVal> & b, double magdiv)
{
  if (a.size() != b.size())
  {
    std::cerr << "Vector size mismatch: " << a.size() << " " << b.size() << std::endl;
    return false;
  }
  for (unsigned i = 0; i < a.size(); ++i)
  {
    if (!CloseEnough(a[i], b[i], magdiv))
    {
      std::cerr << "Value mismatch" << std::endl;
      return false;
    }
  }
  return true;
}


template <typename TVal>
bool
CloseEnough(const itk::VariableLengthVector<TVal> & a, const itk::VariableLengthVector<TVal> & b, double magdiv)
{
  if (a.GetSize() != b.GetSize())
  {
    std::cerr << "Vector size mismatch: " << a.GetSize() << " " << b.GetSize() << std::endl;
    return false;
  }
  for (unsigned i = 0; i < a.GetSize(); ++i)
  {
    if (!CloseEnough(a[i], b[i], magdiv))
    {
      std::cerr << "Value mismatch" << std::endl;
      return false;
    }
  }
  return true;
}


template <typename TVal>
void
PrintVec(const TVal & a)
{
  std::cerr << a;
}

template <typename TVal>
void
PrintVec(const std::vector<TVal> & vec)
{
  std::cerr << "[";

  for (unsigned i = 0; i < vec.size(); ++i)
  {
    PrintVec(vec[i]);
    if (i < vec.size() - 1)
    {
      std::cerr << " ";
    }
  }
  std::cerr << "]" << std::endl;
}

#endif // DWIConvertUtils
