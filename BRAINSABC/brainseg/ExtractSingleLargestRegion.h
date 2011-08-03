#ifndef __ExtractSingleLargestRegion_h
#define __ExtractSingleLargestRegion_h

#include "itkImage.h"

extern itk::Image<unsigned char,
                  3>::Pointer ExtractSingleLargestRegionFromMask(itk::Image<unsigned char,
                                                                            3>::Pointer Mask, const int openingSize,
                                                                 const int closingSize, const int safetySize,
                                                                 itk::Image<unsigned char,
                                                                            3>
                                                                 ::Pointer inputLabelImage);

extern
itk::Image<unsigned char, 3>::Pointer ExtractSingleLargestRegion(const unsigned char threshold_low,
                                                                 const unsigned char threshold_high,
                                                                 const int openingSize, const int closingSize,
                                                                 const int safetySize, itk::Image<unsigned char,
                                                                                                  3>
                                                                 ::Pointer inputLabelImage);

#endif // ExtractSingleLargestRegion_h
