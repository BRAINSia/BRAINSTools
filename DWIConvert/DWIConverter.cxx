//
// Created by Johnson, Hans J on 12/3/16.
//

#include "DWIConverter.h"
#include "itkFlipImageFilter.h"

static bool DirectionNeedsFlipping(const Volume4DType::DirectionType & dir, const size_t ind)
{
  static const double FSLDesiredDirectionFlipsWRTLPS[4] = {1,-1,1,1};
  return ( FSLDesiredDirectionFlipsWRTLPS[ind]*dir(ind,ind) < -0.5 ); // i.e. a negative magnitude greater than 0.5
}

Volume4DType::Pointer DWIConverter::OrientForFSLConventions ()
{
  this->ConvertBVectorsToIdentityMeasurementFrame();
  this->ConvertToMutipleBValuesUnitScaledBVectors();


  Volume4DType::Pointer image4D = ThreeDToFourDImage(this->GetDiffusionVolume());
  Volume4DType::DirectionType direction=image4D->GetDirection();
  direction.GetVnlMatrix().get_row(0).magnitude();
  //LPS to RAI as FSL desires images to be formatted for viewing purposes.
  // This conversion makes FSLView display the images in
  // a way that is most easily interpretable.
  typedef itk::FlipImageFilter<Volume4DType> FlipperType;
  FlipperType::Pointer myFlipper = FlipperType::New();
  myFlipper->SetInput( image4D ) ;
  FlipperType::FlipAxesArrayType arrayAxisFlip;
  for(size_t i=0; i< Volume4DType::ImageDimension; ++i)
  {
    arrayAxisFlip[i] =  DirectionNeedsFlipping(direction,i);
    //This is necesssary to ensure that the BVEC file is consistent with FSL orientation assumptions
    const double second_dir_modifier = (i==2)? -1 : 1;
    for(size_t g =0 ; g < this->m_DiffusionVectors.size(); ++g)
    {
      this->m_DiffusionVectors[g][i]  *= arrayAxisFlip[i] ? second_dir_modifier * -1 : second_dir_modifier * 1 ;
    }
  }
  //
  // FSL wants the second and third dimensions flipped with regards to LPS orientation
  // FSL wants the second and third dimeinsions flipped with regards to LPS orientation
  myFlipper->SetFlipAxes(arrayAxisFlip);
  myFlipper->FlipAboutOriginOff();  //Flip the image and direction cosignes
  // this is similar to a transform of [1 0 0; 0 -1 0; 0 0 -1]
  myFlipper->Update();
  Volume4DType::Pointer temp = myFlipper->GetOutput();
  temp->SetMetaDataDictionary( image4D->GetMetaDataDictionary());
  this->m_Volume = FourDToThreeDImage(temp);
  return temp;
}
