#ifndef __vtkTalairachConversion_h
#define __vtkTalairachConversion_h

#include "vtkDataObjectAlgorithm.h"
#include "vtkTalairachGrid.h"
#include "vtkObjectFactory.h"
#include "itkImage.h"
#include <list>
#include <string>

class vtkStructuredGrid;

class VTK_FILTERING_EXPORT vtkTalairachConversion : public vtkDataObjectAlgorithm
{
public:
  static vtkTalairachConversion * New();

  /* Declare image type */
  typedef unsigned char PixelType;
  static const int Dimension = 3;
  typedef itk::Image<PixelType, Dimension> ImageType;

  /* Define the Hemisphere Types */
  enum { right = 0, left = 1, both = 2 };

  /* Description:
   * Clear out the memory */
  void Initialize();

  /* Description:
   * Write out the appropriate data when the talairachConversion object
   * is added to an IO stream */
  void PrintSelf(ostream & os, vtkIndent indent);

  /* Description:
   * Set the origin, direction, etc. of the mask image from
   * an input image */
  void SetImageInformation(ImageType::Pointer exampleImage);

  /* Description:
   * Create a binary mask image from all talairach grids in the grid list
   * that are currently turned on; size, origin and spacing should be set
   * prior to running this function */
  void GenerateImage();

  /* Description:
   * Returns the binary mask image generated from the talairach grids;
   * returned as an ITK image */
  ImageType::Pointer GetImage();

  /* Description:
   * Set the talairach grid */
  void SetTalairachGrid( vtkStructuredGrid *grid);

  /* Description:
   * Return the talairach grid */
  vtkStructuredGrid * GetTalairachGrid();

  /* Description:
   * Add a talairach grid range to the binary mask image */
  int AddTalairachBox(std::string talairachBox);

  /* Description:
   * Remove a talairach grid range from the binary mask image */
  void RemoveTalairachBox(int index);

  void RemoveTalairachBox(std::string talairachBox);

  /* Description:
   * Clear the Talairach Box List */
  void EraseTalairachBoxList();

  /* Description:
   * Set the Mode for Talairach Box Generation */
  void SetHemisphereMode(int mode);

  void SetHemisphereModeBoth();

  void SetHemisphereModeLeft();

  void SetHemisphereModeRight();

  /* Description:
   * Get the Mode for Talairach Box Generation */
  int GetHemisphereMode();

  /* Description:
   * Set the Mode for Expanded Segmentation Mode Boxes */
  void SetSegmentationMode(bool mode);

  void SetSegmentationModeOn();

  void SetSegmentationModeOff();

  /* Description:
   * Get the Mode for Talairach Box Generation */
  bool GetSegmentationMode();

  /* Description:
   * Returns the number of boxes in the list */
  int GetTalairachBoxLength();

  /* Description:
   * Update the image */
  void Update();

protected:

  vtkTalairachConversion();
  ~vtkTalairachConversion();

  /* Description:
   * Process a box file to calculate the regions of active masking */
  void ProcessBOX(bool _left);

private:

  vtkStructuredGrid *    TalairachGrid;
  std::list<std::string> TalairachBoxList;

  bool SegmentationMode;
  int  HemisphereMode;

  ImageType::Pointer MaskImage;
};

#endif
