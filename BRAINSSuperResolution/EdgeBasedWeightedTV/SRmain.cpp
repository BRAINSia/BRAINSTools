// \author Hans J. Johnson
// \date 2016-07-10
// Test program for evaluating matlab to ITK conversions

#include <iostream>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "FFTWUpsample.h"
#include "OpWeightedL2.h"

#include <itkTimeProbe.h>

int main(int argc, char *argv[]) {
  if (argc != 4) {
    std::cout << "ERROR: Incorrrect number of arguments <Intensity_LR> <edgement_HR> <output>" << std::endl;
  }
  FFTWInit(""); //Just use the default in the home account
  itk::TimeProbe tp;
  tp.Start();
  typedef itk::ImageFileReader<FloatImageType> ReaderType;
  typedef itk::ImageFileWriter<FloatImageType> WriterType;

  const std::string edgeFileName = argv[2];
  ReaderType::Pointer edgeReader = ReaderType::New();
  edgeReader->SetFileName(edgeFileName);
  edgeReader->Update();
  FloatImageType::Pointer highResEdgeImage = edgeReader->GetOutput();

  const std::string lriFileName = argv[1];
  ReaderType::Pointer intensityReader = ReaderType::New();
  intensityReader->SetFileName(lriFileName);
  intensityReader->Update();
  FloatImageType::Pointer lriImage = intensityReader->GetOutput();
  FloatImageType::Pointer X_lr = NormalizeDataComponent(lriImage);
  lriImage=NULL;

  FloatImageType::Pointer SRImage = OpWeightedL2(X_lr, highResEdgeImage);

  const std::string hriFileName = argv[3];
  WriterType::Pointer hriWriter = WriterType::New();
  hriWriter->SetInput(SRImage);
  hriWriter->SetFileName(hriFileName);
  hriWriter->Update();

  tp.Stop();
  std::cout << tp.GetTotal() << tp.GetUnit() << std::endl;
  return EXIT_SUCCESS;
}
