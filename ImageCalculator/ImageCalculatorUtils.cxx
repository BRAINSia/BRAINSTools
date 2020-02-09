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
#include "ImageCalculatorUtils.h"
#include "Imgmath.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vcl_compiler.h>
#include <iostream>
#include <cmath>
#include <iostream>
#include <metaCommand.h>
#include <iostream>

void
PrintDataTypeStrings()
{
  // Prints the Input and output data type strings.
  std::cout << "UCHAR" << std::endl;
  std::cout << "SHORT" << std::endl;
  std::cout << "USHORT" << std::endl;
  std::cout << "INT" << std::endl;
  std::cout << "UINT" << std::endl;
  std::cout << "FLOAT" << std::endl;
  std::cout << "DOUBLE" << std::endl;
}

void
ReplaceSubWithSub(std::string & s, const char * o, const char * n)
{
  if (!s.empty())
  {
    std::string            from(o), to(n);
    std::string::size_type start = 0;
    while ((start = s.find(from, start)) != std::string::npos)
    {
      s.replace(start, from.size(), to);
      start += to.size();
    }
  }
}

int
CompareNoCase(const std::string & s, const std::string & s2)
{
  // Compare strings.
  std::string::const_iterator p = s.begin();
  std::string::const_iterator p2 = s2.begin();

  while (p != s.end() && p2 != s2.end())
  {
    if (toupper(*p) != toupper(*p2))
    {
      return (toupper(*p) < toupper(*p2)) ? -1 : 1;
    }
    ++p;
    ++p2;
  }

  return (s2.size() == s.size()) ? 0 : (s.size() < s2.size()) ? -1 : 1;
}

// Call ImageCalculator process for 2d images.
extern void
ImageCalculatorProcess2D(const std::string & InType, MetaCommand & command);

// Call ImageCalculator process for 3d images.
extern void
ImageCalculatorProcess3D(const std::string & InType, MetaCommand & command);

bool
ValidPixelType(const std::string & PixelType)
{
  const char * s = PixelType.c_str();

  // check to see if valid type
  if ((CompareNoCase(s, std::string("UCHAR"))) && (CompareNoCase(s, std::string("SHORT"))) &&
      (CompareNoCase(s, std::string("USHORT"))) && (CompareNoCase(s, std::string("INT"))) &&
      (CompareNoCase(s, std::string("UINT"))) && (CompareNoCase(s, std::string("FLOAT"))) &&
      (CompareNoCase(s, std::string("DOUBLE"))))
  {
    return false;
  }
  return true;
}

// Ensure that the template code is only compiled once for both the real program and the test programs.
int
PrimaryImageCalculatorRoutine(int argc, char * argv[])
{
  MetaCommand command;

  /*Input Image filenames. Any  number of input images may be given. The input filenames must be preceded by the number
    of inputs given.*/
  command.SetOption("in", "", true, "InputFile names");
  command.SetOptionLongTag("in", "in");
  command.AddOptionField("in", "in", MetaCommand::STRING, true);
  command.SetOptionComplete("in", true);

  // The images will be read using the Input pixel type.All the operations are performed in this pixel type.
  command.SetOption("InputPixelType", "", false, "InputPixel Type");
  command.SetOptionLongTag("InputPixelType", "intype");
  command.AddOptionField("InputPixelType", "PixelType", MetaCommand::STRING, false, "FLOAT");

  // The dimensions of the input images. All the images should be of the same dimension.
  command.SetOption("InputDimensions", "d", false, "Input Dimension 2 or 3");
  command.AddOptionField("InputDimensions", "dims", MetaCommand::INT, false, "3");

  // Output filename.
  command.SetOption("OutputFilename", "", false, "OutputFile name");
  command.SetOptionLongTag("OutputFilename", "out");
  command.AddOptionField("OutputFilename", "filename", MetaCommand::STRING, false);

  // The images will be written in this type. The default is input pixel type.
  command.SetOption("OutputPixelType", "", false, "OutputPixel Type");
  command.SetOptionLongTag("OutputPixelType", "outtype");
  command.AddOptionField("OutputPixelType", "PixelType", MetaCommand::STRING, false, "FLOAT");

  // Add the images.
  command.SetOption("Add", "", false, "Add Images");
  command.SetOptionLongTag("Add", "add");
  command.AddOptionField("Add", "add", MetaCommand::FLAG, false);

  // Subtract the images.
  command.SetOption("Sub", "", false, "Subtract Images");
  command.SetOptionLongTag("Sub", "sub");
  command.AddOptionField("Sub", "sub", MetaCommand::FLAG, false);

  // Divide the images.
  command.SetOption("Div", "", false, "Divide Images");
  command.SetOptionLongTag("Div", "div");
  command.AddOptionField("Div", "div", MetaCommand::FLAG, false);

  // Multiply the images.
  command.SetOption("Mul", "", false, "Multiply Images");
  command.SetOptionLongTag("Mul", "mul");
  command.AddOptionField("Mul", "mul", MetaCommand::FLAG, false);

  // Get the variance image
  command.SetOption("Var", "", false, "Variance of Images");
  command.SetOptionLongTag("Var", "var");
  command.AddOptionField("Var", "var", MetaCommand::FLAG, false);

  // Get the Average image
  command.SetOption("Avg", "", false, "Average Images");
  command.SetOptionLongTag("Avg", "avg");
  command.AddOptionField("Avg", "avg", MetaCommand::FLAG, false);

  // Multiply the output with a constant scalar value.
  command.SetOption("OMulC", "", false, "Multiply Output Image with constant value");
  command.SetOptionLongTag("OMulC", "ofmulc");
  command.AddOptionField("OMulC", "constant", MetaCommand::FLOAT, false);

  // Multiply the inputs with a constant scalar value.
  command.SetOption("IMulC", "", false, "Multiply Accumulator Image with constant value");
  command.SetOptionLongTag("IMulC", "ifmulc");
  command.AddOptionField("IMulC", "constant", MetaCommand::FLOAT, false);

  // Divide the output with a constant scalar value.
  command.SetOption("ODivC", "", false, "Divide Output Image with constant value");
  command.SetOptionLongTag("ODivC", "ofdivc");
  command.AddOptionField("ODivC", "constant", MetaCommand::FLOAT, false);

  // Divide the inputs with a constant scalar value.
  command.SetOption("IDivC", "", false, "Divide Accumulator Image with constant value");
  command.SetOptionLongTag("IDivC", "ifdivc");
  command.AddOptionField("IDivC", "constant", MetaCommand::FLOAT, false);

  // Add a constant scalar value to the output.
  command.SetOption("OAddC", "", false, "Add Output Image with constant value");
  command.SetOptionLongTag("OAddC", "ofaddc");
  command.AddOptionField("OAddC", "constant", MetaCommand::FLOAT, false);

  // Add a constant scalar value to the inputs.
  command.SetOption("IAddC", "", false, "Add Accumulator Image with constant value");
  command.SetOptionLongTag("IAddC", "ifaddc");
  command.AddOptionField("IAddC", "constant", MetaCommand::FLOAT, false);

  // Subtract a constant scalar value from the output.
  command.SetOption("OSubC", "", false, "Subtract Output Image with constant value");
  command.SetOptionLongTag("OSubC", "ofsubc");
  command.AddOptionField("OSubC", "constant", MetaCommand::FLOAT, false);

  // Subtract a constant scalar value from the inputs.
  command.SetOption("ISubC", "", false, "Subtract Accumulator Image with constant value");
  command.SetOptionLongTag("ISubC", "ifsubc");
  command.AddOptionField("ISubC", "constant", MetaCommand::FLOAT, false);

  // Square the output image pixel values.
  command.SetOption("OSqr", "", false, "Square Accumulator Image before writing");
  command.SetOptionLongTag("OSqr", "ofsqr");
  command.AddOptionField("OSqr", "ofsqr", MetaCommand::FLAG, false);

  // Input Binary Image.
  command.SetOption("Ifbin", "", false, "Input Binary Image.");
  command.SetOptionLongTag("Ifbin", "ifbin");
  command.AddOptionField("Ifbin", "ifbin", MetaCommand::FLAG, false);

  // Output Binary Image.
  command.SetOption("Ofbin", "", false, "Output Binary Image.");
  command.SetOptionLongTag("Ofbin", "ofbin");
  command.AddOptionField("Ofbin", "ofbin", MetaCommand::FLAG, false);

  // Gaussian Filter with value of sigma to the output.
  command.SetOption("OGaussianSigma", "", false, "Gaussian smooth output image with sigma value");
  command.SetOptionLongTag("OGaussianSigma", "ofgaussiansigma");
  command.AddOptionField("OGaussianSigma", "constant", MetaCommand::FLOAT, false, "");

  // Gaussian Filter with value of sigma to the inputs.
  command.SetOption("IGaussianSigma", "", false, "Gaussian smooth input image with sigma value");
  command.SetOptionLongTag("IGaussianSigma", "ifgaussiansigma");
  command.AddOptionField("IGaussianSigma", "constant", MetaCommand::FLOAT, false, "");

  // If the accumulator buffer is not empty, then every subsequent image is histogram equalized to the current
  // accumulator buffer.
  command.SetOption(
    "IHisteq", "", false, "Ifhisteq equalizes each image to the current accumulator buffer after reading.");
  command.SetOptionLongTag("IHisteq", "ifhisteq");
  command.AddOptionField("IHisteq", "constant", MetaCommand::INT, false, "");

  // Square the input image pixel values.
  command.SetOption("ISqr", "", false, "Square Accumulator Image After reading");
  command.SetOptionLongTag("ISqr", "ifsqr");
  command.AddOptionField("ISqr", "ifsqr", MetaCommand::FLAG, false);

  // Get the square root of the output image pixel values.
  command.SetOption("OSqrt", "", false, "Square Root  Accumulator Image before writing");
  command.SetOptionLongTag("OSqrt", "ofsqrt");
  command.AddOptionField("OSqrt", "ofsqrt", MetaCommand::FLAG, false);

  // Get the square root of the input image pixel values.
  command.SetOption("ISqrt", "", false, "Sqrt Accumulator Image After reading");
  command.SetOptionLongTag("ISqrt", "ifsqrt");
  command.AddOptionField("ISqrt", "ifsqrt", MetaCommand::FLAG, false);

  // Get the average pixel value of the output image.
  command.SetOption("StatAvg", "", false, "Average Output Image Value");
  command.SetOptionLongTag("StatAvg", "statAVG");
  command.AddOptionField("StatAvg", "statAVG", MetaCommand::FLAG, false);

  // Get the variance of the pixel value of the output image.
  command.SetOption("StatVAR", "", false, "Variance of output Image");
  command.SetOptionLongTag("StatVAR", "statVAR");
  command.AddOptionField("StatVAR", "statVAR", MetaCommand::FLAG, false);

  // Get the sum of the pixel value of the output image.
  command.SetOption("StatSUM", "", false, "Sum of output Image Values");
  command.SetOptionLongTag("StatSUM", "statSUM");
  command.AddOptionField("StatSUM", "statSUM", MetaCommand::FLAG, false);

  // Mask the output image with another image. The statistics are given for the masked portion.
  command.SetOption("Statmask", "", false, "Image to mask against.");
  command.SetOptionLongTag("Statmask", "statmask");
  command.AddOptionField("Statmask", "File Name", MetaCommand::STRING, false);

  // If a mask is given then a pixel value should be entered and statitsics will be calculated for input image under
  // this value in the mask.
  command.SetOption(
    "Statmaskvalue", "", false, "Statistics in the image will be calculated for the pixels masked by this value.");
  command.SetOptionLongTag("Statmaskvalue", "statmaskvalue");
  command.AddOptionField("Statmaskvalue", "constant", MetaCommand::INT, false, "");

  // Get the number of pixels in the output image.
  command.SetOption("StatNPX", "", false, "Number of Pixels.");
  command.SetOptionLongTag("StatNPX", "statNPX");
  command.AddOptionField("StatNPX", "statNPX", MetaCommand::FLAG, false);

  // Get the maximum pixel value of the output image.
  command.SetOption("StatMAX", "", false, "Maximum of  Pixels.");
  command.SetOptionLongTag("StatMAX", "statMAX");
  command.AddOptionField("StatMAX", "statMAX", MetaCommand::FLAG, false);

  // Get the minimum pixel value of the output image.
  command.SetOption("StatMIN", "", false, "Minimum of Pixels.");
  command.SetOptionLongTag("StatMIN", "statMIN");
  command.AddOptionField("StatMIN", "statMIN", MetaCommand::FLAG, false);

  // Get the absolute maximum pixel value of the output image.
  command.SetOption("StatAMX", "", false, "Absolute Maximum of Pixels.");
  command.SetOptionLongTag("StatAMX", "statAMX");
  command.AddOptionField("StatAMX", "statAMX", MetaCommand::FLAG, false);

  // Get the absolute minimum pixel value of the output image.
  command.SetOption("StatAMN", "", false, "Absolute Minimum of Pixels.");
  command.SetOptionLongTag("StatAMN", "statAMN");
  command.AddOptionField("StatAMN", "statAMN", MetaCommand::FLAG, false);

  // Show the description of the stat values which can be calculated.
  command.SetOption("Statallcodes", "", false, "Prints the coding of statistical varibles.");
  command.SetOptionLongTag("Statallcodes", "statallcodes");
  command.AddOptionField("Statallcodes", "statallcodes", MetaCommand::FLAG, false);

  command.SetParseFailureOnUnrecognizedOption(true);
  if (!command.Parse(argc, argv))
  {
    return 1;
  }

  // Test if the input data type is valid
  const std::string PixelType(command.GetValueAsString("InputPixelType", "PixelType"));
  if (!PixelType.empty())
  {
    if (!ValidPixelType(PixelType))
    {
      std::cout << "Error. Invalid data type string specified with -intype!" << std::endl;
      std::cout << "Use one of the following:" << std::endl;
      PrintDataTypeStrings();
      throw;
    }
  }

  const std::string OutPixelType(command.GetValueAsString("OutputPixelType", "PixelType"));

  if (!OutPixelType.empty())
  {
    // check to see if valid type
    if (!ValidPixelType(OutPixelType))
    {
      std::cout << "Error. Invalid data type string specified with -outtype!" << std::endl;
      std::cout << "Use one of the following:" << std::endl;
      PrintDataTypeStrings();
      throw;
    }
  }

  // Test that only one operation is set
  int opcount = 0;
  if (command.GetValueAsBool("Add", "add"))
  {
    ++opcount;
  }
  if (command.GetValueAsBool("Sub", "sub"))
  {
    ++opcount;
  }
  if (command.GetValueAsBool("Mul", "mul"))
  {
    ++opcount;
  }
  if (command.GetValueAsBool("Div", "div"))
  {
    ++opcount;
  }
  if (command.GetValueAsBool("Var", "var"))
  {
    ++opcount;
  }
  if (command.GetValueAsBool("Avg", "avg"))
  {
    ++opcount;
  }
  if (opcount > 1)
  {
    itkGenericExceptionMacro(<< "Can only supply one operation to do [-add|-sub|-mul|-div|-var|-avg]");
  }

  // Call the ImageCalculatorReadWrite function based on the dimension.
  const std::string InType(command.GetValueAsString("InputPixelType", "PixelType"));
  const int         dims = command.GetValueAsInt("InputDimensions", "dims");

  switch (dims)
  {
    case 2:
    {
      ImageCalculatorProcess2D(InType, command);
    }
    break;
    case 3:
    {
      ImageCalculatorProcess3D(InType, command);
    }
    break;
    default:
      return 1;
  }

  return 0;
}
