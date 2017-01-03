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
/*=========================================================================
Computing (NAMIC), funded by the National Institutes of Health
through the NIH Roadmap for Medical Research, Grant U54 EB005149.

See License.txt or http://www.slicer.org/copyright/copyright.txt for details.

 ***
 This program converts Diffusion weighted MR images in Dicom format into
 NRRD format.

Assumptions:

1) Uses left-posterior-superior (Dicom default) as default space for philips and siemens.
This is the default space for NRRD header.
2) For GE data, Dicom data are arranged in volume interleaving order.
3) For Siemens data, images are arranged in mosaic form.
4) For oblique collected Philips data, the measurement frame for the
gradient directions is the same as the ImageOrientationPatient

Reference materials:
DICOM Data Dictionary: http://medical.nema.org/Dicom/2011/11_06pu.pdf
=========================================================================*/


#include "DWIConverter.h"

/** the real computation goes on in DWIConverter classes, of which
 * there is one for each manufacturer we encounter.
 */
#include "DWIConverterFactory.h"
#include "FSLDWIConverter.h"
#include "NRRDDWIConverter.h"

#include <BRAINSCommonLib.h>




#undef HAVE_SSTREAM
#include "DWIConvertCLP.h"
#include "DWIConvertLib.h"

int main(int argc, char *argv[])
{
    PARSE_ARGS;
    //const std::string version = commandLine.getVersion();
    //BRAINSRegisterAlternateIO();

    std::cout << "======= DWI Convert Public Lib Ctest =========" << std::endl;
    DWIConvert dWIConvert;

    dWIConvert.setInputFileType(inputVolume, inputDicomDirectory);
    dWIConvert.setInputBValues (inputBValues);
    dWIConvert.setInputBVectors (inputBVectors);
    dWIConvert.setGradientVectorFile (gradientVectorFile);
    dWIConvert.setSmallGradientThreshold (smallGradientThreshold);

    dWIConvert.setfMRIOutput (fMRIOutput);
    dWIConvert.setTranspose (transpose);
    dWIConvert.setAllowLossyConversion (allowLossyConversion);
    dWIConvert.setUseIdentityMeasurementFrame (useIdentityMeaseurementFrame);
    dWIConvert.setUseBMatrixGradientDirections (useBMatrixGradientDirections);

    dWIConvert.setOutputFileType(outputVolume);
    dWIConvert.setOutputDirectory(outputDirectory);
    dWIConvert.setOutputBValues(outputBValues);
    dWIConvert.setOutputBVectors(outputBVectors);

    int result = dWIConvert.read();
    if (EXIT_SUCCESS == result) {
        return dWIConvert.write(outputVolume);
    }
    else return result;
}
