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
/*
 * Author: Hans J. Johnson
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * Neither the name of the Nathan Kline Institute nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <fstream>
#include <cstring>
#include <cmath>
#include <map>

#include "BRAINSThreadControl.h"
#include <itkIntensityWindowingImageFilter.h>
#include "itkReflectiveCorrelationCenterToImageMetric.h"

#include "BRAINSConstellationModelerCLP.h"
#include "landmarksConstellationCommon.h"
#include "landmarksConstellationTrainingDefinitionIO.h"
#include "landmarkIO.h"
#include "GenericTransformImage.h"

#include "itkOrthogonalize3DRotationMatrix.h"

// //////////////////////////////////////////////////////////////
// Computes the unbiased sample variance of a set of n observations
// {x_1,x_2,...,x_n} from a given distribution.
// //////////////////////////////////////////////////////////////
template <typename DType>
void sample_variance(const std::vector<DType> & x, DType *mean, DType *var)
{
  const DType n = x.size();

  if( n == 1 )  // not possible to compute variance of a single sample
    {
    *mean = x[0];
    *var = 0.0;
    }

  DType       sum_of_sq = 0.0;
  DType       sum = 0.0;
  const typename std::vector<DType>::const_iterator theEnd = x.end();
  for( typename std::vector<DType>::const_iterator it = x.begin(); it != theEnd; ++it )
    {
    const DType & value = *it;
    sum_of_sq += value * value;
    sum += value;
    }
  *mean = sum / n;
  *var = ( sum_of_sq - sum * sum / n ) / ( n - 1.0 ); // var = E[X^2]-(E[X])^2
}

//
//
// ===========================================================================
int main(int argc, char *argv[])
{
  std::cout.precision(20);

  // ///////////////////////////////////////////////////////////////////
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  if( ( inputTrainingList.compare("") == 0 )
      || ( outputModel.compare("") == 0 ) )
    {
    std::cerr << "To run the program please specify the training filename and the "
              << "output outputModel filename." << std::endl;
    std::cerr << "Type " << argv[0] << " -h for more help." << std::endl;
    return EXIT_FAILURE;
    }

  LMC::globalverboseFlag = verbose;

  globalResultsDir = resultsDir;
  globalImagedebugLevel = writedebuggingImagesLevel;
  // /////////////////////////////////////////////////////////////////////////////////////////////
  short BackgroundFillValue;
  if( backgroundFillValueString == std::string("BIGNEG") )
    {
    BackgroundFillValue = -32768;
    }
  else
    {
    BackgroundFillValue = atoi( backgroundFillValueString.c_str() );
    }

  // /////////////////////////////////////////////////////////////////////////////////////////////
  // read information from the setup file, allocate some memories, and
  // initialize some variables

  landmarksConstellationTrainingDefinitionIO mDef(inputTrainingList);
  // //////////////////////////////////////////////////////////////////////////
  landmarksConstellationModelIO myModel;
  myModel.CopyFromModelDefinition(mDef);
  myModel.InitializeModel(true);

  // the following tables store the inputed manually selected RP, AC, PC, etc
  // locations for each volume
  // NOTE:  3 because before projection into Mid-sagital-plane (msp) need 3
  // coords
  // these tables store the RP, AC, PC, etc locations for each volume projected
  // on the MSP
  // NOTE:  only 2 coords need after projectsion into MSP, but using 3 to keep
  // math simple.
  // TODO:  Will need to change all index to get y and z (which is currently
  // coded as coord index 0 and 1
  std::vector<SImageType::PointType> rp_InMSPAlignedSpace( myModel.GetNumDataSets() );
  std::vector<SImageType::PointType> ac_InMSPAlignedSpace( myModel.GetNumDataSets() );
  std::vector<SImageType::PointType> pc_InMSPAlignedSpace( myModel.GetNumDataSets() );
  std::vector<SImageType::PointType> vn4_InMSPAlignedSpace( myModel.GetNumDataSets() );
  std::vector<SImageType::PointType> cec_InMSPAlignedSpace( myModel.GetNumDataSets() );
  std::vector<SImageType::PointType> cm_InMSPAlignedSpace( myModel.GetNumDataSets() );

  std::ofstream MSPOptFile;
  std::string   newOptimizedLandmarksTrainingFile;
  if( saveOptimizedLandmarks )
    {
    newOptimizedLandmarksTrainingFile = inputTrainingList + optimizedLandmarksFilenameExtender;

    // Write out the optimized points to a file so that training can be more
    // accurate.
    MSPOptFile.open( newOptimizedLandmarksTrainingFile.c_str() ); // open setup
    // file for
    // writing
    if( !MSPOptFile.is_open() )
      {
      std::cerr << "Can't write " << newOptimizedLandmarksTrainingFile << std::endl;
      std::cerr.flush();
      }

    // We are copying the first 6 lines rather than quoting them in this program
    // to give the user some control over what they are.
    std::ifstream prevMSPOptFile;
    prevMSPOptFile.open( inputTrainingList.c_str() ); // open previous setup
                                                      // file for
    // reading
    char linebuf[300];
    for( int lineNum = 0; lineNum < 6; ++lineNum )
      {
      prevMSPOptFile.getline(linebuf, 300);  // read to end of line;  line size
      // ends up in gcount()
      MSPOptFile.write( linebuf, prevMSPOptFile.gcount() );
      }
    prevMSPOptFile.close();
    }
  // //////////////////////////////////////////////////////////////////////////
  const unsigned int & mDefNumDataSets = mDef.GetNumDataSets();
  for( unsigned int currentDataset = 0; currentDataset < mDefNumDataSets; ++currentDataset )
    {
    std::cout << "====================================================================================" << std::endl;
    std::cout << "PROCESSING:" <<  mDef[currentDataset].GetImageFilename() << std::endl;
    // Since these are oriented images, the reorientation should not be
    // necessary.
    // //////////////////////////////////////////////////////////////////////////
    SImageType::Pointer volOrig = itkUtil::ReadImage<SImageType>( mDef[currentDataset].GetImageFilename() );
    if( volOrig.IsNull() )
      {
      printf( "\nCould not open image %s, aborting ...\n\n", mDef[currentDataset].GetImageFilename().c_str() );
      return EXIT_FAILURE;
      }
    SImageType::Pointer image;

    if( rescaleIntensities == true )
      {
      itk::StatisticsImageFilter<SImageType>::Pointer stats =
        itk::StatisticsImageFilter<SImageType>::New();
      stats->SetInput(volOrig);
      stats->Update();
      SImageType::PixelType minPixel( stats->GetMinimum() );
      SImageType::PixelType maxPixel( stats->GetMaximum() );

      if( trimRescaledIntensities > 0.0 )
        {
        // REFACTOR: a histogram would be traditional here, but seems
        // over-the-top;
        // I did this because it seemed to me if I knew mean, sigma, max and
        // min,
        // then I know Something about extreme outliers.

        double meanOrig( stats->GetMean() );
        double sigmaOrig( stats->GetSigma() );

        // REFACTOR:  In percentiles, 0.0005 two-tailed has worked in the past.
        // It only makes sense to trim the upper bound since the lower bound
        // would most likely
        // represent a large region of air around the head.  But this is not so
        // when using a mask.
        // For one-tailed, an error of 0.001 corresponds to 3.29052 standard
        // deviations of normal.
        // For one-tailed, an error of 0.0001 corresponds to 3.8906 standard
        // deviations of normal.
        // For one-tailed, an error of 0.00001 corresponds to 4.4172 standard
        // deviations of normal.
        // Naturally, the constant should default at the command line, ...

        double variationBound( ( maxPixel - meanOrig ) / sigmaOrig );
        double trimBound(variationBound - trimRescaledIntensities);
        if( trimBound > 0.0 )
          {
          maxPixel = static_cast<SImageType::PixelType>( maxPixel - trimBound * sigmaOrig );
          }
        }

      itk::IntensityWindowingImageFilter<SImageType, SImageType>::Pointer remapIntensityFilter =
        itk::IntensityWindowingImageFilter<SImageType, SImageType>::New();
      remapIntensityFilter->SetInput(volOrig);
      remapIntensityFilter->SetOutputMaximum(rescaleIntensitiesOutputRange[1]);
      remapIntensityFilter->SetOutputMinimum(rescaleIntensitiesOutputRange[0]);
      remapIntensityFilter->SetWindowMinimum(minPixel);
      remapIntensityFilter->SetWindowMaximum(maxPixel);
      remapIntensityFilter->Update();

      image = remapIntensityFilter->GetOutput();
      }
    else
      {
      image = volOrig;
      }

    SImageType::PointType origin;
    origin.Fill(0);

    // This section assumes that the landmarks are defined as
    // ITK compliant physical space
    // That are consistent with ITK images that they were collected from Dicom
    // coordinate system
    // No conversion is necessary
    SImageType::PointType origRP = mDef[currentDataset].GetNamedPoint("RP");
    SImageType::PointType origAC = mDef[currentDataset].GetNamedPoint("AC");
    SImageType::PointType origPC = mDef[currentDataset].GetNamedPoint("PC");
    SImageType::PointType origVN4 = mDef[currentDataset].GetNamedPoint("VN4");
    SImageType::PointType origLE = mDef[currentDataset].GetNamedPoint("LE");
    SImageType::PointType origRE = mDef[currentDataset].GetNamedPoint("RE");
    SImageType::PointType origCEC;
    origCEC.SetToMidPoint(origLE, origRE);

    SImageType::PointType centerOfHeadMass = GetCenterOfHeadMass(image);

    // original input volume from the training set
    // transforms image to MSP aligned voxel lattice
    RigidTransformType::Pointer Tmsp = RigidTransformType::New();
    SImageType::Pointer         volumeMSP;
    double                      c_c = 0;

    ComputeMSP(image, Tmsp, volumeMSP, centerOfHeadMass, mspQualityLevel, c_c);

    if( globalImagedebugLevel > 2 )
      {
      const std::string MSP_ImagePlane( globalResultsDir + "/MSP_PLANE_"
                                        + itksys::SystemTools::GetFilenameName
                                          ( mDef[currentDataset].GetImageFilename() ) );
      CreatedebugPlaneImage(volumeMSP, MSP_ImagePlane);
      }

    // Compute the transform from original space to the AC-PC aligned space using Reflective Correlation method
    //    RigidTransformType::Pointer invTmsp = RigidTransformType::New();
    //    Tmsp->GetInverse(invTmsp);

    // Instead of RC method, now we compute all the ac-pc aligned transforms by estimating the plane passing through RP,
    // AC and PC points
    RigidTransformType::Pointer ACPC_AlignedTransform = computeTmspFromPoints(origRP, origAC, origPC, origin);

    // We cannot easily compute the Inverse transform by the following to lines, we need to use versor for percise
    // transformation
    //    RigidTransformType::Pointer ACPC_AlignedTransform_INV = RigidTransformType::New();
    //    ACPC_AlignedTransform_INV->GetInverse(ACPC_AlignedTransform);

    // AC-PC aligned TRANSFORM
    VersorTransformType::Pointer finalTransform = VersorTransformType::New();
    finalTransform->SetFixedParameters( ACPC_AlignedTransform->GetFixedParameters() );
    itk::Versor<double>               versorRotation; // was commented before
    const itk::Matrix<double, 3, 3> & NewCleanedOrthogonalized = itk::Orthogonalize3DRotationMatrix(
        ACPC_AlignedTransform->GetMatrix() );
    versorRotation.Set( NewCleanedOrthogonalized );
    finalTransform->SetRotation( versorRotation );
    finalTransform->SetTranslation( ACPC_AlignedTransform->GetTranslation() );
    // inverse transform
    VersorTransformType::Pointer  ACPC_AlignedTransform_INV = VersorTransformType::New();
    const SImageType::PointType & centerPoint = finalTransform->GetCenter();
    ACPC_AlignedTransform_INV->SetCenter( centerPoint );
    ACPC_AlignedTransform_INV->SetIdentity();
    finalTransform->GetInverse( ACPC_AlignedTransform_INV );
    //////////////////////////////////////////////////////////////////

    // Transform points based on the new transform
    cm_InMSPAlignedSpace[currentDataset] = ACPC_AlignedTransform_INV->TransformPoint(centerOfHeadMass);
    rp_InMSPAlignedSpace[currentDataset] = ACPC_AlignedTransform_INV->TransformPoint(origRP);
    ac_InMSPAlignedSpace[currentDataset] = ACPC_AlignedTransform_INV->TransformPoint(origAC);
    pc_InMSPAlignedSpace[currentDataset] = ACPC_AlignedTransform_INV->TransformPoint(origPC);
    vn4_InMSPAlignedSpace[currentDataset] = ACPC_AlignedTransform_INV->TransformPoint(origVN4);
    cec_InMSPAlignedSpace[currentDataset] = ACPC_AlignedTransform_INV->TransformPoint(origCEC);

    if( globalImagedebugLevel > 3 )
      {
      const SImageType::PointType finalRP = ACPC_AlignedTransform_INV->TransformPoint(origRP);
      const SImageType::PointType finalAC = ACPC_AlignedTransform_INV->TransformPoint(origAC);
      const SImageType::PointType finalPC = ACPC_AlignedTransform_INV->TransformPoint(origPC);
      const SImageType::PointType finalVN4 = ACPC_AlignedTransform_INV->TransformPoint(origVN4);

      SImageType::Pointer volumeACPC_Aligned =
        TransformResample<SImageType, SImageType>( image.GetPointer(), image.GetPointer(), BackgroundFillValue,
                                                   GetInterpolatorFromString<SImageType>("Linear").GetPointer(),
                                                   ACPC_AlignedTransform.GetPointer() );

      itkUtil::WriteImage<SImageType>( volumeACPC_Aligned, resultsDir + "/ACPC_Aligned_"
                                       + itksys::SystemTools::GetFilenameName( mDef[currentDataset].GetImageFilename() ) );
      MakeLabelImage( volumeACPC_Aligned, finalRP, finalAC, finalPC, finalVN4, resultsDir + "/Mask_Resampled_"
                      + itksys::SystemTools::GetFilenameName( mDef[currentDataset].GetImageFilename() ) );
      }

    // Build template for each landmark
    for(
      std::map<std::string, SImageType::PointType>::iterator it = mDef[currentDataset].begin();
      it != mDef[currentDataset].end(); ++it )
      {
      std::cout << "Training template for " << it->first << std::endl;
      const SImageType::PointType origPoint = it->second;
      const SImageType::PointType transformedPoint = ACPC_AlignedTransform_INV->TransformPoint(origPoint);
      /* PRINT FOR TEST /////////////////////////////////////////////
          std::cout << "original point: " << it->second << std::endl;
          std::cout << "transformed point: " << transformedPoint << std::endl;
      /////////////////////////////////////////////////////////////*/
      for( unsigned int currentAngle = 0; currentAngle < myModel.GetNumRotationSteps(); currentAngle++ )
        {
        // //////  create a rotation about the center with respect to the
        // current test rotation angle
        const float degree_current_angle = myModel.GetInitialRotationAngle() + myModel.GetInitialRotationStep()
          * currentAngle;
        const float current_angle = degree_current_angle * vnl_math::pi / 180;

        RigidTransformType::Pointer Point_Rotate = RigidTransformType::New();
        Point_Rotate->SetCenter(transformedPoint);
        Point_Rotate->SetRotation(current_angle, 0, 0);

        SImageType::Pointer image_TestRotated =
          CreateTestCenteredRotatedImage2(ACPC_AlignedTransform, origPoint, image, Point_Rotate);
        if( globalImagedebugLevel > 5 )
          {
          std::stringstream s("");
          s << "image_" << it->first << "_TestRotated_" << current_angle << "_";
          const std::string rotatedName = resultsDir + "/" + s.str().c_str()
            + itksys::SystemTools::GetFilenameName( mDef[currentDataset].GetImageFilename() );
          std::cout << "Writing file: " << rotatedName << std::endl;
          itkUtil::WriteImage<SImageType>(image_TestRotated, rotatedName);
          }

        // TODO:  The following 3 function calls may be a performance problem,
        // and it should be  straight forward to refactor this into a single function
        // extractZeroMeanNormalizedVector that has the same signature as extractArray, but has many fewer
        // loop iterations.
        LinearInterpolatorType::Pointer imInterp = LinearInterpolatorType::New();
        imInterp->SetInputImage(image_TestRotated);
        //Resize the model internal vector.
        myModel.AccessTemplate(it->first, currentDataset, currentAngle).resize( myModel.m_VectorIndexLocations[it->first].size() );
        extractArray( imInterp, transformedPoint, myModel.m_VectorIndexLocations[it->first],
                      myModel.AccessTemplate(it->first, currentDataset, currentAngle) );
        removeVectorMean( myModel.AccessTemplate(it->first, currentDataset, currentAngle) );
        normalizeVector( myModel.AccessTemplate(it->first, currentDataset, currentAngle) );
        }
      }

    if( saveOptimizedLandmarks )
      {
      MSPOptFile.close();
      }
    }

  /* PRINT FOR TEST ////////////////////////////////////////////
  std::cout << "\nPROCESSING AC transformed values in MSP aligned space by 'Reflective Correlation' method:" << std::endl;
  for( unsigned int currentDataset = 0; currentDataset < mDefNumDataSets; currentDataset++ )
    {
        std::cout << "====================================================================================" << std::endl;
        std::cout << currentDataset+1 << "#: " << ac_InMSPAlignedSpace[currentDataset] << std::endl;
    }
  //////////////////////////////////////////////////////////////*/

  // -------------------------------
  std::cout << "\nCompute vector means:" << std::endl;

  // Compute MPJtoPCMean
    {
    std::vector<SImageType::PointType::VectorType::ComponentType> xComponents( myModel.GetNumDataSets() );
    std::vector<SImageType::PointType::VectorType::ComponentType> yComponents( myModel.GetNumDataSets() );
    std::vector<SImageType::PointType::VectorType::ComponentType> zComponents( myModel.GetNumDataSets() );
    for( unsigned int currentDataset = 0; currentDataset < myModel.GetNumDataSets(); currentDataset++ )
      {
      xComponents[currentDataset] = pc_InMSPAlignedSpace[currentDataset][0] - rp_InMSPAlignedSpace[currentDataset][0];
      yComponents[currentDataset] = pc_InMSPAlignedSpace[currentDataset][1] - rp_InMSPAlignedSpace[currentDataset][1];
      zComponents[currentDataset] = pc_InMSPAlignedSpace[currentDataset][2] - rp_InMSPAlignedSpace[currentDataset][2];
      }
    SImageType::PointType::VectorType RPtoPCmean;
    SImageType::PointType::VectorType RPtoPCvar;
    sample_variance( xComponents, &( RPtoPCmean[0] ), &( RPtoPCvar[0] ) );
    sample_variance( yComponents, &( RPtoPCmean[1] ), &( RPtoPCvar[1] ) );
    sample_variance( zComponents, &( RPtoPCmean[2] ), &( RPtoPCvar[2] ) );
    myModel.SetRPtoXMean("PC", RPtoPCmean);
    std::cout << "RPtoPC mean = " << RPtoPCmean << "mm" << std::endl;
    std::cout << "RPtoPC var = " << RPtoPCvar << "mm^2" << std::endl;
    }

  // Compute MPJtoVN4Mean
    {
    std::vector<SImageType::PointType::VectorType::ComponentType> xComponents( myModel.GetNumDataSets() );
    std::vector<SImageType::PointType::VectorType::ComponentType> yComponents( myModel.GetNumDataSets() );
    std::vector<SImageType::PointType::VectorType::ComponentType> zComponents( myModel.GetNumDataSets() );
    for( unsigned int currentDataset = 0; currentDataset < myModel.GetNumDataSets(); currentDataset++ )
      {
      xComponents[currentDataset] = vn4_InMSPAlignedSpace[currentDataset][0] - rp_InMSPAlignedSpace[currentDataset][0];
      yComponents[currentDataset] = vn4_InMSPAlignedSpace[currentDataset][1] - rp_InMSPAlignedSpace[currentDataset][1];
      zComponents[currentDataset] = vn4_InMSPAlignedSpace[currentDataset][2] - rp_InMSPAlignedSpace[currentDataset][2];
      }
    SImageType::PointType::VectorType RPtoVN4mean;
    SImageType::PointType::VectorType RPtoVN4var;
    sample_variance( xComponents, &( RPtoVN4mean[0] ), &( RPtoVN4var[0] ) );
    sample_variance( yComponents, &( RPtoVN4mean[1] ), &( RPtoVN4var[1] ) );
    sample_variance( zComponents, &( RPtoVN4mean[2] ), &( RPtoVN4var[2] ) );
    myModel.SetRPtoXMean("VN4", RPtoVN4mean);
    std::cout << "RPtoVN4 mean = " << RPtoVN4mean << "mm" << std::endl;
    std::cout << "RPtoVN4 var = " << RPtoVN4var << "mm^2" << std::endl;
    }

  // Compute MPJtoCECMean
    {
    std::vector<SImageType::PointType::VectorType::ComponentType> xComponents( myModel.GetNumDataSets() );
    std::vector<SImageType::PointType::VectorType::ComponentType> yComponents( myModel.GetNumDataSets() );
    std::vector<SImageType::PointType::VectorType::ComponentType> zComponents( myModel.GetNumDataSets() );
    for( unsigned int currentDataset = 0; currentDataset < myModel.GetNumDataSets(); currentDataset++ )
      {
      xComponents[currentDataset] = cec_InMSPAlignedSpace[currentDataset][0] - rp_InMSPAlignedSpace[currentDataset][0];
      yComponents[currentDataset] = cec_InMSPAlignedSpace[currentDataset][1] - rp_InMSPAlignedSpace[currentDataset][1];
      zComponents[currentDataset] = cec_InMSPAlignedSpace[currentDataset][2] - rp_InMSPAlignedSpace[currentDataset][2];
      }
    SImageType::PointType::VectorType RPtoCECmean;
    SImageType::PointType::VectorType RPtoCECvar;
    sample_variance( xComponents, &( RPtoCECmean[0] ), &( RPtoCECvar[0] ) );
    sample_variance( yComponents, &( RPtoCECmean[1] ), &( RPtoCECvar[1] ) );
    sample_variance( zComponents, &( RPtoCECmean[2] ), &( RPtoCECvar[2] ) );
    myModel.SetRPtoCECMean(RPtoCECmean);
    std::cout << "RPtoCEC mean = " << RPtoCECmean << "mm" << std::endl;
    std::cout << "RPtoCEC var = " << RPtoCECvar << "mm^2" << std::endl;
    }

  // Compute MPJtoACMean
    {
    std::vector<SImageType::PointType::VectorType::ComponentType> xComponents( myModel.GetNumDataSets() );
    std::vector<SImageType::PointType::VectorType::ComponentType> yComponents( myModel.GetNumDataSets() );
    std::vector<SImageType::PointType::VectorType::ComponentType> zComponents( myModel.GetNumDataSets() );
    for( unsigned int currentDataset = 0; currentDataset < myModel.GetNumDataSets(); currentDataset++ )
      {
      xComponents[currentDataset] = ac_InMSPAlignedSpace[currentDataset][0] - rp_InMSPAlignedSpace[currentDataset][0];
      yComponents[currentDataset] = ac_InMSPAlignedSpace[currentDataset][1] - rp_InMSPAlignedSpace[currentDataset][1];
      zComponents[currentDataset] = ac_InMSPAlignedSpace[currentDataset][2] - rp_InMSPAlignedSpace[currentDataset][2];
      }
    SImageType::PointType::VectorType RPtoACmean;
    SImageType::PointType::VectorType RPtoACvar;
    sample_variance( xComponents, &( RPtoACmean[0] ), &( RPtoACvar[0] ) );
    sample_variance( yComponents, &( RPtoACmean[1] ), &( RPtoACvar[1] ) );
    sample_variance( zComponents, &( RPtoACmean[2] ), &( RPtoACvar[2] ) );
    std::cout << "RPtoAC mean = " << RPtoACmean << "mm" << std::endl;
    std::cout << "RPtoAC var = " << RPtoACvar << "mm^2" << std::endl;
    myModel.SetRPtoXMean("AC", RPtoACmean);
    }

  // ??  What does RPPC_to_RPAC_angle and RPAC_over_RPPC mean?
  std::vector<float> RPPC_to_RPAC_angle( myModel.GetNumDataSets() );
  std::vector<float> RPAC_over_RPPC( myModel.GetNumDataSets() );
  // Initializing CMtoRPMean
  SImageType::PointType::VectorType CMtoRPMean;
  CMtoRPMean.Fill(0.0);
  // This for loop does two jobs
  for( unsigned int currentDataset = 0; currentDataset < myModel.GetNumDataSets(); currentDataset++ )
    {
    // JOB1: RPPC_to_RPAC_angle and RPAC_over_RPPC
    double curr_RPPC_to_RPAC_angle;
    double curr_RPAC_over_RPPC;
    decomposeRPAC(rp_InMSPAlignedSpace[currentDataset],
                  pc_InMSPAlignedSpace[currentDataset],
                  ac_InMSPAlignedSpace[currentDataset],
                  &curr_RPPC_to_RPAC_angle,
                  &curr_RPAC_over_RPPC);
    RPPC_to_RPAC_angle[currentDataset] = curr_RPPC_to_RPAC_angle;
    RPAC_over_RPPC[currentDataset] = curr_RPAC_over_RPPC;

    // JOB2: CMtoRPMean
    // // NOTE:  This needs to be the average distance from the center of gravity.
    SImageType::PointType::VectorType RPDistanceFromCenterOfMass =
      rp_InMSPAlignedSpace[currentDataset] - cm_InMSPAlignedSpace[currentDataset];
    CMtoRPMean += RPDistanceFromCenterOfMass;
    }
  CMtoRPMean /= static_cast<float>( myModel.GetNumDataSets() );

  // Now compute the average angle and average ratio
  float RPPC_to_RPAC_angleMean, RPAC_over_RPPCMean;
    {
    float dummy_mean(0.0);
    float dummy_var(0.0);
    // HACK:  Need to compute only the mean, since that is all that is needed.
    sample_variance(RPPC_to_RPAC_angle, &dummy_mean, &dummy_var);
    RPPC_to_RPAC_angleMean = dummy_mean;
    sample_variance(RPAC_over_RPPC, &dummy_mean, &dummy_var);
    RPAC_over_RPPCMean = dummy_mean;
    }

  myModel.SetRPPC_to_RPAC_angleMean(RPPC_to_RPAC_angleMean);
  myModel.SetRPAC_over_RPPCMean(RPAC_over_RPPCMean);

  myModel.SetCMtoRPMean(CMtoRPMean);

  myModel.WriteModelFile(outputModel);
  myModel.PrintHeaderInfo();
  return 0;
}
