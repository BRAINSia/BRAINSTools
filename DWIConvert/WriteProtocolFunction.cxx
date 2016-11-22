//
// Created by Johnson, Hans J on 11/22/16.
//


if( writeProtocolGradientsFile == true )
{
std::cerr << "ERROR: DEPRECATED IMPLEMENTED" << std::endl;
return EXIT_FAILURE;
#if 0
itk::NumberToString<double> DoubleConvert;
    //////////////////////////////////////////////
    // writeProtocolGradientsFile write protocolGradientsFile file
    // This part follows a DWI NRRD file in NRRD format 5.
    // There should be a better way using itkNRRDImageIO.

    std::ofstream protocolGradientsFile;
    // std::string protocolGradientsFileFileName = outputDir + "/" + outputFileName;

    const std::string protocolGradientsFileName = outputVolumeHeaderName + ".txt";
    DWIConverter::RotationMatrixType LPSDirCos = converter->GetLPSDirCos();
    DWIConverter::RotationMatrixType MeasurementFrame = converter->GetMeasurementFrame();

    protocolGradientsFile.open( protocolGradientsFileName.c_str() );
    protocolGradientsFile << "ImageOrientationPatient (0020|0032): "
                          << DoubleConvert(LPSDirCos[0][0]) << "\\"
                          << DoubleConvert(LPSDirCos[1][0])
                          << "\\" << DoubleConvert(LPSDirCos[2][0]) << "\\"
                          << DoubleConvert(LPSDirCos[0][1])
                          << "\\" << DoubleConvert(LPSDirCos[1][1]) << "\\"
                          << DoubleConvert(LPSDirCos[2][1]) << "\\"
                          << std::endl;
    protocolGradientsFile << "==================================" << std::endl;
    protocolGradientsFile << "Direction Cosines: " << std::endl << LPSDirCos << std::endl;
    protocolGradientsFile << "==================================" << std::endl;
    protocolGradientsFile << "MeasurementFrame: " << std::endl << MeasurementFrame << std::endl;
    protocolGradientsFile << "==================================" << std::endl;
    for( unsigned int k = 0; k < converter->GetNVolume(); ++k )
      {
      protocolGradientsFile << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << k << "=["
                            << DoubleConvert(BvalueScaledDiffusionVectors[k][0]) << ";"
                            << DoubleConvert(BvalueScaledDiffusionVectors[k][1]) << ";"
                            << DoubleConvert(BvalueScaledDiffusionVectors[k][2]) << "]" << std::endl;
      }
    protocolGradientsFile << "==================================" << std::endl;
    for( unsigned int k = 0; k < converter->GetNVolume(); ++k )
      {
      const vnl_vector_fixed<double, 3u> ProtocolGradient = InverseMeasurementFrame * BvalueScaledDiffusionVectors[k];
      protocolGradientsFile << "Protocol_gradient_" << std::setw(4) << std::setfill('0') << k << "=["
                            << DoubleConvert(ProtocolGradient[0]) << ";"
                            << DoubleConvert(ProtocolGradient[1]) << ";"
                            << DoubleConvert(ProtocolGradient[2]) << "]" << std::endl;
      }
    protocolGradientsFile << "==================================" << std::endl;
    protocolGradientsFile.close();
#endif
}