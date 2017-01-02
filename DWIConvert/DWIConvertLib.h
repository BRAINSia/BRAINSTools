//
// Created by Hui Xie on 12/26/16.
//

#ifndef BRAINSTOOLS_DWICONVERTLIB_H
#define BRAINSTOOLS_DWICONVERTLIB_H

#include "DWIConverter.h"
#include "DWIConverterFactory.h"


class DWIConvert {
public:
    DWIConvert();
    DWIConvert(const std::string& outputVolume, const std::string& inputVolume, const std::string& inputDicomDirectory = "");
    ~DWIConvert();

    std::string setConversionMode(); //according input params, automatically set conversion mode
    void setConversionMode(const std::string conversionMode);

    std::string getConversionMode();

    //rewrite DWIConvert1, there is a test failure needing to debug
    //Please use readWrite interface instead
    //int DWIConvert1();

    DWIConverter *CreateDicomConverter(
            const std::string inputDicomDirectory,
            const bool useBMatrixGradientDirections,
            const bool transpose,
            const double smallGradientThreshold,
            const bool allowLossyConversion);

    //according Hans's code, encapsulate below readWrite method
    int readWrite();

    int read();

    int write();



//application case use in gtractConcatDWI.cxx
/*  in DWIConverter class:
 *  const DWIMetaDataDictionaryValidator::GradientTableType &GetDiffusionVectors();
 *  const std::vector<double> &GetBValues() const;
 *  Volume3DUnwrappedType::PointType GetOrigin() const;
 *  Volume3DUnwrappedType::Pointer getVolumePointer();
 * */

    // get and set methods for private data members

    const std::string &getInputVolume() const;

    void setInputVolume(const std::string &inputVolume);

    const std::string &getInputDicomDirectory() const;

    void setInputDicomDirectory(const std::string &inputDicomDirectory);

    const std::string &getInputBValues() const;

    void setInputBValues(const std::string &inputBValues);

    const std::string &getInputBVectors() const;

    void setInputBVectors(const std::string &inputBVectors);

    const std::string &getGradientVectorFile() const;

    void setGradientVectorFile(const std::string &gradientVectorFile);

    double getSmallGradientThreshold() const;

    void setSmallGradientThreshold(double smallGradientThreshold);

    bool isfMRIOutput() const;

    void setfMRIOutput(bool fMRIOutput);

    bool isTranspose() const;

    void setTranspose(bool transpose);

    bool isAllowLossyConversion() const;

    void setAllowLossyConversion(bool allowLossyConversion);

    bool isUseIdentityMeasurementFrame() const;

    void setUseIdentityMeasurementFrame(bool useIdentityMeasurementFrame);

    bool isUseBMatrixGradientDirections() const;

    void setUseBMatrixGradientDirections(bool useBMatrixGradientDirections);

    const std::string &getOutputVolume() const;

    void setOutputVolume(const std::string &outputVolume);

    const std::string &getOutputDirectory() const;

    void setOutputDirectory(const std::string &outputDirectory);

    const std::string &getOutputBValues() const;

    void setOutputBValues(const std::string &outputBValues);

    const std::string &getOutputBVectors() const;

    void setOutputBVectors(const std::string &outputBVectors);

private:
    std::string findFilenameExt(const std::string filename);

    //one of ["DicomToNrrd", "DicomToFSL", "NrrdToFSL", "FSLToNrrd","NrrdToNrrd", "FSLToFSL"]
    std::string m_conversionMode;

    std::string m_inputVolume;
    std::string m_inputDicomDirectory;
    std::string m_inputBValues;  //default: ""  ??
    std::string m_inputBVectors; //default: ""  ??
    std::string m_gradientVectorFile; //deprecated
    double m_smallGradientThreshold; //default = 0.2

    bool m_fMRIOutput; //default: false
    bool m_transpose; //default:false
    bool m_allowLossyConversion; //defualt: false
    bool m_useIdentityMeasurementFrame; //default: false
    bool m_useBMatrixGradientDirections; //default: false

    std::string m_outputVolume;
    std::string m_outputDirectory;  //default: "."
    std::string m_outputBValues; //default: ""
    std::string m_outputBVectors;//default: ""

    DWIConverter *m_converter;

};





#endif //BRAINSTOOLS_DWICONVERTLIB_H