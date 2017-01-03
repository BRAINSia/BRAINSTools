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
    DWIConvert(const std::string& outputVolume, const std::string& inputVolume, const std::string& inputDicomDirectory);
    DWIConvert(const std::string& inputVolume, const std::string& inputDicomDirectory);
    ~DWIConvert();

    int read();

    int write(const std::string& outputVolume);
    int write();

    DWIConverter *getConverter() const;

    // get and set methods for private data members
    void setInputFileType(const std::string& inputVolume, const std::string& inputDicomDirectory);

    void setOutputFileType(const std::string& outputVolume);

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

    std::string getInputFileType();
    std::string getOutputFileType();
    std::string findFilenameExt(const std::string filename);

    DWIConverter *CreateDicomConverter(
            const std::string inputDicomDirectory,
            const bool useBMatrixGradientDirections,
            const bool transpose,
            const double smallGradientThreshold,
            const bool allowLossyConversion);

    std::string m_inputFileType;
    std::string m_outputFileType;

    std::string m_inputVolume;
    std::string m_inputDicomDirectory;
    std::string m_inputBValues;  //default: ""  for FSL file
    std::string m_inputBVectors; //default: ""  for FSL file
    std::string m_gradientVectorFile; //deprecated
    double m_smallGradientThreshold; //default = 0.2

    bool m_fMRIOutput; //default: false
    bool m_transpose; //default:false
    bool m_allowLossyConversion; //defualt: false
    bool m_useIdentityMeasurementFrame; //default: false
    bool m_useBMatrixGradientDirections; //default: false

    std::string m_outputVolume;
    std::string m_outputDirectory;  //default: "."
    std::string m_outputBValues; //default: ""  for FSL file
    std::string m_outputBVectors;//default: ""  for FSL file

    DWIConverter *m_converter;

};





#endif //BRAINSTOOLS_DWICONVERTLIB_H