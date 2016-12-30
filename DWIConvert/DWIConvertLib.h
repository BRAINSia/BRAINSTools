//
// Created by Hui Xie on 12/26/16.
//

#ifndef BRAINSTOOLS_DWICONVERTLIB_H
#define BRAINSTOOLS_DWICONVERTLIB_H

#include "DWIConverter.h"
#include "DWIConverterFactory.h"


class DWIConvert {
public:
    //todo: get and set functions later
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

    const std::string &getM_inputVolume() const;

    void setM_inputVolume(const std::string &m_inputVolume);

    const std::string &getM_inputDicomDirectory() const;

    void setM_inputDicomDirectory(const std::string &m_inputDicomDirectory);

    const std::string &getM_inputBValues() const;

    void setM_inputBValues(const std::string &m_inputBValues);

    const std::string &getM_inputBVectors() const;

    void setM_inputBVectors(const std::string &m_inputBVectors);

    const std::string &getM_gradientVectorFile() const;

    void setM_gradientVectorFile(const std::string &m_gradientVectorFile);

    double getM_smallGradientThreshold() const;

    void setM_smallGradientThreshold(double m_smallGradientThreshold);

    bool isM_fMRIOutput() const;

    void setM_fMRIOutput(bool m_fMRIOutput);

    bool isM_transpose() const;

    void setM_transpose(bool m_transpose);

    bool isM_allowLossyConversion() const;

    void setM_allowLossyConversion(bool m_allowLossyConversion);

    bool isM_useIdentityMeasurementFrame() const;

    void setM_useIdentityMeasurementFrame(bool m_useIdentityMeasurementFrame);

    bool isM_useBMatrixGradientDirections() const;

    void setM_useBMatrixGradientDirections(bool m_useBMatrixGradientDirections);

    const std::string &getM_outputVolume() const;

    void setM_outputVolume(const std::string &m_outputVolume);

    const std::string &getM_outputDirectory() const;

    void setM_outputDirectory(const std::string &m_outputDirectory);

    const std::string &getM_outputBValues() const;

    void setM_outputBValues(const std::string &m_outputBValues);

    const std::string &getM_outputBVectors() const;

    void setM_outputBVectors(const std::string &m_outputBVectors);

    std::string m_outputDirectory;  //default: "."
    std::string m_outputBValues; //default: ""
    std::string m_outputBVectors;//default: ""

    DWIConvert();

    std::string setConversionMode(); //according input params, automatically set conversion mode
    void setConversionMode(const std::string conversionMode);

    std::string getConversionMode();

    //rewrite DWIConvert1, there is a test failure needing to debug
    //Please use DWIConvert2 interface instead
    int DWIConvert1();

    DWIConverter *CreateDicomConverter(
            const std::string inputDicomDirectory,
            const bool useBMatrixGradientDirections,
            const bool transpose,
            const double smallGradientThreshold,
            const bool allowLossyConversion);

    //according Hans's code, encapsulate the DWIConvert2
    int DWIConvert2();

    int DWIConvertRead();

    int DWIConvertWrite();



//application case use in gtractConcatDWI.cxx
/*  in DWIConverter class:
 *  const DWIMetaDataDictionaryValidator::GradientTableType &GetDiffusionVectors();
 *  const std::vector<double> &GetBValues() const;
 *  Volume3DUnwrappedType::PointType GetOrigin() const;
 *  Volume3DUnwrappedType::Pointer getVolumePointer();
 * */


private:
    std::string findFilenameExt(const std::string filename);

    std::string m_conversionMode;

};





#endif //BRAINSTOOLS_DWICONVERTLIB_H