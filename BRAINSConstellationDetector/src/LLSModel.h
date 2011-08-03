#ifndef LLSModel_h
#define LLSModel_h

#include <string>
#include <map>
#include <vector>
#include "vnl/vnl_matrix.h"
#include "itkMacro.h"

#if  ITK_VERSION_MAJOR >= 4
#include "itk_hdf5.h"
#include "itk_H5Cpp.h"
#else
#include "hdf5.h"
#include "H5Cpp.h"
#endif

class LLSModel
{
public:
  LLSModel();
  ~LLSModel();
  typedef std::map<std::string, std::vector<double> > LLSMeansType;
  typedef vnl_matrix<double>                          MatrixType;
  typedef std::map<std::string, MatrixType>           LLSMatricesType;
  typedef std::map<std::string, double>               LLSSearchRadiiType;

  void SetFileName(const std::string & fileName);

  //
  // read & write return zero on success -1 otherwise.
  int Read();

  int Write();

  void SetLLSMeans(const LLSMeansType & llsMeans);

  const LLSMeansType & GetLLSMeans();

  void SetLLSMatrices(const LLSMatricesType & llsMatrices);

  const LLSMatricesType & GetLLSMatrices();

  void SetSearchRadii(const LLSSearchRadiiType & llsSearchRadii);

  const LLSSearchRadiiType & GetSearchRadii();

private:
  // private methods
  void WriteVector(const std::string & path, const std::vector<double> & vec);

  void WriteMatrix(const std::string & path, const MatrixType & matrix);

  void WriteScalar(const std::string & path, const double & value);

  double ReadScalar(const std::string & DataSetName);

  std::vector<double> ReadVector(const std::string & DataSetName);

  MatrixType ReadMatrix(const std::string & DataSetName);

private:
  std::string               m_FileName;
  LLSMeansType              m_LLSMeans;
  LLSMatricesType           m_LLSMatrices;
  LLSSearchRadiiType        m_LLSSearchRadii;
  H5::H5File *              m_H5File;
  static const char * const m_LLSMeansGroupName;
  static const char * const m_LLSMatricesGroupName;
  static const char * const m_LLSSearchRadiiGroupName;
};

#endif // LLSModel_h
