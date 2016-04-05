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
#ifndef LLSModel_h
#define LLSModel_h

#include <string>
#include <map>
#include <vector>
#include "vnl/vnl_matrix.h"
#include "itkMacro.h"

#include "itk_hdf5.h"
#include "itk_H5Cpp.h"

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

  const LLSMeansType & GetLLSMeans() const ;

  void SetLLSMatrices(const LLSMatricesType & llsMatrices);

  const LLSMatricesType & GetLLSMatrices() const ;

  void SetSearchRadii(const LLSSearchRadiiType & llsSearchRadii);

  const LLSSearchRadiiType & GetSearchRadii() const;

private:
  // private methods
  void WriteVector(const std::string & path, const std::vector<double> & vec);

  void WriteMatrix(const std::string & path, const MatrixType & matrix);

  void WriteScalar(const std::string & path, const double & value);

  void WriteString(const std::string & path, const std::string & strname);

  double ReadScalar(const std::string & DataSetName);

  std::string ReadString(const std::string & DataSetName);

  std::vector<double> ReadVector(const std::string & DataSetName);

  MatrixType ReadMatrix(const std::string & DataSetName);

private:
  std::string               m_FileName;
  LLSMeansType              m_LLSMeans;
  LLSMatricesType           m_LLSMatrices;
  LLSSearchRadiiType        m_LLSSearchRadii;
  H5::H5File *              m_H5File;
  static const char * const m_LLSVersionGroupName;
  static const char * const m_LLSMeansGroupName;
  static const char * const m_LLSMatricesGroupName;
  static const char * const m_LLSSearchRadiiGroupName;
};

#endif // LLSModel_h
