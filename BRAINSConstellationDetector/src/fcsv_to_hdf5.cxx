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
 %
 % Objective:
 % Train a linear model that can evolutionary estimate new landmarks from
 % already known landmarks.
 %
 % Inputs of the algorithm:
 % baselandmarks     - Base landmarks in training datasets
 % newLandmarks      - EPCA landmarks in training datasets
 % numBaseLandmarks  - Number of base landmarks
 % numNewLandmarks   - Number of EPCA landmarks
 % numDatasets       - Number of training datasets
 % dim               - Dimension of the landmarks
 %
 % Outputs: (variables stored in an output LLSModel files)
 % M (LLSMatrices)                  - Optimal linear combination of principal components
 %                                    in each iteration
 % s (LLSMeans)                     - The 3-by-1 mean of landmark vector space in each
 %                                    iteration
 % search_radius (LLSSearchRadii)   - Search radius of each EPCA landmark
 %
 */

/*
 fcsv_to_hdf5:
A program to automatically generate a linear model estimation for the secondary landmark points.

Arguments:
--landmarkTypesList {arg}         : A list of the landmarks, grouped by landmarkTypes ("AC,base\nPC,base\n\n")
                                    Landmarks should be ordered from the most prominent to the least prominent.

--landmarkGlobPattern {arg}       : A glob of landmark files to include (default: "*.fcsv") (Be sure to use quotes!).

--landmarksInformationFile {arg}  : The name of the output H5 file that saves all landmarks information.
                                    (This output is not needed by BCD)

--modelFile {arg}                 : MAIN OUTPUT OF THIS PROGRAM:
                                    name of HDF5 file containing BRAINSConstellationDetector Model file
                                    (LLSMatrices, LLSMeans and LLSSearchRadii).

Example Usage:
    ~/fcsv_to_hdf5 \
    --landmarkTypesList inputLandmarksNames.list \
    --landmarkGlobPattern ${DATA_PATH}/"\*_est.fcsv" \
    --landmarksInformationFile landmarksInfo.h5 \
    --modelFile outputLLSModel.h5

*/

/*
Development notes: Yes, I do use returns of stl containers instead of
passing the return containers by reference as arguments.  While this is slow
in 03, in 0x it's fast due to move constructors, and anyway, this code
benefits more from readability than speed.
*/

// I N C L U D E S ////////////////////////////////////////////////////////////

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <glob.h>
#include <wordexp.h>

#include "Slicer3LandmarkIO.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_matlab_print.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/vnl_matlab_write.h>

#include "itkPoint.h"
#include "fcsv_to_hdf5CLP.h"
#include "itkNumberToString.h"

#include "LLSModel.h"
#include "BRAINSThreadControl.h"
#include "itk_hdf5.h"
#include "itk_H5Cpp.h"
#include <BRAINSCommonLib.h>

// D E F I N E S //////////////////////////////////////////////////////////////

// will exit on the first error.
// If false, the program will
// continue and use "* * *" as the missing data.
typedef std::vector<std::pair<std::string, std::string> > SubjectFilenameVector;

// C L A S S E S //////////////////////////////////////////////////////////////

// F U N C T I O N S //////////////////////////////////////////////////////////

// ----------
// Remove matlab-invalid characters from a name
static std::string sanitize(const std::string& landmarkname)
{
  std::string ret(landmarkname);

  std::replace(ret.begin(), ret.end(), '-', '_');

  return ret;
}

// ----------
// Read in a csv file
static std::vector<std::vector<std::string> > read_csv(const std::string& file)
{
  std::vector<std::vector<std::string> > ret;

  std::ifstream infile(file.c_str() );
  std::string   line;

  while( getline(infile, line) )
    {
    std::istringstream       linestream(line);
    std::string              item;
    std::vector<std::string> line_vec;
    while( getline(linestream, item, ',') )
      {
      line_vec.push_back(item);
      }

    ret.push_back(line_vec);
    }

  return ret;
}

// ----------
// Extract the subjectID from the path
static std::string get_subjectid(const std::string& content)
{
  int start = content.length() - 1;
  int end = start;

  for( ; start >= 0; start-- )
    {
    if( content[start] == '_' )
      {
      end = start;
      }
    else if( content[start] == '/' )
      {
      start++;
      break;
      }
    }
  if( start == -1 )
    {
    start = 0;
    }
  return content.substr(start, end);
}

// ----------
// Check to see if a field is a comment
static bool is_comment(const std::string& field)
{
  for( unsigned int i = 0; i < field.length(); i++ )
    {
    if( (field[i] != ' ') && (field[i] != '\t') )
      {
      if( i + 1 >= field.length() )
        {
        return false;
        }
      if( field[i] == '#' )
        {
        return true;
        }
      else
        {
        return false;
        }
      }
    }
  return false;
}

// ----------
// Prints the arguments we're using for the user (both what they specified and
// the defaults that they didn't).
static void print_arguments_for_user(const std::string& file_glob_string, const std::string& outfile,
                                     const std::string& landmark_types_list_file)
{
  std::cout << std::endl;
  std::cout << "Landmark files:          " << file_glob_string << std::endl;
  std::cout << "Category list:           " << landmark_types_list_file << std::endl;
  std::cout << "Output file:             " << outfile << std::endl;
  std::cout << std::endl;
}

// ----------
// Build up list of subjectid/filename tuples
SubjectFilenameVector get_subject_filename_tuples(const std::string& file_glob_string)
{
  wordexp_t p;

  wordexp(file_glob_string.c_str(), &p, 0);
  char* *                  w = p.we_wordv;
  std::vector<std::string> files;
  for( unsigned int i = 0; i < p.we_wordc; i++ )
    {
    files.push_back(w[i]);
    // std::cout << "Adding file to processing list: " << w[i] << std::endl;
    }

  std::vector<std::pair<std::string, std::string> > subjects;
  for( std::vector<std::string>::const_iterator file_iter = files.begin(); file_iter != files.end(); ++file_iter )
    {
    subjects.push_back(std::make_pair(get_subjectid(*file_iter), *file_iter) );
    }
  return subjects;
}

// A map between the read in filename and the LandmarksMapType
typedef std::map<std::string, LandmarksMapType> FileToLandmarksMapType;

// ----------
// Get the data for the individual matrices for each landmark type
static FileToLandmarksMapType get_allFileToLandmarkMap(const std::vector<std::pair<std::string,
                                                                                   std::string> >& subjects)
{
  FileToLandmarksMapType allLandmarks;

  // Process all subjects
  for( std::vector<std::pair<std::string, std::string> >::const_iterator subject_iter = subjects.begin();
       subject_iter != subjects.end();
       ++subject_iter )
    {
    std::cout << "Processing " << subject_iter->second << "(" << subject_iter->first << ")" << std::endl;
    allLandmarks[subject_iter->first] = ReadSlicer3toITKLmk( subject_iter->second );
    if( allLandmarks[subject_iter->first].size() == 0 )
      {
      itkGenericExceptionMacro(<< "FAILED TO PROCESS FILE:  " << subject_iter->second);
      }
    }
  return allLandmarks;;
}

// Get the data for what landmark_types each landmark type is.
static std::vector<std::pair<std::string, std::vector<std::string> > > get_landmark_types(
  const std::string& landmark_types_list_file)
{
  static std::vector<std::pair<std::string, std::vector<std::string> > > landmark_types;
  static std::vector<std::string>                                        baseLandmarksNames;
  static std::vector<std::string>                                        newLandmarksNames;
  const std::vector<std::vector<std::string> >                           csv_data = read_csv(landmark_types_list_file);

  for( std::vector<std::vector<std::string> >::const_iterator iter = csv_data.begin(); iter != csv_data.end(); ++iter )
    {
    // Check for comment
    if( is_comment( (*iter)[0]) )
      {
      continue;
      }

    if( iter->size() != 2 )
      {
      std::cout << "Row in " << landmark_types_list_file << " must have precisely two entries:" << std::endl;
      for( unsigned int i = 0; i < iter->size(); i++ )
        {
        std::cout << (*iter)[i] << " ";
        }
      itkGenericExceptionMacro(<< "Row in " << landmark_types_list_file
                               << " must have precisely two entries:");
      }

    // Add it to the list.  We'll assume that there either are no duplicates, or
    // that any duplicates are intentional.
    if( (*iter)[1] == "baseLandmarks" )
      {
      baseLandmarksNames.push_back(sanitize( (*iter)[0]) );
      }
    else if( (*iter)[1] == "newLandmarks" )
      {
      newLandmarksNames.push_back(sanitize( (*iter)[0]) );
      }
    }

  std::pair<std::string, std::vector<std::string> > basePair;
  basePair = make_pair("baseLandmarks", baseLandmarksNames);
  landmark_types.push_back(basePair);

  std::pair<std::string, std::vector<std::string> > newPair;
  newPair = make_pair("newLandmarks", newLandmarksNames);
  landmark_types.push_back(newPair);

  return landmark_types;
}

static
void
WriteHDFStringList(H5::H5File & file,
                   const char * const name,
                   const std::vector<std::string> & stringList)
{
  const hsize_t numStrings(stringList.size() );
  H5::StrType   strType(H5::PredType::C_S1, H5T_VARIABLE);

  H5::DataSpace strSpace(1, &numStrings);
  H5::DataSet   strSet = file.createDataSet(name, strType, strSpace);

  std::vector<char const *> stringListCstr;

  stringListCstr.reserve(numStrings);
  for( std::vector<std::string>::const_iterator it = stringList.begin();
       it != stringList.end(); ++it )
    {
    stringListCstr.push_back( it->c_str() );
    }
  strSet.write(&(stringListCstr[0]), strType);
  strSet.close();
  strSpace.close();
  strType.close();
}

// M A I N ////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  if( outputFile == "" )
    {
    std::cerr << "Missing output file name" << std::endl;
    return EXIT_FAILURE;
    }
  if( landmarkTypesFile == "" )
    {
    std::cerr << "Missing landmark types filename name" << std::endl;
    return EXIT_FAILURE;
    }

  std::string cleanedGlobPattern = landmarkGlobPattern;
  while( cleanedGlobPattern.find("\\") != std::string::npos )
    {
    cleanedGlobPattern.replace(cleanedGlobPattern.find("\\"), 1, ""); // Replace
                                                                      // escape
                                                                      // characters
    }

  print_arguments_for_user(cleanedGlobPattern, outputFile, landmarkTypesFile);

  SubjectFilenameVector subjects = get_subject_filename_tuples(cleanedGlobPattern);

  std::vector<std::string> subjectIDVec;
  for( SubjectFilenameVector::iterator it = subjects.begin();
       it != subjects.end();
       ++it )
    {
    subjectIDVec.push_back(it->first);
    }

  const FileToLandmarksMapType allFileToLandmarkMap =
    get_allFileToLandmarkMap(subjects);

  typedef std::vector<std::pair<std::string, std::vector<std::string> > > LandmarkClassMapsTypes;

  const LandmarkClassMapsTypes landmark_types = get_landmark_types(landmarkTypesFile);

  std::vector<std::string> landmarkNames;

  // "byClassLandmarkMatrix" object contains information of each landmark:
  // - landmark type
  // - landmark name
  // - landmark coordinates in each data set (a 3xN matrix)

  std::map<std::string, std::vector<std::pair<std::string, vnl_matrix<double> > > > byClassLandmarkMatrix;
  for( LandmarkClassMapsTypes::const_iterator cit = landmark_types.begin();
       cit != landmark_types.end();
       ++cit )
    {
    std::vector<std::pair<std::string, vnl_matrix<double> > > perLandmarkMatrix;
    for( std::vector<std::string>::const_iterator lit = cit->second.begin();
         lit != cit->second.end();
         ++lit )
      {
      landmarkNames.push_back(*lit);
      vnl_matrix<double> CurrentLandmarkMatrix(allFileToLandmarkMap.size(), 3);
      int                subjCount = 0;
      for( FileToLandmarksMapType::const_iterator fit = allFileToLandmarkMap.begin();
           fit != allFileToLandmarkMap.end();
           ++fit )
        {
        const itk::Point<double, 3> & currPoint = fit->second.find(*lit)->second;
        CurrentLandmarkMatrix[subjCount][0] = currPoint[0];
        CurrentLandmarkMatrix[subjCount][1] = currPoint[1];
        CurrentLandmarkMatrix[subjCount][2] = currPoint[2];
        subjCount++;
        }
      perLandmarkMatrix.push_back( make_pair( (*lit), CurrentLandmarkMatrix) );
      }
    byClassLandmarkMatrix[cit->first] = perLandmarkMatrix;
    }

  // ================================
  // Write out all matricies to .hd5 format in addition to dumping to screen:
  H5::H5File output(outputFile.c_str(), H5F_ACC_TRUNC);
  WriteHDFStringList(output, "/SubjectIDs", subjectIDVec);
  WriteHDFStringList(output, "/LandmarkNames", landmarkNames);

  for( LandmarkClassMapsTypes::const_iterator cit = landmark_types.begin();
       cit != landmark_types.end();
       ++cit )
    {
    std::string groupName("/");
    groupName += cit->first; // baseLandmarks or newLandmarks
    H5::Group group = output.createGroup(groupName);

    for( std::vector<std::pair<std::string, vnl_matrix<double> > >::const_iterator lit =
           byClassLandmarkMatrix[cit->first].begin();
         lit != byClassLandmarkMatrix[cit->first].end();
         ++lit )
      {
      vnl_matrix<double> CurrentLandmarkMatrix = lit->second;
      hsize_t dims[2];
      dims[0] = CurrentLandmarkMatrix.rows();
      dims[1] = CurrentLandmarkMatrix.columns();

      // Print info on the screen
      std::cout << "=========" << cit->first << " "  << lit->first << "========== "
                << dims[0] << "x" << dims[1]
                << std::endl;
      vnl_matlab_print(std::cout, CurrentLandmarkMatrix, (lit->first).c_str() );

      // Write info to a .hd5 file
      try
        {
        H5::DataSpace dataspace(2, dims);
        H5::FloatType dataType(H5::PredType::NATIVE_DOUBLE );

        std::string matName(groupName);
        matName += "/";
        matName += lit->first;

        H5::DataSet dataset = output.createDataSet(matName, dataType, dataspace);
        dataset.write(CurrentLandmarkMatrix.data_block(), H5::PredType::NATIVE_DOUBLE, dataspace);
        dataspace.close();
        dataType.close();
        dataset.close();
        }
      // catch failure caused by the H5File operations
      catch( H5::FileIException& error )
        {
        error.printErrorStack();
        return -1;
        }

      // catch failure caused by the DataSet operations
      catch( H5::DataSetIException& error )
        {
        error.printErrorStack();
        return -1;
        }

      // catch failure caused by the DataSpace operations
      catch( H5::DataSpaceIException& error )
        {
        error.printErrorStack();
        return -1;
        }

      // catch failure caused by the DataSpace operations
      catch( H5::DataTypeIException& error )
        {
        error.printErrorStack();
        return -1;
        }
      }
    }
  output.close();

  const unsigned int numNewLandmarks = byClassLandmarkMatrix["newLandmarks"].size();
  const unsigned int numBaseLandmarks = byClassLandmarkMatrix["baseLandmarks"].size();
  const unsigned int dim = byClassLandmarkMatrix["newLandmarks"].begin()->second.cols(); //
                                                                                         // dim=3
                                                                                         // means
                                                                                         // 3D
                                                                                         // points
  const unsigned int numDatasets = byClassLandmarkMatrix["newLandmarks"].begin()->second.rows();

  // What the hell is this s matrix?
  vnl_matrix<double> s(dim, numNewLandmarks);
  /*
  //=====print for test====
  std::cout << "\n\n====================== " << s.rows() << "x" << s.cols() << std::endl;
  vnl_matlab_print(std::cout,s,"s");
  //=======================
  */
  std::vector<vnl_matrix<double> > lmk_est(numDatasets, s); // Initialize each
                                                            // element with s

  std::vector<vnl_matrix<double> > err(numDatasets, s); // to find search_radius

  std::vector<vnl_matrix<double> > W;
  std::vector<vnl_matrix<double> > M;

  // Initialize the landmark vector space
  vnl_matrix<double> Xi(dim * (numBaseLandmarks - 1), numDatasets);
  unsigned int       nb = 2;
  for( std::vector<std::pair<std::string, vnl_matrix<double> > >::const_iterator kit =
         ++(byClassLandmarkMatrix["baseLandmarks"].begin() );
       kit != byClassLandmarkMatrix["baseLandmarks"].end();
       ++kit )
    {
    vnl_matrix<double> lmkVec =
      ( (kit->second) - (byClassLandmarkMatrix["baseLandmarks"].begin()->second) ).transpose();
    Xi.update(lmkVec, (nb - 2) * dim, 0);
    nb++;
    }

  // Compute all principal components

  // Test Stuff

  vnl_matrix<double> ratioPC1(numNewLandmarks, 1, 0.0);
  vnl_matrix<double> ratioPC(numNewLandmarks, 1, 0.0);
  for( unsigned int i = 0; i < numNewLandmarks; i++ )
    {
    if( i > 0 )
      {
      vnl_matrix<double> lmkVec =
        ( (byClassLandmarkMatrix["newLandmarks"][i
                                                 - 1].second)
          - (byClassLandmarkMatrix["baseLandmarks"][0].second) ).transpose();

      // instead of  " Xi = [Xi ; lmkVec] " we write:
      vnl_matrix<double> Xi_temp = Xi;
      Xi.set_size(Xi.rows() + lmkVec.rows(), Xi.columns() );
      Xi.update(Xi_temp, 0, 0);
      Xi.update(lmkVec, Xi_temp.rows(), 0);
      }

    // ///remove Xi mean//////////
    vnl_matrix<double> Xi_mean(dim, 1);

    vnl_vector<double> mean_stacked(Xi.rows() ); // "mean_stack" is a vector
                                                 // with the size of Xi rows
                                                 // which is included the means
                                                 // of every rows of Xi
    for( unsigned int l = 0; l != Xi.rows(); l++ )
      {
      // mean_stacked[l] = Xi.extract(1,Xi.columns(),l,0).mean();
      mean_stacked[l] = Xi.get_row(l).mean();
      }
    for( unsigned int d = 0; d < dim; d++ )
      {
      double       sum = 0;
      unsigned int me = d;
      while( me <= (mean_stacked.size() - dim + d) )
        {
        sum += mean_stacked[me];
        me += dim;
        }

      Xi_mean(d, 0) = (sum / (numBaseLandmarks + i - 1) ); // "i" starts from
                                                           // zero instead of
                                                           // one, so we write
                                                           // "i-1" instead of
                                                           // "i-2"
      }

    s.set_column(i, Xi_mean.get_column(0) ); // s(:,i)=Xi_mean;

    vnl_matrix<double> I_si( (Xi_mean.rows() ) * (numBaseLandmarks + i - 1), (Xi_mean.cols() ) * numDatasets);
    for( unsigned int co = 0; co < numDatasets; co++ )
      {
      unsigned int ro = 0;
      while( ro < I_si.rows() )
        {
        I_si.update(Xi_mean, ro, co);
        ro += Xi_mean.rows();
        }
      }

    vnl_matrix<double> Xi_demeaned(Xi.rows(), Xi.cols() );
    Xi_demeaned = Xi - I_si;

    vnl_symmetric_eigensystem<double> eig( Xi_demeaned * Xi_demeaned.transpose() );
    vnl_matrix<double>                V = eig.V;
    vnl_matrix<double>                D = eig.D;

    ratioPC1(i, 0) = D.absolute_value_max() / D.absolute_value_sum();

    // Number of PCs should be chosen so that the following condition is
    // met:
    // sum(sum(D(:, end-numPCs+1))) / sum(D(:)) > 99%;
    // in this case, numPCs >= 1 will meet the requirement
    // or we can use as must PCs as possible
    // tol argument of the rank func is set conservatively
    // Adjust the tol argument for different training datasets
    // tolXi = 100;
    // numPCs = rank(Xi_demeaned, tolXi);

    unsigned int numPCs = 10;
    W.push_back( V.extract(V.rows(), numPCs, 0, V.cols() - numPCs) );

    /*  //////PRINT FOR TEST/////
    for (unsigned int wi=0; wi<W.size(); wi++)
    {
        std::cout << "\n\n====================== W[" << wi << "]==" << W[wi].rows() << "x" << W[wi].cols() << std::endl;//print for test
        vnl_matlab_print(std::cout,W[wi],"W[i]"); //, vnl_matlab_print_format_short_e );
      }
    */

    // Training phase-2
    // Train optimal linear relationship between already known landmark
    // vector space and the EPCA landmark to be estimated in each iteration.
    vnl_matrix<double> Zi = W[i].transpose() * Xi_demeaned; // PCA mapped space
    vnl_matrix<double> Yi = (byClassLandmarkMatrix["newLandmarks"][i].second)
      - (byClassLandmarkMatrix["baseLandmarks"][0].second);

    vnl_matrix<double> Zinv = vnl_matrix_inverse<double>(Zi * Zi.transpose() );
    vnl_matrix<double> Ci = Zinv * (Zi * Yi);
    M.push_back( W[i] * Ci );

    // Compute the estimation errors for training datasets
    vnl_matrix<double> Xi_t( (numBaseLandmarks + i - 1) * dim, 1, 0.0);
    vnl_matrix<double> x1_t = byClassLandmarkMatrix["baseLandmarks"][0].second;

    vnl_matrix<double> I_si_t( (Xi_mean.rows() ) * (numBaseLandmarks + i - 1), (Xi_mean.cols() ) * 1); //
                                                                                                       // I_si_t
                                                                                                       // =
                                                                                                       // repmat(s(:,
                                                                                                       // i),
                                                                                                       // numBaseLandmarks+i-2,
                                                                                                       // 1);
    unsigned int roo = 0;
    while( roo < I_si_t.rows() )
      {
      I_si_t.update(Xi_mean, roo, 0);
      roo += Xi_mean.rows();
      }

    vnl_matrix<double> xk_t;
    vnl_matrix<double> Xi_demeaned_t;
    vnl_matrix<double> lmk_groundTruth;
    for( unsigned int j = 0; j < numDatasets; j++ )
      {
      for( unsigned int k = 1; k < (numBaseLandmarks + i); k++ )
        {
        if( k < numBaseLandmarks )
          {
          xk_t = byClassLandmarkMatrix["baseLandmarks"][k].second;
          }
        else
          {
          xk_t = byClassLandmarkMatrix["newLandmarks"][k - numBaseLandmarks].second;
          }
        Xi_t.update( (xk_t.extract(1, xk_t.columns(), j,
                                   0).transpose() - x1_t.extract(1, x1_t.columns(), j,
                                                                 0).transpose() ), (k - 1) * dim, 0);
        }

      Xi_demeaned_t = Xi_t - I_si_t;
      lmk_est[j].update( (x1_t.extract(1, x1_t.columns(), j,
                                       0).transpose() + (M[i].transpose() * Xi_demeaned_t) ), 0, i);
      lmk_groundTruth = byClassLandmarkMatrix["newLandmarks"][i].second;
      err[j].update( (lmk_est[j].extract(lmk_est[j].rows(), 1, 0,
                                         i) - lmk_groundTruth.extract(1, lmk_groundTruth.columns(), j,
                                                                      0).transpose() ), 0, i);
      }
    }

  vnl_matrix<double> err_dist(numNewLandmarks, numDatasets, 0.0);
  for( unsigned int j = 0; j < numNewLandmarks; j++ )
    {
    for( unsigned int k = 0; k < numDatasets; k++ )
      {
      err_dist(j, k) = err[k].get_column(j).two_norm();
      }
    }

  vnl_matrix<double> err_mean(err_dist.rows(), 1); // "err_mean" is a vector
                                                   // with the size of err_dist
                                                   // rows which is included the
                                                   // means of every rows of Xi
  for( unsigned int l = 0; l != err_dist.rows(); l++ )
    {
    err_mean(l, 0) = err_dist.get_row(l).mean();
    }

  // //colculating "err_std"
  vnl_matrix<double> rep_err_mean( err_mean.rows() * 1, err_mean.cols() * numDatasets); //
                                                                                        // repmat(err_mean,1,numDatasets)
  for( unsigned int co = 0; co < rep_err_mean.cols(); co++ )
    {
    rep_err_mean.update(err_mean, 0, co);
    }

  vnl_matrix<double> sum_err_diff(numNewLandmarks, 1); // sum( (err_dist -
                                                       // repmat(err_mean,1,numDatasets)).^2,
                                                       // 2)
  for( unsigned int l = 0; l != numNewLandmarks; l++ )
    {
    sum_err_diff(l, 0) = (err_dist - rep_err_mean).get_row(l).squared_magnitude();
    }

  vnl_matrix<double> err_std = (sum_err_diff / numDatasets).apply(sqrt);

  vnl_matrix<double> err_max(err_dist.rows(), 1); // err_max =
                                                  // max(err_dist,[],2);
  for( unsigned int l = 0; l != err_dist.rows(); l++ )
    {
    err_max(l, 0) = err_dist.get_row(l).max_value();
    }

  vnl_matrix<double> err_min(err_dist.rows(), 1); // err_min =
                                                  // min(err_dist,[],2);
  for( unsigned int l = 0; l != err_dist.rows(); l++ )
    {
    err_min(l, 0) = err_dist.get_row(l).min_value();
    }
  // multiply matrix "err_std" by 3
  for( unsigned int l = 0; l != err_std.rows(); l++ )
    {
    err_std.scale_row(l, 3);
    }
  // multiply matrix "err_max" by 1.2
  for( unsigned int l = 0; l != err_max.rows(); l++ )
    {
    err_max.scale_row(l, 1.2);
    }

  vnl_matrix<double> search_radius = err_mean + err_std;
  vnl_matrix<double> search_radius_max = err_max;
  vnl_matrix<double> search_radius_min(numNewLandmarks, 1, 1.6);
  // if max<min, take min
  for( unsigned int l = 0; l < numNewLandmarks; l++ )
    {
    if( search_radius(l, 0) > search_radius_max(l, 0) )
      {
      search_radius(l, 0) = search_radius_max(l, 0);
      }
    if( search_radius(l, 0) < search_radius_min(l, 0) )
      {
      search_radius(l, 0) = search_radius_min(l, 0);
      }
    }

  // KENT:  there are 3 files with hard coded file names created by this program
  //       the information from these files is all inter-connected and
  //       combined they make up a single definition of the "model" that is to
  // be used
  //       These three files are read into the BRAINSConstellationDetector by
  // the function
  //       loadLLSModelMat( inputEPCAModelMat, inputEPCAModelTxt, llsMeans,
  // llsMatrices, searchRadii );
  //       at line 175 of BRAINSConstellationDetectorPrimary.cxx
  //       I want the following:
  //       1)
  //       Please write all the information encoded in these separate files
  // single
  //       HDF5 file as a single function called writeLLSModelMat(
  // commandLineSpecifiedModelFileName, byClassLandmarkMatrix, .... )
  //       2)
  //       The name of the model being generated should NOT be hard-coded, but
  // should be supplied on the command line.
  //       3)
  //       Convert loadLLSModelMat to read from the HDF5 file rather than these
  // 3 hard coded files names.
#if 1
  if( modelFile == "" )
    {
    std::cerr << "Missing model file name" << std::endl;
    return 1;
    }
  LLSModel theModel;
  theModel.SetFileName(modelFile);

  LLSModel::LLSMeansType       means;
  LLSModel::LLSMatricesType    matrices;
  LLSModel::LLSSearchRadiiType searchRadii;
  for( unsigned int i = 0; i < numNewLandmarks; i++ )
    {
    std::string         lmName(byClassLandmarkMatrix["newLandmarks"][i].first);
    std::vector<double> curVec(dim);
    for( unsigned int d = 0; d < dim; d++ )
      {
      curVec[d] = s(d, i);
      }
    means[lmName] = curVec;
    matrices[lmName] = M[i].transpose();
    searchRadii[lmName] = search_radius(i, 0);
    }
  theModel.SetLLSMeans(means);
  theModel.SetLLSMatrices(matrices);
  theModel.SetSearchRadii(searchRadii);
  theModel.Write();

#else
  // Write model parameters to file
  // txt version
  std::ofstream fid;
  fid.open("LME_EPCA.txt");
  itk::NumberToString<double> doubleToString;
  for( unsigned int i = 0; i < numNewLandmarks; i++ )
    {
    fid << byClassLandmarkMatrix["newLandmarks"][i].first << std::endl;
    for( unsigned int d = 0; d < dim; d++ )
      {
      fid << doubleToString(s(d, i) ) << " ";
      }
    fid << std::endl;
    fid << doubleToString(search_radius(i, 0) ) << std::endl;
    fid << M[i].rows() << std::endl;
    unsigned int l = 0;
    while( l < dim )
      {
      vnl_vector<double> temp = M[i].get_column(l);
      for( unsigned int ll = 0; ll < M[i].rows(); ll++ )
        {
        fid << doubleToString(temp[ll]) << " ";
        }
      fid << std::endl;
      l++;
      }

    fid << std::endl;
    }
  fid.close();

  // mat version
  fid.open("processingList.txt");
  for( unsigned int i = 0; i < numNewLandmarks; i++ )
    {
    fid << byClassLandmarkMatrix["newLandmarks"][i].first << std::endl;
    fid << doubleToString(search_radius(i, 0) ) << std::endl;
    }
  fid.close();

  fid.open("LME_EPCA.m");
  fid << "clear" << std::endl;
  for( unsigned int i = 0; i < numNewLandmarks; i++ )
    {
    std::string name = byClassLandmarkMatrix["newLandmarks"][i].first;

    // write s
    if( name.size() > 16 )
      {
      fid << name.substr(0, 16) << "__s = [";
      }
    else
      {
      fid << name << "__s = [";
      }
    for( unsigned int d = 0; d < dim; d++ )
      {
      fid << doubleToString(s(d, i) ) << " ";
      }
    fid << "];" << std::endl;

    // write M
    if( name.size() > 16 )
      {
      fid << name.substr(0, 16) << "__M = [";
      }
    else
      {
      fid << name << "__M = [";
      }
    unsigned int l = 0;
    while( l < dim )
      {
      vnl_vector<double> temp = M[i].get_column(l);
      for( unsigned int ll = 0; ll < M[i].rows(); ll++ )
        {
        fid << doubleToString(temp[ll]) << " ";
        }
      fid << ";" << std::endl;
      l++;
      }

    fid << "];" << std::endl;
    }
  fid.close();

  fid.open("LME_EPCA.mat");
  for( unsigned int i = 0; i < numNewLandmarks; i++ )
    {
    std::string name = byClassLandmarkMatrix["newLandmarks"][i].first;
    std::string string_name_s;
    std::string string_name_M;

    // write s
    if( name.size() > 16 )
      {
      string_name_s = name.substr(0, 16) + "__s";
      }
    else
      {
      string_name_s = name + "__s";
      }
    vnl_matrix<double> name_s(1, dim, 0.0);
    for( unsigned int d = 0; d < dim; d++ )
      {
      name_s(0, d) = s(d, i);
      }
    vnl_matlab_write(fid, name_s.data_array(), name_s.rows(), name_s.cols(), string_name_s.c_str() );

    // write M
    if( name.size() > 16 )
      {
      string_name_M = name.substr(0, 16) + "__M";
      }
    else
      {
      string_name_M = name + "__M";
      }
    vnl_matrix<double> name_M(dim, M[i].rows(), 0.0);
    unsigned int       l = 0;
    while( l < dim )
      {
      vnl_vector<double> temp = M[i].get_column(l);
      for( unsigned int ll = 0; ll < M[i].rows(); ll++ )
        {
        name_M(l, ll) = temp[ll];
        }
      l++;
      }

    vnl_matlab_write(fid, name_M.data_array(), name_M.rows(), name_M.cols(), string_name_M.c_str() );
    }
  fid.close();
#endif

  std::cout << "LME-EPCA coefficients have been written to disk." << std::endl;

  return 0;
}
