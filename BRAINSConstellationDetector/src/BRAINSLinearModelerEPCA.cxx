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
 * Author: Wei Lu
 * at Psychiatry Imaging Lab, University of Iowa Health Care, 2010
 */

#include "BRAINSLinearModelerEPCA.h"
#include "BRAINSLinearModelerEPCACLP.h"
#include "BRAINSThreadControl.h"
#include <BRAINSCommonLib.h>

int main( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  LmkDBType                                             baseLmkDB; // in the format of [landmarkID][datasetID]
  LmkDBType                                             EPCALmkDB;
  CreateLmkDB( inputTrainingList, baseLmkDB, EPCALmkDB );

  MatrixMapType MMatrixMap; // Principal components of landmark vector
                            // space in each iteration
  VectorMapType SVectorMap; // s vectors in each iteration
  ComputeEPCAModel( MMatrixMap, SVectorMap, baseLmkDB, EPCALmkDB );

  // Write model to Matlab binary file

  return EXIT_SUCCESS;
}

void
CreateLmkDB( std::string filename, LmkDBType & baseLmkDB, LmkDBType & EPCALmkDB )
{
  // Read in list of landmark list file
  std::ifstream myfile( filename.c_str() );

  if( !myfile.is_open() )
    {
    itkGenericExceptionMacro(<< "Cannot open training landmark list file!"
                             << filename.c_str() );
    }

  // for each file enlisted on the list of landmark list file
  std::string line;
  while( getline( myfile, line ) )
    {
    // for each landmark in the landmark list file
    LandmarksMapType                 lmkMap = ReadSlicer3toITKLmk( line );
    LandmarksMapType::const_iterator itLmk = lmkMap.begin();
    while( itLmk != lmkMap.end() )
      {
      std::string name = itLmk->first;
      // Save base landmarks
      if( ( name.compare( "AC" ) == 0 ) ||
          ( name.compare( "PC" ) == 0 ) ||
          ( name.compare( "RP" ) == 0 ) ||
          ( name.compare( "VN4" ) == 0 ) ||
          ( name.compare( "LE" ) == 0 ) ||
          ( name.compare( "RE" ) == 0 ) )
        {
        baseLmkDB[name][line] = itLmk->second;
        }
      else
        {
        EPCALmkDB[name][line] = itLmk->second;
        }
      ++itLmk;
      }
    }

  // Sanity check
  // Check if the number of landmarks are the same in each landmark list file
  // First check for number of base landmarks
    {
    LmkDBType::const_iterator itDB = baseLmkDB.begin();
    unsigned int              numLmks = itDB->second.size();
    while( itDB != baseLmkDB.end() )
      {
      if( itDB->second.size() != numLmks )
        {
        itkGenericExceptionMacro(<< "Error: number of landmark \"" << itDB->first
                                 << "\" in training list files mismatched!");
        }
      ++itDB;
      }
    }

  // Then check for the number of EPCA landmarks
    {
    LmkDBType::const_iterator itDB = EPCALmkDB.begin();
    unsigned int              numLmks = itDB->second.size();
    while( itDB != EPCALmkDB.end() )
      {
      if( itDB->second.size() != numLmks )
        {
        itkGenericExceptionMacro(<< "Error: number of landmark \"" << itDB->first
                                 << "\" in training list files mismatched!");
        }
      ++itDB;
      }
    }
}

MatrixType
InitializeXi( LmkDBType & baseLmkDB )
{
  const unsigned int numBaseLmks( baseLmkDB.size() );
  const unsigned int numDatasets( baseLmkDB.begin()->second.size() );

  MatrixType X_i;

  X_i.set_size( ( numBaseLmks - 1 ) * PointDim, numDatasets );

  // Assert RP (MPJ) exists
  if( baseLmkDB.find( "RP" ) == baseLmkDB.end() )
    {
    itkGenericExceptionMacro(<< "Error: RP (MPJ) landmark is missing!")
    }

  // for each base landmark
  unsigned int              k = 0; // landmark index
  LmkDBType::const_iterator itDB = baseLmkDB.begin();
  while( itDB != baseLmkDB.end() )
    {
    if( itDB->first.compare( "RP" ) != 0 )
      {
      unsigned int                     j = 0; // dataset index
      DatasetMapType                   datasetMap = itDB->second;
      LandmarksMapType::const_iterator itDataset = datasetMap.begin();
      while( itDataset != datasetMap.end() )
        {
        std::string     datasetId = itDataset->first;
        const PointType lmkRP( baseLmkDB["RP"][datasetId] );
        for( unsigned int dim = 0; dim < PointDim; ++dim )
          {
          X_i( k * PointDim + dim, j ) = itDataset->second[dim] - lmkRP[dim];
          }
        ++j;
        ++itDataset;
        }

      ++k;
      }
    ++itDB;
    }

  return X_i;
}

void
ComputeEPCAModel( MatrixMapType & MMatrixMap, VectorMapType & SVectorMap,
                  LmkDBType & baseLmkDB, LmkDBType & EPCALmkDB )
{
  // Deliberately add the following line to eliminate the "unused warning"
  MatrixType bogusMMatrix;

  MMatrixMap["bogus M_i matrix"] = bogusMMatrix;

  // Initialize the landmark vector space X_i matrix
  MatrixType X_i = InitializeXi( baseLmkDB );

  // Evolutionarily construct X_i, compute s_i, W_i, and M_i in each iteration
  LmkDBType::const_iterator itDB = EPCALmkDB.begin();
  unsigned int              k = 0; // landmark index
  // while ( itDB != EPCALmkDB.end() - 1 ) // NO end() - 1 in map iterator?
  const unsigned int numEPCALmks = EPCALmkDB.size();
  while( k < numEPCALmks )
    {
    // Update X_i
    if( k > 0 )
      {
      MatrixType X_iLast( X_i );
      X_i.set_size(X_iLast.rows() + PointDim, X_iLast.columns() );
      for( unsigned int row = 0; row < X_iLast.rows(); ++row )
        {
        X_i.set_row( row, X_iLast.get_row( row ) );
        }

      const unsigned int               numBaseLmks( baseLmkDB.size() );
      unsigned int                     j = 0; // dataset index
      DatasetMapType                   datasetMap = itDB->second;
      LandmarksMapType::const_iterator itDataset = datasetMap.begin();
      while( itDataset != datasetMap.end() )
        {
        std::string     datasetId = itDataset->first;
        const PointType lmkRP( baseLmkDB["RP"][datasetId] );
        for( unsigned int dim = 0; dim < PointDim; ++dim )
          {
          X_i( ( numBaseLmks + k - 2 ) * PointDim + dim, j ) = itDataset->second[dim] - lmkRP[dim];
          }
        ++j;
        ++itDataset;
        }
      }

    // Compute si
    VectorType s_i( ComputeSVector( X_i ) );
    SVectorMap[EPCALmkDB.begin()->first] = s_i;

    // remove I_si of X_i
    MatrixType I_si( ComputeIsiMatrix( X_i.rows(), X_i.columns(), s_i ) );

    MatrixType X_i0Mean = X_i - I_si;

    // Compute W_i
    MatrixType         W_i; // principal components/eigenvectors of X_i*X_i'
    vnl_vector<double> D;   // eigenvalue of X_i*X_i'
    if( !vnl_symmetric_eigensystem_compute(X_i0Mean * X_i0Mean.transpose(), W_i, D) )
      {
      itkGenericExceptionMacro(<< "Error: vnl_symmetric_eigensystem_compute failed.")
      }

    // TODO -------------------------------------------------------

    // Compute W_ri

    // Construct Y_i

    // Compute C_opt

    // Compute and save M_i

    // DEBUG
    // std::cout << "X_i( i = " << k << " ) = \n" << X_i << std::endl;
    // std::cout << "I_si = \n" << I_si << std::endl;
    // std::cout << "s_i( i = " << k << " ) = " << s_i << std::endl;
    // std::cout << "X_i0Mean( i = " << k << " ) = \n" << X_i0Mean << std::endl;
    // std::cout << "D( i = " << k << " ) = \n" << D << std::endl;
    std::cout << "W_i( i = " << k << " ) = \n" << W_i << std::endl;

    if( k++ > 0 )
      {
      ++itDB;
      }
    }
}

VectorType
ComputeSVector( const MatrixType & X_i )
{
  const unsigned int numDataset( X_i.columns() );

  // Sanity check for X_i
  if( X_i.rows() % PointDim != 0 )
    {
    itkGenericExceptionMacro(<< "Error: Bad X_i!");
    }
  const unsigned int numLmks( X_i.rows() / PointDim );
  VectorType         s_i;
  s_i.fill( 0 );
  for( unsigned int dim = 0; dim < PointDim; ++dim )
    {
    for( unsigned int j = 0; j < numDataset; ++j )
      {
      for( unsigned int k = 0; k < numLmks; ++k )
        {
        s_i[dim] += X_i( k * PointDim + dim, j );
        }
      }
    s_i[dim] /= numDataset * numLmks;
    }
  return s_i;
}

MatrixType
ComputeIsiMatrix( const unsigned int rows, const unsigned int columns, const VectorType & s_i )
{
  const unsigned int numDataset( columns );

  // Sanity check for X_i
  if( rows % PointDim != 0 )
    {
    itkGenericExceptionMacro(<< "Error: Bad X_i!" );
    }
  const unsigned int numLmks( rows / PointDim );

  MatrixType I_si;
  I_si.set_size( rows, columns );
  for( unsigned int j = 0; j < numDataset; ++j )
    {
    for( unsigned int k = 0; k < numLmks; ++k )
      {
      for( unsigned int dim = 0; dim < PointDim; ++dim )
        {
        I_si( k * PointDim + dim, j ) = s_i[dim];
        }
      }
    }
  return I_si;
}
