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

#ifndef __BRAINSLinearModelerEPCA__h
#define __BRAINSLinearModelerEPCA__h

#include "Slicer3LandmarkIO.h"
#include "itkPoint.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include <string>
#include <fstream>
#include <iostream>

// typedef
const unsigned int PointDim = 3;
typedef itk::Point<double, PointDim>          PointType;
typedef std::map<std::string, PointType>      DatasetMapType;
typedef std::map<std::string, DatasetMapType> LmkDBType;
typedef vnl_matrix<double>                    MatrixType;
typedef vnl_vector_fixed<double, PointDim>    VectorType;
typedef std::map<std::string, MatrixType>     MatrixMapType;
typedef std::map<std::string, VectorType>     VectorMapType;

/*
 * Description:
 * Training linear model using EPCA
 * Implementation based on my MS thesis,
 * "A METHOD FOR AUTOMATED LANDMARK CONSTELLATION DETECTION USING
 * EVOLUTIONARY PRINCIPAL COMPONENTS AND STATISTICAL SHAPE MODELS"
 */
int BRAINSLinearModelerEPCAPrimary( int argc, char * argv[] );

/*
 * Build up the landmark database from a list of fcsv files
 * TODO: Explain each argument in the func
 * Input:
 * filename ...
 */
void CreateLmkDB( std::string filename, LmkDBType & baseLmkDB, LmkDBType & EPCALmkDB );

/*
 * Initialize X_i matrix from base landmarks
 */
MatrixType InitializeXi( LmkDBType & baseLmkDB );

/*
 * Compute the principal components of landmark vector space
 * in each iteration
 * TODO: Explain each argument in the func
 * Input:
 */
void ComputeEPCAModel( MatrixMapType & MMatrixMap, VectorMapType & SVectorMap, LmkDBType & baseLmkDB,
                       LmkDBType & EPCALmkDB );

/*
 * Compute the s_i vector from X_i matrix
 */
VectorType ComputeSVector( const MatrixType & X_i );

/*
 * Compute the I_si matrix from X_i and s_i
 */
MatrixType ComputeIsiMatrix( const unsigned int rows, const unsigned int columns, const VectorType & s_i );

#endif
