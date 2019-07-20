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
//
//
// //////////////////////////////////////////////////////////////////////////////
//
// Clustering using Minimum Spanning Tree edge breaking
// MST constructed using Delaunay triangulation followed by Kruskal's algorithm
// (best for sparse graphs)
//
// Designed to be used iteratively
//
// Follows the heuristic described in:
// Cocosco, C.A., Zijdenbos, A.P., Evans, A.C.: A fully automatic and robust
// brain MRI tissue classification method. Medical Image Analysis 7 (2003)
// 513-527
//
//
//
// //////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 5/2008

#ifndef __QHullMSTClusteringProcess_h
#define __QHullMSTClusteringProcess_h

#include <vector>
#include "MSTEdge.h"

#include "vnl/vnl_vector.h"
/**
 * \class QHullMSTClusteringProcess
 */
class QHullMSTClusteringProcess
{
public:
  QHullMSTClusteringProcess();
  ~QHullMSTClusteringProcess();

  using VertexType = vnl_vector< float >;
  using EdgeType = MSTEdge;

  using VertexList = std::vector< VertexType >;

  void
  SetInputVertices( const VertexList & l );

  // Break edges with threshold value T, and then cluster based on MST edge
  // connectivity
  //
  // Returns the number of clusters and fills the cluster map array
  unsigned int
  GetClusters( unsigned int * maps, float T );

  inline void
  SortOn()
  {
    m_SortFlag = true;
  }

  inline void
  SortOff()
  {
    m_SortFlag = false;
  }

private:
  unsigned int m_NumberOfVertices;

  MSTEdge * m_MSTEdges;

  float * m_NodeAverages;

  bool m_SortFlag;
};

#endif
