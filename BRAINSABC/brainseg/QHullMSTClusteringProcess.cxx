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
#include "Heap.h"
#include "QHullMSTClusteringProcess.h"

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <cstring>

extern "C"
{
#include "qhull.h"
#include "mem.h"
#include "qset.h"
}

#define BUF_SIZE 2048
#define QHULL_COMMAND "qhull d QJ Tv "
// #define QHULL_COMMAND "qhull d QJ "

/**
 * \class MSTCluster
 * For sorting cluster maps based on size (descending)
 */
class MSTCluster
{
public:
  unsigned int map;
  unsigned int size;

  MSTCluster()
  {
    this->map = 0;
    this->size = 0;
  }

  ~MSTCluster() {}

  inline MSTCluster &
  operator=(const MSTCluster & c)
  {
    this->map = c.map;
    this->size = c.size;
    return *this;
  }

  inline bool
  operator<(const MSTCluster & c) const
  {
    return this->size > c.size;
  }
};

QHullMSTClusteringProcess ::QHullMSTClusteringProcess()
{
  m_NumberOfVertices = 0;

  m_MSTEdges = 0;

  m_NodeAverages = 0;

  m_SortFlag = false;
}

QHullMSTClusteringProcess ::~QHullMSTClusteringProcess()
{
  delete[] m_MSTEdges;
  delete[] m_NodeAverages;
}

void
QHullMSTClusteringProcess ::SetInputVertices(const VertexList & vlist)
{
  delete[] m_MSTEdges;
  delete[] m_NodeAverages;

  m_NumberOfVertices = vlist.size();
  if (m_NumberOfVertices == 0)
  {
    return;
  }

  const unsigned int dim = vlist[0].size();

  // Allocate memory to store elements of the list
  //  coordT *points = (coordT*)
  // malloc((dim+1)*m_NumberOfVertices*sizeof(coordT));
  coordT * const points = new coordT[(dim + 1) * m_NumberOfVertices];
  // Store each coordinate in an array element
  for (unsigned int k = 0; k < m_NumberOfVertices; k++)
  {
    const VertexType v = vlist[k];
    const long       pos = k * (dim + 1);
#if 1
    double sumS = 0;
    for (unsigned int j = 0; j < dim; j++)
    {
      points[pos + j] = v[j];
      sumS += v[j] * v[j];
    }
    points[pos + dim] = sumS;
#else
    points[pos + dim] = 0.0;
#endif
  }

  // Call out to qhull library to compute DeLaunay triangulation
  qh_init_A(stdin, stdout, stderr, 0, NULL);
  const int exitcode = setjmp(qh errexit);

  // Compute Delaunay triangulation and insert it to heap
  Heap<MSTEdge> delaunayHeap;
  if (!exitcode)
  {
    char options[BUF_SIZE];
    // Add extra options here
    strcpy(options, QHULL_COMMAND);

    qh_initflags(options);
    qh_setdelaunay(dim + 1, m_NumberOfVertices, points);
    qh_init_B(points, m_NumberOfVertices, dim + 1, False);
    qh_qhull();
    qh_check_output();

    facetT * facet;
    long     numfacets = 0;
    FORALLfacets
    {
      if (!facet->upperdelaunay)
      {
        ++numfacets;
      }
    }

    vertexT *  vertex;
    vertexT ** vertexp;
    // Insert DT edges to heap
    FORALLfacets
    {
      if (!facet->upperdelaunay)
      {
        std::vector<unsigned int> ids;
        ids.resize(dim + 1);
        FOREACHvertex_(facet->vertices) { ids.push_back(qh_pointid(vertex->point)); }
        for (unsigned int s = 0; s < ids.size(); s++)
        {
          for (unsigned int t = s + 1; t < ids.size(); t++)
          {
            MSTEdge e;
            e.i = ids[s];
            e.j = ids[t];
            const VertexType dij = vlist[e.i] - vlist[e.j];
            e.dist = dij.squared_magnitude();
            delaunayHeap.Insert(e);
          }
        }
      }
    }
  }
  delete[] points;

  // Free allocated memory
  qh NOerrexit = True;
  qh_freeqhull(False);
  int curlong, totlong;
  qh_memfreeshort(&curlong, &totlong);
  if (curlong || totlong)
  {
    std::cerr << "qhull internal warning (main): did not free " << totlong << "bytes of long memory (" << curlong
              << "pieces)" << std::endl;
    throw "QHull error";
  }

  // Map vertex to set (tree) it belongs to, for cycle test
  unsigned int * const treeMap = new unsigned int[m_NumberOfVertices];
  for (unsigned int i = 0; i < m_NumberOfVertices; i++)
  {
    treeMap[i] = i;
  }

  // Number of edges in MST
  unsigned int edgeCount = 0;

  // Build MST using Kruskal's algorithm
  // Edges added in ascending order
  m_MSTEdges = new MSTEdge[m_NumberOfVertices - 1];
  while (!delaunayHeap.IsEmpty())
  {
    MSTEdge minEdge = delaunayHeap.ExtractMinimum();

    const unsigned int a = minEdge.i;
    const unsigned int b = minEdge.j;

    const unsigned int map1 = treeMap[a];
    const unsigned int map2 = treeMap[b];

    // Skip if they belong to the same tree (will form cycle)
    if (map1 == map2)
    {
      continue;
    }
    // Merge trees
    for (unsigned int k = 0; k < m_NumberOfVertices; k++)
    {
      if (treeMap[k] == map2)
      {
        treeMap[k] = map1;
      }
    }
    m_MSTEdges[edgeCount] = minEdge;
    ++edgeCount;

    // See if a tree is formed already
    if (edgeCount == (m_NumberOfVertices - 1))
    {
      break;
    }
  }

  delete[] treeMap;
  if (edgeCount != (m_NumberOfVertices - 1))
  {
    std::cerr << "MST construction failed, E != (V-1)" << std::endl;
    throw "MST construction failed, E != (V-1)";
  }

  m_NodeAverages = new float[m_NumberOfVertices];
  // Compute node averages
  for (unsigned int k = 0; k < m_NumberOfVertices; k++)
  {
    m_NodeAverages[k] = 0.0;
  }

  unsigned int * const countArray = new unsigned int[m_NumberOfVertices];
  for (unsigned int k = 0; k < m_NumberOfVertices; k++)
  {
    countArray[k] = 0;
  }

  {
    std::vector<std::vector<double>> nodeDistances;
    nodeDistances.resize(m_NumberOfVertices + 1);
    for (unsigned int i = 0; i <= m_NumberOfVertices + 1; i++)
    {
      std::vector<double> dlist;
      nodeDistances.push_back(dlist);
    }
    for (unsigned int k = 0; k < (m_NumberOfVertices - 1); k++)
    {
      const unsigned int a = m_MSTEdges[k].i;
      const unsigned int b = m_MSTEdges[k].j;

      m_NodeAverages[a] += m_MSTEdges[k].dist;
      countArray[a]++;

      m_NodeAverages[b] += m_MSTEdges[k].dist;
      countArray[b]++;

      nodeDistances[a].push_back(m_MSTEdges[k].dist);
      nodeDistances[b].push_back(m_MSTEdges[k].dist);
    }
    for (unsigned int k = 0; k < m_NumberOfVertices; k++)
    {
      // Use mean OR...
      // if (countArray[k] != 0)
      //  m_NodeAverages[k] /= countArray[k];
      // Median instead?
      if (countArray[k] != 0)
      {
        m_NodeAverages[k] = heapMedian(nodeDistances[k], nodeDistances[k].size());
      }
    }
  }
  delete[] countArray;
}

unsigned int
QHullMSTClusteringProcess ::GetClusters(unsigned int * treeMap, float T)
{
  // Get number of vertices and edges
  const unsigned int v = m_NumberOfVertices;
  const unsigned int e = v - 1;

  // Allocate edge break flag array
  unsigned char * const breakArray = new unsigned char[e];

  for (unsigned int k = 0; k < e; k++)
  {
    breakArray[k] = 0;
  }

  // Break edges
  unsigned int numBroken = 0;
  for (unsigned int i = 0; i < v; i++)
  {
    const float thres = T * m_NodeAverages[i];
    // Break the coinciding long edges
    for (unsigned int k = 0; k < e; k++)
    {
      if (breakArray[k] != 0)
      {
        continue;
      }

      const unsigned int a = m_MSTEdges[k].i;
      const unsigned int b = m_MSTEdges[k].j;
      const bool         incident = (i == a) || (i == b);

      // Never break zero length edges
      if (incident && (m_MSTEdges[k].dist > thres))
      {
        breakArray[k] = 1;
        ++numBroken;
      }
    }
  }

  if (numBroken == 0)
  {
    std::cerr << "No edges broken" << std::endl;
    delete[] breakArray;
    // INFO: FIXME:
    // return whole tree with same label
    return 0;
  }
  // Figure out distinct trees, merge connected vertices
  for (unsigned int k = 0; k < v; k++)
  {
    treeMap[k] = k;
  }
  for (unsigned int i = 0; i < v; i++)
  {
    unsigned int map1 = treeMap[i];
    // Check incident edges
    for (unsigned int j = 0; j < e; j++)
    {
      if (breakArray[j] != 0)
      {
        continue;
      }

      const unsigned int a = m_MSTEdges[j].i;
      const unsigned int b = m_MSTEdges[j].j;
      const bool         incident = (i == a) || (i == b);

      if (!incident)
      {
        continue;
      }

      // Get the map of the other id
      unsigned int map2 = treeMap[b];
      if (i == b)
      {
        map2 = treeMap[a];
      }
      if (map1 == map2)
      {
        continue;
      }
      for (unsigned int k = 0; k < v; k++)
      {
        if (treeMap[k] == map2)
        {
          treeMap[k] = map1;
        }
      }
    } // for neighbors
  }   // for i

  delete[] breakArray;

  if (!m_SortFlag)
  {
    return numBroken + 1;
  }

  // Sort the cluster maps based on cluster size (descending)
  // Cluster 0 is the largest cluster

  // Obtain cluster info
  MSTCluster * const clusters = new MSTCluster[v];
  unsigned int       numNonZero = 0;
  for (unsigned int i = 0; i < v; i++)
  {
    clusters[i].map = i;
    unsigned int s = 0;
    for (unsigned int j = 0; j < v; j++)
    {
      if (treeMap[j] == i)
      {
        ++s;
      }
    }
    if (s > 0)
    {
      ++numNonZero;
    }
    clusters[i].size = s;
  }

  Heap<MSTCluster> heap;
  heap.Allocate(numNonZero);
  for (unsigned int i = 0; i < v; i++)
  {
    if (clusters[i].size != 0)
    {
      heap.Insert(clusters[i]);
    }
  }

  delete[] clusters;

  unsigned int * const sortedMap = new unsigned int[v];
  for (unsigned int i = 0; i < numNonZero; i++)
  {
    const MSTCluster   c = heap.ExtractMinimum();
    const unsigned int m = c.map;
    for (unsigned int j = 0; j < v; j++)
    {
      if (treeMap[j] == m)
      {
        sortedMap[j] = i;
      }
    }
  }
  for (unsigned int i = 0; i < v; i++)
  {
    treeMap[i] = sortedMap[i];
  }

  delete[] sortedMap;
  return numBroken + 1;
}
