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
// //////////////////////////////////////////////////////////////////////////////
//
// Binary heap using arrays
//
// Element type must have: default constructor, copy ctor, operator=, operator<
//
//
// //////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 9/2003

#ifndef __Heap_h
#define __Heap_h

#include <vector>

/**
 * \class Heap
 */
template <typename T>
class Heap
{
public:

  Heap();
  Heap(const Heap & h);
  ~Heap();

  Heap & operator=(const Heap & h);

  void Allocate(unsigned int size);

  inline void Clear()
  {
    m_Elements.clear();
  }

  T ExtractMinimum();

  inline unsigned int GetNumberOfElements()
  {
    return m_Elements.size();
  }

  void Insert(const T & e);

  bool IsEmpty();

  T * GetElements()
  {
    return m_Elements.GetRawArray();
  }

  void UpdateElementAt(unsigned int i);

private:

  std::vector<T> m_Elements;

  void PreserveHeapOrder();
};

// Get the first k sorted elements using heap sort
template <typename T>
T * heapFirstK(std::vector<T> & array, unsigned int n, unsigned int k);

// Get the k-th element using heap sort
template <typename T>
T heapKthElement(std::vector<T> & array, unsigned int n, unsigned int k);

// Get median using heap sort
template <typename T>
T heapMedian(std::vector<T> & array, unsigned int n);

#ifndef MU_MANUAL_INSTANTIATION
#include "Heap.hxx"
#endif

#endif
