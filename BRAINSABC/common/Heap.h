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
template <class T>
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
template <class T>
T * heapFirstK(std::vector<T> & array, unsigned int n, unsigned int k);

// Get the k-th element using heap sort
template <class T>
T heapKthElement(std::vector<T> & array, unsigned int n, unsigned int k);

// Get median using heap sort
template <class T>
T heapMedian(std::vector<T> & array, unsigned int n);

#ifndef MU_MANUAL_INSTANTIATION
#include "Heap.hxx"
#endif

#endif
