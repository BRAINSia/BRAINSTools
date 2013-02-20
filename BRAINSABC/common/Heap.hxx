#ifndef __Heap_hxx
#define __Heap_hxx

#include "Heap.h"

#include "muException.h"

// Macros to convert binary tree traversal to array access
// Tree root is 1
#define HEAP_PARENT(x) ( ( x ) / 2 )
#define HEAP_LEFT(x) ( 2 * ( x ) )
#define HEAP_RIGHT(x) ( 2 * ( x ) + 1 )

template <class T>
Heap<T>
::Heap()
{
}

template <class T>
Heap<T>
::Heap(const Heap<T> & h) :
  m_Elements(h.m_Elements)
{
}

template <class T>
Heap<T>
::~Heap()
{
}

template <class T>
Heap<T> &
Heap<T>
::operator =( const Heap & h )
{
  m_Elements = h.m_Elements;
}

template <class T>
void
Heap<T>
::Allocate(unsigned int size)
{
  m_Elements.resize(size);
}

template <class T>
T
Heap<T>
::ExtractMinimum()
{
  if( this->IsEmpty() )
    {
    muExceptionMacro(<< "[Heap::ExtractMinimum] Heap is empty");
    }

  T minElem = m_Elements[0];

  unsigned int last = m_Elements.size() - 1;

  m_Elements[0] = m_Elements[last];
  m_Elements.pop_back();

  this->PreserveHeapOrder();

  return minElem;
}

template <class T>
bool
Heap<T>
::IsEmpty()
{
  return m_Elements.empty();
}

template <class T>
void
Heap<T>
::Insert(const T & e)
{
  m_Elements.push_back(e);

  unsigned int i = m_Elements.size();

  while( ( i > 1 ) && ( e < m_Elements[HEAP_PARENT(i) - 1] ) )
    {
    m_Elements[i - 1] = m_Elements[HEAP_PARENT(i) - 1];
    i = HEAP_PARENT(i);
    }

  m_Elements[i - 1] = e;
}

template <class T>
void
Heap<T>
::PreserveHeapOrder()
{
  unsigned int numElements = m_Elements.size();

  unsigned int current;
  unsigned int left;
  unsigned int right;
  unsigned int smallest;

  // Start at root node
  current = 1;

  do
    {
    left = HEAP_LEFT(current);
    right = HEAP_RIGHT(current);

    smallest = current;

    if( ( left <= numElements )
        &&
        ( m_Elements[left - 1] < m_Elements[smallest - 1] ) )
      {
      smallest = left;
      }

    if( ( right <= numElements )
        &&
        ( m_Elements[right - 1] < m_Elements[smallest - 1] ) )
      {
      smallest = right;
      }

    if( smallest == current )
      {
      break;
      }

    T temp = m_Elements[current - 1];
    m_Elements[current - 1] = m_Elements[smallest - 1];
    m_Elements[smallest - 1] = temp;

    current = smallest;
    }
  while( true );
}

template <class T>
void
Heap<T>
::UpdateElementAt(unsigned int i)
{
  unsigned int loc = i + 1;

  // Root update
  if( loc == 1 )
    {
    this->PreserveHeapOrder();
    }

  //
  // Find new location while swapping values between parent-child
  //

  // Loop until we reach root or until the element is in the right order
  // relative to parent
  while(
    ( loc > 1 )
    &&
    ( m_Elements[loc - 1] < m_Elements[HEAP_PARENT(loc) - 1] ) )
    {
    // Swap parent and current
    T t = m_Elements[HEAP_PARENT(loc) - 1];
    m_Elements[HEAP_PARENT(loc) - 1] = m_Elements[loc - 1];
    m_Elements[loc - 1] = t;

    // Process parent at next iteration
    loc = HEAP_PARENT(loc);
    }
}

template <class T>
T *
heapFirstK(std::vector<T> & array, unsigned int n, unsigned int k)
{
  if( k >= n )
    {
    return 0;
    }
  T *firstk = new T[k];

  Heap<T> heap;
  heap.Allocate(n);
  for( unsigned int i = 0; i < n; i++ )
    {
    heap.Insert(array[i]);
    }
  for( unsigned int i = 0; i < k; i++ )
    {
    firstk[i] = heap.ExtractMinimum();
    }

  return firstk;
}

template <class T>
T
heapKthElement(std::vector<T> & array, unsigned int n, unsigned int k)
{
  if( k >= n )
    {
    return array[n - 1];
    }

  Heap<T> heap;
  heap.Allocate(n);
  for( unsigned int i = 0; i < n; i++ )
    {
    heap.Insert(array[i]);
    }
  // Throw away first k-1 values
  for( unsigned int i = 0; i < k; i++ )
    {
    heap.ExtractMinimum();
    }

  return heap.ExtractMinimum();
}

template <class T>
T
heapMedian(std::vector<T> & array, unsigned int n)
{
  if( n == 0 )
    {
    return 0;
    }

  if( n == 1 )
    {
    return array[0];
    }

  if( n == 2 )
    {
    return ( array[0] + array[1] ) / 2;
    }

  if( ( n % 2 ) == 0 )
    {
    unsigned int k = ( n / 2 ) + 1;

    T *tmp = heapFirstK(array, n, k);
    T  mid = ( tmp[k - 2] + tmp[k - 1] ) / 2;

    delete[] tmp;

    return mid;
    }
  else
    {
    return heapKthElement(array, n, ( n + 1 ) / 2);
    }
}

#endif
