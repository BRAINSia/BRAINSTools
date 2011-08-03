#include "Heap.h"

#include <iostream>

class WrapInt
{
public:
  WrapInt()
  {
    v = 0;
  }

  WrapInt(int i)
  {
    v = i;
  }

  bool operator<(const WrapInt & w) const
  {
    return this->v < w.v;
  }

  WrapInt & operator=(const WrapInt & w)
  {
    this->v = w.v; return *this;
  }

  int v;
};

int main()
{
  Heap<WrapInt> heap;

  heap.Insert( WrapInt(3) );
  heap.Insert( WrapInt(99) );
  heap.Insert( WrapInt(-1) );
  heap.Insert( WrapInt(7) );
  heap.Insert( WrapInt(21) );
  heap.Insert( WrapInt(-9) );
  heap.Insert( WrapInt(30) );
  heap.Insert( WrapInt(5) );

  std::cout << "Heap size = " << heap.GetNumberOfElements() << std::endl;

  unsigned int origSize = heap.GetNumberOfElements();
  for( unsigned int i = 0; i < origSize; i++ )
    {
    WrapInt mm = heap.ExtractMinimum();
    std::cout << mm.v << " :: ";
    }
  std::cout << std::endl;

  Heap<float> heapf;

  heapf.Insert(9.9);
  heapf.Insert(-2.5);
  heapf.Insert(8.5);
  heapf.Insert(3.0);
  heapf.Insert(7.0);
  heapf.Insert(-1.0);
  heapf.Insert(6.0);

  origSize = heapf.GetNumberOfElements();
  for( unsigned int i = 0; i < origSize; i++ )
    {
    std::cout << heapf.ExtractMinimum() << " :: ";
    }
  std::cout << std::endl;

  double x[] = {7.1, 2.4, -1.2, 12.1, 50.1, 1.9};

  std::cout << "x = ";
  for( unsigned int i = 0; i < 6; i++ )
    {
    std::cout << x[i] << ", ";
    }
  std::cout << std::endl;

  std::cout << "x median = " << heapMedian(x, 6) << std::endl;
  std::cout << "x largest = " << heapKthElement(x, 6, 5) << std::endl;

  double *x3 = heapFirstK(x, 6, 3);
  std::cout << "x3 = ";
  for( unsigned int i = 0; i < 3; i++ )
    {
    std::cout << x3[i] << ", ";
    }
  std::cout << std::endl;

  return 0;
}
