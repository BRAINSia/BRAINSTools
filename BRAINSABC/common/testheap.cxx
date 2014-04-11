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
