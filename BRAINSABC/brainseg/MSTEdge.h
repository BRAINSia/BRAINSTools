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
#ifndef __MSTEdge_h
#define __MSTEdge_h

/** \class MSTEdge
 */
class MSTEdge
{
public:
  // Pair of graph vertex indices
  unsigned int i;
  unsigned int j;

  // Distance between the vertices
  float dist;

  MSTEdge()
  {
    this->i = 0;
    this->j = 0;
    this->dist = 0;
  }

  ~MSTEdge() {}

  MSTEdge &
  operator=( const MSTEdge & e )
  {
    this->i = e.i;
    this->j = e.j;
    this->dist = e.dist;
    return *this;
  }

  bool
  operator<( const MSTEdge & e ) const
  {
    return this->dist < e.dist;
  }
};

#endif
