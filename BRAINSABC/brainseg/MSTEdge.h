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
    this->i = 0; this->j = 0; this->dist = 0;
  }

  ~MSTEdge()
  {
  }

  const MSTEdge & operator=(const MSTEdge & e)
  {
    this->i = e.i;
    this->j = e.j;
    this->dist = e.dist;
    return *this;
  }

  bool operator<(const MSTEdge & e) const
  {
    return this->dist < e.dist;
  }
};

#endif
