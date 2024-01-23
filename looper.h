// LOOPER.H

#if     !defined(_LOOPER_H)
#define _LOOPER_H      1       //flags that this file has been included

#include "quads.h"

class Quadlooper {

public:
  // Iterator through Quads (up to units) with norms from nn to nm (or indefinite if nm==0)
  Quadlooper(long nn, long nm, int conj=0)
    :d(Quad::d), disc(Quad::disc), n(nn), nmax(nm), include_conjugates(conj)
    {
      setblims();
    }
  void operator++();
  operator Quad() const {return val;}
  vector<Quad> values_with_current_norm(int sorted=0);
  vector<Quad> values_with_norm_up_to(const INT& m, int sorted=0);
  int ok() const {return nmax==0 || n<=nmax;}

private:
  long d;
  INT disc,n,n4,b,db2,bmin,bmax,nmax;
  Quad val;
  int include_conjugates;
  int testb();
  void nstep();
  void bstep();
  void setblims();
};

// Lists of quads of norm in ranges (up to units, omitting conjugates if conj==0)
vector<Quad> quads_of_norm_between(const INT& n1, const INT& n2, int conj=1, int sorted=0);
vector<Quad> quads_of_norm(const INT& n, int conj=1, int sorted=0);
inline vector<Quad> quads_of_norm_up_to(const INT& n, int conj=1, int sorted=0)
{
  return quads_of_norm_between(1,n,conj,sorted);
}

#endif
