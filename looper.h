// LOOPER.H

#if     !defined(_LOOPER_H)
#define _LOOPER_H      1       //flags that this file has been included

#include "quads.h"

class Quadlooper {

public:
  // Iterator through Quads (up to units) with norms from n1 to n2 (or indefinite if n2==0)
  Quadlooper(long n1=1, long n2=0, int conj=1)
    :d(Quad::d), disc(Quad::disc), n(n1), nmax(n2), include_conjugates(conj)
    {
      setblims();
    }
  void operator++();
  operator Quad() const {return val;}
  vector<Quad> values_with_current_norm(int sorted=0);
  vector<Quad> values_with_norm_up_to(const INT& m, int sorted=0);
  int ok() const {return nmax.is_zero() || n<=nmax;}

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
  return quads_of_norm_between(ONE,n,conj,sorted);
}

#endif
