// LOOPER.H

#if     !defined(_LOOPER_H)
#define _LOOPER_H      1       //flags that this file has been included

#include "quads.h"

class Quadlooper {

public:
  Quadlooper(long nn, long nm, int conj=0)
    :d(Quad::d), disc(Quad::disc), n(nn), nmax(nm), include_conjugates(conj)
    {
      setblims();
    }
  void operator++();
  operator Quad() const {return val;}
  int ok() const {return n<=nmax;}

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

// Lists of elements of norm in ranges (up to units, excluding 0)
vector<Quad> elements_of_norm_between(const INT& n1, const INT& n2);
inline vector<Quad> elements_of_norm(const INT& n) {return elements_of_norm_between(n,n);}
inline vector<Quad> elements_of_norm_up_to(const INT& n) {return elements_of_norm_between(1,n);}

#endif
