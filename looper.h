// LOOPER.H

#if     !defined(_LOOPER_H)
#define _LOOPER_H      1       //flags that this file has been included

#include "quads.h"

class Quadlooper {

public:
  Quadlooper(long nn, long nmax, int conj=0)
    :d(Quad::d), disc(Quad::disc), n(nn), nlim(nmax), include_conjugates(conj)
    { setblims();
      b=bmin; while(!finda()) bstep();
    }
  void operator++();
  operator Quad() const {return Quad(a,b);}
  int ok() const {return n<=nlim;}

private:
  long d;
  INT disc,n,a,b,bmin,bmax,nlim;
  int include_conjugates;
  int finda();
  void nstep();
  void bstep();
  void setblims();
};

#endif
