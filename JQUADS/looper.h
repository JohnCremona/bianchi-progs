// FILE looper.h
// modified by JSB to allow underlying Field to be non-Euclidean,
// and to support options on conjugates and associates

//////////////////////////////////////////////////////////////////////
// Usage:
//     Quad a;
//     for (Quadlooper alpha(firstn,lastn,bothconj,all); alpha.ok(); alpha++)
//         { a=(Quad)alpha;
//           ...
//         }
// for looping through Quads of norm from firstn to lastn.
//
// There are two flags controlling which Quads to return:
//
// bothconj=1, all=1 : absolutely all Quads of the stated norm
// bothconj=0, all=1 : only one from each pair of conjugates
// bothconj=1, all=0 : only one from each set of associates
// bothconj=0, all=0 : only one from each set {associates and conjugates}
//
// The last option is the default, for compatibility with the old
// version of this class, which did not support such options.
//////////////////////////////////////////////////////////////////////

#include "quads.h"

class Quadlooper {

public: 
  Quadlooper(long nn, long nmax, int bothconj=0, int allassocs=0)
    : n(nn), nlim(nmax), both(bothconj), all(allassocs)
    { initbloop();
      while(!finda()) bstep();
      conjflag=0;
      nretd=0;
    }
  void operator++(int);
  operator Quad() const;
  int ok() const {return n<=nlim;}

private: 
  long n,a,b,bsqlim,nlim,nretd;  // nretd is no. of assocs returned so far
  int both, all, conjflag;
  int finda();
  void bstep();
  void initbloop();
};

int issquare(long asq, long& a);

// END OF FILE looper.h
