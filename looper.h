#include "quads.h"

class Quadlooper {

public:
  Quadlooper(long nn, long nmax, int conj=0)
    :d(Quad::d), n(nn), nlim(nmax), include_conjugates(conj)
    { setblims();
      b=bmin; while(!finda()) bstep();
    }
  void operator++();
  operator Quad() const {return Quad(a,b);}
  int ok() const {return n<=nlim;}

private:
  long d,n,a,b,bmin,bmax,nlim;
  int include_conjugates;
  int finda();
  void nstep();
  void bstep();
  void setblims();
};

int issquare(long asq, long& a);
