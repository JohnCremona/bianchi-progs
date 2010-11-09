#include "quads.h"

class Quadlooper {

public: 
  Quadlooper(int dd, int nn, int nmax) :d(dd), n(nn), nlim(nmax)
    { setbsqlim();
      b=0; while(!finda()) bstep();
    }
  void operator++();
  operator Quad() const {return Quad(a,b);}
  int ok() const {return n<=nlim;}

private: 
  int d,n,a,b,bsqlim,nlim;
  int finda();
  void nstep();
  void bstep();
  void setbsqlim();
};

int issquare(int asq, int& a);
