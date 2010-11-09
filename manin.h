// File MANIN.H

#if     !defined(_MANIN_H)
#define _MANIN_H      1       //flags that this file has been included

#include "oldforms.h"
#include "newforms.h"

class manin :public newforms {
private:
  int easy; 
  Quad nq, dq;
  vector<long> qdotlist;
  vec initvec;
  void findq();    //Computes nq, dq, qdotlist
  void getoneap(const Quad& p, int output, ofstream& out, int verbose=0);
public:
  manin(const Quad& n, int useolddata=0, int disp=0);
  void getap(int first, int last, int output, string eigfile, int verbose=0);
};

#endif
