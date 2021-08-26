// FILE CUSP.H

#if     !defined(_CUSP_H)
#define _CUSP_H      1       //flags that this file has been included

#include "mat22.h"
#include "ratquads.h"
#include "qideal.h"

class cusplist {
 private:
  Qideal N;
  int plusflag;
  vector<RatQuad> cusps;
 public:
  cusplist(const Qideal& level, int plus)
    :N(level), plusflag(plus)
  {;}
  int index(const RatQuad& a); // adds to list if new
  RatQuad item(int n) const {return cusps[n];}
  void display() const;
  int count() const {return cusps.size();}
};

#endif
