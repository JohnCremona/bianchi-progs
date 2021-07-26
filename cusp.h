// FILE CUSP.H

#include "moddata.h"
#include "ratquads.h"

class cusplist {
 private:
  Quad N;
  int plusflag;
  vector<RatQuad> cusps;
 public:
  cusplist(const Quad& level, int plus)
    :N(level), plusflag(plus)
  {;}
  int index(const RatQuad& a); // adds to list if new
  RatQuad item(int n) const {return cusps[n];}
  void display() const;
  int count() const {return cusps.size();}
};
