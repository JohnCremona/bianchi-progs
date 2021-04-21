// FILE CUSP.H

#include "moddata.h"
#include "ratquads.h"

// class cusp :public RatQuad {
//  public:
//   cusp(const RatQuad& r) :RatQuad(r) {;}
//   cusp(int n=0, int d=0) :RatQuad(n,d) {;}
//   int eq(const cusp& c) const;
// };

class cusplist {
 private:
    const moddata* N;
    vector<RatQuad> cusps;
 public:
  cusplist(int n=0, const moddata* iN=0)
    :N(iN)
  {;}
  int index(const RatQuad& a); // adds to list if new
  RatQuad item(int n) const {return cusps[n];}
  void display() const;
  int count() const {return cusps.size();}
};
