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
    RatQuad *list;
    int number,maxnumber;
    int cuspeq(const RatQuad& c1, const RatQuad& c2) const;
 public:
  cusplist(int n=0, const moddata* iN=0)
    :N(iN), number(0), maxnumber(n)
  {list=new RatQuad[n];}
    ~cusplist() {delete[] list;}
    int index(const RatQuad& a);
    RatQuad item(int n) const {return list[n];}  //should check n really
    void display() const {for(int i=0; i<number; i++) cout<<i<<"\t"<<list[i]<<endl;}
    int count() const {return number;}
};
