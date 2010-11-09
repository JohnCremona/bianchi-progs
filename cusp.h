// FILE CUSP.H

#include "moddata.h"
#include "ratquads.h"

class cusp :public RatQuad {
 public:
  cusp(const RatQuad& r) :RatQuad(r) {;}
  cusp(int n=0, int d=0) :RatQuad(n,d) {;}
  int eq(const cusp& c) const;
};

class cusplist {
 private:
    cusp *list;
    int num,maxnum;
 public:
    cusplist(int n=0) {maxnum=n; num=0; list=new cusp[n];}
    ~cusplist() {delete[] list;}
    int index(const cusp& a);
    cusp item(int n) const {return list[n];}  //should check n really
    void display() const {for(int i=0; i<num; i++) cout<<i<<"\t"<<list[i]<<endl;}
    int count() const {return num;}
};   
