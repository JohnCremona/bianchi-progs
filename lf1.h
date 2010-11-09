// lf1.h : class for integrating newforms

#if     !defined(_LF1_H)
#define _LF1_H      1       //flags that this file has been included

#include "newforms.h"

class period_via_lf1chi {
private:
  Quad N; int debug;
  long limitnorm;
  double factor, sum;
  Quad lambda; Quadlist lambdares;
  int *chitable; // created and deleted in the constructor
  int chi(const Quad& n) {return chitable[lambdares.locate(n%lambda)];}
  int nap;  vector<long> aplist;
  rational ratio;
  double lf1chivalue, period;

  void use(const Quad& n, int an);
  void add(const Quad& n, int pindex, int y, int z);

public:
  period_via_lf1chi (newform* f, int db=0); 
  double getlf1chivalue() {return lf1chivalue;}
  double getperiod() {return period;}
  rational getratio(){/* cout<<"<getratio: "<<ratio<<">";*/ return ratio;}
};

/*
class periods_direct {
private:
  int N,limit;
  double theta1,theta2,efactor;
  Complex sum;
  int nap;  longlist aplist;  longlist primelist;
  Complex *periods;

  void use(int n, int an);
  void add(int n, int pindex, int y, int z);

public:
  periods_direct (h1newform* f); 
  ~periods_direct () {delete periods;}

  Complex* getperiods() {return periods;}
};

Complex epi(double);
*/

#endif
