// lf1.h : class for integrating newforms

#if     !defined(_LF1_H)
#define _LF1_H      1       //flags that this file has been included

#include "newforms.h"

// Class to compute L(F,1) from a newform

class period_via_lf1chi {
private:
  Qideal N;
  int debug;
  long limitnorm;
  double factor, sum;
  Quad lambda;
  vector<Quad> lambdares;
  vector<int> chitable;
  int chi(const Quad& n)
  {
    return chitable[std::distance(lambdares.begin(), std::find(lambdares.begin(), lambdares.end(), n%lambda))];
  }
  int nap;
  vector<long> aplist;
  rational ratio;
  double lf1chivalue, period;

  void use(const Quad& n, int an);
  void add(const Quad& n, int pindex, int y, int z);

public:
  explicit period_via_lf1chi (newform* f, int db=0);
  double get_lf1chivalue()
  {return lf1chivalue;}
  double get_period()
  {return period;}
  rational get_ratio()
  {return ratio;}
};

// Class to compute the period of F, by computing the integral of F
// from a to M(a) where M=[a,b;nu*c,d] in Gamma_0(nu) is a matrix
// chosen in the newform constructor to be such that this integral is
// a nonzero multiple of the period

class period_direct {
private:
  Qideal N;
  int debug;
  long limitnorm;
  bigcomplex theta1, theta2;
  double factor, sum, period;
  int nap;
  vector<long> aplist;
  double psi_factor(const Quad& n);
  void use(const Quad& n, int an);
  void add(const Quad& n, int pindex, int y, int z);

public:
  period_direct (newform* f, int db=0);

  double get_period() {return period;}
};

#endif
