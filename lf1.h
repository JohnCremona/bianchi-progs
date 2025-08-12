// lf1.h : class for integrating newforms

#if     !defined(_LF1_H)
#define _LF1_H      1       //flags that this file has been included

#include "newforms.h"

// Class to compute L(F,1) or L(F,chi,1) or L'(F,1) from a newform

class lf1 {
private:
  Qideal N;
  int debug;
  long limitnorm;
  double factor, sum;
  Quad lambda;
  vector<Quad> lambdares;
  vector<int> chitable;
  int chi_is_trivial;
  vector<long> aplist;
  rational loverp;
  int sfe, ar;
  double lf1chi, ldash1, period;

  int chi(const Quad& n);
  double K(double x);
  void use(const Quad& n, int an);
  void add(const Quad& n, int pindex, int y, int z);

public:
  explicit lf1(newform* f, int r, int db=0);
  double get_lf1chi()  {return lf1chi;}
  double get_period()  {return period;}
  double get_ldash1()  {return ldash1;}
};

// Class to compute the period of F, by computing the integral of F
// from a to M(a) where M=[a,b;nu*c,d] in Gamma_0(nu) is a matrix
// chosen in the newform constructor to be such that this integral is
// a nonzero multiple of the period

class period_direct {
private:
  Qideal N;
  Quad nu;
  int debug;
  long limitnorm, maxnormp;
  bigcomplex eta, z1, z2;
  double rootdisc, factor, sum;
  vector<long> aplist;
  Quad fa, fb, fc, fd; // matrix stored with f
  long period_multiple;
  double psi_factor(const Quad& n);
  void use(const Quad& n, int an);
  void add(const Quad& n, int pindex, int y, int z);

public:
  period_direct (newform* f, int db=0);
  // Compute the period along {.,g(.)} for g in Gamma_0(N)
  double compute_period(const Quad& a, const Quad& b, const Quad& c, const Quad& d);
  // Compute the period along {.,g(.)} for f's stored matrix g and divide by period_multiple
  double compute_base_period();
};

// absolute value of Quad x
double realnorm(const Quad& z);

// convert to a complex number
bigcomplex to_bigcomplex(const Quad& a);

#endif
