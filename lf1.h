// lf1.h : class for integrating newforms

#if     !defined(_LF1_H)
#define _LF1_H      1       //flags that this file has been included

#include "newforms.h"
#include "real.h"

// Class to compute L(F,1) or L(F,chi,1) or L'(F,1) from a newform

class lf1 {
private:
  Qideal N;
  int debug;
  long limitnorm;
  REAL factor, sum;
  Quad lambda;
  vector<Quad> lambdares;
  vector<int> chitable;
  int chi_is_trivial;
  vector<long> aplist;
  rational loverp;
  int sfe, ar;
  REAL lf1chi, ldash1, period;

  int chi(const Quad& n);
  REAL K(REAL x);
  void use(const Quad& n, int an);
  void add(const Quad& n, int pindex, int y, int z);

public:
  explicit lf1(newform* f, int r, int db=0);
  REAL get_lf1chi()  {return lf1chi;}
  REAL get_period()  {return period;}
  REAL get_ldash1()  {return ldash1;}
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
  REAL factor, sum;
  vector<long> aplist;
  Quad fa, fb, fc, fd; // matrix stored with f
  Quad Ma, Mb, Mc, Md; // matrix of current integration
  long period_multiple;
  REAL psi_factor(const Quad& alpha);
  void use(const Quad& n, int an);
  void add(const Quad& n, int pindex, int y, int z);

public:
  period_direct (newform* f, int db=0);
  // Compute the period along {.,g(.)} for g in Gamma_0(N)
  REAL compute_period(const Quad& a, const Quad& b, const Quad& c, const Quad& d);
  // Compute the period along {.,g(.)} for f's stored matrix g and divide by period_multiple
  REAL compute_base_period();
};

// sqrt(Quad::d)
inline REAL rootd()
{
  static REAL rd = sqrt(REAL(Quad::d));
  return rd;
}
// sqrt(Quad::absdisc)
inline REAL rootdisc()
{
  static REAL rD = sqrt(Quad::absdisc);
  return rD;
}
// real part of Quad x
REAL real_part(const Quad& z);
// imag part of Quad x
REAL imag_part(const Quad& z, int norootd=0);
// absolute value of Quad x
REAL realnorm(const Quad& z);

#endif
