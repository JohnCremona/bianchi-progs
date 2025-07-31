// lf1.h : class for integrating newforms

#if     !defined(_LF1_H)
#define _LF1_H      1       //flags that this file has been included

#include "newforms.h"

inline double realnorm(const Quad& z) {  return sqrt(to_double(z.norm()));}
// inline bigfloat psif(bigcomplex z) {  return cos(4*PI*real(z));}
// inline bigfloat psig(bigcomplex z) {  return sin(4*PI*real(z));}


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
  double getlf1chivalue()
  {return lf1chivalue;}
  double getperiod()
  {return period;}
  rational getratio()
  {return ratio;}
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
