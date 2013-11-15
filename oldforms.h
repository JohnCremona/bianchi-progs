// File OLDFORMS.H

#if     !defined(_OLDFORMS_H)
#define _OLDFORMS_H      1       //flags that this file has been included

#include "moddata.h"

class eigdata {

// Contains eigs info about a sub-level.  The constructor can be used
// to get eigs from a lower level to construct oldforms, and also to
// retrieve data for this level.

public:
  const level *N;
  Quad sublevel;
  int nforms,nforms2;
  int nap;
  vector<vector<long> > aqs, aps, eigs;
  vector<vector<int> > intdata;  // sfe, pdot, dp0, cuspidalfactor,
                                 // lambdadot, matdot
  vector<vector<Quad> > Quaddata; // lambda, a, b, c, d
  eigdata(const level *iN, const Quad& m, int neigs=-1, int verbose=0);
};


class oldforms {
 public:
  int noldclasses, nap, ntp;
  int olddim1, olddim2, olddimall;  // total dim of rational oldforms, resp. non-rational, all oldform
 private:
  const level* N;
  vector< vector<long> >  oldformap;
  vector<long> oldclassdims;
  vector<Quad> oldlevels;
  void getoldclasses(const Quad& d, int verbose);
 public:
  oldforms(const level* iN, int verbose=0);
  long dimoldpart(const vector<long> aplist);
  void display(void) const;
};

#endif
