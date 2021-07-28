// File OLDFORMS.H

#if     !defined(_OLDFORMS_H)
#define _OLDFORMS_H      1       //flags that this file has been included

#include "primes.h"

string eigfile(const Quad& N);    //returns filename for eigs at level N
string eigfile(Qideal& N);    //returns filename for eigs at level N

class eigdata {

// Contains eigs info about a sub-level.  The constructor can be used
// to get eigs from a lower level to construct oldforms, and also to
// retrieve data for this level.

public:
  Qideal N;
  Qideal M;
  int nforms,nforms2;
  int nap;
  vector<vector<long> > aqs, aps, eigs;
  vector<vector<int> > intdata;  // sfe, pdot, dp0, cuspidalfactor,
                                 // lambdadot, matdot
  vector<vector<Quad> > Quaddata; // lambda, a, b, c, d
  eigdata(Qideal& iN, Qideal& iM, int neigs=-1, int verbose=0);
};


class oldforms {
 public:
  int noldclasses, nap;
  int olddim1, olddim2, olddimall;  // total dim of rational oldforms, resp. non-rational, all oldform
 private:
  Qideal N;
  vector<Quadprime> plist;
  vector< vector<long> >  oldformap;
  vector<long> oldclassdims;
  vector<Qideal> oldlevels;
  void getoldclasses(Qideal& D, int verbose);
 public:
  oldforms(Qideal& iN, const vector<Quadprime>& pr, int verbose=0);
  long dimoldpart(const vector<long> aplist);
  void display(void) const;
};

#endif
