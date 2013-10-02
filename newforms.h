// File NEWFORMS.H

#if     !defined(_NEWFORMS_H)
#define _NEWFORMS_H      1       //flags that this file has been included

#include <eclib/xsplit.h>   // which includes method.h
#include <eclib/rat.h>
#include "oldforms.h"
#include "homspace.h"

class newforms;

class newform {
friend class newforms;
public:
  vec basis; 
  vector<long> aplist; 
  Quad lambda;  // twisting prime
//Quad a,b,c,d,dot;
  rational loverp;
  int dp0, pdot, sfe;            // sign of F.E.
  newform(void) :basis(0), aplist(0) {;}
  newform(const newforms* nfs, const vec& v, const vector<long>& ap);
  newform(const newform& nf) 
    :basis(nf.basis),aplist(nf.aplist),
     lambda(nf.lambda), loverp(nf.loverp),
     dp0(nf.dp0),pdot(nf.pdot),sfe(nf.sfe) {;}
  void operator=(const newform& nf) 
    {
      basis =nf.basis;    aplist   =nf.aplist; 
      lambda=nf.lambda;
      loverp=nf.loverp;   sfe      =nf.sfe;
      pdot  =nf.pdot;     dp0      =nf.dp0;
    }
  void display(void) const;
};

class newforms :public level, splitter_base {
friend class newform;
private:
  int dimsplit, maxdepth, upperbound;
  mat opmat(int i, int d, int v=0) {return h1->opmat(i,d,v);}
  mat opmat_restricted(int i, const subspace& s, int d, int v=0) 
  {return h1->opmat_restricted(i,s,d,v);}
  smat s_opmat(int i, int d, int v=0) 
  {return h1->s_opmat(i,d,v);}
  smat s_opmat_restricted(int i, const ssubspace& s, int d, int v=0) 
  {return h1->s_opmat_restricted(i,s,d,v);}
  long matdim(void) {return h1->dimension;} 
  long matden(void) {return h1->denom3;}
  vector<long> eigrange(int i) {return h1->eigrange(i);}
  long dimoldpart(const vector<long> l) {return of->dimoldpart(l);}
protected:
  oldforms* of;
  Quad p0; vec mvp;
public:
  int verbose, n1ds,n2ds, nnflist, nap, ntp, nwq;
  homspace* h1;
  long hmod;
  vector<newform> nflist;
  newforms(const Quad& n, int useolddata=0, int disp=0);
  ~newforms(void) {
                   if(h1)delete h1;
                  }
  //NB tpmats,tpknown,wmats are created and deleted in the constructor
  void display(void) const;
  vec proj(const vec&);   //returns vec of components in each eig-space
  void allproj(void);  //Replaces "coord" member of homspace with projections
                       //onto eigenspaces, to save time
  void makeh1plus(void);
  void use(const vec& b1, const vec& b2, const vector<long> l); 
// add newform with basis b1, eiglist l to current list (b2 not used)
 private:
  void createfromscratch();
  void createfromeigs();
};

#endif
