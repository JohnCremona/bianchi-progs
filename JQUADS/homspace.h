// FILE homspace.h   ---   Declaration of class homspace

#if     !defined(_HOMSPACE_H)
#define _HOMSPACE_H      1       //flags that this file has been included

#include "subspace.h"
#include "symb.h"

class rel_matrix : public mat {
  long numrel, maxnumrel, nc;
public:
  rel_matrix(long r, long c) : mat(r,c), numrel(0), maxnumrel(r), nc(c) {;}
  void userel(vec& rel);
  void trim() { mat* mp=this; *mp = mp->slice(numrel,nc); }
  long get_numrel() const {return numrel;}
};

class homspace :public symbdata {
//friend class newforms;
private:
  void relate_two_term();
  void relate_faces();
  void restrict_to_kerdelta();

  int verbose;
  int GL_vs_SL_flag;
  long *coordindex;
  long* gens;        // needed only to share data between relate_two_term...
  long ngens,rk,denom1,denom2;
  subspace kern;
  n_modsym *freemods;
  int *needed;
//   mat opmat(long,long=0);
//   long matdim(void) {return dimension;} 
//   long matden(void) {return denom3;}
//   longlist eigrange(long i);
protected:
  mat coord, projcoord;   // needed in derived class manin
  long dimension, denom3;
  long ncusps;                // so far, never needed
 
public:
  homspace(const Qideal& n, int hp, int verb);
  ~homspace() { delete[] needed; delete[] freemods; delete[] coordindex;}
  long h1dim() const {return dimension;}  // No confusion with subspace::dim
  long h1denom() const {return denom3;}
  long h1ncusps() const {return ncusps;}

//  int getGL_vs_SL_flag() const {return GL_vs_SL_flag;}
 public:

  vec chain_signed_symb(long k, long eps) const;
  vec reconvert(const n_cusp&) const;
  vec reconvert(const n_modsym&) const;

  mat calc_conj() const;           // only makes sense if level is self-conj
  mat calcop(const matop&) const;

  mat chi_p(const Qideal&p) const;
  mat hecke_p(const Qideal&p) const;
  mat hecke_u_princ(const Qideal&p) const;
  mat hecke_npnp(const Qideal&p) const;
  mat hecke_npnpchi(const Qideal&p, const Qideal&r) const;
  mat hecke_npnq(const Qideal&p, const Qideal&q) const;
  mat hecke_npnqchi(const Qideal&p, const Qideal&q, const Qideal&r) const;

/*
   vec chain(const symb& s) const;
   vec chaincd(const Quad& c, const Quad& d) const;
   vec projchaincd(const Quad& c, const Quad& d) const;
   vec chain(const Quad& nn, const Quad& dd) const;
   vec chain(const RatQuad& r) const 
        {return chain(num(r),den(r));}
   vec kernelpart(const vec& v) const 
        {return v[pivots(kern)];}
   vec cycle(const symb& s) const 
        {return kernelpart(chain(s));}
   vec cycle(const Quad& n, const Quad& d) const 
        {return kernelpart(chain(n,d));}
   vec cycle(const RatQuad& r) const 
        {return kernelpart(chain(num(r),den(r)));}
   vec cycle(const modsym& m) const
     {return cycle(m.beta())-cycle(m.alpha());}
   vec projcycle(const Quad& n, const Quad& d) const;
   vec applyop(const matop& mlist, const RatQuad& q) const;
   vec applyop(const matop& mlist, const modsym& m) const 
                   {return applyop(mlist,m.beta())-applyop(mlist,m.alpha());} 
   mat calcop(char* opname, const Quad& p, const matop& mlist, int display=0) const;
   mat heckeop(const Quad& p, int display=0) const;
   mat wop(const Quad& q, int display=0) const;
   mat fricke(int display=0) const;
//   mat conj(int display=0) const;
   vec maninvector(const Quad& p) const;
   vec projmaninvector(const Quad& p) const;
   vec manintwist(const Quad& lambda, const Quadlist& res, int* chitable) const;
   vec newhecke(const Quad& p, const Quad& n, const Quad& d) const;
*/
};

#endif

// END OF FILE homspace.h
