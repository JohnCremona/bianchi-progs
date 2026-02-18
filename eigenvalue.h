// File EIGENVALUE.H: classes for working with for Hecke eigenvalues
////////////////////////////////////////////////////////////////////

#ifndef _EIGENVALUE_H
#define _EIGENVALUE_H      1

#include "field.h"

class FieldModSq;   // finite subgroups of (F^*)/(F^*)^2
class Eigenvalue;   // a FieldElement a and an index i into FieldModSq representing a*sqrt(elt(i))

// Class to handle finite subgroups of (F^*)/(F^*)^2

// Contains a list gens of length r of nonzero elements of F which are
// independent modulo squares.  Elements of the group are represented
// by an unsigned int between 0 and 2^r-1; if the i'th bit is b_i then
// the element is the product of those gens whose index i has b_i=1.

class FieldModSq {
private:
  Field* F; // this will be totally real
  unsigned int r;
  vector<FieldElement> gens;
  vector<FieldElement> elements;
  int real_flag; // 1 unless we have included -1 as first generator
  // from r gens make the list of 2^r elements
  void make_elements();
public:
  FieldModSq() :F(FieldQQ), r(0), elements({F->one()}), real_flag(1) {;}
  explicit FieldModSq(Field* HF) :F(HF), r(0), elements({F->one()}), real_flag(1) {;}
  FieldModSq(Field* HF, vector<FieldElement>& g)
    :F(HF), r(g.size()), gens(g)
  {
    real_flag = !(r>0 && gens[0]==F->minus_one());
    make_elements();
  }
  Field* field() {return F;}
  FieldElement gen(unsigned int i) const {return gens.at(i);}
  FieldElement elt(unsigned int i) const {return elements.at(i);}
  vector<FieldElement> elts() const {return elements;}
  int is_real() const {return real_flag;}
  int is_complex() const {return !real_flag;}

  // Compute the index of a nonzero element.

  // If a belongs to the current group return i and set s, where a =
  // elements[i]*s^2.

  // If a does not belong to the subgroup (mod squares):
  //   if update (default):
  //      append a to gens, increment r, set s=1 return the new r;
  //   else:
  //      do not change the group, return -1.
  unsigned int get_index(const FieldElement& a, FieldElement& s, int update=1);
  string elt_str(unsigned int i) const;
  unsigned int rank() const {return r;}
  int order() const {return elements.size();}
  // string representation, either raw (suitable for reinput using >>)
  // or more readable for display:
  string str(int raw=0) const;
  // x must be initialised with field F, this inputs rank and gens and
  // sets real_flag and elements:
  friend istream& operator>>(istream& s, FieldModSq& x);

  // Return an embedding into an absolute field (optionally
  // polredabs'ed) together with a list of images of the gens.  If the
  // rank is 0 return the identity.
  FieldIso absolute_field_embedding(vector<FieldElement>& im_gens, string newvar, int reduce=0) const;
};

inline ostream& operator<<(ostream& s, const FieldModSq& x)
{ s << x.str(); return s;}


// An instance of class Eigenvalue consists of:
// a: FieldElement;
// root_index: index i into a SqCl;
// SqCl: a subgroup of F^*/(F^*)^2
// xf: 0,+1,-1
// Where the value is a*sqrt(r)*(1+xf*i) where r is the SqCl->elt(root_index)
class Eigenvalue {
private:
  FieldElement a;
  unsigned int root_index;
  FieldModSq* SqCl;
  int xf;
public:
  Eigenvalue() {;}
  Eigenvalue(const FieldElement& x, FieldModSq* S, unsigned int i=0, int f=0)
    : a(x), root_index(i), SqCl(S), xf(f)
  {;}
  FieldModSq* parent() const {return SqCl;}
  FieldElement base() const {return a;}
  unsigned int index() const {return root_index;}
  FieldElement root_part() const  { return SqCl->elt(root_index); }
  string extra_factor() const {return (xf>0? "(1+i)" : (xf<0? "(1-i)" : ""));}
  int xfac() const {return xf;}
  // When i=sqrt(-1) is the first element of SqCl normalise using sqrt(-r)*(1+i)=-sqrt(r)*(1-i) and similar
  void normalise();
  Eigenvalue operator*(const Eigenvalue& b) const;
  Eigenvalue operator/(const Eigenvalue& b) const;
  Eigenvalue operator-() const {return Eigenvalue(-a, SqCl, root_index, xf);}
  Eigenvalue inverse() const; // raise error if zero      // inverse
  // Eigenvalue times_i() const;
  // Eigenvalue times_minus_i() const;
  Eigenvalue conj() const; // swap 1+i and 1-i factors (complex conjugation)
  int operator==(const Eigenvalue& b) const;
  int operator!=(const Eigenvalue& b) const;
  int is_zero() const {return a.is_zero();}
  int is_one() const {return a.is_one() && root_index==0 && xf==0;}
  int is_minus_one() const {return a.is_minus_one() && root_index==0 && xf==0;}
  bigrational norm() const;
  bigrational trace() const;
  void negate() {a.negate();} // negate in place

  // as a pretty string, or (if raw) a raw string suitable for re-input:
  string str(int raw=0) const;

  // Input (from a raw string format)
  friend istream& operator>>(istream& s, Eigenvalue& x);
};

inline ostream& operator<<(ostream& s, const Eigenvalue& x)
{ s << x.str(); return s;}

// integer multiple of i, assuming not real
Eigenvalue eye(FieldModSq* S, const ZZ& n = ZZ(1));

// embed an Eigenvalue into the absolute field Fabs, given an
// embedding of F into Fabs and images of the MieldModSq gens in Fabs
FieldElement embed_eigenvalue(const Eigenvalue& ap, const FieldIso& emb, const vector<FieldElement>& im_gens);

#endif
