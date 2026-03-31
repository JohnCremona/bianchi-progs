// File EIGENVALUE.H: classes for working with for Hecke eigenvalues
////////////////////////////////////////////////////////////////////

#ifndef _EIGENVALUE_H
#define _EIGENVALUE_H      1

#include "eclib.h"

class FieldMQExt;   // finite subgroups of (F^*)/(F^*)^2
class FieldMQElement;   // a FieldElement a and an index i into FieldMQExt representing a*sqrt(elt(i))

// Class for a multiquadratic extensions of a base Field F, defined by
// a pointer to F and a finite subgroups of (F^*)/(F^*)^2

// Contains a list gens of r nonzero elements of F which are
// independent modulo squares.  Elements of the group are represented
// by an unsigned int between 0 and 2^r-1; if the i'th bit is b_i then
// the element is the product of those gens whose index i has b_i=1.

// In the application, F is totally real and the extensions is either
// also t.r. (with real_flag=1) or totally complex (with real_flag=0),
// in which case r is at least 1 and the first gen is -1.

class FieldMQExt {
private:
  const Field* F; // this will be totally real
  unsigned int r;
  vector<FieldElement> gens;
  vector<FieldElement> elements;
  int real_flag;         // 0 iff we have included -1 as first generator, else 1
  void make_elements();  // from r gens make the list of 2^r elements
public:
  // default constructor, a trivial extension of Q
  FieldMQExt() :F(FieldQQ), r(0), elements({F->one()}), real_flag(1) {;}
  // constructor from any field, defining the trivial extension
  explicit FieldMQExt(Field* F0) :F(F0), r(0), elements({F0->one()}), real_flag(1) {;}
  // constructor from any field, given a list of gens (assumed to be
  // multiplicatively independent modulo squares)
  FieldMQExt(Field* F0, vector<FieldElement>& g)
    :F(F0), r(g.size()), gens(g)
  {
    real_flag = !(r>0 && gens[0]==F->minus_one());
    make_elements();
  }
  const Field* base() {return F;}
  int base_degree() const {return F->degree();}
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
  friend istream& operator>>(istream& s, FieldMQExt& x);

  // Return an embedding into an absolute field (optionally
  // polredabs'ed) together with a list of images of the gens.  If the
  // rank is 0 return the identity.
  FieldIso absolute_field_embedding(vector<FieldElement>& im_gens, string newvar, int reduce=0) const;
};

inline ostream& operator<<(ostream& s, const FieldMQExt& x)
{ s << x.str(); return s;}


// An instance of class FieldMQElement consists of:
// a: FieldElement;
// root_index: index i into a SqCl;
// SqCl: a subgroup of F^*/(F^*)^2
// xf: 0,+1,-1
// Where the value is a*sqrt(r)*(1+xf*i) where r is the SqCl->elt(root_index)
class FieldMQElement {
private:
  FieldElement a;
  unsigned int root_index;
  FieldMQExt* SqCl;
  int xf;
public:
  FieldMQElement() {;}
  FieldMQElement(const FieldElement& x, FieldMQExt* S, unsigned int i=0, int f=0)
    : a(x), root_index(i), SqCl(S), xf(f)
  {;}
  FieldMQExt* parent() const {return SqCl;}
  FieldElement base() const {return a;}
  unsigned int index() const {return root_index;}
  FieldElement root_part() const  { return SqCl->elt(root_index); }
  string extra_factor() const {return (xf>0? "(1+i)" : (xf<0? "(1-i)" : ""));}
  int xfac() const {return xf;}
  // When i=sqrt(-1) is the first element of SqCl normalise using sqrt(-r)*(1+i)=-sqrt(r)*(1-i) and similar
  void normalise();
  FieldMQElement operator*(const FieldMQElement& b) const;
  FieldMQElement operator/(const FieldMQElement& b) const;
  FieldMQElement operator-() const {return FieldMQElement(-a, SqCl, root_index, xf);}
  FieldMQElement inverse() const; // raise error if zero      // inverse
  // FieldMQElement times_i() const;
  // FieldMQElement times_minus_i() const;
  FieldMQElement conj() const; // swap 1+i and 1-i factors (complex conjugation)
  int operator==(const FieldMQElement& b) const;
  int operator!=(const FieldMQElement& b) const;
  int is_zero() const {return a.is_zero();}
  int is_one() const {return a.is_one() && root_index==0 && xf==0;}
  int is_minus_one() const {return a.is_minus_one() && root_index==0 && xf==0;}
  bigrational norm() const;
  bigrational trace() const;
  void negate() {a.negate();} // negate in place

  // as a pretty string, or (if raw) a raw string suitable for re-input:
  string str(int raw=0) const;

  // Input (from a raw string format)
  friend istream& operator>>(istream& s, FieldMQElement& x);
};

inline ostream& operator<<(ostream& s, const FieldMQElement& x)
{ s << x.str(); return s;}

// integer multiple of i, assuming not real
FieldMQElement eye(FieldMQExt* S, const ZZ& n = ZZ(1));

// embed an FieldMQElement into the absolute field HFabs, given an
// embedding of F into HFabs and images of the MieldModSq gens in HFabs
FieldElement embed_eigenvalue(const FieldMQElement& ap, const FieldIso& emb, const vector<FieldElement>& im_gens);

#endif
