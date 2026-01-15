// File FIELD.H: class for working with number fields for Hecke eigenvalues
//////////////////////////////////////////////////////////////////////////

#ifndef _FIELD_H
#define _FIELD_H      1

#include "eclib.h"

class Field;
class FieldElement;
class FieldModSq;   // finite subgroups of (F^*)/(F^*)^2
class Eigenvalue;   // a FieldElement a and an index i into FieldModSq representing a*sqrt(elt(i))

extern Field* FieldQQ;

class Field {
  friend class FieldElement;
  friend class Newform;
private:
  string var;   // name of generator
  int d;        // degree
  ZZX minpoly;  // irreducible poly of degree d
  ZZX abspoly;  // irreducible poly of degree d (polredabs of minpoly)
  ZZ denom; // minpoly is the (integral) min poly of A/denom
  mat_m A;        // dxd matrix with scaled min.poly. minpoly
  mat_m B, Binv;  // Binv*B = Bdet*I
  // B has first column Bdet1*e1, Binv has first column Bdet2*e1,
  // Bdet=Bdet1*Bdet2, Bdet3 = denom*Bdet2
  ZZ Bdet, Bdet1, Bdet2, Bdet3, Bcontent;
  // Binv*A*B = Bdet*denom*C
  // cols of Binv are Bdet * coeffs of basis w.r.t. a-powers
  mat_m C;        // dxd companion matrix with min.poly. minpoly
  vector<mat_m> Cpowers;  // C^i for i=0,1,...,d-1
public:
  //Field(const ZZX& p);
  Field(); // defaults to Q
  Field(const mat_m& m, const ZZ& den = to_ZZ(1), string a="a", int verb=0);
  Field(const ZZX& p, string a="a", int verb=0);
  FieldElement rational(const bigrational& x);
  FieldElement rational(const ZZ& x); // {return rational(bigrational(x));}
  FieldElement rational(long x); // {return rational(bigrational(ZZ(x)));}
  FieldElement rational(int x); // {return rational(bigrational(ZZ(x)));}
  FieldElement one(); // {return rational(1);}
  FieldElement minus_one(); // {return rational(-1);}
  FieldElement two(); // {return rational(2);}
  FieldElement minus_two(); // {return rational(-2);}
  FieldElement zero(); // {return rational(0);}
  FieldElement gen();
  FieldElement element(const vec_m& c, const ZZ& d=to_ZZ(1), int raw=0);
  int degree() const {return d;}
  int isQ() const {return this==FieldQQ;}
  ZZX poly() const {return minpoly;}
  ZZX reduced_poly() const {return abspoly;}
  mat_m basis() const {return Binv;} // columns are Bfactor * coeffs of basis w.r.t. a-powers
  ZZ basis_factor() const {return Bdet;}
  mat_m inv_basis() const {return B;} // columns are coeffs of a-powers w.r.t. basis
  void display(ostream&s = cout, int raw=0); // if raw, also display raw basis
  void display_bases(ostream&s = cout) const; // display powers of A and C and bases in both embeddings
  string get_var() const {return var;}
  void set_var(const string& v)  {var = v;}
  // String for pretty output, like "Q" or "Q(i) = Q[X]/(X^2+1)", or
  // (if raw) raw output, suitable for re-input, like "Q" or "i [1 0 1]":
  string str(int raw=0) const;
  friend ostream& operator<<(ostream& s, const Field& F);
  friend istream& operator>>(istream& s, Field** F);
};

class FieldElement {
  friend class Field;
  friend class FieldModSq;
  friend class Eigenvalue;
  friend FieldElement evaluate(const ZZX& f, const FieldElement a);
private:
  Field* F;

  // In general the field element is (1/denom)*coords-combination of power basis of F
  // NB On construction every element will be reduced using cancel()
  vec_m coords; // length F->d
  ZZ denom;     // >=1
  void cancel(); // divides through by gcd(content(coords, denom))

  // When F is Q this is just a wrapper round eclib's bigrational
  // class and coords and denom are ignored
  bigrational val;
public:
  FieldElement()
    :F(FieldQQ) {;}
  FieldElement( Field* HF)
    :F(HF), coords(vec_m(HF->d)), denom(to_ZZ(1))  {if (HF==FieldQQ) val = bigrational(0);}
  // raw means the given coords are w.r.t. the B-basis
  FieldElement( Field* HF, const vec_m& c, const ZZ& d=to_ZZ(1), int raw=0);
  // creation from a rational (general F)
  FieldElement( Field* HF, const ZZ& a, const ZZ& d=to_ZZ(1))
    :F(HF), coords(a*vec_m::unit_vector(HF->d, 1)), denom(d), val(bigrational(a,d)) { cancel();}
  // creation from a rational (F=Q)
  FieldElement( const bigrational& r)
    :F(FieldQQ), val(r) {;}
  // creation from a rational (general F)
  FieldElement( Field* HF, const bigrational& r)
    :F(HF), coords(r.num()*vec_m::unit_vector(HF->d, 1)), denom(r.den()), val(r) {;}


  // String for pretty printing, used in default <<, or (if raw) raw
  // output, suitable for re-input:
  string str(int raw=0) const;

  Field* field() {return F;}
  mat_m matrix() const; // ignores denom
  ZZX charpoly() const;
  ZZX minpoly() const;
  int degree() const {return deg(minpoly());}
  int is_zero() const;
  int is_one() const;
  int is_minus_one() const;
  int is_generator() const {return degree()==F->d;}
  bigrational get_val() const {return val;}
  int operator==(const FieldElement& b) const;
  int operator!=(const FieldElement& b) const;

  FieldElement operator+(const FieldElement& b) const; // add
  FieldElement operator+(const ZZ& b) const {return operator+(FieldElement(F,b));} // add
  void operator+=(const FieldElement& b); // add b to this
  void operator+=(const ZZ& b) { operator+=(FieldElement(F,b));} // add b

  FieldElement operator-(const FieldElement& b) const; // subtract
  FieldElement operator-(const ZZ& b) const {return operator-(FieldElement(F,b));} // subtract
  void operator-=(const FieldElement& b); // subtract b from this
  void operator-=(const ZZ& b) { operator-=(FieldElement(F,b));} // subtract b
  FieldElement operator-() const;                           // unary minus

  FieldElement operator*(const FieldElement& b) const; // product
  FieldElement operator*(const ZZ& b) const {return operator*(FieldElement(F,b));} // product
  void operator*=(const FieldElement& b); // multiply by b
  void operator*=(const ZZ& b) { operator*=(FieldElement(F,b));} // multiply by b

  FieldElement inverse() const; // raise error if zero      // inverse
  FieldElement operator/(const FieldElement& b) const; // divide (raise error if b is zero)
  FieldElement operator/(const ZZ& b) const {return operator/(FieldElement(F,b));} // divide
  void operator/=(const FieldElement& b);                        // divide by b
  void operator/=(const ZZ& b) { operator/=(FieldElement(F,b));} // divide by b

  // NB for a in F, either [Q(sqrt(a))=Q(a)] or [Q(sqrt(a)):Q(a)]=2.
  // The first function only applies when a has maximal degree:
  // return 1 and r s.t. r^2=this, with deg(r)=degree(), else 0
  int is_absolute_square(FieldElement& r)  const;
  // Same as above if the denom is 1
  int is_absolute_integral_square(FieldElement& r)  const;
  // The second function applies in general: return 1 and r
  // s.t. r^2=this, with deg(r)=degree(), else 0. Here, ntries is the
  // number of squares this is multiplied by to get odd co-degree.
  int is_square(FieldElement& r, int ntries=100) const;

  // x must be initialised with a Field before input to x
  friend istream& operator>>(istream& s, FieldElement& x);
};

FieldElement evaluate(const ZZX& f, const FieldElement a);

ostream& operator<<(ostream& s, const FieldElement& x);

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
public:
  FieldModSq() :F(FieldQQ), r(0), elements({F->one()}), real_flag(1) {;}
  FieldModSq(Field* HF) :F(HF), r(0), elements({F->one()}), real_flag(1) {;}
  FieldModSq(Field* HF, vector<FieldElement>& g, int is_real) :F(HF), r(g.size()), gens(g), real_flag(is_real)
  {
    elements = {F->one()};
    for (auto r: gens)
      {
        vector<FieldElement> new_elements(elements.size(), FieldElement(F));
        std::transform(elements.begin(), elements.end(), new_elements.begin(),
                       [r](const FieldElement& x){return r*x;});
        elements.insert(elements.end(), new_elements.begin(), new_elements.end());
      }
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
  void display();
  string str() const;
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
  FieldElement coeff() const {return a;}
  FieldElement root_part() const  { return SqCl->elt(root_index); }
  string extra_factor() const {return (xf>0? "(1+i)" : (xf<0? "(1-i)" : ""));}
  int xfac() const {return xf;}
  // When i=sqrt(-1) is the first element of SqCl normalise using sqrt(-r)*(1+i)=-sqrt(r)*(1-i) and similar
  void normalise();
  Eigenvalue operator*(Eigenvalue b) const;
  Eigenvalue operator/(Eigenvalue b) const;
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

  // as a pretty string, or (if raw) a raw string suitable for re-input:
  string str(int raw=0) const;
};

inline ostream& operator<<(ostream& s, const Eigenvalue& x)
{ s << x.str(); return s;}

// integer multiple of i, assuming not real
Eigenvalue eye(FieldModSq* S, const ZZ& n = ZZ(1));


#endif
