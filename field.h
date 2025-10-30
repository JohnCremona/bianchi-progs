// File FIELD.H: class for working with number fields for Hecke eigenvalues
//////////////////////////////////////////////////////////////////////////

#ifndef _FIELD_H
#define _FIELD_H      1

#include "eclib.h"

class Field;
class FieldElement;
class FieldModSq;   // finite subgroups of (F^*)/(F^*)^2
class Eigenvalue;   // a FieldElement a and an index i into FieldModSq representing a*sqrt(elt(i))

class Field {
  friend class FieldElement;
  friend class Newform;
private:
  string var;   // name of generator
  int d;        // degree
  ZZX minpoly;  // irredicible poly of degree d
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
  FieldElement one();
  FieldElement zero();
  FieldElement gen();
  FieldElement element(const vec_m& c, const ZZ& d=to_ZZ(1), int raw=0);
  int degree() const {return d;}
  ZZX poly() const {return minpoly;}
  mat_m basis() const {return Binv;} // columns are Bfactor * coeffs of basis w.r.t. a-powers
  ZZ basis_factor() const {return Bdet;}
  mat_m inv_basis() const {return B;} // columns are coeffs of a-powers w.r.t. basis
  void display(ostream&s = cout, int raw=0); // if raw, also display raw basis
  void display_bases(ostream&s = cout) const; // display powers of A and C and bases in both embeddings
};

extern Field* FieldQQ;

class FieldElement {
  friend class Field;
  friend class FieldModSq;
  friend class Eigenvalue;
  friend FieldElement evaluate(const ZZX& f, const FieldElement a);
private:
  Field* F;
  vec_m coords; // length F->d
  ZZ denom;     // >=1
  // the field element is (1/denom)*coords-combination of power basis of F
  // NB On construction every element will be reduced using cancel()
  void cancel(); // divides through by gcd(content(coords, denom))
public:
  FieldElement()
    :F(FieldQQ) {;}
  FieldElement( Field* HF)
    :F(HF), coords(vec_m(HF->d)), denom(to_ZZ(1))  {;}
  // raw means the given coords are w.r.t. the B-basis
  FieldElement( Field* HF, const vec_m& c, const ZZ& d=to_ZZ(1), int raw=0);
  // creation from a rational
  FieldElement( Field* HF, const ZZ& a, const ZZ& d=to_ZZ(1));

  string str() const;
  Field* field() {return F;}
  mat_m matrix() const; // ignores denom
  ZZX charpoly() const {return ::charpoly(matrix());}
  ZZX minpoly() const;
  int degree() const {return deg(minpoly());}
  int is_zero() const;
  int is_one() const;
  int is_minus_one() const;
  int is_generator() const {return degree()==F->d;}
  int operator==(const FieldElement& b) const;

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
  // Same as above if the min poly is known
  int is_absolute_square(FieldElement& r, const ZZX& minpol)  const;
  // The second function applies in general: return 1 and r
  // s.t. r^2=this, with deg(r)=degree(), else 0. Here, ntries is the
  // number of squares this is multiplied by to get odd co-degree.
  int is_square(FieldElement& r, int ntries=100) const;
};

FieldElement evaluate(const ZZX& f, const FieldElement a);

inline ostream& operator<<(ostream& s, const FieldElement& x)
{ s << x.str(); return s;}

// Class to handle finite subgroups of (F^*)/(F^*)^2

// Contains a list gens of length r of nonzero elements of F which are
// independent modulo squares.  Elements of the group are represented
// by an unsigned int between 0 and 2^r-1; if the i'th bit is b_i then
// the element is the product of those gens whose index i has b_i=1.

class FieldModSq {
private:
  Field* F;
  unsigned int r;
  vector<FieldElement> gens;
  vector<FieldElement> elements;
public:
  FieldModSq(){;}
  FieldModSq(Field* HF) :F(HF), r(0), elements({F->one()}) {;}
  FieldElement gen(unsigned int i) const {return gens.at(i);}
  FieldElement elt(unsigned int i) const {return elements.at(i);}
  vector<FieldElement> elts() const {return elements;}
  // Compute the index of a nonzero element. If a belongs to the
  // current group return i and set s, where a = elements[i]*s^2. If
  // a does not belong to the subgroup (mod squares) it is appended to
  // gens, r is incremented, s=1 and the new r is returned.
  unsigned int get_index(const FieldElement& a, FieldElement& s);
  string elt_str(unsigned int i) const;
  unsigned int rank() const {return r;}
  int order() const {return elements.size();}
  void display();
};

// an instance of class Eigenvalue consist of FieldElement a and an
// index i into a FieldModSq representing a*sqrt(elt(i))
class Eigenvalue {
private:
  FieldElement a;
  unsigned int root_index;
  FieldModSq* SqCl;
public:
  Eigenvalue() {;}
  Eigenvalue(const FieldElement& x, FieldModSq* S, unsigned int i=0)
    : a(x), root_index(i), SqCl(S)
  {;}
  FieldElement coeff() const {return a;}
  FieldElement root_part() const  { return SqCl->elt(root_index); }
  Eigenvalue operator*(Eigenvalue b) const;
  Eigenvalue operator/(Eigenvalue b) const;
  int is_zero() const {return a.is_zero();}
  string str() const;
};

inline ostream& operator<<(ostream& s, const Eigenvalue& x)
{ s << x.str(); return s;}

#endif
