// File HECKEFIELD.H: class for working with number fields for Hecke eigenvalues
//////////////////////////////////////////////////////////////////////////

#ifndef _HECKEFIELD_H
#define _HECKEFIELD_H      1

#include "matprocs.h"

class HeckeField;
class HeckeFieldElement;

class HeckeField {
  friend class HeckeFieldElement;
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
  //HeckeField(const ZZX& p);
  HeckeField(); // defaults to Q
  HeckeField(const mat_m& m, const ZZ& den = to_ZZ(1), string a="a", int verb=0);
  HeckeField(const ZZX& p, string a="a", int verb=0);
  HeckeFieldElement one();
  HeckeFieldElement zero();
  HeckeFieldElement gen();
  HeckeFieldElement element(const vec_m& c, const ZZ& d=to_ZZ(1), int raw=0);
  int degree() const {return d;}
  ZZX poly() const {return minpoly;}
  mat_m basis() const {return Binv;} // columns are Bfactor * coeffs of basis w.r.t. a-powers
  ZZ basis_factor() const {return Bdet;}
  mat_m inv_basis() const {return B;} // columns are coeffs of a-powers w.r.t. basis
  void display(ostream&s = cout, int raw=0); // if raw, also display raw basis
  void display_bases(ostream&s = cout) const; // display powers of A and C and bases in both embeddings
};

class HeckeFieldElement {
  friend class HeckeField;
private:
  HeckeField* F;
  vec_m coords; // length F->d
  ZZ denom;     // >=1
  // the field element is (1/denom)*coords-combination of power basis of F
  // NB On construction every element will be reduced using cancel()
  void cancel(); // divides through by gcd(content(coords, denom))
public:
  HeckeFieldElement( HeckeField* HF)
    :F(HF), coords(vec_m(HF->d)), denom(to_ZZ(1))  {};
  // raw means the given coords are w.r.t. the B-basis
  HeckeFieldElement( HeckeField* HF, const vec_m& c, const ZZ& d=to_ZZ(1), int raw=0);
  // creation from a rational
  HeckeFieldElement( HeckeField* HF, const ZZ& a, const ZZ& d=to_ZZ(1));

  string str() const;
  const HeckeField* field() const {return F;}
  mat_m matrix() const; // ignores denom
  ZZX charpoly() const {return ::charpoly(matrix());}
  ZZX minpoly() const;
  int degree() const {return deg(minpoly());}
  int is_zero() const;
  int is_one() const;
  int is_minus_one() const;
  int is_generator() const {return degree()==F->d;}
  int operator==(const HeckeFieldElement& b) const;

  HeckeFieldElement operator+(const HeckeFieldElement& b) const; // add
  void operator+=(const HeckeFieldElement& b); // add b to this
  void operator+=(const ZZ& b) { operator+=(HeckeFieldElement(F,b));} // add b

  HeckeFieldElement operator-(const HeckeFieldElement& b) const; // subtract
  void operator-=(const HeckeFieldElement& b); // subtract b from this
  void operator-=(const ZZ& b) { operator-=(HeckeFieldElement(F,b));} // subtract b
  HeckeFieldElement operator-() const;                           // unary minus

  HeckeFieldElement operator*(const HeckeFieldElement& b) const; // product
  void operator*=(const HeckeFieldElement& b); // multiply by b
  void operator*=(const ZZ& b) { operator*=(HeckeFieldElement(F,b));} // multiply by b

  HeckeFieldElement inverse() const; // raise error if zero      // inverse
  HeckeFieldElement operator/(const HeckeFieldElement& b) const; // divide (raise error if b is zero)
  void operator/=(const HeckeFieldElement& b);                        // divide by b
  void operator/=(const ZZ& b) { operator/=(HeckeFieldElement(F,b));} // divide by b

  // NB for a in F, either [Q(sqrt(a))=Q(a)] or [Q(sqrt(a)):Q(a)]=2.
  // The first function only applies when a has maximal degree:
  // return 1 and r s.t. r^2=this, with deg(r)=degree(), else 0
  int is_absolute_square(HeckeFieldElement& r)  const;
  // Same as above if the min poly is known
  int is_absolute_square(HeckeFieldElement& r, const ZZX& minpol)  const;
  // The second function applies in general: return 1 and r
  // s.t. r^2=this, with deg(r)=degree(), else 0. Here, ntries is the
  // number of squares this is multiplied by to get odd co-degree.
  int is_square(HeckeFieldElement& r, int ntries=100) const;
};

HeckeFieldElement evaluate(const ZZX& f, const HeckeFieldElement a);

inline ostream& operator<<(ostream& s, const HeckeFieldElement& x)
{ s << x.str(); return s;}

#endif
