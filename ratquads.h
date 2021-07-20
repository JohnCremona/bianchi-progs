// ratquads.h

#if     !defined(_RATQUADS_H)
#define _RATQUADS_H      1       //flags that this file has been included

#include "quads.h"

class RatQuad {

public:
  // constructors
  RatQuad(Quad nn=0, Quad d=1, int reduce=0);
  RatQuad(long a, long b, long dd, int reduce=0); // (a+b*w)/dd

  // RatQuad manipulations

  // reduce to lowest terms: when non-principal this method only divides
  // n and d by the gcd of the content:
  int reduce();
  // reduce, ensuring that resulting ideal is coprime to N
  void reduce(const Qideal& N);
  int is_integral() const {return d.nm==1;} // assumes reduced
  Quad num() const {return n;}
  Quad den() const {return d;}
  RatQuad recip() const {return RatQuad(d, n);}  // no reduction needed
  RatQuad translation_reduce() const             // reduce mod Quads
  {
    return (d==0? RatQuad(1,0): RatQuad(n%d,d));
  }
  Quad round() const {return n/d;}               // nearest Quad, using rounded division of Quads
  Qideal ideal() const;                          // ideal (n,d)
  Qideal denominator_ideal() const;              // (d)/(n,d)
  int is_principal() const;

  // Binary Operator Functions
  friend RatQuad operator+(const RatQuad&, const RatQuad&);
  friend RatQuad operator+(Quad, const RatQuad&);
  friend RatQuad operator+(const RatQuad&, Quad);
  friend RatQuad operator-(const RatQuad&, const RatQuad&);
  friend RatQuad operator-(Quad, const RatQuad&);
  friend RatQuad operator-(const RatQuad&, Quad);
  friend RatQuad operator*(const RatQuad&, const RatQuad&);
  friend RatQuad operator*(const RatQuad&, Quad);
  friend RatQuad operator*(Quad, const RatQuad&);
  friend RatQuad operator/(const RatQuad&, const RatQuad&);
  friend RatQuad operator/(const RatQuad&, Quad);
  friend RatQuad operator/(Quad, const RatQuad&);
  friend int operator==(const RatQuad&, const RatQuad&);
  friend int operator!=(const RatQuad&, const RatQuad&);
  friend ostream& operator<< (ostream&s, const RatQuad&);
  RatQuad& operator+=(const RatQuad&);
  RatQuad& operator+=(Quad);
  RatQuad& operator-=(const RatQuad&);
  RatQuad& operator-=(Quad);
  RatQuad& operator*=(const RatQuad&);
  RatQuad& operator*=(Quad);
  RatQuad& operator/=(const RatQuad&);
  RatQuad& operator/=(Quad);
  RatQuad operator+();
  RatQuad operator-();

  friend int cuspeq(const RatQuad& c1, const RatQuad& c2, const Quad& N, int plusflag);
  friend int cuspeq(const RatQuad& c1, const RatQuad& c2, const Qideal& N, int plusflag);

  // Implementation
private:
  Quad d, n;
};


// Inline RatQuad functions

inline RatQuad::RatQuad(Quad nn, Quad dd, int reduce)
{
  n=nn; d=dd;
  if (reduce)
    (*this).reduce();
}

inline RatQuad::RatQuad(long a, long b, long dd, int reduce) // (a+b*w)/dd
{
  n=Quad(a,b);
  d=Quad(dd);
  if (reduce)
    (*this).reduce();
}

inline RatQuad RatQuad::operator+()
{
  return *this;
}

inline RatQuad RatQuad::operator-()
{
  return RatQuad(-n, d);
}


// Definitions of compound-assignment operator member functions

inline RatQuad& RatQuad::operator+=(const RatQuad& r)
{
  n = n*r.d + d*r.n;
  d *= r.d;
  (*this).reduce();
  return *this;
}

inline RatQuad& RatQuad::operator+=(Quad q)
{
  n += d*q;
  (*this).reduce();
  return *this;
}

inline RatQuad& RatQuad::operator-=(const RatQuad& r)
{
  n = n*r.d - d*r.n;
  d *= r.d;
  (*this).reduce();
  return *this;
}

inline RatQuad& RatQuad::operator-=(Quad q)
{
  n -= d*q;
  (*this).reduce();
  return *this;
}

inline RatQuad& RatQuad::operator*=(Quad q)
{
  n*=q;
  (*this).reduce();
  return *this;
}

inline RatQuad& RatQuad::operator/=(Quad q)
{
  d*=q;
  (*this).reduce();
  return *this;
}

// Definitions of non-member binary operator functions

inline RatQuad operator+(const RatQuad& r1, const RatQuad& r2)
{
  return RatQuad(r1.n*r2.d + r2.n*r1.d, r1.d*r2.d, 1);
}

inline RatQuad operator+(Quad q, const RatQuad& r)
{
  Quad n3 = q*r.d;
  n3 += r.n;
  return RatQuad(r.n + q*r.d, r.d, 1);
}

inline RatQuad operator+(const RatQuad& r, Quad q)
{
  return RatQuad(r.n + q*r.d, r.d, 1);
}

inline RatQuad operator-(const RatQuad& r1, const RatQuad& r2)
{
  return RatQuad(r1.n*r2.d - r2.n*r1.d, r1.d*r2.d, 1);
}

inline RatQuad operator-(Quad q, const RatQuad& r)
{
  return RatQuad(q*r.d - r.n, r.d, 1);
}

inline RatQuad operator-(const RatQuad& r, Quad q)
{
  return RatQuad(r.n - q*r.d, r.d, 1);
}

inline RatQuad operator*(const RatQuad& r, Quad q)
{
  return RatQuad(q*r.n, r.d, 1);
}

inline RatQuad operator*(Quad q, const RatQuad& r)
{
  return RatQuad(q*r.n, r.d, 1);
}

inline RatQuad operator/(const RatQuad& r, Quad q)
{
  return RatQuad(r.n, q*r.d, 1);
}

inline RatQuad operator/(Quad q, const RatQuad& r)
{
  return RatQuad(q*r.d, r.n, 1);
}

inline RatQuad operator/(const RatQuad& r1, const RatQuad& r2)
{
  return RatQuad(r1.n*r2.d, r2.n*r1.d, 1);
}

inline int operator==(const RatQuad& r1, const RatQuad& r2)
{
  return r1.n*r2.d == r2.n*r1.d;
}

inline int operator!=(const RatQuad& r1, const RatQuad& r2)
{
  return r1.n*r2.d != r2.n*r1.d;
}

inline ostream& operator<<(ostream& s, const RatQuad& r)
{
   if (r.d==Quad(0)) s<<"oo";
   else if (r.d==Quad(1)) s<<r.n;
   else s << "(" << r.n << ")/(" << r.d << ")";
   return s;
}

// Cusp equivalence mod Gamma_0(N)

// Special case: N is a Quad and c1, c2 assumed principal and reduced:

int cuspeq(const RatQuad& c1, const RatQuad& c2, const Quad& N, int plusflag);

// General case: N is a Qideal

int cuspeq(const RatQuad& c1, const RatQuad& c2, const Qideal& N, int plusflag);

#endif
