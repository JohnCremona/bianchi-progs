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
  void cancel();                           // cancel *this in situ
  friend Quad num(const RatQuad&);         // the numerator
  friend Quad den(const RatQuad&);         // the denominator
  friend RatQuad recip(const RatQuad&);    // reciprocal
  friend RatQuad reduce(const RatQuad& q); // reduce mod Quads
  friend Quad round(const RatQuad&);       // nearest Quad

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

  // Implementation
private:
  Quad d, n;
};


// Inline RatQuad functions

inline void RatQuad::cancel()                     // cancel *this in situ
{
 Quad g = quadgcd(n,d);
 if (quadnorm(g)>1) {n/=g; d/=g;}
 while (!pos(d)) {n*=fundunit; d*=fundunit;}
}

inline RatQuad::RatQuad(Quad nn, Quad dd, int reduce)
{
  n=nn; d=dd;
  if (reduce)
    (*this).cancel();
}

inline RatQuad::RatQuad(long a, long b, long dd, int reduce) // (a+b*w)/dd
{
  n=Quad(a,b);
  d=Quad(dd);
  if (reduce)
    (*this).cancel();
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
  (*this).cancel();
  return *this;
}

inline RatQuad& RatQuad::operator+=(Quad q)
{
  n += d*q;
  (*this).cancel();
  return *this;
}

inline RatQuad& RatQuad::operator-=(const RatQuad& r)
{
  n = n*r.d - d*r.n;
  d *= r.d;
  (*this).cancel();
  return *this;
}

inline RatQuad& RatQuad::operator-=(Quad q)
{
  n -= d*q;
  (*this).cancel();
  return *this;
}

inline RatQuad& RatQuad::operator*=(Quad q)
{
  n*=q;
  (*this).cancel();
  return *this;
}

inline RatQuad& RatQuad::operator/=(Quad q)
{
  d*=q;
  (*this).cancel();
  return *this;
}

// Definitions of non-member RatQuad functions

inline Quad num(const RatQuad& r)
{
  return r.n;
}

inline Quad den(const RatQuad& r)
{
  return r.d;
}

inline RatQuad recip(const RatQuad& r)
{
  return RatQuad(r.d, r.n); // no cancelling needed
}

inline Quad round(const RatQuad& r)
{
  return r.n / r.d;   // uses rounded division of Quads
}

inline RatQuad reduce(const RatQuad& r)
{
  if (r.d==Quad(0)) return RatQuad(1,0);
  return RatQuad(r.n%r.d,r.d);
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
   if (r.d==Quad(1)) s<<r.n;
   else s << "(" << r.n << ")/(" << r.d << ")";
   return s;
}

#endif
