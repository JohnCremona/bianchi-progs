// ratquads.h

#if     !defined(_RATQUADS_H)
#define _RATQUADS_H      1       //flags that this file has been included

#include "quads.h"

class RatQuad {

public:
  // constructors
  RatQuad(Quad nn=0, Quad d=1);
  RatQuad(long a, long b, long dd); // (a+b*w)/dd

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

inline RatQuad::RatQuad(Quad nn, Quad dd)
{
  n=nn; d=dd;
  (*this).cancel();
}

inline RatQuad::RatQuad(long a, long b, long dd) // (a+b*w)/dd
{
  n=Quad(a,b);
  d=Quad(dd);
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

inline RatQuad& RatQuad::operator+=(const RatQuad& q2)
{
  Quad n1=n*q2.d, n2=d*q2.n;
  n = n1+n2;
  d *= q2.d;
  (*this).cancel();
  return *this;
}

inline RatQuad& RatQuad::operator+=(Quad n2)
{
  n += d*n2;
  return *this;
}

inline RatQuad& RatQuad::operator-=(const RatQuad& q2)
{
  Quad n1=n*q2.d, n2=d*q2.n;
  n = n1-n2;
  d *= q2.d;
  (*this).cancel();
  return *this;
}

inline RatQuad& RatQuad::operator-=(Quad n2)
{
  n -= d*n2;
  return *this;
}

inline RatQuad& RatQuad::operator*=(Quad n2)
{
  n*=n2;
  (*this).cancel();
  return *this;
}

inline RatQuad& RatQuad::operator/=(Quad n2)
{
  d*=n2;
  (*this).cancel();
  return *this;
}

// Definitions of non-member RatQuad functions

inline Quad num(const RatQuad& q)
{
  return q.n;
}

inline Quad den(const RatQuad& q)
{
  return q.d;
}

inline RatQuad recip(const RatQuad& q)
{
  return RatQuad(q.d, q.n);
}

inline Quad round(const RatQuad& q)
{
  return q.n / q.d;   // uses rounded division of Quads
}

inline RatQuad reduce(const RatQuad& q)
{
  if (q.d==Quad(0)) return RatQuad(1,0);
  return RatQuad(q.n%q.d,q.d);
}


// Definitions of non-member binary operator functions

inline RatQuad operator+(const RatQuad& q1, const RatQuad& q2)
{
  Quad n1 = q1.n*q2.d;
  Quad n2 = q2.n*q1.d;
  Quad n3 = n1+n2;
  Quad d3 = q1.d * q2.d;
  return RatQuad(n3, d3);
}

inline RatQuad operator+(Quad n1, const RatQuad& q2)
{
  Quad n3 = n1*q2.d;
  n3 += q2.n;
  return RatQuad(n3, q2.d);
}

inline RatQuad operator+(const RatQuad& q1, Quad n2)
{
  return RatQuad(q1.n + n2*q1.d, q1.d);
}

inline RatQuad operator-(const RatQuad& q1, const RatQuad& q2)
{
  Quad n1 = q1.n*q2.d;
  Quad n2 = q2.n*q1.d;
  Quad n3 = n1-n2;
  Quad d3 = q1.d * q2.d;
  return RatQuad(n3, d3);
}

inline RatQuad operator-(Quad n1, const RatQuad& q2)
{
  Quad n3 = n1*q2.d;
  n3 -= q2.n;
  return RatQuad(n3, q2.d);
}

inline RatQuad operator-(const RatQuad& q1, Quad n2)
{
        return RatQuad(q1.n - n2*q1.d, q1.d);
}

inline RatQuad operator*(const RatQuad& q1, Quad n2)
{
        return RatQuad(q1.n*n2, q1.d);
}

inline RatQuad operator*(Quad n1, const RatQuad& q2)
{
        return RatQuad(q2.n*n1, q2.d);
}

inline RatQuad operator/(const RatQuad& q1, Quad n2)
{
        return RatQuad(q1.n, q1.d*n2);
}

inline int operator==(const RatQuad& q1, const RatQuad& q2)
{
        return q1.n*q2.d == q2.n*q1.d;
}

inline int operator!=(const RatQuad& q1, const RatQuad& q2)
{
        return q1.n*q2.d != q2.n*q1.d;
}

inline ostream& operator<<(ostream& s, const RatQuad& q)
{
   if (q.d==Quad(1)) s<<q.n;
   else s << "(" << q.n << ")/(" << q.d << ")";
   return s;
}

#endif
