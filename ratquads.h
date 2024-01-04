// ratquads.h

#if     !defined(_RATQUADS_H)
#define _RATQUADS_H      1       //flags that this file has been included

#include "quads.h"
#include "rat.h"

class RatQuad {

public:
  // constructors
  RatQuad(const Quad& nn=Quad::zero, const Quad& d=Quad::one, int reduce=0);
  RatQuad(INT a, INT b, INT dd, int reduce=0); // (a+b*w)/dd
  RatQuad(const RAT& a);

  // RatQuad manipulations

  // reduce to lowest terms: when non-principal this method only divides
  // n and d by the gcd of the content:
  int reduce();
  // reduce, ensuring that resulting ideal is coprime to N
  void reduce(const Qideal& N);
  void reduce(long n);
  int is_integral() const {return div(d,n);}
  int is_integral(Quad& a) const {return div(d,n,a);}
  int is_infinity() const {return d.is_zero();}
  int is_finite() const {return !d.is_zero();}
  int is_zero() const {return n.is_zero();}
  Quad num() const {return n;}
  Quad den() const {return d;}
  RatQuad recip() const {return RatQuad(d, n);}  // no reduction needed
  RatQuad conj() const
  {
    Quad nc = n.conj(), dc = d.conj();
    while (!pos(dc)) {nc*=fundunit; dc*=fundunit;}
    return RatQuad(nc, dc);   // no reduction needed
  }
  // reduce mod Quads -- see also reduce_to_rectangle()
  RatQuad translation_reduce() const
  {
    if (d.is_zero())
      return RatQuad(n,d);
    else
      return RatQuad(n%d,d);
  }
  void normalise();                              // scale so ideal is a standard class rep
  Quad round() const {return n/d;}               // nearest Quad, using rounded division of Quads
  Qideal ideal() const;                          // ideal (n,d)
  Qideal denominator_ideal() const;              // (d)/(n,d)
  int is_principal() const;

  RAT norm() const {return RAT(n.norm(), d.norm());}

  // Rational {x,y} such that this=x+y*w (or x+y*sqrt(-d) if rectangle=1)
  vector<RAT> coords(int rectangle=0) const;
  RAT x_coord(int rectangle=0) const {return coords(rectangle)[0];}
  RAT y_coord(int rectangle=0) const {return coords(rectangle)[1];}
  int in_rectangle() const;            // x in (-1/2,1/2] and y in (-1/2,1/2] (even d) or (-1/4,1/4] (odd d)
  int in_quarter_rectangle() const;    // x in [0,1/2] and y in [0,1/2] (even d) or [0,1/4] (odd d)
  friend RatQuad reduce_to_rectangle(const RatQuad&, Quad&);   // subtract Quad to put z into rectangle
  friend vector<Quad> nearest_quads(const RatQuad&, int just_one);    // list of Quad(s) a s.t. N(a-z)<1
  // Binary Operator Functions
  friend RatQuad operator+(const RatQuad&, const RatQuad&);
  friend RatQuad operator+(const Quad&, const RatQuad&);
  friend RatQuad operator+(const RatQuad&, const Quad&);
  friend RatQuad operator+(const RatQuad&, long);
  friend RatQuad operator-(const RatQuad&, const RatQuad&);
  friend RatQuad operator-(const Quad&, const RatQuad&);
  friend RatQuad operator-(const RatQuad&, const Quad&);
  friend RatQuad operator-(const RatQuad&, long);
  friend RatQuad operator*(const RatQuad&, const RatQuad&);
  friend RatQuad operator*(const RatQuad&, const Quad&);
  friend RatQuad operator*(const Quad&, const RatQuad&);
  friend RatQuad operator*(const RatQuad&, long);
  friend RatQuad operator/(const RatQuad&, const RatQuad&);
  friend RatQuad operator/(const RatQuad&, const Quad&);
  friend RatQuad operator/(const Quad&, const RatQuad&);
  friend RatQuad operator/(const RatQuad&, long);
  friend int operator==(const RatQuad&, const RatQuad&);
  friend int operator!=(const RatQuad&, const RatQuad&);
  friend ostream& operator<< (ostream&s, const RatQuad&);
  void operator+=(const RatQuad&);
  void operator+=(const Quad&);
  void operator-=(const RatQuad&);
  void operator-=(const Quad&);
  void operator*=(const RatQuad&);
  void operator*=(const Quad&);
  void operator/=(const RatQuad&);
  void operator/=(const Quad&);
  RatQuad operator+() const;
  RatQuad operator-() const;

  friend int cuspeq(const RatQuad& c1, const RatQuad& c2, const Quad& N, int plusflag);
  friend int cuspeq(const RatQuad& c1, const RatQuad& c2, const Qideal& N, int plusflag);

  // Implementation
private:
  Quad d, n;
};

class modsym {
 private:
    RatQuad a,b;
 public:
  modsym() :a(), b() {}
    modsym(const RatQuad& ra, const RatQuad& rb) :a(ra),b(rb) {}
    explicit modsym(const mat22& M, int type=0);              //conversion from (c:d)
    RatQuad alpha() const {return a;}
    RatQuad  beta() const {return b;}
    modsym reverse() const {return modsym(b,a);}
    modsym conj() const {return modsym(a.conj(), b.conj());}
    friend int operator==(const modsym& m1, const modsym& m2);
    friend int operator!=(const modsym& m1, const modsym& m2);
  friend ostream& operator<< (ostream& s, const modsym& m); //inline below
};

// Inline RatQuad functions

inline RatQuad::RatQuad(const Quad& nn, const Quad& dd, int reduce)
  :d(dd), n(nn)
{
  if (reduce)
    this->reduce();
}

inline RatQuad::RatQuad(INT a, INT b, INT dd, int reduce) // (a+b*w)/dd
  :d(a,b), n(dd,1)
{
  if (reduce)
    this->reduce();
}

inline RatQuad::RatQuad(const RAT& a)
  :d(a.den()), n(a.num())
{};


inline RatQuad RatQuad::operator+() const
{
  return *this;
}

inline RatQuad RatQuad::operator-() const
{
  return RatQuad(-n, d);
}


// Definitions of compound-assignment operator member functions

inline void RatQuad::operator+=(const RatQuad& r)
{
  n = mma(n,r.d,d,r.n); // n*r.d + d*r.n
  d *= r.d;
}

inline void RatQuad::operator+=(const Quad& q)
{
  n.addprod(d,q);
}

inline void RatQuad::operator-=(const RatQuad& r)
{
  n = mms(n,r.d,d,r.n); // n*r.d - d*r.n
  d *= r.d;
}

inline void RatQuad::operator-=(const Quad& q)
{
  n.subprod(d,q);
}

inline void RatQuad::operator*=(const Quad& q)
{
  n*=q;
}

inline void RatQuad::operator/=(const Quad& q)
{
  d*=q;
}

// Definitions of non-member binary operator functions

inline RatQuad operator+(const RatQuad& r1, const RatQuad& r2)
{
  return RatQuad(mma(r1.n,r2.d,r2.n,r1.d), r1.d*r2.d, 1);
}

inline RatQuad operator+(const Quad& q, const RatQuad& r)
{
  Quad n = r.n;
  n.addprod(q,r.d); // r.n + q*r.d
  return RatQuad(n, r.d, 1);
}

inline RatQuad operator+(const RatQuad& r, const Quad& q)
{
  Quad n = r.n;
  n.addprod(q,r.d); // r.n + q*r.d
  return RatQuad(n, r.d, 1);
}

inline RatQuad operator+(const RatQuad& r, long q)
{
  Quad n = r.n;
  n.addprod(q,r.d); // r.n + q*r.d
  return RatQuad(n, r.d, 1);
}

inline RatQuad operator-(const RatQuad& r1, const RatQuad& r2)
{
  return RatQuad(mms(r1.n,r2.d, r2.n,r1.d), r1.d*r2.d, 1);
}

inline RatQuad operator-(const Quad& q, const RatQuad& r)
{
  Quad n = -(r.n);
  n.addprod(q,r.d); // -r.n + q*r.d
  return RatQuad(n, r.d, 1);
}

inline RatQuad operator-(const RatQuad& r, const Quad& q)
{
  Quad n = r.n;
  n.subprod(q,r.d); // r.n - q*r.d
  return RatQuad(n, r.d, 1);
}

inline RatQuad operator-(const RatQuad& r, long q)
{
  Quad n = r.n;
  n.subprod(q,r.d); // r.n - q*r.d
  return RatQuad(n, r.d, 1);
}

inline RatQuad operator*(const RatQuad& q, const RatQuad& r)
{
  return RatQuad(q.n*r.n, q.d*r.d, 1);
}

inline RatQuad operator*(const RatQuad& r, const Quad& q)
{
  return RatQuad(q*r.n, r.d, 1);
}

inline RatQuad operator*(const RatQuad& r, long q)
{
  return RatQuad(q*r.n, r.d, 1);
}

inline RatQuad operator*(const Quad& q, const RatQuad& r)
{
  return RatQuad(q*r.n, r.d, 1);
}

inline RatQuad operator/(const RatQuad& r, const Quad& q)
{
  return RatQuad(r.n, q*r.d, 1);
}

inline RatQuad operator/(const RatQuad& r, long q)
{
  return RatQuad(r.n, q*r.d, 1);
}

inline RatQuad operator/(const Quad& q, const RatQuad& r)
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
  if ((r.d).is_zero()) s<<"oo";
  else if (r.d==Quad::one) s<<r.n;
   else s << "(" << r.n << ")/(" << r.d << ")";
   return s;
}

inline int operator==(const modsym& m1, const modsym& m2)
{
  return m1.a==m2.a && m1.b==m2.b;
}

inline int operator!=(const modsym& m1, const modsym& m2)
{
  return m1.a!=m2.a || m1.b!=m2.b;
}

// Cusp equivalence mod Gamma_0(N)

// Special case: N is a Quad and c1, c2 assumed principal and reduced:

int cuspeq(const RatQuad& c1, const RatQuad& c2, const Quad& N, int plusflag);

// General case: N is a Qideal

int cuspeq(const RatQuad& c1, const RatQuad& c2, const Qideal& N, int plusflag);

// Debugging version, does a second conjugate test and compares
int cuspeq_conj(const RatQuad& c1, const RatQuad& c2, const Qideal& N, int plusflag);

inline ostream& operator<< (ostream& s, const modsym& m)
{
   s << "{" << (m.a) << "," << (m.b) << "}";
   return s;
}

// Finding cusp in list, with or without translation

// Return index of c in clist, or -1 if not in list
int cusp_index(const RatQuad& c, const vector<RatQuad>& clist);

// Return index i of c mod O_K in clist, with t=c-clist[i], or -1 if not in list
int cusp_index_with_translation(const RatQuad& c, const vector<RatQuad>& clist, Quad& t);

// Comparison function (based only on norm of denominator)
struct Cusp_comparison {
  bool operator() (const RatQuad& lhs, const RatQuad& rhs) const
  {return lhs.den().norm()<rhs.den().norm();}
};

extern Cusp_comparison Cusp_cmp;

#endif
