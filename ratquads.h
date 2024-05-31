// ratquads.h

#if     !defined(_RATQUADS_H)
#define _RATQUADS_H      1       //flags that this file has been included

#include "quads.h"
#include "rat.h"

class RatQuad {

public:
  // constructors
  RatQuad(const Quad& nn=Quad::zero, const Quad& d=Quad::one, int red=0);
  RatQuad(INT a, INT b, INT dd, int red=0); // (a+b*w)/dd
  RatQuad(long a, long b, long dd, int red=0); // (a+b*w)/dd
  explicit RatQuad(long a);
  explicit RatQuad(const RAT& a);
  static RatQuad infinity() {return RatQuad(Quad::one, Quad::zero, 0);}
  // RatQuad manipulations

  // reduce to "lowest terms", so the ideal is one of the class representatives
  void reduce();
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
  int ideal_class() const;                       // in {0,1,...,h-1}

  RAT norm() const {return RAT(n.norm(), d.norm());}

  // Rational {x,y} such that this=x+y*w (or x+y*sqrt(-d) if rectangle=1)
  vector<RAT> coords(int rectangle=0) const;
  RAT x_coord(int rectangle=0) const {return coords(rectangle)[0];}
  RAT y_coord(int rectangle=0) const {return coords(rectangle)[1];}
  int in_rectangle() const;            // x in (-1/2,1/2] and y in (-1/2,1/2] (even d) or (-1/4,1/4] (odd d)
  int in_quarter_rectangle() const;    // x in [0,1/2] and y in [0,1/2] (even d) or [0,1/4] (odd d)
  friend RatQuad reduce_to_rectangle(const RatQuad&, Quad&);   // subtract Quad to put z into rectangle
  friend RatQuad reduce_to_rectangle(const RatQuad&);   // subtract Quad to put z into rectangle
  // list of Quad(s) q s.t. N(q-z)<1:
  friend vector<Quad> nearest_quads(const RatQuad&, int just_one);
  // list of Quad(s) q s.t. N(q-a/b)<1, i.e. N(a-b*q)<N(b):
  friend vector<Quad> nearest_quads_to_quotient(const Quad&, const Quad&, int just_one);
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
  friend RatQuad operator*(const RatQuad&, int);
  friend RatQuad operator*(const RatQuad&, const INT&);
  friend RatQuad operator*(const RatQuad&, const RAT&);
  friend RatQuad operator/(const RatQuad&, const RatQuad&);
  friend RatQuad operator/(const RatQuad&, const Quad&);
  friend RatQuad operator/(const Quad&, const RatQuad&);
  friend RatQuad operator/(const RatQuad&, const INT&);
  friend RatQuad operator/(const RatQuad&, const RAT&);
  friend RatQuad operator/(const RatQuad&, long);
  friend RatQuad operator/(const RatQuad&, int);
  friend int operator==(const RatQuad&, const RatQuad&);
  friend int operator!=(const RatQuad&, const RatQuad&);
  friend int operator<(const RatQuad&, const RatQuad&); // used for sorting, not a mathematical order
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

  // [a,b,c,d]=(a-b)(c-d)/(a-d)(c-b)
  friend RatQuad crossratio(const RatQuad& a, const RatQuad& b, const RatQuad& c, const RatQuad& d);
  // [oo,b,c,d]=(c-d)/(c-b)
  friend RatQuad crossratio3(const RatQuad& b, const RatQuad& c, const RatQuad& d);
  // Sign of imagainary part of c.r.:
  friend int sign_im_cr(const RatQuad& a, const RatQuad& b, const RatQuad& c, const RatQuad& d);
  friend int sign_im_cr(const RatQuad& b, const RatQuad& c, const RatQuad& d);

  friend int cuspeq(const RatQuad& c1, const RatQuad& c2, const Quad& N, int plusflag);
  friend int cuspeq(const RatQuad& c1, const RatQuad& c2, const Qideal& N, int plusflag);
  friend int integral_difference(const RatQuad& a, const RatQuad& b, Quad& t);

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
  RatQuad alpha() const {return a;}
  RatQuad  beta() const {return b;}
  modsym reverse() const {return modsym(b,a);}
  modsym conj() const {return modsym(a.conj(), b.conj());}
  int operator==(const modsym& m2) const {return a==m2.a && b==m2.b; }
  int operator!=(const modsym& m2) const {return a!=m2.a || b!=m2.b; }
  inline friend ostream& operator<< (ostream& s, const modsym& m)
  {
    s << "{" << (m.a) << "," << (m.b) << "}";
    return s;
  }
};

// Inline RatQuad functions

inline RatQuad::RatQuad(const Quad& nn, const Quad& dd, int red)
  :d(dd), n(nn)
{
  if (red)
    reduce();
}

inline RatQuad::RatQuad(INT a, INT b, INT dd, int red) // (a+b*w)/dd
  :d(dd,INT(0)), n(a,b)
{
  if (red)
    reduce();
}

inline RatQuad::RatQuad(long a, long b, long dd, int red) // (a+b*w)/dd
  :d(dd), n(a,b)
{
  if (red)
    reduce();
}

inline RatQuad::RatQuad(long a)
  :d(1), n(a)
{};

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
  return q+r;
}

inline RatQuad operator+(const RatQuad& r, long q)
{
  return Quad(q)+r;
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

inline RatQuad operator*(long q, const RatQuad& r)
{
  return r*q;
}

inline RatQuad operator*(const RatQuad& r, int q)
{
  return r*long(q);
}

inline RatQuad operator*(int q, const RatQuad& r)
{
  return r*q;
}

inline RatQuad operator*(const RatQuad& r, const INT& q)
{
  return RatQuad(r.n*q, r.d, 1);
}

inline RatQuad operator*(const INT& q, const RatQuad& r)
{
  return r*q;
}

inline RatQuad operator*(const RatQuad& r, const RAT& q)
{
  return RatQuad(r.n*q.num(), r.d*q.den(), 1);
}

inline RatQuad operator*(const RAT& q, const RatQuad& r)
{
  return r*q;
}

inline RatQuad operator*(const Quad& q, const RatQuad& r)
{
  return RatQuad(q*r.n, r.d, 1);
}

inline RatQuad operator/(const RatQuad& r, const Quad& q)
{
  return RatQuad(r.n, q*r.d, 1);
}

inline RatQuad operator/(const RatQuad& r, const INT& q)
{
  return RatQuad(r.n, q*r.d, 1);
}

inline RatQuad operator/(const RatQuad& r, const RAT& q)
{
  return RatQuad(r.n*q.den(), r.d*q.num(), 1);
}

inline RatQuad operator/(const RatQuad& r, long q)
{
  return RatQuad(r.n, q*r.d, 1);
}

inline RatQuad operator/(const RatQuad& r, int q)
{
  return r/long(q);
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

inline int operator<(const RatQuad& r1, const RatQuad& r2) // used for sorting, not a mathematical order
{
  return (r1.d<r2.d) || ((r1.d==r2.d) && (r1.n<r2.n));
}

inline ostream& operator<<(ostream& s, const RatQuad& r)
{
  if ((r.d).is_zero()) s<<"oo";
  else if (r.d==Quad::one) s<<r.n;
   else s << "(" << r.n << ")/(" << r.d << ")";
   return s;
}

// Cusp equivalence mod Gamma_0(N)

// Special case: N is a Quad and c1, c2 assumed principal and reduced:

int cuspeq(const RatQuad& c1, const RatQuad& c2, const Quad& N, int plusflag);

// General case: N is a Qideal

int cuspeq(const RatQuad& c1, const RatQuad& c2, const Qideal& N, int plusflag);

// Debugging version, does a second conjugate test and compares
int cuspeq_conj(const RatQuad& c1, const RatQuad& c2, const Qideal& N, int plusflag);

// Finding cusp in list, with or without translation

// Return index of c in clist, or -1 if not in list
int cusp_index(const RatQuad& c, const vector<RatQuad>& clist);

// Return index i of c mod O_K in clist, with t=c-clist[i], or -1 if not in list
int cusp_index_with_translation(const RatQuad& c, const vector<RatQuad>& clist, Quad& t);

// list of Quad(s) q s.t. N(q-z)<1:
vector<Quad> nearest_quads(const RatQuad&, int just_one);
// list of Quad(s) q s.t. N(q-a/b)<1, i.e. N(a-b*q)<N(b):
vector<Quad> nearest_quads_to_quotient(const Quad&, const Quad&, int just_one);

// Return list of 0, 1 or 2 sqrts of a rational r in k
vector<RatQuad> sqrts_in_k(const RAT& r);

typedef vector<RatQuad> CuspList;  // may refactor using sets later
typedef std::set<RatQuad> CuspPair;

#endif
