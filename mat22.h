// FILE MAT22.H: 2x2 matrices over the ring of integers

#if     !defined(_MAT22_H)
#define _MAT22_H      1       //flags that this file has been included

#include <eclib/arith.h>
#include <eclib/unimod.h>
#include <assert.h>
#include "ratquads.h"
#include "primes.h"

class mat22 {
private:
   Quad a,b,c,d;
public:
  mat22() :a(Quad::zero),b(Quad::zero),c(Quad::zero),d(Quad::zero) {}
  mat22(const Quad& ia, const Quad& ib, const Quad& ic, const Quad& id)
    :a(ia),b(ib),c(ic),d(id) {}
  static mat22 identity;
  static mat22 J;
  static mat22 S;
  static mat22 Tmat(const Quad& x) {return mat22(Quad::one,x,Quad::zero,Quad::one);}
  static mat22 TS;
  static mat22 TiS;
  static mat22 R;

  // access to entries
  Quad entry(int i, int j) const
  {
    return (i==0? (j==0?a:b): (j==0?c:d));
  }

  // matrix multiplcation
  mat22 operator*(const mat22& M) const
  {
    Quad a1 = a*M.a+b*M.c, b1 = a*M.b+b*M.d, c1 = c*M.a+d*M.c, d1 = c*M.b+d*M.d;
    return mat22(a1, b1, c1, d1);
  }

  void operator*=(const mat22& M)
  {
    Quad a1 = a*M.a+b*M.c, b1 = a*M.b+b*M.d, c1 = c*M.a+d*M.c, d1 = c*M.b+d*M.d;
    a = a1;
    b = b1;
    c = c1;
    d = d1;
  }

  // matrix inverse (mod scalars)
  mat22 inverse() const
  {
    return mat22(d,-b,-c,a);
  }

  // left action on r/s as column vector, changing in place:
  void apply_left(Quad& r, Quad& s) const
  {
    Quad t = a*r+b*s;
    s = c*r+d*s;
    r = t;
  }
  RatQuad operator()(const RatQuad& q) const;
  RatQuad image_oo() const;
  RatQuad preimage_oo() const;
  RatQuad image_0() const;
  RatQuad preimage_0() const;

  // right action on (c:d) symbols as row vectors, changing in place
  void apply_right(Quad& sc, Quad& sd) const
  {
    Quad t = a*sc + c*sd;
    sd = b*sc + d*sd;
    sc = t;
  }
  // right action on (c:d) symbols as row vectors, changing in place
  void apply_right_inverse(Quad& sc, Quad& sd) const
  {
    Quad t = d*sc - c*sd;
    sd = -b*sc + a*sd;
    sc = t;
  }
  Quad det() const {return a*d-b*c;}
  Quad trace() const {return a+d;}

  int is_scalar() const {return (b.is_zero() && c.is_zero() && (a==d));}
  int is_unimodular(int strict=1) const {Quad dt=det(); return ((dt==Quad::one) || ((!strict) && (dt.norm()==1)));}

  friend ostream& operator<< (ostream&, const mat22&); // inline below
  friend void pseudo_euclidean_step(Quad&, Quad&, int&, Quad&, Quad&, Quad&, Quad&);
  friend class modsym;
};

inline ostream& operator<< (ostream& s, const mat22& m)
{
   s << "[" << (m.a) << "," << (m.b) << "; " << (m.c) << "," << (m.d) << "]";
   return s;
}

// Atkin-Lehner and Hecke operators

inline string opname(const Quad& p, const Quad& n)
{
  ostringstream ans;
  ans << (div(p,n) ? "W" : "T") << "(" << p << ")";
  return ans.str();
}

inline string opname(const Quadprime& P, const Qideal& N)
{
  ostringstream ans;
  ans << (P.divides(N) ? "W" : "T") << "(" << P << ")";
  return ans.str();
}

inline string opname(const Qideal& N)
{
  ostringstream ans;
  ans << "W(" << N << ")";
  return ans.str();
}

mat22 AtkinLehner(const Quad& p, const Quad& n); // P=(p) principal prime dividing N=(n)
mat22 AtkinLehner(const Quad& p, Qideal& N); // P=(p) principal prime dividing N
mat22 AtkinLehner(Qideal& M1, Qideal& M2); // assume [M1] square and M1,M2 coprime
mat22 AtkinLehnerP(Quadprime& P, const Qideal& N); // =AL(P^e,N/P^e) where P^e||N

vector<mat22> Hecke(const Quad& p);  // P=(p) principal prime
vector<mat22> Hecke(const Quad& p, Qideal& N); //  P=(p) principal prime not dividing N
vector<mat22> Hecke(Quadprime& P, Qideal& N); // assume [P] square
vector<mat22> HeckeSq(Quadprime& P, Qideal& N); // T_{P^2}, when P^2 is principal (and P not)
vector<mat22> HeckePQ(Quadprime& P, Quadprime& Q, Qideal& N); // assume P*Q principal, P,Q not dividing N

inline mat22 Fricke(const Quad& n)
{
  return mat22(Quad::zero,-Quad::one, n,Quad::zero);
}

inline mat22 Fricke(Qideal& N) // assumes [N] square
{
  Qideal One(Quad::one);
  return AtkinLehner(N, One);
}

// Matrix inducing T_{A,A} at level N, when A^2 is principal and A+N=1

mat22 Char(Qideal& A, const Qideal& N);

class matop {  // formal sum of 2x2 matrices
 public:
  vector<mat22> mats;
  string the_name;
  matop() {;}
  explicit matop(const mat22& m, const string& n="") :mats({m}), the_name(n) {;}
  explicit matop(const vector<mat22>& mlist, const string& n="") :mats(mlist), the_name(n) {;}
  mat22 operator[](int i) const {return mats[i];}
  int length() const {return mats.size();}
  string name() const {return the_name;}
};

inline matop AtkinLehnerOp(const Quad& p, const Quad& n)
{
  return matop(AtkinLehner(p,n), opname(p,n));
}

inline matop AtkinLehnerOp(Quadprime& P, const Qideal& N)
{
  return matop(AtkinLehnerP(P,N), opname(P,N));
}

inline matop HeckeOp(Quadprime& P, Qideal& N)
{
  return matop(Hecke(P,N), opname(P,N));
}

inline matop HeckeSqOp(Quadprime& P, Qideal& N)
{
  ostringstream s;
  s << "T(" << P << "^2)";
  return matop(HeckeSq(P,N), s.str());
}

inline matop HeckePQOp(Quadprime& P, Quadprime& Q, Qideal& N)
{
  ostringstream s;
  s << "T(" << P << "*" << Q << ")";
  return matop(HeckePQ(P,Q,N), s.str());
}

inline matop AtkinLehnerOrHeckeOp(Quadprime& P, Qideal& N)
{
  return (P.divides(N)? AtkinLehnerOp(P,N): HeckeOp(P,N));
}

inline matop FrickeOp(Qideal& N)
{
  return matop(Fricke(N), opname(N));
}

inline matop CharOp(Qideal& A, const Qideal& N)
{
  ostringstream s;
  s << "nu";
  return matop(Char(A,N), s.str());
}


#endif

// END OF FILE MAT22.H
