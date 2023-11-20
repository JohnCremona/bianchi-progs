// FILE MAT22.H: 2x2 matrices over the ring of integers

#if     !defined(_MAT22_H)
#define _MAT22_H      1       //flags that this file has been included

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
  mat22(long ia, long ib, long ic, long id)
    :a(QUINT(ia)),b(QUINT(ib)),c(QUINT(ic)),d(QUINT(id)) {}
  static mat22 identity;
  static mat22 J;
  static mat22 S;
  static mat22 Tmat(const Quad& x) {return mat22(Quad::one,x,Quad::zero,Quad::one);}
  static mat22 diag(const Quad& x, const Quad& y) {return mat22(x,Quad::zero,Quad::zero,y);}
  static mat22 scalar(const Quad& x) {return mat22(x,Quad::zero,Quad::zero,x);}
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


// return a matrix [a, b; c, d] with det=1 and (c:d)=(cc:dd) in P^1(N)
mat22 lift_to_SL2(Qideal& N, const Quad& cc, const Quad& dd);

// return a matrix [a, b; c, d] with det=1 and c in M and (c:d)=(cc:dd) in P^1(N)
// If (u,v)!=(0,0) they should satisfy u+v=1 with u in N, v in M,
// otherwise such u,v will be computed.
mat22 lift_to_Gamma_0(Qideal& M, Qideal& N, const Quad& cc, const Quad& dd, const Quad& u = Quad::zero, const Quad& v = Quad::zero);

#endif

// END OF FILE MAT22.H
