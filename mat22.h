// FILE MAT22.H: 2x2 matrices over the ring of integers

#if     !defined(_MAT22_H)
#define _MAT22_H      1       //flags that this file has been included

#include <eclib/arith.h>
#include <eclib/unimod.h>
#include <assert.h>
#include "ratquads.h"

class mat22 {
private:
   Quad a,b,c,d;
public:
  mat22() :a(0),b(0),c(0),d(0) {}
  mat22(const Quad ia, const Quad ib, const Quad ic, const Quad id)
    :a(ia),b(ib),c(ic),d(id) {}
  static mat22 identity;
  static mat22 J;
  static mat22 S;
  static mat22 Tmat(const Quad& x) {return mat22(1,x,0,1);}
  static mat22 TS;
  static mat22 TiS;
  static mat22 R;

  // access to entries
  Quad entry(int i, int j) const
  {
    return (i==0? (j==0?a:b): (j==0?c:d));
  }

  // matrix multiplcation
  mat22 operator*(const mat22 M) const
  {
    return mat22(a*M.a+b*M.c, a*M.b+b*M.d, c*M.a+d*M.c, c*M.b+d*M.d);
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
  RatQuad operator()(const RatQuad& q) const; // implemented in cusp.cc

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

  friend ostream& operator<< (ostream&, const mat22&); // inline below
  friend void pseudo_euclidean_step(Quad&, Quad&, int&, Quad&, Quad&, Quad&, Quad&);
  friend class modsym;
};

inline ostream& operator<< (ostream& s, const mat22& m)
{
   s << "[" << (m.a) << "," << (m.b) << "; " << (m.c) << "," << (m.d) << "]";
   return s;
}

class matop {  // formal sum of 2x2 matrices
 public:
  vector<mat22> mats;
  matop(const Quad& p, const Quad& n);  // constructor for hecke ops (non-ideal version)
  matop(Qideal& P, Qideal& N);          // constructor for hecke ops (ideal version)
  mat22 operator[](int i) const {return mats[i];}
  int length() const {return mats.size();}
};

#endif

// END OF FILE MAT22.H
