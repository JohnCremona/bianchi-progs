// RAT.h: interface to FLINT fmpq rational numbers

#if     !defined(_RAT_H)
#define _RAT_H      1       //flags that this file has been included

#include "arith_extras.h"
#include "int.h"
#include <flint/fmpq.h>

class RAT {

private:
  fmpq_t q;

public:
  // constructors
  RAT() {
    fmpq_init(q); // sets to 0
  }
  RAT(const INT& n) {
    fmpz_t den; fmpz_init(den); fmpq_set_fmpz_frac(q, n.z, den);
  }
  RAT(long n) {
    fmpz_t num, den; fmpz_init(den); fmpz_set_si(num, n); fmpq_set_fmpz_frac(q, num, den);
  }
  RAT(const INT& n, const INT& d) {
    fmpq_set_fmpz_frac(q, n.z, d.z);
  }
  RAT(long n, long d) {
    fmpz_t num, den; fmpz_set_si(den, d); fmpz_set_si(num, n); fmpq_set_fmpz_frac(q, num, den);
  }
  RAT(const RAT& x) {
    fmpq_set(q, x.q);
  }
  ~RAT() {
    fmpq_clear(q);
  }

  RAT& operator=(const RAT& x) {
    fmpq_set(q, x.q); return *this;
  }
  RAT operator=(long n) {
    RAT x(n); return x;
  }

  // RAT manipulations
  void cancel() {
    fmpq_canonicalise(q); // cancel *this in situ
  }
  INT num() const       // the numerator
  {
    INT n; fmpz_init_set(n.z, fmpq_numref(q)); return n;
  }
  INT den() const        // the denominator
  {
    INT n; fmpz_init_set(n.z, fmpq_denref(q)); return n;
  }
  RAT recip() const  // the reciprocal
  {
    RAT x;  fmpq_inv(x.q, q); return x;
  }

  INT round() const;      // nearest integer

  // Binary Operator Functions
  RAT operator+(const RAT& b) const {
    RAT c; fmpq_add(c.q, q, b.q); return c;
  }
  RAT operator+(const INT& b) const {
    RAT c; fmpq_add_fmpz(c.q, q, b.z); return c;
  }
  friend inline RAT operator+(const INT& a, const RAT& b) {return b+a;}
  RAT operator-(const RAT& b) const {
    RAT c; fmpq_sub(c.q, q, b.q); return c;
  }
  RAT operator-(INT b) const {
    RAT c; fmpq_sub_fmpz(c.q, q, b.z); return c;
  }
  friend inline RAT operator-(const INT& a, const RAT& b) {return -(b-a);}

  RAT operator*(const RAT& b) const {
    RAT c; fmpq_mul(c.q, q, b.q); return c;
  }
  RAT operator*(const INT& b) const {
    RAT c; fmpq_mul_fmpz(c.q, q, b.z); return c;
  }
  friend inline RAT operator*(const INT& a, const RAT& b) {return b*a;}

  RAT operator/(const RAT& b) const {
    RAT c; fmpq_div(c.q, q, b.q); return c;
  }
  RAT operator/(const INT& b) const {
    RAT c; fmpq_div_fmpz(c.q, q, b.z); return c;
  }
  friend inline RAT operator/(const INT& a, const RAT& b) {RAT c = b/a; return c.recip();}

  int operator==(const RAT& b) const {
    return fmpq_equal(q, b.q);
  }
  int operator==(const INT& b) const {
    return fmpq_cmp_fmpz(q, b.z);
  }
  int operator==(long b) const {
    return fmpq_cmp_si(q, b);
  }
  friend inline int operator==(const INT& a, const RAT& b) {return b==a;}
  friend inline int operator==(long a, const RAT& b) {return b==a;}
  friend inline int operator!=(const RAT& a, const RAT& b) {return !(a==b);}
  friend inline int operator!=(const RAT& a, const INT& b) {return !(a==b);}
  friend inline int operator!=(const INT& a, const RAT& b) {return b!=a;}
  friend inline int operator!=(long a, const RAT& b) {return b!=INT(a);}

  int operator<(const RAT& b) const {return fmpq_cmp(q, b.q)<0;}
  int operator<(const INT& b) const {return fmpq_cmp_fmpz(q, b.z)<0;}
  int operator<(long b) const {return fmpq_cmp_si(q, b)<0;}
  int operator>(const RAT& b) const {return fmpq_cmp(q, b.q)>0;}
  int operator>(const INT& b) const {return fmpq_cmp_fmpz(q, b.z)>0;}
  int operator>(long b) const {return fmpq_cmp_si(q, b)>0;}
  friend inline int operator<(const INT& a, const RAT& b) {return b>a;}
  friend inline int operator<(long a, const RAT& b) {return b>a;}
  friend inline int operator>(const INT& a, const RAT& b) {return b<a;}
  friend inline int operator>(long a, const RAT& b) {return b<a;}

  friend ostream& operator<< (ostream&, const RAT&);
  friend istream& operator>> (istream&, RAT&);

  void operator+=(const RAT& b) {fmpq_add(q, q, b.q);}
  void operator+=(const INT& b) {fmpq_add_fmpz(q, q, b.z);}
  void operator-=(const RAT& b) {fmpq_sub(q, q, b.q);}
  void operator-=(const INT& b) {fmpq_sub_fmpz(q, q, b.z);}
  void operator*=(const RAT& b) {fmpq_mul(q, q, b.q);}
  void operator*=(const INT& b) {fmpq_mul_fmpz(q, q, b.z);}
  void operator/=(const RAT& b) {fmpq_div(q, q, b.q);}
  void operator/=(const INT& b) {fmpq_div_fmpz(q, q, b.z);}

  RAT operator+() const {
    RAT x(*this); return x;
  }
  RAT operator-() const {
    RAT x(*this); fmpq_neg(x.q, x.q); return x;
  }
  RAT abs() const {
    RAT x(*this); fmpq_abs(x.q, x.q); return x;
  }
  INT floor() const;
  INT ceil() const;
};


// Definitions of non-member RAT functions


inline INT RAT::round() const {
  return rounded_division(num(), den());
}

inline ostream& operator<<(ostream& s, const RAT& q)
{
  INT n(q.num()), d(q.den());
  if (is_zero(d))
    s<<"oo";
  else
    {
      s << n;
      if (d!=1)
        s << "/" << d;
    }
  return s;
}

inline istream& operator>> (istream& is, RAT& r)
{
  std::string n;
  is>>n;
  fmpq_set_str(r.q, n.c_str(), 10);
  fmpq_canonicalise(r.q);
  return is;
}

inline INT RAT::floor() const
{
  INT q, r;
  divrem(num(), den(), q, r);
  return q;
}

inline INT RAT::ceil() const
{
  INT q, r;
  divrem(num(), den(), q, r);
  return (is_zero(r)? q : q+1);
}

#endif
