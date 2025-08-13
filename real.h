// FILE REAL.H: minimal wrapper for FLINT's fmpz type

#if     !defined(_real_H)
#define _REAL_H      1       //flags that this file has been included

#include <iostream>
#include <vector>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/arb.h>
#include <flint/arb_hypgeom.h>

#define ARB_PREC 50
#define ARB_OUTPUT_DIGITS 15

class REAL {
private:
  arb_t x;
public:
  // Constructors:
  REAL() {arb_init(x);} // sets to 0
  explicit REAL(int a) {arb_init(x); arb_set_si(x, a);}
  explicit REAL(long a) {arb_init(x); arb_set_si(x, a); }
  explicit REAL(const INT& a) {arb_init(x); arb_set_fmpz(x, a.z); }
  REAL(const REAL& a) {arb_init(x); arb_set(x, a.x);}
  explicit REAL(arb_t a) {arb_init(x); arb_set(x, a);}
  // Destructor:
  ~REAL() {arb_clear(x);}
  REAL& operator=(int a) {arb_init(x); arb_set_si(x, a); return *this;}
  REAL& operator=(long a) {arb_init(x); arb_set_si(x, a); return *this;}
  REAL& operator=(const REAL& a) {arb_init(x); arb_set(x, a.x); return *this;}
  int operator==(const REAL& a) const {return arb_equal(x, a.x);}
  int operator==(int a) const {return arb_equal_si(x, a);}
  int operator==(long a) const {return arb_equal_si(x, a);}
  int operator!=(const REAL& a) const {return !arb_equal(x, a.x);}
  int operator!=(int a) const {return !arb_equal_si(x, a);}
  int operator!=(long a) const {return !arb_equal_si(x, a);}
  REAL operator-() const {REAL b; arb_neg(b.x,x); return b;}
  REAL abs() const {REAL b; arb_abs(b.x,x); return b;}
  REAL inv() const {REAL b; arb_inv(b.x,x,ARB_PREC); return b;}
  REAL operator+(const REAL& a) const {REAL b; arb_add(b.x,x,a.x,ARB_PREC); return b;}
  REAL operator+(int a) const {REAL b; arb_add_si(b.x,x,a,ARB_PREC); return b;}
  REAL operator+(long a) const {REAL b; arb_add_si(b.x,x,a,ARB_PREC); return b;}
  REAL operator-(const REAL& a) const {REAL b; arb_sub(b.x,x,a.x,ARB_PREC); return b;}
  REAL operator-(int a) const {REAL b; arb_sub_si(b.x,x,a,ARB_PREC); return b;}
  REAL operator-(long a) const {REAL b; arb_sub_si(b.x,x,a,ARB_PREC); return b;}
  REAL operator*(int a) const {REAL b; arb_mul_si(b.x,x,a,ARB_PREC); return b;}
  REAL operator*(long a) const {REAL b; arb_mul_si(b.x,x,a,ARB_PREC); return b;}
  REAL operator*(const INT& a) const {REAL b; arb_mul_fmpz(b.x,x,a.z,ARB_PREC); return b;}
  REAL operator*(const REAL& a) const {REAL b; arb_mul(b.x,x,a.x,ARB_PREC); return b;}
  REAL operator/(const REAL& a) const {REAL b; arb_div(b.x,x,a.x,ARB_PREC); return b;}
  REAL operator/(int a) const {REAL b; arb_div_si(b.x,x,a,ARB_PREC); return b;}
  REAL operator/(long a) const {REAL b; arb_div_si(b.x,x,a,ARB_PREC); return b;}
  void operator +=(const REAL& a) {arb_add(x,x,a.x,ARB_PREC);}
  void operator +=(int a) {arb_add_si(x,x,a,ARB_PREC);}
  void operator +=(long a) {arb_add_si(x,x,a,ARB_PREC);}
  void operator -=(const REAL& a) {arb_sub(x,x,a.x,ARB_PREC);}
  void operator -=(int a) {arb_sub_si(x,x,a,ARB_PREC);}
  void operator -=(long a) {arb_sub_si(x,x,a,ARB_PREC);}
  void operator *=(const REAL& a) {arb_mul(x,x,a.x,ARB_PREC);}
  void operator *=(int a) {arb_mul_si(x,x,a,ARB_PREC);}
  void operator *=(long a) {arb_mul_si(x,x,a,ARB_PREC);}
  void operator /=(const REAL& a) {arb_div(x,x,a.x,ARB_PREC);}
  void operator /=(int a) {arb_div_si(x,x,a,ARB_PREC);}
  void operator /=(long a) {arb_div_si(x,x,a,ARB_PREC);}
  int is_zero() const {return arb_is_zero(x);}
  int is_nonzero() const {return !arb_is_zero(x);}
  REAL sqrt() const {REAL b; arb_sqrt(b.x,x,ARB_PREC); return b;}
  friend std::ostream& operator<<(std::ostream& s, const REAL& a);
  friend std::istream& operator>>(std::istream& s, REAL& x);
  friend void swap(REAL& a, REAL& b);
  friend REAL hypot(const REAL& x, const REAL& y);
  friend REAL KBessel(int nu, const REAL& x);
  friend REAL cos(const REAL& x);
  friend REAL sin(const REAL& x);

  friend class INT;
  friend class RAT;

  static REAL Pi() {REAL b; arb_const_pi(b.x,ARB_PREC); return b;}
};

inline REAL operator+(int a, const REAL& b) {return b+a;}
inline REAL operator+(long a, const REAL& b) {return b+a;}
inline REAL operator-(int a, const REAL& b) {return -(b-a);}
inline REAL operator-(long a, const REAL& b) {return -(b-a);}
inline REAL operator*(int a, const REAL& b) {return b*a;}
inline REAL operator*(long a, const REAL& b) {return b*a;}
inline REAL abs(const REAL& a) {return a.abs();}
inline int operator==(int a, const REAL& b) {return b==a;}
inline int operator==(long a, const REAL& b) {return b==a;}
inline void swap(REAL& a, REAL& b) {arb_swap(a.x, b.x);}
inline int is_zero(const REAL& a) {return a.is_zero();}
inline int is_nonzero(const REAL& a) {return a.is_nonzero();}
inline REAL sqrt(const REAL& a) {return a.sqrt();}
inline REAL sqrt(const INT& a) {return REAL(a).sqrt();}
inline REAL hypot(const REAL& x, const REAL& y) {REAL a; arb_hypot(a.x, x.x, y.x, ARB_PREC); return a;}
inline REAL KBessel(int nu, const REAL& x) {REAL a; arb_hypgeom_bessel_k(a.x, REAL(nu).x, x.x, ARB_PREC); return a;}
inline REAL cos(const REAL& x) {REAL b; arb_cos(b.x,x.x,ARB_PREC); return b;}
inline REAL sin(const REAL& x) {REAL b; arb_sin(b.x,x.x,ARB_PREC); return b;}
inline REAL K0(const REAL& x) {return KBessel(0,x);}
inline REAL K1(const REAL& x) {return KBessel(1,x);}

inline std::ostream& operator<<(std::ostream& s, const REAL& a)
{
  char* st = arb_get_str(a.x, ARB_OUTPUT_DIGITS, 0);
  s << std::string(st);
  flint_free(st);
  return s;
}

inline std::istream& operator>>(std::istream& s, REAL& x)
{
  std::string n;
  s>>n;
  arb_set_str(x.x, n.c_str(), ARB_PREC);
  return s;
}

#endif
