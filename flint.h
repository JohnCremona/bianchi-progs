// FILE FLINT.H: wrapper for FLINT's fmpz type

#if     !defined(_FLINT_H)
#define _FLINT_H      1       //flags that this file has been included

// for convenience we put all the includes from eclib here, since all
// other files include this one directly or indirectly:

#include <string>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>

class INT {
private:
  fmpz_t z;
public:
  INT() {fmpz_init(z);} // sets to 0
  INT(int a) {fmpz_init_set_si(z, a);}
  INT(long a) {fmpz_init_set_si(z, a); }
  INT(const INT& a) {fmpz_init_set(z, a.z);}
  INT& operator=(int a) {fmpz_set_si(z, a); return *this;}
  INT& operator=(long a) {fmpz_set_si(z, a); return *this;}
  INT& operator=(const INT& a) {fmpz_set(z, a.z); return *this;}

  int operator==(const INT& a) const {return fmpz_equal(z, a.z);}
  int operator==(int a) const {return fmpz_equal_si(z, a);}
  int operator==(long a) const {return fmpz_equal_si(z, a);}
  int operator!=(const INT& a) const {return !fmpz_equal(z, a.z);}
  int operator!=(int a) const {return !fmpz_equal_si(z, a);}
  int operator!=(long a) const {return !fmpz_equal_si(z, a);}

  int sign() const {return fmpz_sgn(z);}
  INT operator-() const {INT b; fmpz_set(b.z,z); fmpz_neg(b.z,b.z); return b;}
  INT abs() const {INT b; fmpz_abs(b.z,z); return b;}
  int is_long() const {return fmpz_fits_si(z);}

  INT operator+(int a) const {INT b; fmpz_set(b.z,z); fmpz_add_si(b.z,b.z,a); return b;}
  INT operator+(long a) const {INT b; fmpz_set(b.z,z); fmpz_add_si(b.z,b.z,a); return b;}
  INT operator+(const INT& a) const {INT b; fmpz_set(b.z,z); fmpz_add(b.z,b.z,a.z); return b;}
  INT operator-(int a) const {INT b; fmpz_set(b.z,z); fmpz_sub_si(b.z,b.z,a); return b;}
  INT operator-(long a) const {INT b; fmpz_set(b.z,z); fmpz_sub_si(b.z,b.z,a); return b;}
  INT operator-(const INT& a) const {INT b; fmpz_set(b.z,z); fmpz_sub(b.z,b.z,a.z); return b;}
  INT operator*(int a) const {INT b; fmpz_set(b.z,z); fmpz_mul_si(b.z,b.z,a); return b;}
  INT operator*(long a) const {INT b; fmpz_set(b.z,z); fmpz_mul_si(b.z,b.z,a); return b;}
  INT operator*(const INT& a) const {INT b; fmpz_set(b.z,z); fmpz_mul(b.z,b.z,a.z); return b;}
  INT operator/(const INT& a) const {INT b; fmpz_set(b.z,z); fmpz_divexact(b.z,b.z,a.z); return b;}
  INT operator/(int a) const {INT b; fmpz_set(b.z,z); fmpz_divexact_si(b.z,b.z,a); return b;}
  INT operator/(long a) const {INT b; fmpz_set(b.z,z); fmpz_divexact_si(b.z,b.z,a); return b;}
  INT operator%(const INT& a) const {INT b; fmpz_mod(b.z,z,a.z); return b;}
  int operator%(int a) const {INT b; fmpz_mod_ui(b.z,z,a); return (int)I2long(b);}  // a must be >0
  long operator%(long a) const {INT b; fmpz_mod_ui(b.z,z,a); return I2long(b);} // a must be >0
  INT operator^(int e) const {INT b; fmpz_pow_ui(b.z,z,e); return b;}
  INT operator^(long e) const {INT b; fmpz_pow_ui(b.z,z,e); return b;}
  void operator +=(const INT& a) {fmpz_add(z,z,a.z);}
  void operator +=(int a) {fmpz_add_si(z,z,a);}
  void operator +=(long a) {fmpz_add_si(z,z,a);}
  void operator -=(const INT& a) {fmpz_sub(z,z,a.z);}
  void operator -=(int a) {fmpz_sub_si(z,z,a);}
  void operator -=(long a) {fmpz_sub_si(z,z,a);}
  void operator *=(const INT& a) {fmpz_mul(z,z,a.z);}
  void operator *=(int a) {fmpz_mul_si(z,z,a);}
  void operator *=(long a) {fmpz_mul_si(z,z,a);}
  void operator /=(const INT& a) {fmpz_divexact(z,z,a.z);}
  void operator /=(int a) {fmpz_divexact_si(z,z,a);}
  void operator /=(long a) {fmpz_divexact_si(z,z,a);}
  int is_square() const {return fmpz_is_square(z);}
  int is_square(INT& a) const
  {
    if (!fmpz_is_square(z)) return 0;
    else {fmpz_sqrt(a.z,z); return 1;}
  }
  INT isqrt() const {INT b; fmpz_sqrt(b.z,z); return b;}
  long valuation(const INT& p) {INT q;   return fmpz_remove(q.z, z, p.z);}

  string to_string() const
  {
    return string(fmpz_get_str(NULL, 10, z));
  }
  friend ostream& operator<<(ostream& s, const INT& a)
  {
    s << a.to_string();
    return s;
  }
  friend istream& operator>>(istream& s, INT& x)
  {
    string n;
    s>>n;
    fmpz_set_str(x.z, n.c_str(), 10);
    return s;
  }
  friend INT fmma(const INT& a, const INT& b, const INT& c, const INT& d); // a*b+c*d
  friend INT fmms(const INT& a, const INT& b, const INT& c, const INT& d); // a*b-c*d
  friend INT mod(const INT& a, const INT& b); // a mod b in range (-b/2,b/2]
  friend INT gcd(const INT& a, const INT& b);
  friend INT bezout(const INT& a, const INT& b, INT& x, INT& y);
  friend int compare(const INT& a, const INT& b);
  friend int compare(const INT& a, int b);
  friend int compare(const INT& a, long b);
  friend void divrem(const INT& a, const INT& b, INT& quo, INT& rem);
  friend vector<INT> pdivs(const INT& a);
  friend vector<INT> sqdivs(const INT& a);
  friend int legendre(const INT& a, const INT& p) {return kronecker(a,p);}
  friend int kronecker(const INT& a, const INT& n);
  friend INT invmod(const INT&a, const INT& p);
  friend long invmod(const INT&a, long p);
  friend long I2long(const INT& a);
// return e and set q, where a=q*f^e, e maximal
  friend long divide_out(const INT&a, const INT& f, INT& q);
  friend void sqrt_mod_p(INT& b, const INT& a, const INT& p);
  friend void swap(INT& a, INT& b);
  };

inline int sign(const INT& a) {return a.sign();}
inline INT operator+(int a, const INT& b) {return b+a;}
inline INT operator+(long a, const INT& b) {return b+a;}
inline INT operator-(int a, const INT& b) {return -(b-a);}
inline INT operator-(long a, const INT& b) {return -(b-a);}
inline INT operator*(int a, const INT& b) {return b*a;}
inline INT operator*(long a, const INT& b) {return b*a;}

inline INT fmma(const INT& a, const INT& b, const INT& c, const INT& d) {INT e; fmpz_fmma(e.z, a.z, b.z, c.z, d.z); return e;}
inline INT fmms(const INT& a, const INT& b, const INT& c, const INT& d) {INT e; fmpz_fmms(e.z, a.z, b.z, c.z, d.z); return e;}

inline INT abs(const INT& a) {return a.abs();}
inline INT posmod(const INT& a, const INT& b) {return a%b;} // in range [0,|b|)
inline INT posmod(const INT& a, long b) {return a%b;} // in range [0,|b|)
inline INT mod(const INT& a, const INT& b) {INT c; fmpz_smod(c.z,a.z,b.z); return c;} // in range (-|b|/2,|b|/2-1]
inline INT gcd(const INT& a, const INT& b) {INT c; fmpz_gcd(c.z,a.z,b.z); return c;}
inline INT bezout(const INT& a, const INT& b, INT& x, INT& y) {INT c; fmpz_xgcd_canonical_bezout(c.z, x.z, y.z, a.z, b.z); return c;}
inline int compare(const INT& a, const INT& b) {return fmpz_cmp(a.z,b.z);}
inline int compare(const INT& a, int b) {return fmpz_cmp_si(a.z,b);}
inline int compare(const INT& a, long b) {return fmpz_cmp_si(a.z,b);}
inline int operator==(int a, const INT& b) {return b==a;}
inline int operator==(long a, const INT& b) {return b==a;}
inline int is_zero(const INT& a) {return a.sign()==0;}
inline int operator<(const INT& a, const INT& b) {return compare(a,b)<0;}
inline int operator>(const INT& a, const INT& b) {return compare(a,b)>0;}
inline int operator<=(const INT& a, const INT& b) {return compare(a,b)<=0;}
inline int operator>=(const INT& a, const INT& b) {return compare(a,b)>=0;}

// ndiv rounds the quotient such that the remainder has the smallest
// absolute value. In case of ties, it rounds the quotient towards zero.
inline void divrem(const INT& a, const INT& b, INT& quo, INT& rem) {fmpz_ndiv_qr(quo.z, rem.z, a.z, b.z);}

inline int divides(const INT& a, const INT& b) {return (b%a)==0;}
inline int divides(int a, const INT& b) {return (b%a)==0;}
inline int divides(long a, const INT& b) {return (b%a)==0;}

inline INT rounded_division(const INT& a, const INT& b)
{
  INT q, r;
  divrem(a,b,q,r);
  INT r2 = 2*r;
  // We want -b <= r2 < +b
  if (r2<-b)
    q-=1;
  else
    if (r2>=b)
      q+=1;
  return q;
}

inline vector<INT> pdivs(const INT& a)
{
  fmpz_factor_t f;
  fmpz_factor_init(f);
  fmpz_factor(f, a.z);
  vector<INT> ans;
  for (int i =0; i< f->num; i++)
    {
      INT p;
      fmpz_init_set(p.z, f->p + i);
      ans.push_back(p);
    }
  return ans;
}

inline vector<INT> sqdivs(const INT& a)
{
  fmpz_factor_t f;
  fmpz_factor_init(f);
  fmpz_factor(f, a.z);
  vector<INT> plist;
  vector<int> elist;
  int nd = 1;
  for (slong i =0; i< f->num; i++)
    {
      INT p;
      fmpz_init_set(p.z, f->p + i);
      plist.push_back(p);
      int e = (f->exp[i])/2;
      elist.push_back(e);
      nd *= (1+e);
    }
  vector<INT> dlist(1, INT(1));
  dlist.resize(nd);
  nd = 1;
  auto pr = plist.begin();
  auto ei = elist.begin();
  while(pr!=plist.end())
    {
      INT p=*pr++;
      int e=*ei++;
      for (int j=0; j<e; j++)
        for (int k=0; k<nd; k++)
          dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
      nd*=(e+1);
    }
  return dlist;
}

inline int kronecker(const INT& a, const INT& p)
{
  return fmpz_kronecker(a.z, p.z);
}

inline int legendre(const INT& a, int p)
{
  return legendre(a,INT(p));
}

inline int legendre(const INT& a, long p)
{
  return legendre(a,INT(p));
}

inline INT invmod(const INT&a, const INT& p)
{
  INT b;
  fmpz_invmod(b.z, a.z, p.z);
  return b;
}

inline long invmod(const INT&a, long p)
{
  fmpz_t b, P;
  fmpz_init_set_si(P,p);
  fmpz_invmod(b, a.z, P);
  return (long)fmpz_get_si(b);
}

inline long I2long(const INT& a)
{
  if (a.is_long())
    return (long)fmpz_get_si(a.z);
  assert(0 && "INT does not fint into an slong");
}

// return e and set q, where a=q*f^e, e maximal
inline long divide_out(const INT&a, const INT& f, INT& q)
{
  return fmpz_remove(q.z, a.z, f.z);
}

// return valuation of a at f
inline long val(const INT&f, const INT& a)
{
  INT q; // to be discarded
  return divide_out(a, f, q);
}

inline void sqrt_mod_p(INT& b, const INT& a, const INT& p)
{
  fmpz_sqrtmod(b.z, a.z, p.z);
}

inline INT sqrt_mod_p(const INT& a, const INT& p)
{
  INT b;
  sqrt_mod_p(b, a, p);
  return b;
}

inline void swap(INT& a, INT& b)
{
  fmpz_swap(a.z, b.z);
}

inline INT pow(const INT& a, int e) {return a^e;}
inline INT pow(const INT& a, long e) {return a^e;}

#endif
