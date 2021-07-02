// FILE QUADS.H

#if     !defined(_QUADS_H)
#define _QUADS_H      1       //flags that this file has been included

#include <eclib/arith.h>
#include <assert.h>

#define PI M_PI

class Quad;
class RatQuad;

// Valid fields
extern vector<int> valid_fields;
int check_field(int d, vector<int> fields=valid_fields);

//functions assigned by Quad::field initializer
extern Quad (*quadconj)(const Quad& a);  //Can't have same names as Complex functions
extern Quad (*mult)(const Quad& a, const Quad& b);
extern Quad (*qdivi)(const Quad& a, long c);
extern int (*pos)(const Quad& a);

//GCD-related functions.
extern Quad (*quadgcd)(const Quad& aa, const Quad& bb);
extern Quad (*quadbezout)(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy);
Quad quadgcd1(const Quad& aa, const Quad& bb);   // Euclidean only
Quad quadgcd2(const Quad& aa, const Quad& bb);   // General
Quad quadgcd_psea(const Quad& aa, const Quad& bb);   // Using (pseudo-)EA
Quad quadbezout1(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy); //Euclidean
Quad quadbezout2(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy); //General
Quad quadbezout_psea(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy); // Using (pseudo-)EA
Quad invmod(const Quad& a, const Quad& p);
int coprime(const Quad& a, const Quad& b);
int invertible(const Quad& a, const Quad& b, Quad& inverse);

//functions defined in quads.cc
int ndiv(const Quad& a, const Quad& b);
vector<Quad> residues(const Quad& a);

//Primes
extern vector<Quad> quadprimes;
extern long nquadprimes;     //The number of them.
extern vector<Quad> quadunits, squareunits;
extern Quad fundunit;

class Quad {
/***********************************************************************
First the static members which only depend on the field.  To initialize,
include the line
        Quad::field(d,maxnorm);
at the top of the program, where d is a suitable positive integer
and maxnorm (default 1000) is the upper bound for the norms of primes.
***********************************************************************/
 public:
  static int     d;          // square-free >0
  static int     t;          // trace of w
  static int     n;          // norm of w
  static char    name;       // name of w for printing
  static int  maxnorm;       // largest norm of primes
  static int     disc;       // discriminant
  static int   nunits;       // number of units
  static int is_Euclidean;   // 1 for Euclidean fields, else 0
  static int is_class_number_one;   // 1 for fields of class number 1, else 0
  static Quad zero;
  static Quad one;
  static Quad w;

  static void field(int dd, int max=1000);
  static void displayfield(ostream& s = cout);
  static int chi(long p); // quadratic character associated to the field
  static void initquadprimes();
  static vector<Quad> primes_above(long p, int& sig);
  static void setup_geometry();

// Now the actual data elements:
 private:
  long r,i, nm; // nm is the cached norm
 public:
//constructors:
  void setnorm()
  {
    nm = r*r + n*i*i;
    if (t) {nm += r*i;};
  }
  Quad(long x=0, long y=0, long nrm=-1) :r(x),i(y), nm(nrm)
  {
    if (nm<0) setnorm();
  }
  Quad(const bigcomplex& z);   //rounds to nearest
  Quad(const Quad& a) :r(a.r), i(a.i), nm(a.nm) {;}
  Quad conj() const {return quadconj(*this);}
  long re() const {return r;}
  long im() const {return i;}
  long norm() const {return nm;}
  Quad pos_assoc() const {return makepos(*this);}

//operators and related functions (friends are inlined below):

  void operator=(const Quad& a) {r=a.r; i=a.i; nm=a.nm;}
  friend long real(const Quad& a) {return a.r;}
  friend long imag(const Quad& a) {return a.i;}
  friend Quad makepos(const Quad& a);
  friend long quadnorm(const Quad& a) {return a.nm;} // used when t=0
  friend Quad quadconj0(const Quad& a);
  friend Quad quadconj1(const Quad& a);
  friend Quad mult0(const Quad& a, const Quad& b);
  friend Quad mult1(const Quad& a, const Quad& b);
  friend Quad qdivi0(const Quad& a, long c);      // used when t=0
  friend Quad qdivi1(const Quad& a, long c);      // used when t=1
  friend int pos13(const Quad& a);
  friend int pos2(const Quad& a);
  friend int div(const Quad& a, const Quad& b);           // implemented in quads.cc
  friend int val(const Quad& factor, const Quad& number); // implemented in quads.cc
  friend Quad quadgcd1(const Quad&, const Quad&);   //Euclidean only
  friend Quad quadgcd2(const Quad&, const Quad&);   //General
  friend Quad quadgcd_psea(const Quad&, const Quad&);   // Using (pseudo-)EA
  friend Quad quadbezout_psea(const Quad&, const Quad&, Quad&, Quad&);   // Using (pseudo-)EA

  int operator== (const Quad& b) const {return (r==b.r) && (i==b.i);}
  int operator== (long b) const {return (r==b) && (i==0);}
  int operator!= (const Quad& b) const {return (r!=b.r) || (i!=b.i);}
  int operator!= (long b) const {return (r!=b) || (i!=0);}
  Quad operator* (const Quad& b) const {return mult(*this,b);}
  void operator*=(const Quad& b) {*this=mult(*this,b);}
  Quad operator* (long m) const {return Quad(m*r,m*i, m*m*nm);}
  void operator*=(long m) {r*=m;i*=m; nm*=(m*m);}
  friend Quad operator*(long m, const Quad& a);
  Quad operator+ (const Quad& b) const {return Quad(r+b.r,i+b.i);}
  Quad operator+ (long b) const {return Quad(r+b,i);}
  friend Quad operator+(long m, const Quad& a);
  void operator+=(const Quad& b) {r+=b.r; i+=b.i; setnorm();}
  void operator+=(long b) {r+=b;}
  Quad operator- (const Quad& b) const {return Quad(r-b.r,i-b.i);}
  Quad operator- (long b) const {return Quad(r-b,i);}
  friend Quad operator-(long m, const Quad& a);
  void operator-=(const Quad& b) {r-=b.r; i-=b.i; setnorm();}
  void operator-=(long b) {r-=b; setnorm();}
  Quad operator- () const {return Quad(-r,-i);}
  Quad operator/ (const Quad& b) const
  {
    if (b.i)
      return qdivi(mult(*this,quadconj(b)), b.nm);
    else
      return qdivi(*this,b.r);
  }
  Quad operator/ (long b) const {return qdivi(*this,b);}
  void operator/=(const Quad& b)
  {
    if (b.i)
      *this=qdivi(mult(*this,quadconj(b)), b.nm);
    else
      *this=qdivi(*this,b.r);
  }
  void operator/=(long b) {*this=qdivi(*this,b);}
  operator bigcomplex() const;

// iostream functions

  friend ostream& operator<<(ostream& s, const Quad& x);
  friend istream& operator>>(istream& s, Quad& x);

  friend class level;
  friend class RatQuad;
  friend void pseudo_euclidean_step(Quad&, Quad&, int&, Quad&, Quad&, Quad&, Quad&);
};

char* to_string(const Quad& a);  // outputs to a (new) string

// Inline definitions of friend functions of class Quad (those not
// here are not inline, and are in quads.cc):

inline Quad operator% (const Quad& a, const Quad& b)
{ return a-(b*(a/b));}
inline Quad makepos(const Quad& a)
{Quad ans=a;
  while(!pos(ans)) {ans*=fundunit; }
 return ans;
}
inline Quad quadconj0(const Quad& a)  {return Quad(a.r , -a.i, a.nm);}
inline Quad quadconj1(const Quad& a)  {return Quad(a.r + a.i, -a.i, a.nm);}
inline Quad mult0(const Quad& a, const Quad& b)
{
  if (b.i==0)
    return b.r * a;
  if (a.i==0)
    return a.r * b;
  return Quad(a.r*b.r-Quad::n*a.i*b.i, a.r*b.i+a.i*b.r, a.nm*b.nm);
}
inline Quad mult1(const Quad& a, const Quad& b)
{
  if (b.i==0)
    return b.r * a;
  if (a.i==0)
    return a.r * b;
  return Quad(a.r*b.r-Quad::n*a.i*b.i, a.r*b.i+a.i*b.r+a.i*b.i, a.nm*b.nm);
}
inline int pos13(const Quad& a)
   {return (((a.i>=0)&&(a.r>0))||((a.r==0)&&(a.i==0)));}
inline int pos2(const Quad& a)
   {return ((a.i>0)||((a.i==0)&&(a.r>=0)));}
inline Quad operator*(long m, const Quad& a) {return Quad(m*a.r,m*a.i, m*m*a.nm);}
inline Quad operator+(long m, const Quad& a) {return Quad(m+a.r,a.i);}
inline Quad operator-(long m, const Quad& a) {return Quad(m-a.r,-a.i);}
inline istream& operator>>(istream& s, Quad& x)
{
  s >> x.r >> x.i;
  x.setnorm();
  return s;
}

vector<long> findminquadcoeffs(const Quad&, const Quad&, Quad&, Quad&);
vector<long> findminquadcoeffs(const Quad&, const Quad&, Quad&);
void findminquad(const Quad&, const Quad&, Quad&, Quad&);
void findminquad(const Quad&, const Quad&, Quad&);


vector<long> HNF(const Quad& alpha);  // returns HNF of ideal (alpha)
string ideal_label(const Quad& alpha);  // returns label of ideal (alpha)
string field_label(); // returns field label, e.g. '2.0.4.1'

Quad primdiv(const Quad&);           // "First" prime divisor
vector<Quad> pdivs(const Quad&);         // list of prime divisors
vector<Quad> posdivs(const Quad&);       // all "positive" divisors (up to units)
vector<Quad> alldivs(const Quad&);       // absolutely all divisors
vector<Quad> sqdivs(const Quad&);        // returns divisors whose square divides
                                 // (up to +/- sign)
vector<Quad> sqfreedivs(const Quad&);    // returns square-free divisors

int are_associate(const Quad&, const Quad&);
int is_ideal_Galois_stable(const Quad&);

// some other general-purpose functions

#include <values.h>

inline double realnorm(const Quad& z) {  return sqrt(double(quadnorm(z)));}
inline double psif(bigcomplex z) {  return to_double(cos(4*PI*real(z)));}
inline double psig(bigcomplex z) {  return to_double(sin(4*PI*real(z)));}
int squaremod(const Quad& a, const Quad& m, const vector<Quad>& reslist);
vector<int> makechitable(const Quad& lambda, const vector<Quad>& reslist);
double gauss(const Quad& m, const vector<Quad>& reslist);

class mat22 {  //2x2 matrix for linear fractional transformations
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
  matop(const Quad& p, const Quad& n);  // constructor for hecke ops
  mat22 operator[](int i) const {return mats[i];}
  int length() const {return mats.size();}
};

#endif

// END OF FILE QUADS.H
