// FILE QUADS.H

#if     !defined(_QUADS_H)
#define _QUADS_H      1       //flags that this file has been included

#include <assert.h>
#include "intprocs.h" // typedefs INT

#define PI M_PI

class Quad;
class mat22;
class RatQuad;
class Qideal;

// Valid fields
extern vector<long> valid_fields;
vector<long> valid_field_discs(long max_disc=0);
int check_field(long d, vector<long> fields=valid_fields);

//functions assigned by Quad::field initializer
extern Quad (*quadconj)(const Quad& a);  //Can't have same names as Complex functions
extern Quad (*mult)(const Quad& a, const Quad& b);
extern Quad (*qdivi)(const Quad& a, const INT c);
extern int (*pos)(const Quad& a);

//GCD-related functions.
extern Quad (*quadgcd)(const Quad& aa, const Quad& bb);
extern Quad (*quadbezout)(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy);
Quad quadgcd_default(const Quad& aa, const Quad& bb);   // Using ideals
Quad quadbezout_default(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy);   // Using ideals
Quad quadgcd_psea(const Quad& aa, const Quad& bb);   // Using (pseudo-)EA
Quad quadbezout_psea(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy); // Using (pseudo-)EA
mat22 generalised_extended_euclid(const Quad& aa, const Quad& bb, int& s);

Quad invmod(const Quad& a, const Quad& p);
int coprime(const Quad& a, const Quad& b);
int principal_gcd(const Quad& a, const Quad& b, Quad& g);
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
  static long      d;          // square-free >0
  static long      t;          // trace of w (= 0 or 1)
  static INT     n;          // norm of w
  static char      name;       // name of w for printing
  static long      maxnorm;       // largest norm of primes
  static INT     disc, absdisc;       // discriminant (<0) and its absolute value
  static vector<INT> prime_disc_factors;   // n2r prime discriminant factors of discriminant
  static vector<INT> all_disc_factors;     // all 2^n2r discriminant factors of discriminant
  static int   nunits;       // number of units
  static int is_Euclidean;   // 1 for Euclidean fields, else 0
  static int class_number;   // the class number
  static int class_group_2_rank; // 2-rank of class group (= 1 less than length of prime_disc_factors)
  static vector<Qideal> class_group;                // list of ideals representing ideal classes (no structure)
  static vector<Qideal> class_group_2_torsion;      // list of ideals representing 2-torsion in class group
  static vector<Qideal> class_group_2_torsion_gens; // list of ideals generating 2-torsion in class group
  static vector<Qideal> class_group_2_cotorsion;      // list of ideals representing class group mod squares
  static vector<Qideal> class_group_2_cotorsion_gens; // list of ideals generating class group mod squares
  static int ideal_class_mod_squares(const Qideal& I); // binary index of [I] from 0 to 2^{2-rank}-1
  static int unramified_character(int i, const Qideal& I); // image (+1/-1) of I under i'th quadratic character
  static Quad zero;
  static Quad one;
  static Quad w;
  static void field(long dd, long max=1000);
  static void displayfield(ostream& s = cout, int info2=0); // if info2, also output info about 2-part of class group
  static int chi(INT p); // quadratic character associated to the field
  static void initquadprimes();
  static vector<Quad> primes_above(long p, int& sig);
  static void fill_class_group(); // implemented in qidloop.cc
  static int geometry_initialised; // set to 0 on init and to 1 when the following function has been called
  static void setup_geometry();   // implemented in geometry.cc

// Now the actual data elements:
 private:
  INT r,i, nm; // nm is the cached norm
 public:
//constructors:
  void setnorm();
  Quad() :r((long)0), i((long)0), nm((long)0)  {}
  Quad(INT x) :r(INT(x)), i(INT(0)), nm(INT(x)*x)  {}
  Quad(INT x, INT y) :r(INT(x)),i(INT(y))
  {
    nm = r*r + n*i*i;
    if (t) {nm += r*i;};
    assert (nm>=0);
  }
  Quad(INT x, INT y, INT nrm) :r(x), i(y), nm(nrm)  {}
  explicit Quad(const bigcomplex& z);   //rounds to nearest
  Quad(const Quad& a) :r(a.r), i(a.i), nm(a.nm) {}
  Quad conj() const {return quadconj(*this);}
  INT re() const {return r;}
  INT im() const {return i;}
  INT norm() const {return nm;}
  INT content() const {return gcd(r,i);}
  Quad pos_assoc() const {return makepos(*this);}
  int is_zero() const {return ::is_zero(nm);}

//operators and related functions (friends are inlined below):

  Quad& operator=(const Quad& a) {r=a.r; i=a.i; nm=a.nm; return *this;}
  friend INT real(const Quad& a) {return a.r;}
  friend INT imag(const Quad& a) {return a.i;}
  friend Quad makepos(const Quad& a);
  friend INT quadnorm(const Quad& a) {return a.nm;}
  friend Quad quadconj0(const Quad& a);
  friend Quad quadconj1(const Quad& a);
  friend Quad mult0(const Quad& a, const Quad& b);
  friend Quad mult1(const Quad& a, const Quad& b);
  friend Quad qdivi0(const Quad& a, INT c);      // used when t=0
  friend Quad qdivi1(const Quad& a, INT c);      // used when t=1
  friend int pos13(const Quad& a);
  friend int pos2(const Quad& a);
  friend int div(const Quad& a, const Quad& b);           // implemented in quads.cc
  friend int val(const Quad& factor, const Quad& number); // implemented in quads.cc
  friend void pseudo_euclidean_step(Quad&, Quad&, int&, Quad&, Quad&, Quad&, Quad&);
  friend Quad quadgcd_psea(const Quad&, const Quad&);   // Using (pseudo-)EA
  friend Quad quadbezout_psea(const Quad&, const Quad&, Quad&, Quad&);   // Using (pseudo-)EA
  friend mat22 generalised_extended_euclid(const Quad& aa, const Quad& bb, int& s);
// Replace alpha, beta by an SL2Z-equivalent pair with the same Z-span.
// The new alpha has the smallest norm.
// In the first form (not used), U holds the unimodular transform
  friend void sl2z_reduce(Quad& alpha, Quad& beta, unimod&U);
  friend void sl2z_reduce(Quad& alpha, Quad& beta);

  int operator== (const Quad& b) const {return (r==b.r) && (i==b.i);}
  int operator== (const INT b) const {return (r==b) && (i==0);}
  int operator!= (const Quad& b) const {return (r!=b.r) || (i!=b.i);}
  int operator!= (const INT b) const {return (r!=b) || (i!=0);}
  Quad operator* (const Quad& b) const {return mult(*this,b);}
  void operator*=(const Quad& b) {*this=mult(*this,b);}
  Quad operator* (INT m) const {return Quad(m*r,m*i, m*m*nm);}
  void operator*=(INT m) {r*=m;i*=m; nm*=(m*m);}
  friend Quad operator*(INT m, const Quad& a);
  Quad operator+ (const Quad& b) const {return Quad(r+b.r,i+b.i);}
  Quad operator+ (INT b) const {return Quad(r+b,i);}
  friend Quad operator+(INT m, const Quad& a);
  void operator+=(const Quad& b) {r+=b.r; i+=b.i; setnorm();}
  void operator+=(INT b) {r+=b;}
  Quad operator- (const Quad& b) const {return Quad(r-b.r,i-b.i);}
  Quad operator- (INT b) const {return Quad(r-b,i);}
  friend Quad operator-(INT m, const Quad& a);
  void operator-=(const Quad& b) {r-=b.r; i-=b.i; setnorm();}
  void operator-=(INT b) {r-=b; setnorm();}
  Quad operator- () const {return Quad(-r,-i);}
  Quad operator/ (const Quad& b) const;
  Quad operator/ (INT b) const {return qdivi(*this,b);}
  void operator/=(const Quad& b);
  void operator/=(INT b) {*this=qdivi(*this,b);}
  Quad operator% (long b) { return Quad(INT(r%b), INT(i%b));}
  operator bigcomplex() const;

// iostream functions

  operator string() const;
  friend ostream& operator<<(ostream& s, const Quad& x);
  friend istream& operator>>(istream& s, Quad& x);

  friend class level;
  friend class RatQuad;
  friend class Qideal;
};

// Inline definitions of friend functions of class Quad (those not
// here are not inline, and are in quads.cc):

inline Quad operator% (const Quad& a, const Quad& b)
{ return a-(b*(a/b));}
inline Quad quadconj0(const Quad& a)  {return Quad(a.r , -a.i, a.nm);}
inline Quad quadconj1(const Quad& a)  {return Quad(a.r + a.i, -a.i, a.nm);}
inline int pos13(const Quad& a)
{return ((sign(a.i)>=0&&sign(a.r)>0)||((sign(a.r)==0)&&(sign(a.i)==0)));}
inline int pos2(const Quad& a)
{return (sign(a.i)>0||(sign(a.i)==0&&sign(a.r)>=0));}
inline Quad operator*(INT m, const Quad& a) {return Quad(m*a.r,m*a.i, m*m*a.nm);}
inline Quad operator+(INT m, const Quad& a) {return Quad(m+a.r,a.i);}
inline Quad operator-(INT m, const Quad& a) {return Quad(m-a.r,-a.i);}

// reduction of gamma modulo Z<alpha,beta>
Quad reduce_mod_zbasis(const Quad& gamma, const Quad& alpha, const Quad& beta);

// vector<INT> findminquadcoeffs(const Quad&, const Quad&, Quad&, Quad&);
// vector<INT> findminquadcoeffs(const Quad&, const Quad&, Quad&);
// void findminquad(const Quad&, const Quad&, Quad&, Quad&);
// void findminquad(const Quad&, const Quad&, Quad&);


vector<INT> HNF(const Quad& alpha);  // returns HNF of ideal (alpha)
string old_ideal_label(const Quad& alpha);  // returns old HNF label of ideal (alpha)
string ideal_label(const Quad& alpha);  // returns new N.i label of ideal (alpha)
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

// convert a binary vector into a discriminant dividing Quad::disc
INT discchar(vector<int> c);
// convert a discriminant dividing Quad::disc into a binary vector
vector<int> chardisc(INT D);

// some other general-purpose functions

#include <values.h>

//inline bigfloat realnorm(const Quad& z) {  return sqrt(to_bigfloat(quadnorm(z)));}
inline bigfloat psif(bigcomplex z) {  return cos(4*PI*real(z));}
inline bigfloat psig(bigcomplex z) {  return sin(4*PI*real(z));}

int squaremod(const Quad& a, const Quad& m, const vector<Quad>& reslist);
vector<int> makechitable(const Quad& lambda, const vector<Quad>& reslist);
bigfloat gauss(const Quad& m, const vector<Quad>& reslist);
string ideal_code(const Quad& N); // string code for a (principal)  ideal

// all_disc_factors modulo D mod squares, i.e. factoring out D.  D
// should be in all_disc_factors.  Returns a list of half the length
// unless D=1.
vector<INT> disc_factors_mod_D(const INT& D);

#endif

// END OF FILE QUADS.H
