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
class SwanData;

// Valid fields
extern vector<long> valid_fields;
vector<long> valid_field_discs(long max_disc=0);
int check_field(long d, vector<long> fields=valid_fields);

//functions assigned by Quad::field initializer
extern Quad (*quadconj)(const Quad& a);  //Can't have same names as Complex functions
extern Quad (*mult)(const Quad& a, const Quad& b);
extern Quad (*qdivi)(const Quad& a, const INT& c);
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
vector<Quad> invertible_residues(const Quad& a);

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
  static vector<Quad> shifts_by_one;
  static void field(long dd, long max=1000);
  static void displayfield(ostream& s = cout, int info2=0); // if info2, also output info about 2-part of class group
  static int chi(const INT& p); // quadratic character associated to the field
  static void initquadprimes();
  static vector<Quad> primes_above(long p, int& sig);
  static void fill_class_group(); // implemented in qidloop.cc
  static SwanData SD;
  static int geometry_initialised; // set to 0 on init and to 1 when the following function has been called
  static void setup_geometry(string subdir="geodata", int debug=0);   // reads or creates SwanData object SD

// Now the actual data elements:
 private:
  INT r,i, nm; // nm is the cached norm
 public:
//constructors:
  void setnorm();
  Quad() :r(0), i(0), nm(0) {}
  explicit Quad(const INT& x) :r(x), i(0), nm(x*x)  {}
  explicit Quad(long x) :r(x), i(0), nm(x*x)  {}
  explicit Quad(int x) :r(x), i(0), nm(x*x)  {}
  explicit Quad(const INT& x, const INT& y) :r(x),i(y)  {setnorm();}
  explicit Quad(const INT& x, const INT& y, const INT& nrm) :r(x), i(y), nm(nrm)  {}
  explicit Quad(int x, int y) :r(x), i(y)  {setnorm();}
  explicit Quad(long x, long y) :r(x), i(y)  {setnorm();}
  //explicit Quad(const bigcomplex& z);   //rounds to nearest
  Quad(const Quad& a) :r(a.r), i(a.i), nm(a.nm) {}
  // racb = real part of (first times conjugate of second)
  friend INT racb(const Quad& a, const Quad& b)
  {
    INT c = fmma(a.r, b.r , a.i, b.i*n); // a.r*b.r + a.i*b.i*n;
    if (Quad::t) c += a.r*b.i;
    return c;
  }
  // iacb = imag part of (first times conjugate of second)
  friend INT iacb(const Quad& a, const Quad& b)
  {
    return fmms(a.i, b.r, b.i, a.r); // a.i*b.r - b.i*a.r
  }
  // a * conj(b)
  friend Quad mult_conj(const Quad& a, const Quad& b)
  {
    INT x = fmma(a.r, b.r , a.i, b.i*n), y = a.r*b.i;
    if (Quad::t) x += y;
    y = a.i*b.r - y;
    return Quad(x, y, a.nm*b.nm);
  }
  Quad conj() const {return quadconj(*this);}
  INT re() const {return r;}
  INT im() const {return i;}
  INT norm() const {
    // INT N = r*r+t*r*i+n*i*i;
    // if (N!=nm) cerr<<"In norm(): Quad "<<(*this)<<" with r = "<<r<<", i = "<<i
    //                <<" and norm "<<N<<" has nm = "<<nm<<endl;
    // assert(N==nm);
    return nm;
  }
  INT content() const {return gcd(r,i);}
  Quad pos_assoc() const {return makepos(*this);}
  int is_zero() const {return nm.is_zero();}
  int is_unit() const {return nm.is_one();}
  int is_one() const {return r.is_one() && i.is_zero();}

//operators and related functions (friends are inlined below):

  Quad& operator=(const Quad& a) {r=a.r; i=a.i; nm=a.nm;  return *this;}
  Quad& operator=(const INT& a) {r=a; i=0; nm=a*a;  return *this;}
  Quad& operator=(long a) {r=a; i=0; nm=a*a;  return *this;}
  Quad& operator=(int a) {r=a; i=0; nm=a*a;  return *this;}
  friend INT real(const Quad& a) {return a.r;}
  friend INT imag(const Quad& a) {return a.i;}
  friend Quad makepos(const Quad& a);
  friend Quad quadconj0(const Quad& a);
  friend Quad quadconj1(const Quad& a);
  friend Quad mult0(const Quad& a, const Quad& b);
  friend Quad mult1(const Quad& a, const Quad& b);
  friend Quad qdivi0(const Quad& a, const INT& c);      // used when t=0
  friend Quad qdivi1(const Quad& a, const INT& c);      // used when t=1
  friend int pos13(const Quad& a);
  friend int pos2(const Quad& a);
  friend int div(const Quad& a, const Quad& b);
  friend int div(const Quad& a, const Quad& b, Quad& quo);
  friend int val(const Quad& factor, const Quad& number);
  friend Quad quadgcd_psea(const Quad&, const Quad&);   // Using (pseudo-)EA
  friend Quad quadbezout_psea(const Quad&, const Quad&, Quad&, Quad&);   // Using (pseudo-)EA
  friend mat22 generalised_extended_euclid(const Quad& aa, const Quad& bb, int& s);

  // Replace alpha, beta by an SL2Z-equivalent pair with the same Z-span.
  // The new alpha has the smallest norm.
  // In the first form (not used), U holds the unimodular transform
  friend void sl2z_reduce(Quad& alpha, Quad& beta, unimod&U);
  friend void sl2z_reduce(Quad& alpha, Quad& beta);

  int operator== (const Quad& b) const {return (r==b.r) && (i==b.i);}
  int operator== (const INT& b) const {return (r==b) && i.is_zero();}
  int operator<(const Quad& b) const // just for sorting
  {
    return (nm<b.nm) || ((nm==b.nm) && ((i<b.i) || ((i==b.i) && (r<b.r))));
  }
  int operator!= (const Quad& b) const {return (r!=b.r) || (i!=b.i);}
  int operator!= (const INT& b) const {return (r!=b) || i.is_nonzero();}
  Quad operator* (const Quad& b) const {return mult(*this,b);}
  void operator*=(const Quad& b) {*this=mult(*this,b);}
  Quad operator* (const INT& m) const {return Quad(m*r,m*i, m*m*nm);}
  Quad operator* (long m) const {return Quad(m*r,m*i, m*m*nm);}
  Quad operator* (int m) const {return Quad(m*r,m*i, m*m*nm);}
  void operator*=(const INT& m) {r*=m; i*=m; nm*=(m*m); }
  void times_w() {INT u = -i*Quad::n; i = (Quad::t? r+i: r); r = u; nm*=Quad::n; }
  friend Quad operator*(const INT& m, const Quad& a);
  friend Quad operator*(long m, const Quad& a);
  friend Quad operator*(int m, const Quad& a);
  Quad operator+ (const Quad& b) const {return Quad(r+b.r,i+b.i);}
  Quad operator+ (const INT& b) const {return Quad(r+b,i);}
  Quad operator+ (long b) const {return Quad(r+b,i);}
  Quad operator+ (int b) const {return Quad(r+b,i);}
  friend Quad operator+(const INT& m, const Quad& a);
  friend Quad operator+(long m, const Quad& a);
  friend Quad operator+(int m, const Quad& a);
  void operator+=(const Quad& b) {r+=b.r; i+=b.i; setnorm();}
  void operator+=(const INT& b) {r+=b; setnorm();}
  void operator+=(long b) {r+=b; setnorm();}
  void operator+=(int b) {r+=b; setnorm();}
  Quad operator- (const Quad& b) const {return Quad(r-b.r,i-b.i);}
  Quad operator- (const INT& b) const {return Quad(r-b,i);}
  Quad operator- (long b) const {return Quad(r-b,i);}
  Quad operator- (int b) const {return Quad(r-b,i);}
  friend Quad operator-(const INT& m, const Quad& a);
  friend Quad operator-(long m, const Quad& a);
  friend Quad operator-(int m, const Quad& a);
  void operator-=(const Quad& b) {r-=b.r; i-=b.i; setnorm();}
  void operator-=(const INT& b) {r-=b; setnorm();}
  void operator-=(long b) {r-=b; setnorm();}
  void operator-=(int b) {r-=b; setnorm();}
  Quad operator- () const {return Quad(-r,-i);}
  Quad operator/ (const Quad& b) const;
  Quad operator/ (const INT& b) const {return qdivi(*this,b);}
  void operator/=(const Quad& b);
  void operator/=(const INT& b) {*this=qdivi(*this,b);}
  void addprod(const Quad& a, const Quad& b); // this +=a*b
  void addprod(long a, const Quad& b); // this +=a*b
  void subprod(const Quad& a, const Quad& b); // this -=a*b
  void subprod(long a, const Quad& b); // this -=a*b
  Quad operator% (long b) const { return Quad(INT(r%b), INT(i%b));}
  Quad operator% (const INT& b) const { return Quad(r%b, i%b);}
  //operator bigcomplex() const;
  // compute a*b+c*d
  friend Quad mma(const Quad& a, const Quad& b, const Quad& c, const Quad& d);
  // compute a*b-c*d
  friend Quad mms(const Quad& a, const Quad& b, const Quad& c, const Quad& d);

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

inline Quad operator% (const Quad& a, const Quad& b) { return a-(b*(a/b));}
inline Quad quadconj0(const Quad& a)  {return Quad(a.r , -a.i, a.nm);}
inline Quad quadconj1(const Quad& a)  {return Quad(a.r + a.i, -a.i, a.nm);}
inline int pos13(const Quad& a) {return ((sign(a.i)>=0&&sign(a.r)>0)||((sign(a.r)==0)&&(sign(a.i)==0)));}
inline int pos2(const Quad& a) {return (sign(a.i)>0||(sign(a.i)==0&&sign(a.r)>=0));}
inline Quad operator*(const INT& m, const Quad& a) {return Quad(m*a.r,m*a.i, m*m*a.nm);}
inline Quad operator*(long m, const Quad& a) {return Quad(m*a.r,m*a.i, m*m*a.nm);}
inline Quad operator*(int m, const Quad& a) {return Quad(m*a.r,m*a.i, m*m*a.nm);}
inline Quad operator+(const INT& m, const Quad& a) {return Quad(m+a.r,a.i);}
inline Quad operator+(long m, const Quad& a) {return Quad(m+a.r,a.i);}
inline Quad operator+(int m, const Quad& a) {return Quad(m+a.r,a.i);}
inline Quad operator-(const INT& m, const Quad& a) {return Quad(m-a.r,-a.i);}
inline Quad operator-(long m, const Quad& a) {return Quad(m-a.r,-a.i);}
inline Quad operator-(int m, const Quad& a) {return Quad(m-a.r,-a.i);}

// reduction modulo Z<alpha,beta>
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
vector<int> chardisc(const INT& D);

// some other general-purpose functions

#include <values.h>

int squaremod(const Quad& a, const Quad& m, const vector<Quad>& reslist);
vector<int> makechitable(const Quad& lambda, const vector<Quad>& reslist);
bigfloat gauss(const Quad& m, const vector<Quad>& reslist);
string ideal_code(const Quad& N); // string code for a (principal)  ideal

// all_disc_factors modulo D mod squares, i.e. factoring out D.  D
// should be in all_disc_factors.  Returns a list of half the length
// unless D=1.
vector<INT> disc_factors_mod_D(const INT& D);

// reduce r mod s so that r/s is in the rectangle
Quad rectify(const Quad& r, const Quad& s);

string eigfile(const Quad& N, long p=0);    //returns directory/filename for eigs at level N, characteristic p
string eigfile(Qideal& N, long p=0);        //returns directory/filename for eigs at level N, characteristic p

// return name of newforsm directory for this field; if
// create_if_necessary, creates the directory if it does not yet exist
string newforms_directory(int create_if_necessary=1);

#endif

// END OF FILE QUADS.H
