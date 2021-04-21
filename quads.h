// FILE QUADS.H

#if     !defined(_QUADS_H)
#define _QUADS_H      1       //flags that this file has been included

#include <eclib/interface.h>
#define PI M_PI

long roundover(long aa, long bb);

class Quad;

//functions assigned by Quad::field initializer
extern Quad (*quadconj)(const Quad& a);  //Can't have same names as Complex functions
extern long (*quadnorm)(const Quad& a);  //
extern Quad (*mult)(const Quad& a, const Quad& b);
extern Quad (*qdivi)(const Quad& a, long c);
extern int (*pos)(const Quad& a);
//GCD-related functions.
extern Quad (*quadgcd)(const Quad& aa, const Quad& bb);
Quad quadgcd1(const Quad& aa, const Quad& bb); //Euclidean only
Quad quadgcd2(const Quad& aa, const Quad& bb); //General
extern Quad (*quadbezout)(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy);
Quad quadbezout1(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy); //Euclidean
Quad quadbezout2(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy); //General
Quad invmod(const Quad& a, const Quad& p);
int coprime(const Quad& a, const Quad& b);
int invertible(const Quad& a, const Quad& b, Quad& inverse);
 

//functions defined in quads.cc
int div(const Quad& a, const Quad& b);
int ndiv(const Quad& a, const Quad& b);
int val(const Quad& factor, const Quad& number);
vector<Quad> residues(const Quad& a);

//Primes
extern vector<Quad> quadprimes;  
extern long nquadprimes;     //The number of them.
extern vector<Quad> quadunits;
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
 public:
  static int     disc;       // discriminant
  static int   nunits;       // number of units

public:
  static void field(int dd, int max=1000);
  static void displayfield(ostream& s = cout);
  static void initquadprimes();
  static vector<Quad> primes_above(long p, int& sig);

// Now the actual data elements:
 private:
  long r,i;
 public:
//constructors:
  Quad(long x=0, long y=0) :r(x),i(y) {;}
  Quad(const bigcomplex& z);   //rounds to nearest
  Quad(const Quad& a) :r(a.r), i(a.i) {;}
//  Quad(const Quadvar& qv) {*this=*(qv.values);}

//operators and related functions (friends are inlined below):

  void operator=(const Quad& a) {r=a.r; i=a.i;}
//  void operator=(const Quadvar& qv) {*this=*(qv.values);}
  friend long real(const Quad& a);
  friend long imag(const Quad& a);
  friend long quadnorm0(const Quad& a);      // used when t=0
  friend long quadnorm1(const Quad& a);      // used when t=1
  friend Quad quadconj0(const Quad& a);
  friend Quad quadconj1(const Quad& a);
  friend Quad mult0(const Quad& a, const Quad& b);
  friend Quad mult1(const Quad& a, const Quad& b);
  friend Quad qdivi0(const Quad& a, long c);
  friend Quad qdivi1(const Quad& a, long c);      // used when t=1
  friend int pos13(const Quad& a);
  friend int pos2(const Quad& a);

  int operator== (const Quad& b) const {return (r==b.r) && (i==b.i);}
  int operator== (long b) const {return (r==b) && (i==0);}
  int operator!= (const Quad& b) const {return (r!=b.r) || (i!=b.i);}
  int operator!= (long b) const {return (r!=b) || (i!=0);}
  Quad operator* (const Quad& b) const {return mult(*this,b);}
  void operator*=(const Quad& b) {*this=mult(*this,b);}
  Quad operator* (long m) const {return Quad(m*r,m*i);}
  void operator*=(long m) {r*=m;i*=m;}
  friend Quad operator*(long m, const Quad& a);
  Quad operator+ (const Quad& b) const {return Quad(r+b.r,i+b.i);}
  Quad operator+ (long b) const {return Quad(r+b,i);}
  friend Quad operator+(long m, const Quad& a);
  void operator+=(const Quad& b) {r+=b.r; i+=b.i;}
  void operator+=(long b) {r+=b;}
  Quad operator- (const Quad& b) const {return Quad(r-b.r,i-b.i);}
  Quad operator- (long b) const {return Quad(r-b,i);}
  friend Quad operator-(long m, const Quad& a);
  void operator-=(const Quad& b) {r-=b.r; i-=b.i;}
  void operator-=(long b) {r-=b;}
  Quad operator- () const {return Quad(-r,-i);}
  Quad operator/ (const Quad& b) const {return qdivi(mult(*this,quadconj(b)),quadnorm(b));}
  Quad operator/ (long b) const {return qdivi(*this,b);}
  void operator/=(const Quad& b) {*this=qdivi(mult(*this,quadconj(b)),quadnorm(b));}
  void operator/=(long b) {*this=qdivi(*this,b);}
  operator bigcomplex() const;   

// iostream functions

  friend ostream& operator<<(ostream& s, const Quad& x);
  friend istream& operator>>(istream& s, Quad& x);
};

char* to_string(const Quad& a);  // outputs to a (new) string

// Inline definitions of friend functions of class Quad (those not here are 
// not inline, and are in quads.cc):

inline Quad operator% (const Quad& a, const Quad& b) 
{ return a-(b*(a/b));}
inline Quad makepos(const Quad& a) 
{Quad ans=a; 
 while(!pos(ans)) ans=ans*fundunit; 
 return ans;
}
inline long real(const Quad& a) {return a.r;}
inline long imag(const Quad& a) {return a.i;}
inline long quadnorm0(const Quad& a)      // used when t=0
   {return a.r*a.r + Quad::n*a.i*a.i;}  
inline long quadnorm1(const Quad& a)      // used when t=1
   {return a.r*a.r + Quad::n*a.i*a.i + a.r*a.i;}
inline Quad quadconj0(const Quad& a)  {return Quad(a.r , -a.i);}
inline Quad quadconj1(const Quad& a)  {return Quad(a.r + a.i, -a.i);}
inline Quad mult0(const Quad& a, const Quad& b)
   {return Quad(a.r*b.r-Quad::n*a.i*b.i, a.r*b.i+a.i*b.r);}
inline Quad mult1(const Quad& a, const Quad& b)
   {return Quad(a.r*b.r-Quad::n*a.i*b.i, a.r*b.i+a.i*b.r+a.i*b.i);}
inline Quad qdivi0(const Quad& a, long c)
   {return Quad(roundover(a.r,c),roundover(a.i,c));}     // used when t=0
inline Quad qdivi1(const Quad& a, long c)
   {Quad b; b.i=roundover(a.i,c); 
            b.r=roundover(2*a.r+a.i-c*b.i,2*c); 
    return b;}                                          // used when t=1
inline int pos13(const Quad& a)
   {return (((a.i>=0)&&(a.r>0))||((a.r==0)&&(a.i==0)));}
inline int pos2(const Quad& a)
   {return ((a.i>0)||((a.i==0)&&(a.r>=0)));}
inline Quad operator*(long m, const Quad& a) {return Quad(m*a.r,m*a.i);}
inline Quad operator+(long m, const Quad& a) {return Quad(m+a.r,a.i);}
inline Quad operator-(long m, const Quad& a) {return Quad(m-a.r,-a.i);}
inline istream& operator>>(istream& s, Quad& x) {return s >> x.r >> x.i;}

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

#endif

// END OF FILE QUADS.H                              
