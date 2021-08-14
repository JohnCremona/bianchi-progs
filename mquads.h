// FILE MQUADS.H -- Multilength mquads

#if     !defined(_MQUADS_H)
#define _MQUADS_H      1       //flags that this file has been included

#include "eclib/marith.h"
#define complex bigcomplex

class mQuad;
class mQuadlist;
class mQuadvar;

//functions assigned by mQuad::field initializer
extern mQuad (*mquadconj)(const mQuad& a);  //Can't have same names as Complex functions
extern bigint (*mquadnorm)(const mQuad& a);  //
extern mQuad (*mult)(const mQuad& a, const mQuad& b);
extern mQuad (*qdivi)(const mQuad& a, const bigint& c);
extern int (*pos)(const mQuad& a);
//GCD-related functions.
extern mQuad (*mquadgcd)(const mQuad& aa, const mQuad& bb);
mQuad mquadgcd1(const mQuad& aa, const mQuad& bb); //Euclidean only
mQuad mquadgcd2(const mQuad& aa, const mQuad& bb); //General
extern mQuad (*mquadbezout)(const mQuad& aa, const mQuad& bb, mQuad& xx, mQuad& yy);
mQuad mquadbezout1(const mQuad& aa, const mQuad& bb, mQuad& xx, mQuad& yy); //Euclidean
mQuad mquadbezout2(const mQuad& aa, const mQuad& bb, mQuad& xx, mQuad& yy); //General
mQuad invmod(const mQuad& a, const mQuad& p);
int coprime(const mQuad& a, const mQuad& b);
int invertible(const mQuad& a, const mQuad& b, mQuad& inverse);
 

//functions defined in mmquads.cc
int div(const mQuad& a, const mQuad& b);
int ndiv(const mQuad& a, const mQuad& b);
int val(const mQuad& factor, const mQuad& number);
mQuadlist residues(const mQuad& a);

//Primes
extern mQuadlist mquadprimes;  
extern long nmquadprimes;     //The number of them.
extern mQuadlist mquadunits;
extern mQuad fundunit;

class mQuad {
/***********************************************************************
First the static members which only depend on the field.  To initialize,
include the line
        mQuad::field(d,maxnorm);
at the top of the program, where d is a suitable positive integer
and maxnorm (default 1000) is the upper bound for the norms of primes.
***********************************************************************/
 public:
  static int     d;          // square-free >0
  static bigint     t;          // trace of w
  static bigint     n;          // norm of w
  static char    name;       // name of w for printing
  static bigint  maxnorm;       // largest norm of primes
 public:
  static int     disc;       // discriminant
  static int   nunits;       // number of units

public:
  static void field(int dd, const bigint& max=BIGINT(1000));
  static void displayfield(ostream& s = cout);
  static void initmquadprimes();

// Now the actual data elements:
 private:
  bigint r,i;
 public:
//constructors:
  mQuad(const bigint& x=BIGINT(0), const bigint& y=BIGINT(0)) :r(x),i(y) {;}
  explicit mQuad(long x) :r((bigint)(x)),i(BIGINT(0)) {;}
  explicit mQuad(const complex& z);   //rounds to nearest
  mQuad(const mQuad& a) :r(a.r), i(a.i) {;}

//operators and related functions:

  mQuad& operator=(const mQuad& a) {r=a.r; i=a.i; return *this;}
  void operator=(const bigint& a) {r=a; i=0; return *this;}
  void operator=(long a) {r=a; i=0; return *this;}
  friend inline bigint real(const mQuad& a) {return a.r;}
  friend inline bigint imag(const mQuad& a) {return a.i;}
  friend inline bigint mquadnorm0(const mQuad& a)      // used when t=0
   {return a.r*a.r + mQuad::n*a.i*a.i;}  
  friend inline bigint mquadnorm1(const mQuad& a)      // used when t=1
   {return a.r*a.r + mQuad::n*a.i*a.i + a.r*a.i;}
  friend inline mQuad mquadconj0(const mQuad& a)  {return mQuad(a.r , -a.i);}
  friend inline mQuad mquadconj1(const mQuad& a)  {return mQuad(a.r + a.i, -a.i);}
  friend inline mQuad mult0(const mQuad& a, const mQuad& b)
   {return mQuad(a.r*b.r-mQuad::n*a.i*b.i, a.r*b.i+a.i*b.r);}
  friend inline mQuad mult1(const mQuad& a, const mQuad& b)
   {return mQuad(a.r*b.r-mQuad::n*a.i*b.i, a.r*b.i+a.i*b.r+a.i*b.i);}
  friend inline mQuad qdivi0(const mQuad& a, const bigint& c)
   {return mQuad(roundover(a.r,c),roundover(a.i,c));}     // used when t=0
  friend inline mQuad qdivi1(const mQuad& a, const bigint& c)
   {mQuad b; b.i=roundover(a.i,c); 
            b.r=roundover(bigint(2)*a.r+a.i-c*b.i,bigint(2)*c); 
    return b;}                                          // used when t=1
  friend inline int pos13(const mQuad& a)
   {return (((a.i>=0)&&(a.r>0))||((a.r==0)&&(a.i==0)));}
  friend inline int pos2(const mQuad& a)
   {return ((a.i>0)||((a.i==0)&&(a.r>=0)));}

  int operator== (const mQuad& b) const {return (r==b.r) && (i==b.i);}
  int operator== (const bigint& b) const {return (r==b) && (i==0);}
  int operator== (long b) const {return (r==(bigint)b) && (i==0);}
  int operator!= (const mQuad& b) const {return (r!=b.r) || (i!=b.i);}
  int operator!= (const bigint& b) const {return (r!=b) || (i!=0);}
  int operator!= (long b) const {return (r!=(bigint)b) || (i!=0);}
  mQuad operator* (const mQuad& b) const {return mult(*this,b);}
  void operator*=(const mQuad& b) {*this=mult(*this,b);}
  mQuad operator* (const bigint& m) const {return mQuad(m*r,m*i);}
  void operator*=(const bigint& m) {r*=m;i*=m;}
  mQuad operator* (long m) const {bigint mm(m); return mQuad(mm*r,mm*i);}
  void operator*=(long m) {bigint mm(m); r*=mm;i*=mm;}
  friend inline mQuad operator*(const bigint& m, const mQuad& a) {return mQuad(m*a.r,m*a.i);}
  friend inline mQuad operator*(long m, const mQuad& a) {bigint mm(m); return mQuad(mm*a.r,mm*a.i);}
  mQuad operator+ (const mQuad& b) const {return mQuad(r+b.r,i+b.i);}
  mQuad operator+ (const bigint& n) const {return mQuad(r+n,i);}
  mQuad operator+ (long n) const {return mQuad(r+bigint(n),i);}
  void operator+=(const mQuad& b) {r+=b.r; i+=b.i;}
  void operator+=(const bigint& b) {r+=b;}
  void operator+=(long b) {r+=b;}
  friend inline mQuad operator+(const bigint& m, const mQuad& a) {return mQuad(m+a.r,a.i);}
  friend inline mQuad operator+(long m, const mQuad& a) {return mQuad(bigint(m)+a.r,a.i);}
  mQuad operator- (const mQuad& b) const {return mQuad(r-b.r,i-b.i);}
  mQuad operator- (const bigint& n) const {return mQuad(r-bigint(n),i);}
  mQuad operator- (long n) const {return mQuad(r-bigint(n),i);}
  void operator-=(const mQuad& b) {r-=b.r; i-=b.i;}
  void operator-=(const bigint& b) {r-=b;}
  void operator-=(long b) {r-=b;}
  friend inline mQuad operator-(const bigint& m, const mQuad& a) {return mQuad(bigint(m)-a.r,-a.i);}
  friend inline mQuad operator-(long m, const mQuad& a) {return mQuad(bigint(m)-a.r,-a.i);}
  mQuad operator- () const {return mQuad(-r,-i);}
  mQuad operator/ (const mQuad& b) const {return qdivi(mult(*this,mquadconj(b)),mquadnorm(b));}
  mQuad operator/ (const bigint& b) const {return qdivi(*this,b);}
  mQuad operator/ (long b) const {return qdivi(*this,b);}
  void operator/=(const mQuad& b) {*this=qdivi(mult(*this,mquadconj(b)),mquadnorm(b));}
  void operator/=(const bigint& b) {*this=qdivi(*this,b);}
  void operator/=(long b) {*this=qdivi(*this,b);}
  operator complex();   

// friend functions

  friend ostream& operator<<(ostream& s, const mQuad& x);
  friend inline istream& operator>>(istream& s, const mQuad& x)
    {s.flags(s.flags()|ios::dec);  // to avoid bug in library, force decimal
     return s >> x.r >> x.i;
    }
};

char* to_string(const mQuad& a);  // outputs to a (new) string

inline mQuad operator% (const mQuad& a, const mQuad& b) 
{ return a-(b*(a/b));}
inline mQuad makepos(const mQuad& a) 
{mQuad ans=a; 
 while(!pos(ans)) ans=ans*fundunit; 
 return ans;
}

class mQuadlist {
public:
  long length;        
private:
  mQuad *items;  
public:
  mQuadlist(long n=0) {length=n; items=new mQuad[n];}
  ~mQuadlist() {length=0; delete[] items;}
  mQuadlist(long n, mQuad* list) {length=n; items=list;} //NB Does NOT copy!
  mQuadlist(const mQuadlist&l) 
    {items=new mQuad[length=l.length]; 
     long n=length; mQuad* x=items, *y=l.items;
     while(n--) {*x++=*y++;}}
  mQuadlist(const mQuadlist&l, long n) 
    {length= n<=l.length ? n : l.length;
     items=new mQuad[length]; 
     n=length; mQuad* x=items, *y=l.items;
     while(n--) {*x++=*y++;}}
  mQuadlist& operator=(const mQuadlist&l)
    {if (this==&l) return *this;
     delete[] items; items=new mQuad[length=l.length]; 
     long n=length; mQuad* x=items, *y=l.items;
     while(n--) {*x++=*y++;}
     return *this;}
  mQuad get(long i) const {
    if ((i>=0)&&(i<length)) {return items[i];}
    else {cerr<<"bad index "<<i<<" in mQuadlist!\n"; return 0;}
  }
  mQuad& operator[](long i) {
    if ((i>=0)&&(i<length)) {return items[i];}
    else {cerr<<"bad index "<<i<<" in mQuadlist!\n"; return items[0];}
  }
  mQuad operator()(long i) const {
    if ((i>=0)&&(i<length)) {return items[i];}
    else {cerr<<"bad index "<<i<<" in mQuadlist!\n"; return items[0];}
  }
  void set(long i, mQuad val) {
    if ((i>=0)&&(i<length)) items[i]=val; 
    else {cerr<<"bad index "<<i<<" in mQuadlist!\n";}
  }
  void truncate(long newlength) 
    { if (newlength<length) length = newlength;}
  // no need to deallocate memory here because copy constructors do so
  // almost always it is the temporary of a function that calls truncate
  void output(ostream& os) const
    {long i;mQuad *it;
     for(i=0,it=items;i<length;i++,it++)
       {if(i)os<<", ";os<<(*it);}
   }
  friend inline ostream& operator<<(ostream& os, const mQuadlist& l)
    {l.output(os);return os;}
  friend class mQuadvar;
};

class mQuadvar {
public:
        long index;        /* current index */
private:
        mQuad* values;
        long maxindex;     /* max index */
public:
        explicit mQuadvar(const mQuadlist& l) 
          {maxindex=l.length-1; index=0; values=l.items;}
        void init(const mQuadlist& l) 
          {maxindex=l.length-1; index=0; values=l.items;}
        long next() {if ((index++)<maxindex) {values++; return 1;}
                    else return 0;
                   }
        mQuad operator++() {if ((index++)<maxindex) {values++; return *values;}
                           else return 0;
                          }
        int ok() const {return index<=maxindex;}
        int more() const {return index<maxindex;}
        mQuad value() const {return *values;}
        operator mQuad() const {return *values;}
};

/* Usage of mQuadvar: to loop through a mQuadlist "alist" (which might be empty):
   mQuad a; either of

   for(mQuadvar avar(alist); avar.ok(); avar++) {a=avar; ... ;}
   for(avar.init(alist); avar.ok(); avar++) {a=avar; ... ;}  
   // second form IFF avar already exists

   Either way, the argument, alist, must be an instance of an existing
   mQuadlist, and must not be a function returning a mQuadlist (else
   unflagged error will follow)
   
*/

mQuad primdiv(const mQuad&);           // "First" prime divisor
mQuadlist pdivs(const mQuad&);         // list of prime divisors
mQuadlist posdivs(const mQuad&);       // all "positive" divisors (up to units)
mQuadlist alldivs(const mQuad&);       // absolutely all divisors
mQuadlist sqdivs(const mQuad&);        // returns divisors whose square divides
                              // (up to +/- sign)
mQuadlist sqfreedivs(const mQuad&);    // returns square-free divisors


#endif

// END OF FILE QUADS.H                              
