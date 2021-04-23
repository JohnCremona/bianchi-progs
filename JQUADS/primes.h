// FILE primes.h

//////////////////////////////////////////////////////////////////////
// All data elements of Class Quadprimes are static.  To initialize,
// include the line 
//    Quadprimes::init(d);
// at the top of the program, where maxnorm (default 1000) is the upper
// bound for the norms of primes.
//////////////////////////////////////////////////////////////////////

// The files primes* extend the capabilities of quadarith.* by providing
// a list of primes and various functions which depend on it, e.g for listing
// the prime divisors, or all divisors, of a Quad.
//
// Two versions can be compiled.  This is controlled by the macro MAX_CLASSNUM.
// If it is 1 (the default) primes are just (irreducible) Quads.
// If it is 2 (higher values are not supported), primes are Qideals, and the
// file qideal.h is included to support their use.

#ifndef __PRIMES_H__
#define __PRIMES_H__

#ifndef MAX_CLASSNUM
#define MAX_CLASSNUM 1
#endif

#include "quadarith.h"

#if MAX_CLASSNUM<2
#define PRIMETYPE Quad
#else
#include "qideal.h"
#define PRIMETYPE Qideal
#endif

#include <iostream>

class Quadprime : public PRIMETYPE {
  long p;
  long flag;
public:
  // constructor
#if MAX_CLASSNUM<2
  compile time error - not implemented!!!;
#else
  Quadprime(long a, long b, long c, long pp, long fflag)
    : Qideal(a,b,c) { ((Qideal)*this)=this->pos_assoc(); p=pp; flag=fflag; }

  Quadprime(const Qideal&x) : Qideal(x) { p=0; flag=-1; }
#endif

  Quadprime() {;}
  Quadprime(const Quadprime&x) : PRIMETYPE(x) { p=x.p; flag=x.flag; }
  friend inline ostream& operator<<(ostream& s, const Quadprime& x);
};

inline ostream& operator<<(ostream& s, const Quadprime& x)
{
  if (x.flag==-1)
    s << (PRIMETYPE)x;
  else
    {
      s << "P" << x.p;
      if (x.flag==1) s << "a";
      if (x.flag==2) s << "b";
    }
  return s;
}

typedef Tlist<Quadprime> Primelist;
typedef Tvar<Quadprime> Primevar;

class Quadprimes {
public:
  static long maxnorm;           // largest norm of primes
  static Primelist list;  // the list of primes
  static void init(long max=1000);     //  sets the list up
  static void display(ostream& s = cout);
};

class prime_factn {
  Primelist    plist;
  Tlist<long>  elist;
public:
  prime_factn(const PRIMETYPE &);           // constructor
  void display(ostream& s = cout) const;

  long num_primes() const { return plist.getlength(); }
  Quadprime prime(long i) const { return plist(i); }
  long expo(long i) const { return elist(i); }
};

// implemented uniformly
Quadprime primdiv(const PRIMETYPE&);        // "First" prime divisor
Primelist pdivs(const PRIMETYPE&);          // list of prime divisors
Tlist<PRIMETYPE> posdivs(const PRIMETYPE&); // all "positive" d.s (up to units)

// implemented separately
Tlist<PRIMETYPE> alldivs(const PRIMETYPE&); // absolutely all divisors
Tlist<PRIMETYPE> sqdivs(const PRIMETYPE&);  // divisors whose square divides
                              // (up to +/- sign)
Tlist<PRIMETYPE> sqfreedivs(const PRIMETYPE&); // returns square-free divisors

#endif

// END OF FILE primes.h
