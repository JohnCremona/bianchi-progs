// FILE primes.h

//////////////////////////////////////////////////////////////////////
// All data elements of Class Quadprimes are static.  To initialize,
// include the line
//    Quadprimes::init(d);
// at the top of the program, where maxnorm (default 1000) is the upper
// bound for the norms of primes.
//////////////////////////////////////////////////////////////////////

// The files primes* extend the capabilities of quadarith.* by
// providing a list of prime ideals, and various functions which
// depend on it, e.g for listing the prime divisors, or all divisors,
// of a Quad.

#ifndef __PRIMES_H__
#define __PRIMES_H__

#include "quads.h"
#include "qideal.h"

#include <iostream>

// The Quadprime class is an enhanced version of the Qideal class

class Quadprime : public Qideal {
  long p;      // underlying prime of Z
  long flag;   // = 1 or 2 for a pair of split primes
public:
  // constructor
  Quadprime(long a, long b, long c, long pp, long fflag)
    : Qideal(a,b,c) { p=pp; flag=fflag; }

  Quadprime(const Qideal&x) : Qideal(x) { p=0; flag=-1; }
  Quadprime() {;}
  Quadprime(const Quadprime&x) : Qideal(x) { p=x.p; flag=x.flag; }
  friend inline ostream& operator<<(ostream& s, const Quadprime& x);
};

inline ostream& operator<<(ostream& s, const Quadprime& x)
{
  if (x.flag==-1)
    s << (Qideal)x;
  else
    {
      s << "P" << x.p;
      if (x.flag==1) s << "a";
      if (x.flag==2) s << "b";
    }
  return s;
}

vector<Quadprime> Quadprimes_above(long p); // p should be an integer prime

class Quadprimes {
public:
  static long maxnorm;           // largest norm of primes
  static vector<Quadprime> list;  // the list of primes
  static void init(long maxn=1000);     //  sets the list up
  static void display(ostream& s = cout);
};

class prime_factn {
  vector<Quadprime> plist;
  vector<long> elist;
public:
  prime_factn(const Qideal &);           // constructor
  void display(ostream& s = cout) const;

  long num_primes() const { return plist.size(); }
  Quadprime prime(long i) const { return plist[i]; }
  long expo(long i) const { return elist[i]; }
  friend vector<Quadprime> pdivs(const Qideal& n);
  friend vector<Qideal> alldivs(const Qideal& a);
  };

Quadprime primdiv(const Qideal&);         // "First" prime divisor
vector<Quadprime> pdivs(const Qideal&);   // list of all prime divisors
vector<Qideal> alldivs(const Qideal&);    // list of all ideal divisors
vector<Qideal> sqdivs(const Qideal&);     // list of ideal divisors whose square divides
vector<Qideal> sqfreedivs(const Qideal&); // list of square-free ideal divisors

#endif

// END OF FILE primes.h
