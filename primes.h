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
  int character; // 0 (ramified), +1 (split), -1 (inert).
public:
  // constructor
  Quadprime(long a, long b, long c, long pp, long ind=1)
    : Qideal(a,b,c) { p=pp; index=ind; character=Quad::chi(p); fill();}

  Quadprime(const Quadprime&x) : Qideal(x) { p=x.p; character=x.character;}
  int is_ramified() const {return character==0;}
  int is_split() const {return character==1;}
  int is_inert() const {return character==-1;}
  long prime() const {return p;}
  int residue_degree() const {return (character==-1? 2: 1);}
  int ramification_degree() const {return (character==0? 2: 1);}
  friend inline ostream& operator<<(ostream& s, const Quadprime& x);
};

inline ostream& operator<<(ostream& s, const Quadprime& x)
{
  s << "P" << x.p;
  if (x.is_split()) s << (x.index==1? "a": "b") ;
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

vector<Quadprime> pdivs(const Qideal&);   // list of all prime divisors
vector<Qideal> alldivs(const Qideal&);    // list of all ideal divisors
vector<Qideal> sqdivs(const Qideal&);     // list of ideal divisors whose square divides
vector<Qideal> sqfreedivs(const Qideal&); // list of square-free ideal divisors

#endif

// END OF FILE primes.h
