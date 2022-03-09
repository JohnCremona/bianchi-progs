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

#include "qideal.h"

// The Quadprime class is an enhanced version of the Qideal class

class Quadprime : public Qideal {
  long p;      // underlying prime of Z
  int character; // 0 (ramified), +1 (split), -1 (inert).
public:
  // constructor
  Quadprime(QUINT a, QUINT b, QUINT c, long pp, long ind=1)
    : Qideal(a,b,c) { p=pp; index=ind; character=Quad::chi(p); fill();}
  Quadprime(const Quadprime& x) : Qideal(x) { p=x.p; character=x.character;}
  Quadprime(Qideal& I); // constructor from an ideal (which should be a nonzero prime ideal)
  Quadprime() : Qideal() { p=0; character=0;}
  int is_ramified() const {return character==0;}
  int is_split() const {return character==1;}
  int is_inert() const {return character==-1;}
  long prime() const {return p;}
  int residue_degree() const {return (character==-1? 2: 1);}
  int ramification_degree() const {return (character==0? 2: 1);}

  vector<int> genus_character(); // vector of unram char values, one per prime discriminant
  int genus_character(const QUINT& D); // one unram char value
  long genus_class() {return from_bits(genus_character());}
  int has_square_class() {return (genus_class()==0);}

  friend inline ostream& operator<<(ostream& s, const Quadprime& x);
  friend istream& operator>>(istream& s, Quadprime& P);
};

inline ostream& operator<<(ostream& s, const Quadprime& x)
{
  s << "P" << x.p;
  if (x.is_split()) s << (x.index==1? "a": "b") ;
  return s;
}

typedef pair<Quadprime, int> QuadprimePower;

inline ostream& operator<<(ostream& s, const QuadprimePower& Q)
{
  s << Q.first;
  if (Q.second>1) s << "^" << Q.second;
  return s;
}

vector<Quadprime> Quadprimes_above(long p); // p should be an integer prime

class Quadprimes {
public:
  static QUINT maxnorm;           // largest norm of primes
  static vector<Quadprime> list;  // the list of primes
  static void init(long maxn=1000);     //  sets the list up
  static void display(ostream& s = cout, long maxn=0, int show_genus=0);
};

class QuadprimeLooper {
public:
  // Constructor: the iterator will skip primes dividing N
  QuadprimeLooper(Qideal level=Qideal())
    :Pi(Quadprimes::list.begin()), N(level)
  {
    P = *Pi;
    while (P.divides(N) && ok())
      {
        ++Pi;
        P = *Pi;
      }
  }
  // increment, if possible, skipping primes dividing N
  void operator++()
  {
    ++Pi;
    if (ok())
      {
        P = *Pi;
        while (P.divides(N) && ok())
          {
            ++Pi;
            if (ok())
              P = *Pi;
          }
      }
  }
  // return the current P
  operator Quadprime()
  {
    return P;
  }
  // test whether we have reached the end of Quadprimes::list
  int ok()
  {
    return Pi!=Quadprimes::list.end();
  }
  int at_end()
  {
    return Pi==Quadprimes::list.end();
  }
  void reset()
  {
    Pi = Quadprimes::list.begin();
    P = *Pi;
    while (P.divides(N)) P = *Pi++;
  }
private:
  vector<Quadprime>::iterator Pi;
  Qideal N;
  Quadprime P;
};


class Factorization {
  Qideal I;                     // the ideal whose factorization this is
// The order of the prime powers in the Factorization is given by the order of the underlying rational primes
  vector<QuadprimePower> Qlist; // prime powers Q (as (P,e) pairs)
  vector<Qideal> QIlist;        // prime powers Q (as ideals)
  vector<Quad> CRT_vector;      // list of Quads =1 mod each Q and =0 mod the others (set when first needed)
  void init_CRT();              // compute the CRT vector
public:
  explicit Factorization(const Qideal &);           // constructor

  long size() const { return Qlist.size(); }

  // In the following methods, the order of the
  // primes/prime-powers/exponents is unchanged, hence is the order of
  // the underlying rational primes.  In particular, primes() returns
  // a list of the primes which is not necessarily sorted by norm.
  Quadprime prime(int i) const { return Qlist[i].first; }
  vector<Quadprime> primes() const;
  int exponent(int i) const { return Qlist[i].second; }
  vector<int> exponents() const;
  Qideal prime_power(int i) const {return QIlist[i];}
  vector<QuadprimePower> prime_powers() const {return Qlist; }
  // The next method sorts the output of primes() into norm order
  vector<Quadprime> sorted_primes() const;

  // Test whether this factorization represents a prime, or a prime power:
  int is_prime() {return Qlist.size()==1 && Qlist[0].second==1;}
  int is_prime_power() {return Qlist.size()==1;}
  Quad solve_CRT(const vector<Quad>& v); // solution to x=v[i] mod Qlist[i]
  friend vector<Quadprime> pdivs(Qideal& n);
  friend vector<Qideal> alldivs(Qideal& a);
  friend inline ostream& operator<<(ostream& s, const Factorization& F);
  };

inline ostream& operator<<(ostream& s, const Factorization& F)
{
  if(F.size())
    {
      for (vector<QuadprimePower>::const_iterator Qi=F.Qlist.begin(); Qi!=F.Qlist.end(); ++Qi)
        {
          if (Qi!=F.Qlist.begin()) s<<"*";
          s << (*Qi);
        }
    }
  else
    {
      s << "()";
    }
  return s;
}

vector<Quadprime> pdivs(Qideal&);   // list of all prime divisors
vector<Qideal> alldivs(Qideal&);    // list of all ideal divisors
vector<Qideal> sqdivs(Qideal&);     // list of ideal divisors whose square divides
vector<Qideal> sqfreedivs(Qideal&); // list of square-free ideal divisors

#endif

// END OF FILE primes.h
