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
  Quadprime(INT a, INT b, INT c, long pp, long ind=1)
    : Qideal(a,b,c) { p=pp; index=ind; character=Quad::chi(INT(p)); fill();}
  Quadprime(const Quadprime& x) : Qideal(x) { p=x.p; character=x.character;}
  explicit Quadprime(Qideal& I); // constructor from an ideal (which should be a nonzero prime ideal)
  Quadprime() : Qideal() { p=0; character=0;}
  Quadprime& operator=(const Quadprime& x) {
    Qideal::operator=(x);
    p=x.p;
    character=x.character;
    return *this;
  }
  int is_ramified() const {return character==0;}
  int is_split() const {return character==1;}
  int is_inert() const {return character==-1;}
  long prime() const {return p;}
  int residue_degree() const {return (character==-1? 2: 1);}
  int ramification_degree() const {return (character==0? 2: 1);}

  vector<int> genus_character(); // vector of unram char values in {0,1}, one per prime discriminant
  int genus_character(const INT& D); // one unram char value in {-1,1}
  long genus_class(int contract=0); // integer in [0,2^r-1] whose bits are genus_character() when contract=0
                                    // or same reduced mod 2^{r-1} when contract=1
  int has_square_class() {return (genus_class()==0);}

  friend inline string prime_label(const Quadprime& x);
  friend istream& operator>>(istream& s, Quadprime& P);
};

inline string prime_label(const Quadprime& x)
{
  stringstream s;
  s << "P" << x.p;
  if (x.is_split()) s << (x.index==1? "a": "b") ;
  return s.str();
}

inline ostream& operator<<(ostream& s, const Quadprime& x)
{
  s << prime_label(x);
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
  static INT maxnorm;           // largest norm of primes
  static vector<Quadprime> list;  // the list of primes
  static void init(long maxn=1000);     //  sets the list up
  static void display(ostream& s = cout, long maxn=0, int show_genus=0);
  Quadprime operator[](int i) {return list[i];}
  static int index(const Quadprime& P)
  {
    return std::find(list.begin(), list.end(), P) - list.begin();
  }
  static int conjugate_index(int i)
  {
    return std::find(list.begin(), list.end(), list[i].conj()) - list.begin();
  }
};

class QuadprimeLooper {
public:
  // Constructor: the iterator will skip primes dividing level
  QuadprimeLooper();
  explicit QuadprimeLooper(const Qideal& level);
  // increment, if possible, skipping primes dividing N
  void operator++();
  // return the current P
  operator Quadprime() {return P;}
  // test whether we have reached the end of Quadprimes::list
  int ok() {return Pi!=Quadprimes::list.end();}
  int at_end() {return Pi==Quadprimes::list.end();}
  void reset();
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
      for (auto Qi=F.Qlist.begin(); Qi!=F.Qlist.end(); ++Qi)
        {
          if (Qi!=F.Qlist.begin())
            s<<"*";
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
int ndivs(Qideal&); // number of ideal divisors
int npdivs(Qideal&);  // number of prime ideal divisors

// Return {-m,m} where m is the largest integer <= +2*sqrt(N(P)), the bounds on a(P)
pair<long,long> eigenvalue_range(const Quadprime& P);

// Return {-m,3*m} where m = N(P), the bounds on a(P^2)=a(P)^2-N(P)
pair<long,long> eigenvalue_sq_range(const Quadprime& P);

// Return list of integers between -2*sqrt(N(P)) and +2*sqrt(N(P)) if [P] is square, else
// list of possible eigs for T(P^2) = T(P)^2-N(P), assuming T(P,P) trivial
vector<long> good_eigrange(Quadprime& P);

// compute a list of primes Q dividing N with Q^e||N such that [Q^e] is square
vector<Quadprime> make_squarebadprimes(const Qideal& N, const vector<Quadprime>& badprimes);

// compute a list of at least nap good primes (excluding those
// dividing characteristic if >0), to include at least one principal
// one which has index iP0;
vector<Quadprime> make_goodprimes(const Qideal& N,  int np, int& iP0, long p);

inline long prime_index(const Quadprime& P)
{
  return find(Quadprimes::list.begin(), Quadprimes::list.end(), P) - Quadprimes::list.begin();
}



#endif

// END OF FILE primes.h
