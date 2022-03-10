// FILE primes.cc

#include "primes.h"
#include "qidloop.h"
#include "intprocs.h"
#include <list>

//Definition of static data members of class Quadprimes:
QUINT Quadprimes::maxnorm;
vector<Quadprime> Quadprimes::list;

// Quadprime constructor from an ideal (which should be a nonzero prime ideal)

Quadprime::Quadprime(Qideal& I) :Qideal(I)
{
  if (I.is_prime())
    {
      *this = I.factorization().prime(0);
    }
  else
    {
      cerr << "***cannot construct a prime ideal from the ideal "<<I<<" = "<<F<<", which is not prime!"<<endl;
    }
}

vector<int> Quadprime::genus_character()
{
  // cout<<"Finding genus character of P="<<(*this)<<endl;

  // NB this function is indirectly called in fill_class_group() at
  // which point Quad::class_group_2_rank has not been set

  int r = Quad::discfactors.size()-1; // = Quad::class_group_2_rank
  vector<int> ans(1+r, 0);
  if (r==0 || is_inert() || is_principal())
    return ans;
  // cout<<" -- not principal, even class number..."<<endl;
  int i=0, i0=-1, tot=0;
  for (auto di = Quad::discfactors.begin(); di!=Quad::discfactors.end(); ++di, ++i)
    {
      int k;
      long d = I2long(*di);
      if (p==2)
        k = (d%2? (posmod(d,8)==1? +1 : -1) : 0);
      else
        k = legendre(d,p);
      if (k)
        {
          ans[i] = int(-1==k); // we want additive characters, values 0,1 mod 2
          tot ^= ans[i];       // binary + (exclusive or)
        }
      else
        {
          i0 = i;
        }
    }
  if (i0!=-1) // then we have a ramified prime
    {
      ans[i0] = tot;
      tot = 0;
    }
  // cout<<"P="<<(*this)<<": genus character "<<ans<<endl;
  assert (tot==0);
  return ans;
}

int Quadprime::genus_character(const QUINT& D)
{
  if (!div_disc(D, Quad::disc))
    return 0;
  int r = Quad::discfactors.size();
  // cout<<"In P.genus_character(D) with P="<<(*this)<<", D="<<D;
  // cout<<" with character values "<< chardisc(D) <<endl;
  // cout<<"genus_class of P is "<<genus_class()<<", bits "<<bits(genus_class(), r)<<endl;
  int dot = dotbits(from_bits(chardisc(D)), genus_class(), r);
  // cout<<"dot product = "<<dot<<" --> "<<(dot? -1: +1)<<endl;
  return (dot? -1: +1);
}

istream& operator>>(istream& s, Quadprime& P)
{
  Qideal I;
  s >> I;
  P = Quadprime(I);
  return s;
}

void Quadprimes::display(ostream& s, long maxn, int show_genus)
{
  s << list.size() << " prime ideals initialised, ";
  s << "max norm = " << maxnorm << endl;
  if (maxn==0) return;
  for (vector<Quadprime>::iterator Pi = Quadprimes::list.begin(); Pi != Quadprimes::list.end(); ++Pi)
    {
      Quadprime P = *Pi;
      if (P.norm()>maxn) break;
      vector<Quad> gg = P.gens();
      Quad pi;
      cout << P << " = " << ideal_label(P) << " = " << (Qideal)P << " = (" << gg[0] <<","<<gg[1] << ")";
      if (P.is_principal(pi))
        cout << " = ("<< pi <<") (principal";
      else
        cout << " (not principal";
      if (show_genus && Quad::class_group_2_rank>0)
        {
          vector<int> c = P.genus_character();
          cout << ": ";
          for (auto ci = c.begin(); ci!=c.end(); ++ci)
            cout << (*ci? "-": "+");
        }
      cout<<")"<<endl;
    }
}

vector<Quadprime> Quadprimes_above(long p) // p should be an integer prime
{
  long d=Quad::d; QUINT disc=Quad::disc;
  int t=Quad::t;
  vector<Quadprime> Plist;
  QUINT zero(0), one(1), P(p);

  if (p==2) // treat as special case
    {
      switch (d%4)
        {
        case 1: // ramified, (2) = (2,1+w)^2
          Plist.push_back(Quadprime(P,one,one, p));
          break;
        case 2: // ramified, (2) = (2,w)^2
          Plist.push_back(Quadprime(P,zero,one, p));
          break;
        case 3: // split, (2) = (2,w)*(2,1+w) or inert, (2) = (2)
          if (d%8==3) // inert
            {
              Plist.push_back(Quadprime(one,zero,P, p));
            }
          else // split
            {
              Plist.push_back(Quadprime(P,zero,one, p, 1));
              Plist.push_back(Quadprime(P,one,one, p, 2));
            }
        }
      return Plist;
    }

  // odd p

  if (d%p==0) // ramified, (p) = (p,b+w)^2 where b=t/2 mod p
    {
      //      cout << "odd ramified case" << endl;
      if (t)
        Plist.push_back(Quadprime(P,-(P+one)/2,one, p));
      else
        Plist.push_back(Quadprime(P,zero,one, p));
      return Plist;
    }

  if (legendre(disc,p) == -1) //inert
    {
      //      cout << "odd inert case" << endl;
      Plist.push_back(Quadprime(one,zero,P, p));
    }
  else
    // split, (p) = (p,b1+w)*(p,b2+w) where b1, b2 are roots of x^2+t*x+n=0 (modp)
    // We order so the HNFs are [p,b,1], [p,b',1] with b<b'
    {
      //      cout << "odd split case" << endl;
      QUINT r, b1, b2, MD(posmod(-d,p));
      // NB NTL's SqrRootMod requires the argument to be positive!
      sqrt_mod_p(r, MD, P);
      if (t)
        {
          b1 = posmod(((r%2)? (r-1)/2: (p+r-1)/2), p);
          b2 = posmod(-1-b1, p);
        }
      else
        {
          b1 = posmod(r, p);
          b2 = posmod(-b1, p);
        }
      if (b1>b2) swap(b1,b2);
      Plist.push_back(Quadprime(P,b1,one, p, 1));
      Plist.push_back(Quadprime(P,b2,one, p, 2));
    }
  return Plist;
}

// The order of the prime powers in the Factorization is given by the order of the underlying rational primes

Factorization::Factorization(const Qideal& II)
  : I(II)
{
  Qideal J(I);    // copy for dividing out primes while leaving I unchanged
  // cout<<" J = "<< J <<endl;
  Qideal Q; int e;
  vector<QUINT> pdivs_norm = pdivs(I.norm());
  //  cout<<"Finding prime factors of "<<I<<" with norm "<<I.norm()<<", primes dividing norm are "<<pdivs_norm<<endl;
  for(vector<QUINT>::const_iterator pi = pdivs_norm.begin(); pi!=pdivs_norm.end(); ++pi)
    {
      vector<Quadprime> PP = Quadprimes_above(I2long(*pi));
      //      cout<<"primes above "<<(*pi)<<" are "<<PP<<endl;
      // at least one, but possibly not both when p splits, divides I
      for(vector<Quadprime>::iterator Pi = PP.begin(); Pi!=PP.end(); ++Pi)
        {
          Quadprime P = *Pi;
          if (P.divides(J))
            {
              e = 1;
              J /= P;
              // cout<<" After dividing by "<<P<<", J = "<< J <<endl;
              Q = P;
              while (P.divides(J))
                {
                  e++;
                  J /= P;
                  // cout<<" After dividing by "<<P<<", J = "<< J <<endl;
                  Q *= P;
                }
              Qlist.push_back({P,e});
              QIlist.push_back(Q);
            }
        }
    }
  // just for testing, do this on construction:
  init_CRT();
}

vector<Quadprime> Factorization::primes() const
{
  vector<Quadprime> plist;
  plist.reserve(size());
  for (vector<QuadprimePower>::const_iterator Qi = Qlist.begin(); Qi!=Qlist.end(); ++Qi)
    plist.push_back(Qi->first);
  return plist;
}

vector<Quadprime> Factorization::sorted_primes() const
{
  vector<Quadprime> plist = primes();
  std::sort(plist.begin(), plist.end(), Qideal_cmp);
  return plist;
}

vector<int> Factorization::exponents() const
{
  vector<int> elist;
  elist.reserve(size());
  for (vector<QuadprimePower>::const_iterator Qi = Qlist.begin(); Qi!=Qlist.end(); ++Qi)
    elist.push_back(Qi->second);
  return elist;
}

vector<Quadprime> pdivs(Qideal& I)  // list of prime divisors
{
  return I.factorization().sorted_primes();
}

vector<Qideal> alldivs(Qideal& a)    // list of all ideal divisors
{
  Factorization F = a.factorization();
  int np = F.size();
  long nu = 1; long nd=nu;
  for (long i=0; i<np; i++) {nd*=(1+F.exponent(i));}
  vector<Qideal> dlist(nd);
  dlist[0]=Qideal(Quad::one);
  nd=nu;
  vector<QuadprimePower>::const_iterator Qi;
  for(Qi = F.Qlist.begin();  Qi != F.Qlist.end();  ++Qi)
    {
      Qideal P = Qi->first;
      int e = Qi->second;
      for (int j=0; j<e; j++)
	for (int k=0; k<nd; k++)
	  dlist[nd*(j+1)+k] = (P*dlist[nd*j+k]);
      nd*=(e+1);
    }
  return dlist;
}

void Factorization::init_CRT()              // compute the CRT vector
{
  if (CRT_vector.size()!=0) return; // already done
  CRT_vector.reserve(size());
  if (size()==1)
    {
      CRT_vector.push_back(1);
      return;
    }
  Qideal Q, J;
  Quad r, s;
  for (int i=0; i<size(); i++)
    {
      Q = prime_power(i);
      J = I/Q;
      int t = J.is_coprime_to(Q, r, s);
      CRT_vector.push_back(r);
      // r+s=1 with r in J, s in Q
      // so r=1 mod Q, r=0 mod Q' for other Q'
      assert (t);
      assert (Q.contains(r-Quad::one));
      assert (J.contains(r));
    }
  // check
  for (int i=0; i<size(); i++)
    {
      Q = prime_power(i);
      for (int j=0; j<size(); j++)
        assert(Q.divides(CRT_vector[j]-int(i==j)));
    }
}

Quad Factorization::solve_CRT(const vector<Quad>& v) // solution to x=v[i] mod Qlist[i]
{
  if (CRT_vector.size()==0) init_CRT();
  Quad a(0);
  int i;
  for (i=0; i<size(); i++)
    a = I.reduce(a + v[i]*CRT_vector[i]);
  for (i=0; i<size(); i++)
    assert (prime_power(i).contains(a-v[i]));
  return a;
}

void Quadprimes::init(long maxn)
{
  maxnorm = maxn;
  vector<Quadprime> list1, list2;

  // First fill up lists of degree 1 and degree 2 primes
  //  cout<<"Computing list of prime ideals of norm up to "<<maxnorm<<endl;

  for (primevar pr; pr.ok()&&pr<=maxnorm; ++pr)
    { long p=pr;
      //            cout<<"p = "<<p<<endl;
      vector<Quadprime> PP = Quadprimes_above(p);
      //            cout<<" primes above: "<< PP<<endl;
      for(vector<Quadprime>::iterator Pi = PP.begin(); Pi!=PP.end(); ++Pi)
        {
          Quadprime P = *Pi;
          QUINT q = P.norm();
          //                    cout<<"P = "<<P<<" with norm "<<q<<endl;
          if(q==p) // degree 1 prime
            list1.push_back(P);
          else
            if(q<=maxnorm)
              list2.push_back(P);
        }
    }

  //    cout<<" - found "<<list1.size() << " degree 1 primes and "<<list2.size()<<" degree 2 primes"<<endl;

  // Now merge these into a single list sorted by norm

  //    cout<<" - merging into a single list" <<endl;
  vector<Quadprime>::iterator Pi = list1.begin(), Qi = list2.begin();
  while(Pi!=list1.end() && Qi!=list2.end())
    {
      Quadprime P = *Pi, Q = *Qi;
      if(P.norm()<Q.norm())
        {
          list.push_back(P);
          ++Pi;
        }
      else
        {
          list.push_back(Q);
          ++Qi;
        }
    }
  // only one of the following will do anything, probably the first:
  while(Pi!=list1.end())
    list.push_back(*Pi++);
  while(Qi!=list2.end())
    list.push_back(*Qi++);
}

// need divisors functions etc

vector<Qideal> sqdivs(Qideal& a) // all divisors whose square divides a, up to +/-
{
  Factorization F = a.factorization();
  int np = F.size();
  long nd=1;
  for(long i=0; i<np; i++) { nd *= ( 1+ F.exponent(i)/2 ) ;}

  vector<Qideal> dlist(nd);
  dlist[0]=Qideal(Quad::one);
  nd=1;
  for(long i=0; i<np; i++)
    {
      Qideal P = F.prime(i);
      long e = F.exponent(i)/2;
      for(long j=0; j<e; j++)
	for(long k=0; k<nd; k++)
	  dlist[nd*(j+1)+k] = P*dlist[nd*j+k];
      nd*=(e+1);
    }
  return dlist;
}

vector<Qideal> sqfreedivs(Qideal& a)       // all square-free divisors
{
  vector<Quadprime> plist=pdivs(a);
  Qideal p;
  long np = plist.size();
  long nd = 1;
  while (np-->0) nd*=2;
  vector<Qideal> dlist(nd);
  dlist[0]=Qideal(Quad::one);
  nd=1;
  for(vector<Quadprime>::const_iterator pr = plist.begin(); pr != plist.end(); ++pr)
    {
      p = *pr;
      for (long k=0; k<nd; k++)
	dlist[nd+k] = p*dlist[k];
      nd*=2;
    }
  return dlist;
}

//#define DEBUG_SORT

//////////////////////////////////////////
// implementation of class Qideal_lists //
//////////////////////////////////////////

map<QUINT, vector<Qideal>> Qideal_lists::N_to_Ilist;

// return a sublist including exactly one of each conjugate pair (the
// one with smaller index).  This is a local function and we assume
// that the list is sorted so the i'th ideal (i>=0) has index i+1.
vector<Qideal> remove_conjugates(const vector<Qideal>& Ilist)
{
  if (!Ilist.size())
    return Ilist;
  vector<Qideal> Ilist_no_conj;
  vector<int> seen(Ilist.size(), 0);
  for (int i=0; i<(int)Ilist.size(); i++)
    {
      if (seen[i])
        continue; // the conjugate is already in the list
      Qideal I = Ilist[i];
      Ilist_no_conj.push_back(I);
      seen[I.conj().get_index()-1] = 1;
    }
  return Ilist_no_conj;
}

vector<Qideal> Qideal_lists::ideals_with_norm(QUINT N, int both_conj)
{
  if (N<1) return {};
  if (N==1) return {Qideal(Quad::one)};
  //  cout<<" looking for ideals of norm "<<N<<endl;
  map<QUINT, vector<Qideal>>::iterator I_N = N_to_Ilist.find(N);
  if (I_N!=N_to_Ilist.end())
    {
      vector<Qideal> Ilist = I_N->second;
      //      cout<<" found "<<Ilist.size()<<" ideals: "<<Ilist<<endl;
      if (both_conj)
        return Ilist;
      else
        return remove_conjugates(Ilist);
    }

  // now we compute and cache the ideals of norm N

#ifdef DEBUG_SORT
  cout<<"Computing sorted list of ideals with norm "<<N<<endl<<flush;
#endif
  if (N==1)
    {
      return (N_to_Ilist[N] = {Qideal()});
    }

  vector<Qideal> ans;
  Qideal I;
  vector<QUINT> pp = pdivs(N);
  int np = pp.size();

  if (np==1) // prime power case
    {
      QUINT p = pp[0];
      long e = val(p,N);
#ifdef DEBUG_SORT
      cout<<"Constructing ideals of norm "<<N<<"="<<p<<"^"<<e<<", chi(p) = "<<Quad::chi(p)<<endl;
#endif
      switch (Quad::chi(p))
        {
        case -1:
          if (e%2)
            return (N_to_Ilist[N] = {});
          I = Qideal(pow(p,(e/2)));
          I.set_index(1);
          return (N_to_Ilist[N] = {I});
        case 0:
          if (e%2==0)
            {
              I = Qideal(pow(p,(e/2)));
            }
          else
            {
              I = Quadprimes_above(I2long(p))[0];
              if (e>1)
                I *= pow(p,((e-1)/2));
            }
          I.set_index(1);
          ans.push_back(I);
          return (N_to_Ilist[N] = {I});
        case +1:
          vector<Quadprime> PP = Quadprimes_above(I2long(p));
          Qideal P = PP[0], Q=PP[1];
          ans = {P, Q};
          for (long i=1; i<e; i++)
            {
              I = Q*ans[i]; // = Q^(i+1)
              for (vector<Qideal>::iterator PQ=ans.begin(); PQ!=ans.end(); PQ++)
                (*PQ) *= P;
              ans.push_back(I);
            }
#ifdef DEBUG_SORT
            {
              cerr<<"**********************Constructed list of "<<ans.size()<<" ideals: ";
              for (vector<Qideal>::iterator J=ans.begin(); J!=ans.end(); ++ans)
                cerr<<*J<<" ";
              cerr<<endl;
            }
#endif
          for (long i=0; i<=e; i++)
            ans[i].set_index(i+1);
#ifdef DEBUG_SORT
          cout<<"(split case) sorted list of ideals with norm "<<N<<": "<<ans<<endl;
#endif
          N_to_Ilist[N] = ans;
          if (!both_conj)
            ans = remove_conjugates(ans);;
          return ans;
        }
    } // end of prime power case

  // General case, use recursion

  // II is a list of sorted lists of ideals of prime power norm, the
  // outer list sorted by size of the underlying prime (not the prime
  // power):
  vector<vector<Qideal>> II;
  for (vector<QUINT>::const_iterator pi=pp.begin(); pi!=pp.end(); ++pi)
    II.push_back(ideals_with_norm(pow(*pi, val(*pi, N))));

  ans = {Qideal()}; // unit ideal

  // "merge" lexicographically
  for (vector<vector<Qideal> >::const_reverse_iterator QQ = II.crbegin(); QQ!=II.crend(); ++QQ)
    {
      vector<Qideal> ans2;
      for (vector<Qideal>::const_iterator Qi = QQ->begin(); Qi!=QQ->end(); ++Qi)
        {
          Qideal Q = *Qi;
          for (vector<Qideal>::const_iterator I0 = ans.begin(); I0!=ans.end(); ++I0)
            {
              I = *I0;
              ans2.push_back(I*Q);
            }
        }
      ans = ans2;
    }
  long i=1;
  for (vector<Qideal>::iterator Ii=ans.begin(); Ii!=ans.end(); ++Ii, ++i)
    Ii->set_index(i);

#ifdef DEBUG_SORT
  cout<<"Sorted list of ideals with norm "<<N<<": "<<ans<<endl;
#endif
  N_to_Ilist[N] = ans;
  if (!both_conj)
    ans = remove_conjugates(ans);;
  return ans;
}

vector<Qideal> Qideal_lists::ideals_with_bounded_norm(QUINT maxnorm, int both_conj)
{
  vector<Qideal> ans;
  for (QUINT N(1); N<=maxnorm; N++)
    {
      vector<Qideal> I_N = ideals_with_norm(N);
      ans.insert(ans.end(), I_N.begin(), I_N.end());
    }
  return ans;
}

// return J coprime to N such that d*J=c*this, or (if anti=1) J such that J*this=(c) and d=1
Qideal Qideal::equivalent_coprime_to(const Qideal& N, Quad& c, Quad& d, int anti)
{
  if (is_principal()) // this = (g0) so return (1)
    {
      if (anti)
        {
          c = g0;
          d = Quad::one;
        }
      else
        {
          c = Quad::one;
          d = g0;
        }
      return Qideal(Quad::one);
    }
  Quad g;
  // cout<<"looking for a prime equivalent to "<<(*this)<<" which is coprime to "<<N<<endl;
  for (vector<Quadprime>::iterator Pi = Quadprimes::list.begin(); Pi != Quadprimes::list.end(); ++Pi)
    {
      Quadprime P = *Pi;
      if (P.residue_degree()==2) continue;  // inert primes are principal so no use
      if (!P.is_coprime_to(N)) continue;    // skip P unless it is coprime to N
      Qideal I = (anti? P: P.conj());
      I *= (*this);
      if (I.is_principal(g))
        {
          if (anti)
            {
              c = g; // P*this = I = (g)
              d = Quad::one;
              assert ((*this)*P==Qideal(c));
            }
          else
            {
              d = g;
              c = P.norm(); // c*this = d*P
            }
          return P;
        }
    }
  cerr << "Unable to find an ideal equivalent to "<<(*this)<<" coprime to "<<N<<endl;
  return Qideal();
}

// return J coprime to N such that J^2*this is principal; if no such
// J exists (i.e., if the ideal class is not a square, return the
// zero ideal.  (implemented in primes.cc)
Qideal Qideal::sqrt_coprime_to(const Qideal& N)
{
  if (is_principal()) // this = (g0) so return (1)
    {
      return Qideal(Quad::one);
    }
  Qideal A = sqrt_class(1); // so A^2*this is principal
  if (A.nm==0) return A;
  Qidealooper looper(2, 1000);
  while (looper.not_finished())
    {
      Qideal J = looper.next();
      while (!A.is_equivalent(J))
        J = looper.next();
      if (J.is_coprime_to(N))
        {
          assert ((J*J*(*this)).is_principal());
          assert (J.is_coprime_to(N));
          return J;
        }
    }
  cerr << "Unable to find an ideal I such that I^2 * "<<(*this)<<" is principal and coprime to "<<N<<endl;
  return Qideal();
}

Factorization Qideal::factorization() // sets F if necessary then returns F
{
  if (F==0)
    {
      F = new Factorization(*this);
    }
  return *F;
}

int Qideal::is_square()
{
  vector<int>ee = factorization().exponents();
  for (vector<int>::const_iterator ei = ee.begin(); ei!=ee.end(); ++ei)
    if ((*ei)%2==1) return 0;
  return 1;
}

long Qideal::genus_class()
{
  long c = 0;
  vector<QuadprimePower> PP = factorization().prime_powers();
  for (auto PPi = PP.begin(); PPi!=PP.end(); ++PPi)
    {
      if ((PPi->second)%2==1)
        c ^= (PPi->first).genus_class();
    }
  return c;
}

vector<int> Qideal::genus_character()
{
  return bits(genus_class(), Quad::discfactors.size());
}

// Test whether an ideal is a prime, or a prime power:
int Qideal::is_prime() {return factorization().is_prime();}
int Qideal::is_prime_power() {return factorization().is_prime_power();}

// END OF FILE primes.cc
