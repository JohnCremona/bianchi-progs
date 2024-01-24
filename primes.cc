// FILE primes.cc

#include "primes.h"
#include "qidloop.h"
#include "intprocs.h"
#include <list>

//Definition of static data members of class Quadprimes:
INT Quadprimes::maxnorm;
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

  // It returns a list of values in {0,1}, one per prime discriminant
  // factor, which add to 0 mod 2

  int r = Quad::prime_disc_factors.size()-1; // = Quad::class_group_2_rank
  vector<int> ans(1+r, 0);
  if (r==0 || is_inert() || is_principal())
    return ans;
  // cout<<" -- not principal, even class number..."<<endl;
  int i=0, i0=-1, tot=0, k;
  for (const auto& di : Quad::prime_disc_factors)
    {
      long d = I2long(di);
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
      i++;
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

int Quadprime::genus_character(const INT& D)
{
  if (!div_disc(D, Quad::disc))
    return 0;
  //cout<<"In P.genus_character(D) with P="<<(*this)<<", D="<<D;
  vector<int> v1 = chardisc(D), v2 = genus_character();
  //cout<<" with character values "<< v1 << " for D and "<< v2 << " for P"<<endl;
  int dot = std::inner_product(v1.begin(), v1.end(), v2.begin(), 0) % 2;
  dot = (dot? -1: +1);
  //cout<<"dot product (mod 2) = "<<dot<<" --> "<< dot <<endl;
  return dot;
}

long Quadprime::genus_class(int contract)
{
  vector<int> v = genus_character();
  long c = from_bits(v);
  if (contract) c %= (1<<Quad::class_group_2_rank);
  return c;
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
  for ( auto& P : Quadprimes::list)
    {
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
          vector<int> cc = P.genus_character();
          cout << ": ";
          for (const auto& c : cc)
            cout << (c? "-": "+");
        }
      cout<<")"<<endl;
    }
}

vector<Quadprime> Quadprimes_above(long p) // p should be an integer prime
{
  long d=Quad::d; INT disc=Quad::disc;
  int t=Quad::t;
  vector<Quadprime> Plist;
  INT zero(0), one(1), P(p);

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
      INT r, b1, b2, MD(posmod(-d,p));
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
  vector<INT> pdivs_norm = pdivs(I.norm());
  //  cout<<"Finding prime factors of "<<I<<" with norm "<<I.norm()<<", primes dividing norm are "<<pdivs_norm<<endl;
  for( const auto& p : pdivs_norm)
    {
      vector<Quadprime> PP = Quadprimes_above(I2long(p));
      //      cout<<"primes above "<<p<<" are "<<PP<<endl;
      // at least one, but possibly not both when p splits, divides I
      for( const auto& P : PP)
        {
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
  for (const auto& Qi : Qlist)
    plist.push_back(Qi.first);
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
  for (const auto& Qi : Qlist)
    elist.push_back(Qi.second);
  return elist;
}

vector<Quadprime> pdivs(Qideal& I)  // list of prime divisors
{
  return I.factorization().sorted_primes();
}

int npdivs(Qideal& I)  // number of prime ideal divisors
{
  return I.factorization().size();
}

int ndivs(Qideal& I) // number of ideal divisors
{
  vector<int> ee = I.factorization().exponents();
  long nd = 1;
  for (auto e : ee)
    nd*=(1+e);
  return nd;
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
  for(const auto& Q: F.Qlist)
    {
      Qideal P = Q.first;
      int e = Q.second;
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
      CRT_vector.push_back(Quad::one);
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
        assert(Q.divides(CRT_vector[j]-INT(i==j)));
    }
}

Quad Factorization::solve_CRT(const vector<Quad>& v) // solution to x=v[i] mod Qlist[i]
{
  if (CRT_vector.size()==0) init_CRT();
  Quad a = Quad::zero;
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
  list.clear();
  vector<Quadprime> list1, list2;

  // First fill up lists of degree 1 and degree 2 primes
  //  cout<<"Computing list of prime ideals of norm up to "<<maxnorm<<endl;

  for (primevar pr; pr.ok()&&pr<=maxn; ++pr)
    { long p=pr;
      //            cout<<"p = "<<p<<endl;
      vector<Quadprime> PP = Quadprimes_above(p);
      //            cout<<" primes above: "<< PP<<endl;
      for(const auto& P : PP)
        {
          INT q = P.norm();
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
  for(const auto& p : plist)
    {
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

map<INT, vector<Qideal>> Qideal_lists::N_to_Ilist;

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

vector<Qideal> Qideal_lists::ideals_with_norm(INT N, int both_conj)
{
  if (N<1) return {};
  if (N==1) return {Qideal(Quad::one)};
  //  cout<<" looking for ideals of norm "<<N<<endl;
  map<INT, vector<Qideal>>::iterator I_N = N_to_Ilist.find(N);
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
  vector<INT> pp = pdivs(N);
  int np = pp.size();

  if (np==1) // prime power case
    {
      INT p = pp[0];
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
              for ( auto& PQ : ans)
                PQ *= P;
              ans.push_back(I);
            }
#ifdef DEBUG_SORT
            {
              cerr<<"**********************Constructed list of "<<ans.size()<<" ideals: ";
              cerr<<ans<<endl;
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
  for ( const auto& p : pp)
    II.push_back(ideals_with_norm(pow(p, val(p, N))));

  ans = {Qideal()}; // unit ideal

  // "merge" lexicographically
  for (vector<vector<Qideal> >::const_reverse_iterator QQ = II.crbegin(); QQ!=II.crend(); ++QQ)
    {
      vector<Qideal> ans2;
      for ( auto Q : *QQ)
        for ( const auto& I : ans)
          ans2.push_back(I*Q);
      ans = ans2;
    }
  long i=1;
  for ( auto& I : ans)
    I.set_index(i++);

#ifdef DEBUG_SORT
  cout<<"Sorted list of ideals with norm "<<N<<": "<<ans<<endl;
#endif
  N_to_Ilist[N] = ans;
  if (!both_conj)
    ans = remove_conjugates(ans);;
  return ans;
}

vector<Qideal> Qideal_lists::ideals_with_bounded_norm(INT maxnorm, int both_conj)
{
  vector<Qideal> ans;
  for (INT N(1); N<=maxnorm; N+=1)
    {
      vector<Qideal> I_N = ideals_with_norm(N, both_conj);
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
  //cout<<"looking for a prime equivalent to "<<(*this)<<" which is coprime to "<<N<<endl;
  for ( auto& P : Quadprimes::list)
    {
      //cout<<" -testing P="<<P<<endl;
      if (P.residue_degree()==2) continue;  // inert primes are principal so no use
      if (!P.is_coprime_to(N)) continue;    // skip P unless it is coprime to N
      //cout<<" -passes first tests, now checking its class"<<endl;
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
          //cout<<" ++OK, returning "<<P<<endl;
          return P;
        }
      else
        {
          //cout<<" --no good, continuing"<<endl;
        }
    }
  cerr << "\nUnable to find an ideal equivalent to "<<(*this)<<" coprime to "<<N<<endl;
  cerr << "More primes need to be initialised! Currently max norm is "<<Quad::maxnorm<<endl;
  assert (0 && "failure in Qideal::equivalent_coprime_to()");
  return Qideal();
}

// return J coprime to N such that (J/this)^2, or (J*this)^2 if anti, is principal
Qideal Qideal::equivalent_mod_2_coprime_to(const Qideal& N, int anti)
{
  Qideal I = (*this) * (*this);
  if (I.is_principal())
    {
      return Qideal(Quad::one);
    }
  Quad g;
  //cout<<"looking for a prime equivalent to "<<(*this)<<" mod 2-torsion which is coprime to "<<N<<endl;
  for ( auto& P : Quadprimes::list)
    {
      //cout<<" -testing P="<<P<<endl;
      if (P.residue_degree()==2) continue;  // inert primes are principal so no use
      if (!P.is_coprime_to(N)) continue;    // skip P unless it is coprime to N
      //cout<<" -passes first tests, now checking its class"<<endl;
      Qideal J = (anti? P: P.conj());
      J = I*J*J;
      if (J.is_principal(g))
        return P;
    }
  cerr << "\nUnable to find an ideal equivalent mod 2-torsion to "<<(*this)<<" coprime to "<<N<<endl;
  cerr << "More primes need to be initialised! Currently max norm is "<<Quad::maxnorm<<endl;
  assert (0 && "failure in Qideal::equivalent_mod_2_coprime_to()");
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

QuadprimeLooper::QuadprimeLooper(Qideal level)
  :Pi(Quadprimes::list.begin()), N(level)
{
  P = *Pi;
  while (P.divides(N) && ok())
    {
      ++Pi;
      P = *Pi;
    }
}

void QuadprimeLooper::operator++()
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

void QuadprimeLooper::reset()
{
  Pi = Quadprimes::list.begin();
  P = *Pi;
  while (P.divides(N)) P = *Pi++;
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
  for (auto e : ee)
    if (e%2==1) return 0;
  return 1;
}

long Qideal::genus_class(int contract)
{
  long c = 0;
  vector<QuadprimePower> PP = factorization().prime_powers();
  for ( auto& P : PP)
    {
      if ((P.second)%2==1)
        c ^= (P.first).genus_class(); // binary + (exclusive or)
    }
  if (contract) c %= (1<<Quad::class_group_2_rank);
  return c;
}

int Qideal::genus_character(const INT& D) // one unram char value
{
  if (!div_disc(D, Quad::disc))
    return 0;
  int r = Quad::prime_disc_factors.size();
  int dot = dotbits(from_bits(chardisc(D)), genus_class(), r);
  return (dot? -1: +1);
}

vector<int> Qideal::genus_character()
{
  return bits(genus_class(), Quad::prime_disc_factors.size());
}

vector<INT> Qideal::possible_unramified_twists()  // sublist of Quad::all_disc_factors() consisting of those D not 1
                                                    // for which chi_D(Q)=+1 for all prime powers Q||N
{
  vector<INT> discs;
  if (Quad::class_group_2_rank == 0)
    return discs;
  vector<QuadprimePower> Qlist = Factorization(*this).prime_powers();
  for( const auto& D : Quad::all_disc_factors)
    {
      if (D==1)
        continue;
      int ok = 1;
      for ( auto& Q : Qlist)
        ok = ok && (((Q.second)%2==0) || ((Q.first).genus_character(D)==+1));
      if (ok)
        discs.push_back(D);
    }
  return discs;
}

// Test whether an ideal is a prime, or a prime power:
int Qideal::is_prime() {return factorization().is_prime();}
int Qideal::is_prime_power() {return factorization().is_prime_power();}


// Return {-m,m} where m is the largest integer <= +2*sqrt(N(P)), the bounds on a(P)
pair<long,long> eigenvalue_range(Quadprime& P)
{
  long normp = I2long(P.norm());
  long aplim=2;
  while (aplim*aplim<=4*normp) aplim++;
  aplim--;
  return {-aplim, aplim};
}

// Return {-m,3*m} whereor m = N(P), the bounds on a(P^2)=a(P)^2-N(P)
pair<long,long> eigenvalue_sq_range(Quadprime& P)
{
  long normp = I2long(P.norm());
  return {-normp, 3*normp};
}

// Return list of integers between -2*sqrt(N(P)) and +2*sqrt(N(P))
vector<long> good_eigrange(Quadprime& P)
{
  long normp = I2long(P.norm());
  if (P.has_square_class())
    {
      pair<long,long> apbounds = eigenvalue_range(P);
      return range(apbounds.first, apbounds.second);
    }
  else // want eigs of T(P^2)=T(P)^2-N(P) such that T(P) has integral eig
    {
      pair<long,long> apbounds = eigenvalue_range(P);
      vector<long> ans = range(0,apbounds.second);
      for_each(ans.begin(), ans.end(), [normp](long& a){a = a*a-normp;});
      return ans;
    }
}

// compute a list of primes Q dividing N with Q^e||N such that [Q^e] is square
vector<Quadprime> make_squarebadprimes(Qideal& N, const vector<Quadprime>& badprimes)
{
  if (Quad::class_group_2_rank==0)
    return badprimes;
  vector<Quadprime> squarebadprimes;
  for ( auto Q : badprimes)
    {
      if (val(Q,N)%2==0 || Q.has_square_class()) // then [Q^e] is a square
        squarebadprimes.push_back(Q);
    }
  return squarebadprimes;
}

// compute a list of at least nap good primes (excluding those
// dividing characteristic if >0), to include at least on principal
// one which has index iP0;
vector<Quadprime> make_goodprimes(Qideal& N,  int np, int& iP0, int p)
{
  vector<Quadprime> goodprimes;
  QuadprimeLooper L(p==0? N : INT(long(p))*N);
  iP0=-1;
  for (int i=0; (i<np) || (iP0<0); i++, ++L)
    {
      Quadprime P = L;
      goodprimes.push_back(P);
      if (P.is_principal() && iP0==-1)
        iP0 = goodprimes.size()-1;
    }
  return goodprimes;
}

// END OF FILE primes.cc
