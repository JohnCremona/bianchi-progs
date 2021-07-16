// FILE primes.cc

#include "primes.h"
#include "intprocs.h"
#include <list>

//Definition of static data members of class Quadprimes:
long Quadprimes::maxnorm;
vector<Quadprime> Quadprimes::list;

void Quadprimes::display(ostream& s)
{
  s << list.size() << " prime ideals initialised, ";
  s << "max norm = " << maxnorm << endl;
}

vector<Quadprime> Quadprimes_above(long p) // p should be an integer prime
{
  long d=Quad::d, disc=-Quad::disc, t=Quad::t;
  vector<Quadprime> Plist;

  if (p==2) // treat as special case
    {
      switch (d%4)
        {
        case 1: // ramified, (2) = (2,1+w)^2
          Plist.push_back(Quadprime(2,1,1, 2));
          break;
        case 2: // ramified, (2) = (2,w)^2
          Plist.push_back(Quadprime(2,0,1, 2));
          break;
        case 3: // split, (2) = (2,w)*(2,1+w) or inert, (2) = (2)
          if (d%8==3) // inert
            {
              Plist.push_back(Quadprime(1,0,2, 2));
            }
          else // split
            {
              Plist.push_back(Quadprime(2,0,1, 2, 1));
              Plist.push_back(Quadprime(2,1,1, 2, 2));
            }
        }
      return Plist;
    }

  // odd p

  if (d%p==0) // ramified, (p) = (p,b+w)^2 where b=t/2 mod p
    {
      //      cout << "odd ramified case" << endl;
      if (t)
        Plist.push_back(Quadprime(p,-(p+1)/2,1, p));
      else
        Plist.push_back(Quadprime(p,0,1, p));
      return Plist;
    }

  if (kronecker(disc,p) == -1) //inert
    {
      //      cout << "odd inert case" << endl;
      Plist.push_back(Quadprime(1,0,p, p));
    }
  else
    // split, (p) = (p,b1+w)*(p,b2+w) where b1, b2 are roots of x^2+t*x+n=0 (modp)
    // We order so the HNFs are [p,b,1], [p,b',1] with b<b'
    {
      //      cout << "odd split case" << endl;
      long r, b1, b2;
      if (t)
        {
          r = sqrt_mod_p(disc, p);
          b1 = posmod((r%2? (r-1)/2: (p+r-1)/2), p);
          b2 = posmod(-1-b1, p);
        }
      else
        {
          r = sqrt_mod_p(-d, p);
          b1 = posmod(r, p);
          b2 = posmod(-b1, p);
        }
      if (b1>b2) swap(b1,b2);
      Plist.push_back(Quadprime(p,b1,1, p, 1));
      Plist.push_back(Quadprime(p,b2,1, p, 2));
    }
  return Plist;
}

Factorization::Factorization(const Qideal& II)
{
  I = Qideal(II);
  long np = 0;
  int e;
  vector<long> pdivs_norm = pdivs(I.norm());
  //  cout<<"Finding prime factors of "<<I<<" with norm "<<I.norm()<<", primes dividing norm are "<<pdivs_norm<<endl;
  Qideal J = I; // so we can divide out primes while leaving I unchanged
  for(vector<long>::const_iterator pi = pdivs_norm.begin(); pi!=pdivs_norm.end(); pi++)
    {
      vector<Quadprime> PP = Quadprimes_above(*pi);
      //      cout<<"primes above "<<(*pi)<<" are "<<PP<<endl;
      // at least one, but possibly not both when p splits, divides I
      for(vector<Quadprime>::const_iterator Pi = PP.begin(); Pi!=PP.end(); Pi++)
        {
          Quadprime P = *Pi;
          if (P.divides(J))
            {
              e = 1;
              J /= P;
              while (P.divides(J))
                {
                  e++;
                  J /= P;
                }
              Qlist.push_back({P,e});
              np +=1;
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
  std::transform(Qlist.begin(), Qlist.end(), back_inserter(plist),
                 [](const QuadprimePower& Q) -> Quadprime { return Q.first; });
  return plist;
}

vector<int> Factorization::exponents() const
{
  vector<int> elist;
  elist.reserve(size());
  std::transform(Qlist.begin(), Qlist.end(), back_inserter(elist),
                 [](const QuadprimePower& Q) -> int { return Q.second; });
  return elist;
}

vector<Quadprime> pdivs(const Qideal& I)  // list of prime divisors
{
  Factorization F(I);
  return F.primes();
}

vector<Qideal> alldivs(const Qideal& a)    // list of all ideal divisors
{
  Factorization F(a);
  int np = F.size();
  long nu = 1; long nd=nu;
  for (long i=0; i<np; i++) {nd*=(1+F.exponent(i));}
  vector<Qideal> dlist(nd);
  dlist[0]=1;
  nd=nu;
  Qideal P;
  int e, j, k;
  vector<QuadprimePower>::const_iterator Qi;
  for(Qi = F.Qlist.begin();  Qi != F.Qlist.end();  Qi++)
    {
      P = Qi->first;
      e = Qi->second;
      for (j=0; j<e; j++)
	for (k=0; k<nd; k++)
	  dlist[nd*(j+1)+k] = (P*dlist[nd*j+k]);
      nd*=(e+1);
    }
  return dlist;
}

Qideal Factorization::prime_power(int i) const
{
  Qideal P = Qlist[i].first, Q;
  int j;
  for (j=1, Q=P; j<Qlist[i].second; j++) Q *= P;
  return Q;
}

void Factorization::init_CRT()              // compute the CRT vector
{
  if (CRT_vector.size()!=0) return; // already done
  CRT_vector.reserve(size());
  if (size()==1)
    {
      CRT_vector[0] = 1;
      return;
    }
  Qideal Q, J;
  Quad r, s;
  for (int i=0; i<size(); i++)
    {
      Q = prime_power(i);
      J = I/Q;
      int t = J.is_coprime_to(Q, r, s); // r+s=1 with r in J, s in Q
                                        // so r=1 mod Q, r=0 mod Q' for other Q'
      assert (t==1);
      CRT_vector[i] = r;
      assert (Q.contains(r-1));
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
  Quad a = 0;
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

  for (primevar pr; pr.ok()&&pr<=maxnorm; pr++)
    { long p=pr;
      //      cout<<"p = "<<p<<endl;
      vector<Quadprime> PP = Quadprimes_above(p);
      //      cout<<" primes above: "<< PP<<endl;
      for(vector<Quadprime>::const_iterator Pi = PP.begin(); Pi!=PP.end(); Pi++)
        {
          Quadprime P = *Pi;
          long q = P.norm();
          //          cout<<"P = "<<P<<" with norm "<<q<<endl;
          if(q==p) // degree 1 prime
            list1.push_back(P);
          else
            if(q<=maxnorm)
              list2.push_back(P);
        }
    }

  //  cout<<" - found "<<list1.size() << " degree 1 primes and "<<list2.size()<<" degree 2 primes"<<endl;

  // Now merge these into a single list sorted by norm

  //  cout<<" - merging into a single list" <<endl;
  vector<Quadprime>::const_iterator Pi = list1.begin(), Qi = list2.begin();
  while(Pi!=list1.end() && Qi!=list2.end())
    {
      Quadprime P = *Pi, Q = *Qi;
      if(P.norm()<Q.norm())
        {
          list.push_back(P);
          Pi++;
        }
      else
        {
          list.push_back(Q);
          Qi++;
        }
    }
  // only one of the following will do anything, probably the first:
  while(Pi!=list1.end())
    list.push_back(*Pi++);
  while(Qi!=list2.end())
    list.push_back(*Qi++);
}

// need divisors functions etc

vector<Qideal> sqdivs(const Qideal& a) // all divisors whose square divides a, up to +/-
{
  Factorization F(a);
  Qideal P;
  int np = F.size();

  long nd=1;
  for(long i=0; i<np; i++) { nd *= ( 1+ F.exponent(i)/2 ) ;}

  vector<Qideal> dlist(nd);
  dlist[0]=1;
  nd=1;
  long e;
  for(long i=0; i<np; i++)
    {
      P = F.prime(i);
      e = F.exponent(i)/2;
      for(long j=0; j<e; j++)
	for(long k=0; k<nd; k++)
	  dlist[nd*(j+1)+k] = P*dlist[nd*j+k];
      nd*=(e+1);
    }
  return dlist;
}

vector<Qideal> sqfreedivs(const Qideal& a)       // all square-free divisors
{
  vector<Quadprime> plist=pdivs(a);
  Qideal p;
  long np = plist.size();
  long nd = 1;
  while (np-->0) nd*=2;
  vector<Qideal> dlist(nd);
  dlist[0]=1;
  nd=1;
  for(vector<Quadprime>::const_iterator pr = plist.begin(); pr != plist.end(); pr++)
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

map<long, vector<Qideal>> Qideal_lists::N_to_Ilist;


vector<Qideal> Qideal_lists::ideals_with_norm(long N)
{
  if (N<1) return vector<Qideal>();
  map<long, vector<Qideal>>::iterator I_N = N_to_Ilist.find(N);
  if (I_N!=N_to_Ilist.end()) return I_N->second;

  // now we compute and cache the ideals of norm N

#ifdef DEBUG_SORT
  cout<<"Computing sorted list of ideals with norm "<<N<<endl<<flush;
#endif
  vector<Qideal> ans;
  Qideal I;

  if (N==1)
    {
      return (N_to_Ilist[N] = {Qideal()});
    }

  vector<long> pp = pdivs(N);
  int np = pp.size();

  if (np==1) // prime power case
    {
      long p = pp[0];
      long e = val(p,N);
#ifdef DEBUG_SORT
      cout<<"Constructing ideals of norm "<<N<<"="<<p<<"^"<<e<<endl;
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
              I = Quadprimes_above(p)[0];
              if (e>1)
                I *= pow(p,((e-1)/2));
            }
          I.set_index(1);
          ans.push_back(I);
          return (N_to_Ilist[N] = {I});
        case +1:
          Qideal P = Quadprimes_above(p)[0];
          std::list<Qideal> II; long k;
          if (e%2)
            {
              k = (e-1)/2;
              I = P;
            }
          else
            {
              k = (e-2)/2;
              I = P*P;
              II.push_back(Qideal(pow(p,k+1)));
            }
          P *= P;
          while (k>=0)
            {
              Qideal J = I*pow(p,k);
              II.push_front(J);
              II.push_back(J.conj());
              I *= P;
              k -= 1;
            }
#ifdef DEBUG_SORT
          cout<<"Constructed list of "<<II.size()<<" ideals"<<endl;
#endif
          ans.insert(ans.end(), II.begin(), II.end());
          for (long i=0; i<=e; i++)
            ans[i].set_index(i+1);
#ifdef DEBUG_SORT
          cout<<"(split case) sorted list of ideals with norm "<<N<<": "<<ans<<endl;
#endif
          return (N_to_Ilist[N] = ans);
        }
    } // end of prime power case

  // General case, use recursion

  vector<long> pplist; // list of prime power factors of N, will be sorted by size
  for (vector<long>::const_iterator pi=pp.begin(); pi!=pp.end(); pi++)
    {
      long p = *pi;
      long e = val(p,N);
      pplist.push_back(pow(p,e));
    }
  std::sort(pplist.begin(), pplist.end());

  vector<vector<Qideal>> II;
  for (vector<long>::const_iterator pi=pplist.begin(); pi!=pplist.end(); pi++)
    II.push_back(ideals_with_norm(pow(*pi, val(*pi, N))));

  ans = {Qideal()}; // unit ideal

  // "merge" lexicographically
  while (!II.empty())
    {
      vector<Qideal> ans0 = ans;
      vector<Qideal> ans1 = II.back();
      II.pop_back();
      vector<Qideal> ans2;
      for (vector<Qideal>::const_iterator I0 = ans0.begin(); I0!=ans0.end(); I0++)
        for (vector<Qideal>::const_iterator I1 = ans1.begin(); I1!=ans1.end(); I1++)
          ans2.push_back((*I0)*(*I1));
      ans = ans2;
    }
  long i=1;
  for (vector<Qideal>::iterator Ii=ans.begin(); Ii!=ans.end(); Ii++, i++)
    Ii->set_index(i);

#ifdef DEBUG_SORT
  cout<<"Sorted list of ideals with norm "<<N<<": "<<ans<<endl;
#endif
  return (N_to_Ilist[N] = ans);
}

vector<Qideal> Qideal_lists::ideals_with_bounded_norm(long maxnorm)
{
  vector<Qideal> ans;
  for (long N=1; N<=maxnorm; N++)
    {
      vector<Qideal> I_N = ideals_with_norm(N);
      ans.insert(ans.end(), I_N.begin(), I_N.end());
    }
  return ans;
}




// END OF FILE primes.cc
