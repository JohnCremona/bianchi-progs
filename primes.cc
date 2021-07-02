// FILE primes.cc

#include "primes.h"
#include "intprocs.h"

//Definition of static data members of class Quadprimes:
long Quadprimes::maxnorm;
vector<Quadprime> Quadprimes::list;

void Quadprimes::display(ostream& s)
{
  s << list.size() << " primes initialised, ";
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
          Plist.push_back(Quadprime(2,1,1, 2, 0));
          break;
        case 2: // ramified, (2) = (2,w)^2
          Plist.push_back(Quadprime(2,0,1, 2, 0));
          break;
        case 3: // split, (2) = (2,w)*(2,1+w) or inert, (2) = (2)
          if (d%8==3) // inert
            Plist.push_back(Quadprime(1,0,2, 2, 0));
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
        Plist.push_back(Quadprime(p,(p+1)/2,1, p,0));
      else
        Plist.push_back(Quadprime(p,0,1, p,0));
      return Plist;
    }

  if (kronecker(disc,p) == -1) //inert
    {
      //      cout << "odd inert case" << endl;
      Plist.push_back(Quadprime(1,0,p, p,0));
    }
  else
    // split, (p) = (p,b1+w)*(p,b2+w) where b1, b2 are roots of x^2+t*x+n=0 (modp)
    // We order so the HNFs are [p,b,1], [p,b',1] with b<b'
    {
      //      cout << "odd split case" << endl;
      long r, b1, b2;
      if (t)
        {
          r = sqrt_mod_p(-disc, p);
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
      Plist.push_back(Quadprime(p,b2,1, p, 1));
    }
  return Plist;
}


Quadprime primdiv(const Qideal& a)    // returns one prime divisor
{
  vector<long> pdivs_norm = pdivs(a.norm());
  for(vector<long>::const_iterator pi = pdivs_norm.begin(); pi!=pdivs_norm.end(); pi++)
    {
      vector<Quadprime> PP = Quadprimes_above(*pi);
      // at least one, but possibly not both when p splits, divides a
      if ((PP.size()==1) || (PP[0].divides(a)))
        return PP[0];
      else
        return PP[1];
    }
  cerr<<"Error: primdiv() called for a="<<a<<" which has no prime ideal divisors"<<endl;
  return Quadprime();
}

prime_factn::prime_factn(const Qideal& n)
{
  Qideal a=n;
  long np = 0;
  vector<long> pdivs_norm = pdivs(a.norm());
  //  cout<<"Finding prime factors of "<<a<<" with norm "<<a.norm()<<", primes dividing norm are "<<pdivs_norm<<endl;
  for(vector<long>::const_iterator pi = pdivs_norm.begin(); pi!=pdivs_norm.end(); pi++)
    {
      vector<Quadprime> PP = Quadprimes_above(*pi);
      //      cout<<"primes above "<<(*pi)<<" are "<<PP<<endl;
      // at least one, but possibly not both when p splits, divides a
      for(vector<Quadprime>::const_iterator Pi = PP.begin(); Pi!=PP.end(); Pi++)
        {
          Quadprime P = *Pi;
          if (P.divides(a))
            {
              plist.push_back(P);
              elist.push_back(1);
              a /= P;
              while (P.divides(a))
                {
                  a/=P;
                  elist[np]+=1;
                }
              np +=1;
            }
        }
    }
}

void prime_factn::display(ostream& s) const
{
  if (plist.size()==0)
    {
      s << "No prime factors." << endl;
    }
  else
    {
      s << "Number of distinct prime factors: "<< plist.size()<<endl;
      s << "Exponents   Primes:" << endl;
      for (long i=0; i<(long)plist.size(); i++)
	s << elist[i] << "\t" << plist[i] << endl;
    }
}

vector<Quadprime> pdivs(const Qideal& n)  // list of prime divisors
{
  return prime_factn(n).plist;
}

vector<Qideal> alldivs(const Qideal& a)    // list of all ideal divisors
{
  prime_factn pp(a);
  long np = pp.num_primes();
  long nu = 1; long nd=nu;
  for (long i=0; i<np; i++) {nd*=(1+pp.expo(i));}
  vector<Qideal> dlist(nd);
  dlist[0]=1;
  nd=nu;
  Qideal p;
  long e, j, k;
  vector<Quadprime>::const_iterator pi;
  vector<long>::const_iterator ei;
  for(pi = pp.plist.begin(), ei = pp.elist.begin();
      pi != pp.plist.end();
      pi++, ei++)
    {
      p = *pi;
      e = *ei;
      for (j=0; j<e; j++)
	for (k=0; k<nd; k++)
	  dlist[nd*(j+1)+k] = (p*dlist[nd*j+k]);
      nd*=(e+1);
    }
  return dlist;
}

void Quadprimes::init(long maxn)
{
  maxnorm = maxn;
  vector<Quadprime> list1, list2;
  Quadprime P, Q;

  // First fill up lists of degree 1 and degree 2 primes
  //  cout<<"Computing list of prime ideals of norm up to "<<maxnorm<<endl;

  for (primevar pr; pr.ok()&&pr<=maxnorm; pr++)
    { long p=pr;
      //      cout<<"p = "<<p<<endl;
      vector<Quadprime> PP = Quadprimes_above(p);
      //      cout<<" primes above: "<< PP<<endl;
      for(vector<Quadprime>::const_iterator Pi = PP.begin(); Pi!=PP.end(); Pi++)
        {
          P = *Pi;
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
      P = *Pi; Q = *Qi;
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
  prime_factn pp(a);
  Qideal p;
  long np = pp.num_primes();

  long nd=1;
  for(long i=0; i<np; i++) { nd *= ( 1+ pp.expo(i)/2 ) ;}

  vector<Qideal> dlist(nd);
  dlist[0]=1;
  nd=1;
  long e;
  for(long i=0; i<np; i++)
    {
      p = pp.prime(i);
      e = pp.expo(i)/2;
      for(long j=0; j<e; j++)
	for(long k=0; k<nd; k++)
	  dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
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

// END OF FILE primes.cc
