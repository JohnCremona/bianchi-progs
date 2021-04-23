// FILE primes.cc

#include "primes.h"

//Definition of static data members of class Quadprimes:
long Quadprimes::maxnorm;
Primelist Quadprimes::list;

void Quadprimes::display(ostream& s)
{
  s << list.getlength() << " primes initialised,";
  s << "max norm = " << maxnorm << endl;
}

Quadprime primdiv(const PRIMETYPE& a)    // returns "first" prime divisor
{
    Quadprime ans,q;
    int found=0;
    long normlimit = a.norm();
    if (normlimit<2)                // return 0;   // must return something!
      {
	cerr << "Error: requesting primdiv where none exists!" << endl;
	exit(1);
      }
    long nq;
    for (Primevar pr(Quadprimes::list);
	 found==0 && pr.ok(); pr++)
      {
	q=(Quadprime)pr;
	if (q.div(a)) { found=1; ans = q; }
	else
	  if (nq=q.norm(), nq*nq>normlimit)
	    { found=1; ans=a.pos_assoc(); }  // a *known* to be prime
      }
    if (found==0) {
      ans=a.pos_assoc(); 
      cerr<<"No prime divisor found for "<<a<<" so assuming prime!\n";
    }
    return ans;
}

prime_factn::prime_factn(const PRIMETYPE& n)
{
  Primelist temp_plist(10);                      //Hope this is enough
  Tlist<long> temp_elist(10);                    //Hope this is enough
  PRIMETYPE a=n;
  long np = 0;
  long normp;
  for (Primevar pr(Quadprimes::list);
       (a.norm()>1) && pr.ok() && (np < 10); pr++)
    { 
      Quadprime p = (Quadprime)pr;
      if (p.div(a))
	{
	  temp_plist[np] = p;
	  temp_elist[np]=0;
	  do 
	    { 
	      a/=p; 
	      temp_elist[np]+=1; 
	    } 
	  while (p.div(a));
	  np++;
	}
      else
	if (normp=p.norm(), normp*normp>a.norm())
	  {
	    temp_plist[np] = a.pos_assoc();
	    a=1;
	    temp_elist[np++]=1;
	  }
    }
  if( np == 10 ) 
    cerr << "problems calculating prime divisors in function pdivs";
  //In case of p-factors outside range, assume the cofactor is prime:
  if (a.norm()>1) {
    temp_plist[np]=a.pos_assoc();
    temp_elist[np++]=1;
    cerr << "Warning: assuming PRIMETYPE "<< a <<"is prime; maybe the stored";
    cerr << " Quadprimes::list is too short to factor this PRIMETYPE." << endl;
  }

  temp_plist.truncate(np);
  temp_elist.truncate(np);
  plist=temp_plist;
  elist=temp_elist;
}

void prime_factn::display(ostream& s) const
{
  if (plist.getlength()==0)
    {
      s << "No prime factors." << endl;
    }
  else
    {
      s << "Number of distinct prime factors: "<< plist.getlength()<<endl;
      s << "Exponents   Primes:" << endl;
      for (long i=0; i<plist.getlength(); i++)
	s << elist(i) << "\t" << plist(i) << endl;
    }
}

Primelist pdivs(const PRIMETYPE& aa)  // rtns list of prime divisors
{
  Primelist plist(10);              //Hope this is enough
  PRIMETYPE a=aa;
  long np = 0;
  long normp;
  for (Primevar pr(Quadprimes::list);
       (a.norm()>1) && pr.ok() && (np < 10); pr++)
    { 
      Quadprime p = (Quadprime)pr;
      if (p.div(a))
	{
	  plist[np++] = p;
	  do {a/=p;} while (p.div(a));
	}
      else
	if (normp=p.norm(), normp*normp>a.norm())
	  {
	    plist[np++] = a.pos_assoc();   a=1;
	  }
    }
  if( np == 10 ) 
    cerr << "problems calculating prime divisors in function pdivs";
  //In case of p-factors outside range, assume the cofactor is prime:
  if (a.norm()>1) {
    plist[np++] = a.pos_assoc();  
    cerr << "Warning: assuming PRIMETYPE "<< a <<"is prime, maybe maxnorm";
    cerr << " is too small to factor this PRIMETYPE.";
  }
  plist.truncate(np);
  return plist;
}

Tlist<PRIMETYPE> posdivs(const PRIMETYPE& a)    // all "positive" divisors
{
  prime_factn pp(a);
  PRIMETYPE p; 
  long np = pp.num_primes();
  long nu = 1; long nd=nu;
  for (long i=0; i<np; i++) {nd*=(1+pp.expo(i));}
  Tlist<PRIMETYPE> dlist(nd);
  dlist[0]=1;
  nd=nu;
  long e;
  for(long i=0; i<np; i++)
    {
      p = pp.prime(i);
      e = pp.expo(i);
      for (long j=0; j<e; j++)
	for (long k=0; k<nd; k++)
	  dlist[nd*(j+1)+k] = (p*dlist[nd*j+k]).pos_assoc();
      nd*=(e+1);
    } 
  return dlist;
}


#if MAX_CLASSNUM<2
#include "primes1.cc"
#else
#include "primes2.cc"
#endif

// END OF FILE primes.cc
