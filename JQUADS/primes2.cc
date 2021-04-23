// FILE primes2.cc

// Included by  primes.cc  if class number up to two is required.
// (If h=1 this may be inefficient, since ideal methods will be used.)

void Quadprimes::init(long max)
{
  maxnorm=max;
// JEC rewrote his library in May'98, making the next line obsolete
//  initprimes(maxnorm);
  long d=Field::d, disc=Field::disc, t=Field::t;
  Tlist<Quadprime> list1(maxnorm),list2(maxnorm); //waste:2*pi(maxnorm) enough!
  int ramflag=0;
  long maxinertp=long(sqrt(0.5+maxnorm));  // truncate towards zero
                                           // (+0.5 to survive rounding errors)
  long count1=0, count2=0;
  long p; long b,dd,d1;

  for (primevar pr; pr.ok()&&pr<=maxnorm; pr++)
    { p=pr;
      long sig = kronecker(disc,p);
      switch (sig) {
      case  0:
	if (t)
	  { list1[count1++] = Quadprime(p,(p-1)/2,1, p, 0);}
	else 
	  { list1[count1++] = (ramflag++) ? Quadprime(p,0,1, p, 0) :
	                                    Quadprime(2,(Field::n%2),1, 2, 0);}
	break; 
      case -1: if (p<=maxinertp) list2[count2++]=Quadprime(1,0,p, p, 0);
	break;
      case  1:
	// solve b^2+bt+n \equiv 0 mod p
	b=0; dd=Field::n; d1=1+t; while (dd%p != 0) { b++; dd+=d1; d1+=2; }
        list1[count1++] = Quadprime(p,b,1, p, 1);
        list1[count1++] = Quadprime(p,p-b-t,1, p, 2);
      }
    }
  
  long npr = count1+count2;  //Total number of prime ideals found
  list=Tlist<Quadprime>(npr);
  long ipr=0, point1=0, point2=0;
  Quadprime alpha=list1[0], beta=list2[0];
  while ((point2<count2) && (point1<count1))
    { if (alpha.norm() < beta.norm())
	{list[ipr++]=alpha;  alpha=list1[++point1];}
    else 
      {list[ipr++]=beta;   beta=list2[++point2];}
    }
  
//Only one of the following (almost certainly the first) will run:

  while (point1<count1) list[ipr++]=list1[point1++];
  while (point2<count2) list[ipr++]=list2[point2++];
}

// need divisors functions etc

Tlist<Qideal> alldivs(const Qideal& a)      // all divisors
{
  return posdivs(a);
}

Tlist<Qideal> sqdivs(const Qideal& a)     
{                         // all divisors whose square divides a, up to +/-
  prime_factn pp(a);
  Qideal p;
  long np = pp.num_primes();

  long nd=1;
  for(long i=0; i<np; i++) { nd *= ( 1+ pp.expo(i)/2 ) ;}

  Tlist<Qideal> dlist(nd);
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

Tlist<Qideal> sqfreedivs(const Qideal& a)       // all square-free divisors
{
  Primelist plist=pdivs(a);
  Qideal p;
  long np = plist.getlength();
  long nd = 1;
  while (np-->0) nd*=2;
  Tlist<Qideal> dlist(nd);
  dlist[0]=1;
  nd=1;
  for(Primevar pr(plist); pr.ok(); pr++)
    {
      p = (Quadprime)pr;
      for (long k=0; k<nd; k++)
	dlist[nd+k] = p*dlist[k];
      nd*=2;
    } 
  return dlist;
}


// END OF FILE primes2.cc
