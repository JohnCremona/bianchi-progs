// FILE PRIMES1.CC

// Included by PRIMES.CC if only class number one is required.

//Functions for computing Quad-primes, initializing the Quadlist
//Quadprimes.  NB all primes are "pos" i.e. normalized w.r.t. units

void factorp0(long p, long& a, long& b, long d);
void factorp1(long p, long& a, long& b, long d);

void Quadprimes::init(long max)
{
  maxnorm=max;
  initprimes(maxnorm);
  long d=Field::d, disc=Field::disc, t=Field::t;
  Quadlist list1(maxnorm),list2(maxnorm);
  long count1=0, count2=0;
  long p ; long a,b;
  
  for (primevar pr; pr.ok()&&pr<=maxnorm; pr++)
    {
      p=pr;
      long sig = kronecker(disc,p);
      switch (sig) {
      case  0: 
	list1[count1++]= d==1 ? Quad(1,1) :
                         d==2 ? Quad(0,1) : 
                         d==3 ? Quad(1,1) :
                                Quad(-1,2); 
	break;
      case -1: 
	if(p*p<=maxnorm) list2[count2++]=Quad(p,0);
	break;
      case +1:
	if(t==0) factorp0(p,a,b,d); else factorp1(p,a,b,d);
	list1[count1++] = Quad(a,b).pos_assoc();
	list1[count1++] = Quad(a,b).conj().pos_assoc();
      }
    }
  
  long npr = count1+count2;  //Total number of Quad-primes found
  list=Quadlist(npr);
  long ipr=0, point1=0, point2=0;
  Quad alpha=list1[0], beta=list2[0];
  while ((point2<count2) && (point1<count1))
    {
      if (alpha.norm() < beta.norm())
	{ list[ipr++]=alpha; alpha=list1[++point1]; }
      else
	{ list[ipr++]=beta;   beta=list2[++point2]; }
    }
  
  //Only one of the following (almost certainly the first) will run:
  
  while (point1<count1) list[ipr++]=list1[point1++];
  while (point2<count2) list[ipr++]=list2[point2++];
}

long round(double x)
{
 return long(x+0.5);  // Should be OK for x>0
}

void factorp0(long p, long& a, long& b, long d)
{ int found=0;
  for (b=1; !found; b++)
  { long a2 = p - d*b*b;
    a = round(sqrt(a2));
    found = (a*a == a2);
  }
  b--;  //Not quite sure why, but works
}
 
void factorp1(long p, long& a, long& b, long d)
{ int found=0; long fourp = 4*p;
  for (b=1; !found; b++)
  { long a2 = fourp -d*b*b;
    a = round(sqrt(a2));
    found = (a*a == a2);
  }
  b--;
  a=(a-b)/2;
}

//Functions for computing Quad-primes, initializing the Quadlist
//Quadprimes.  NB all primes are "pos" i.e. normalized w.r.t. units

Quadlist alldivs(const Quad& a)       // all divisors
{
    Quadlist plist=pdivs(a); Quad p; 
    long np = plist.length;
    long e, nu = Field::nunits; long nd=nu;
    long* elist = new long[np];
    for(Quadvar pr(plist); pr.ok(); pr++) 
	{p=(Quad)pr; elist[pr.index]=e=val(p,a); nd*=(1+e);}
    Quadlist dlist(nd);
    for(long i=0; i<nu; i++) dlist[i]=Quadunits::u(i);
    nd=nu;
    for(pr.init(plist); pr.ok(); pr++) {
	p = (Quad)pr; long e = elist[pr.index];
	for (long j=0; j<e; j++)
	    for (long k=0; k<nd; k++)
		dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
	nd*=(e+1);
    } 
    delete[] elist;
    return dlist;
}

Quadlist sqdivs(const Quad& a)     
{                         // all divisors whose square divides a, up to +/-
    Quadlist plist=pdivs(a); Quad p;
    long np = plist.length;
    long e, nu = Field::nunits/2; long nd=nu;
    long* elist = new long[np];
    for(Quadvar pr(plist); pr.ok(); pr++) 
	{p=(Quad)pr; elist[pr.index]=e=val(p,a)/2; nd*=(1+e);}
    Quadlist dlist(nd);
    for(long i=0; i<nu; i++) dlist[i]=Quadunits::u(i);
    nd=nu;
    for(pr.init(plist); pr.ok(); pr++) {
	p = (Quad)pr; long e = elist[pr.index];
	for (long j=0; j<e; j++)
	    for (long k=0; k<nd; k++)
		dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
	nd*=(e+1);
    } 
    delete[] elist;
    return dlist;
}

Quadlist sqfreedivs(const Quad& a)       // all square-free divisors
{
    Quadlist plist=pdivs(a); Quad p;
    long np = plist.length;
    long e, nu = 2; long nd=nu;
    long* elist = new long[np];
    for(Quadvar pr(plist); pr.ok(); pr++) 
	{elist[pr.index]=e=1; nd*=(1+e);}
    Quadlist dlist(nd);
    for(long i=0; i<nu; i++) dlist[i]=Quadunits::u(i);
    nd=nu;
    for(pr.init(plist); pr.ok(); pr++) {
	p = (Quad)pr; long e = elist[pr.index];
	for (long j=0; j<e; j++)
	    for (long k=0; k<nd; k++)
		dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
	nd*=(e+1);
    } 
    delete[] elist;
    return dlist;
}

// obsolete versions follow

/*
Quad primdiv(const Quad& a)
{
    Quad p=0,q;
    if (a.norm()<2) return 0;   // must return something!
    for (Quadvar pr(Quadprimes); pr.ok() && p==0; pr++) {
	q=(Quad)pr; long nq;
	if (q.div(a)) p = q;
	else if (nq=q.norm(), nq*nq > a.norm() ) p=a.pos_assoc();
    }
    if (p==0) {
	p=a.pos_assoc(); 
	cerr<<"No prime divisor found for "<<a<<" so assuming prime!\n";
    }
    return p;
}

Quadlist pdivs(const Quad& aa)
{
    Quad a=aa; long normp;
    long np = 0; Quadlist plist(10);              //Hope this is enough
    for (Quadvar pr(Quadprimes); (a.norm()>1) && pr.ok() && (np < 10); pr++)
	{ 
	    Quad p = (Quad)pr;
	    if (p.div(a)) {
		plist[np++] = p;
		while (p.div(a)) a/=p; 
	    }
	    else
	      if (normp=p.norm(),normp*normp>a.norm())
		{
		  plist[np++] = a.pos_assoc();   a=1;
		}
	}
    if( np == 10 ) 
	cerr << "problems calculating prime divisors in function pdivs";
    //In case of p-factors outside range, assume the cofactor is prime:
    if (a.norm()>1) {
	plist[np++] = a.pos_assoc();  
	cerr << "Warning: assuming Quad "<< a <<"is prime, maybe maxnorm";
	cerr << " is too small to factor this quad.";
    }
    plist.truncate(np);
    return plist;
}

Quadlist posdivs(const Quad& a)       // all "positive" divisors (up to units)
{
    Quadlist plist=pdivs(a); Quad p; 
    long np = plist.length;
    long e, nu = 1; long nd=nu;
    long* elist = new long[np];
    for(Quadvar pr(plist); pr.ok(); pr++) 
	{p=(Quad)pr; elist[pr.index]=e=val(p,a); nd*=(1+e);}
    Quadlist dlist(nd);
    dlist[0]=1;
    nd=nu;
    for(pr.init(plist); pr.ok(); pr++) {
	p = (Quad)pr; 
	long e = elist[pr.index];
	for (long j=0; j<e; j++)
	    for (long k=0; k<nd; k++)
	      dlist[nd*(j+1)+k] = (p*dlist[nd*j+k]).pos_assoc();
	nd*=(e+1);
    } 
    delete[] elist;
    return dlist;
}
*/



// END OF FILE PRIMES1.CC
