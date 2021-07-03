// FILE QUADS.CC

#define testbezout    // define this to turn on self-verification of bezout

#include <iostream>

#include <eclib/arith.h>
#include "intprocs.h"
#include "quads.h"
#include "geometry.h"

//Declare static data members of class Quad:

int Quad::d;
int Quad::disc;
int Quad::t;
int Quad::n;
char Quad::name;
int Quad::maxnorm;
int Quad::nunits;
int Quad::is_Euclidean;
int Quad::is_class_number_one;
Quad Quad::w;
Quad Quad::zero;
Quad Quad::one;

//Primes
vector<Quad> quadprimes;  //Initialised by initquadprimes, see below
long nquadprimes;         //The number of them.
vector<Quad> quadunits, squareunits;
Quad fundunit;

vector<int> euclidean_fields = {1,2,3,7,11};
vector<int> class_number_one_fields = {1,2,3,7,11,19,43,67,163};
vector<int> valid_fields = {1,2,3,7,11,19,43,67,163, 5};

int check_field(int d, vector<int> fields)
{
  return (std::find(fields.begin(), fields.end(), d) != fields.end());
}

// declaration of "extern" functions declared in quads.h:
Quad (*mult)(const Quad& a, const Quad& b);
Quad (*qdivi)(const Quad& a, long c);
int (*pos)(const Quad& a);
Quad (*quadgcd)(const Quad& aa, const Quad& bb);
Quad (*quadbezout)(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy);
Quad (*quadconj)(const Quad& a);

Quad qdivi0(const Quad& a, long c) // c>0,    // used when t=0
{
  Quad ans;
  if (c>0)
    {
      ans.r = roundover(a.r,c);
      ans.i = roundover(a.i,c);
    }
  else
    {
      ans.r = roundover(-a.r,-c);
      ans.i = roundover(-a.i,-c);
    }
  ans.setnorm();
  return ans;
}

Quad qdivi1(const Quad& a, long c) // used when t=1
{
  Quad ans;
  if (c>0)
    {
      ans.i = roundover(a.i,c);
      ans.r = roundover(2*a.r+a.i-c*ans.i,2*c);
    }
  else
    {
      ans.i = roundover(-a.i,-c);
      ans.r = roundover(-2*a.r-a.i+c*ans.i,-2*c);
    }
  ans.setnorm();
  return ans;
}

void Quad::field(int dd, int max)
{
  // if (!check_field(dd))
  //   {
  //     cerr<<"field "<<dd<<" is not implemented: it must be one of: "<<valid_fields<<endl;
  //     exit(1);
  //   }
  d = dd;
  is_Euclidean = check_field(d, euclidean_fields);
  is_class_number_one = check_field(d, class_number_one_fields);

  if ((d+1)%4) {t=0; disc=4*d; n=d;
               quadconj=&quadconj0;
               mult=&mult0; qdivi=&qdivi0;
              }
 else         {t=1; disc=d;   n=(d+1)/4;
               quadconj=&quadconj1;
               mult=&mult1; qdivi=&qdivi1;
              }
  w = Quad(0,1, n);
  zero = Quad(0,0, 0);
  one = Quad(1,0, 1);

  switch (d) {
  case 1:  pos=&pos13; name='i'; nunits=4; fundunit=w; break;
  case 2:  pos=&pos2;  name='t'; nunits=2; fundunit=-one; break;
  case 3:  pos=&pos13; name='w'; nunits=6; fundunit=w; break;
  default: pos=&pos2;  name='a'; nunits=2; fundunit=-one;
  }

  quadgcd=&quadgcd_psea;
  quadbezout=&quadbezout_psea;

  int i;
  quadunits.push_back(1);
  quadunits.push_back(fundunit);
  for(i=2; i<nunits; i++)
    quadunits.push_back(fundunit*quadunits[i-1]);
  for(i=0; 2*i<nunits; i++)
    squareunits.push_back(quadunits[2*i]);
  maxnorm=max;
  if(is_class_number_one)
    {
      initquadprimes();
      setup_geometry();
    }
}

void Quad::displayfield(ostream& s)
{s<<"Q(sqrt("<<-d<<"))\tdiscriminant = "<<disc;
 s<<"\tmin poly("<<name<<") = "<<name<<"^2"; if(t) s<<"-"<<name; 
 s<<"+"<<n<<".\n";
 switch (d) {
 case 1: case 2: case 3: case 7: case 11:       // Euclidean
   cout << "Euclidean" << endl;
   break;
 case 19: case 43: case 67: case 163:           // Non-Euclidean
   cout << "Non-Euclidean" << endl;
   break;
 default:
   cout << "Class number > 1" << endl;
 }
 s<<nquadprimes<<" primes initialised, max norm = " << maxnorm << endl;
}

int Quad::chi(long p)
{
  return (p==2? (d%4==3? (d%8==3? -1: +1): 0):  legendre(-d,p));
}

Quad::Quad(const bigcomplex& z)
{bigfloat x=real(z), y=imag(z);
 if(d>1) y/=sqrt(to_bigfloat(d));
 if(d>2) {x-=y; y*=2.0;}
 longify(x, r); longify(y, i);    //Rounded
 setnorm();
}

Quad::operator bigcomplex() const
{bigfloat x=to_bigfloat(r), y=to_bigfloat(i);
 if(d>2) {y/=2.0; x+=y;}
 if(d>1) y*=sqrt(to_bigfloat(d));
 return bigcomplex(x,y);
}

int div(const Quad& a, const Quad& b)
{
 if (a.nm==0) return (b.nm==0);
 if (b.nm==0) return 1;
 if (b.nm%a.nm!=0) return 0;
 Quad c = b*quadconj(a);
 return (((c.r)%a.nm)==0) && (((c.i)%a.nm)==0);
}

int ndiv(const Quad& a, const Quad& b)
{
 return !div(a,b);
}

int val(const Quad& factor, const Quad& number)
{
  if ((number.nm==0) || (factor.nm<=1))
    {
      cout << "Error in val(): factor = "<<factor<< " should not be a unit"<<endl;
      exit(1);
    }
  int e = 0; Quad n=number, f=factor, nf;
  while (nf=n/f, f*nf==n) {e++; n=nf;}
  return e;
}

vector<Quad> residues(const Quad& a)
{
  long norma = a.norm(), m = gcd(a.re(), a.im());
  long rednorma = (norma/m)/m;
  vector<Quad> ans;
  for(int j=0; j<m*rednorma; j++)
    for(int k=0; k<m; k++) 
      ans.push_back(Quad(j,k)%a) ;
  return ans;
}

ostream& operator<<(ostream& s, const Quad& a)
{
long  i = a.i, r = a.r;
if (i==0) s<<r;
else {if (r==0) 
       {
        if (i==1) ;
        else if (i==-1) s << "-";
             else s << i;
       }
       else 
       {s<<r; 
        if(i>0) s<<"+"; else s<<"-";
        if (abs(i)>1) s<<abs(i);
       }
      s<<Quad::name;
     }
return s;
}


//Functions for computing quad-primes, initializing the vector<Quad>
//quadprimes.  NB all primes are "pos" i.e. normalized w.r.t. units

void factorp0(long p, long& a, long& b, int d)
// finds a,b s.t. a^2+d*b^2=0 (mod p)
{ int found=0;
  for (b=1; !found; b++)
  { long a2 = p - d*b*b;
    a = (long)round(sqrt(double(a2)));
    found = (a*a == a2);
  }
  b--;
}
 
void factorp1(long p, long& a, long& b, int d)
// finds a,b s.t. a^2+a*b+((d+1)/4)*b^2=0 (mod p)
{ int found=0; long fourp = 4*p;
  for (b=1; !found; b++)
  { long a2 = fourp -d*b*b;
    a = (long)round(sqrt(double(a2)));
    found = (a*a == a2);
  }
  b--;
  a=(a-b)/2;
}

vector<Quad> Quad::primes_above(long p, int& sig)
{
  int d=Quad::d, disc=-Quad::disc, t=Quad::t;
  long a,b;  Quad pi, piconj;
  vector<Quad> list;
  sig = kronecker(disc,p);
  switch (sig) {
  case  0: // ramified
    pi =  (d==1 ? Quad(1,1, 2) :
           d==2 ? Quad(0,1, 2) :
           d==3 ? Quad(1,1, 3) :
           Quad(-1,2, d));
    list.push_back(pi);
    break;
  case -1: // inert
    pi = Quad(p,0, p*p);
    list.push_back(pi);
    break;
  case 1: // split
    if(t==0) factorp0(p,a,b,d); else factorp1(p,a,b,d);
    pi = makepos(Quad(a,b, p));
    piconj = makepos(quadconj(pi));
    // We choose the ordering so the HNFs are [p,c,1], [p,c',1] with c<c'
    int c = posmod(a*invmod(b,p),p);
    if (2*c<p)
      {
        list.push_back(pi);
        list.push_back(piconj);
      }
    else
      {
        list.push_back(piconj);
        list.push_back(pi);
      }
  }
  //cout << "primes_above("<<p<<") = " << list << endl;
  return list;
}

void Quad::initquadprimes()
{
  long p; int sig;
  vector<Quad> list, list1, list2;
  vector<Quad>::iterator pi, alpha, beta;
  for (primevar pr; pr.ok()&&pr<maxnorm; pr++)
    { p=pr;
      list = Quad::primes_above(p, sig);
      for (pi = list.begin(); pi!=list.end(); )
        switch (sig) {
        case  0:
          list1.push_back(*pi++);
          break;
        case -1:
          if(p*p<=maxnorm)
            list2.push_back(*pi);
          pi++;
          break;
        case 1:
          {
            list1.push_back(*pi++);
            list1.push_back(*pi++);
          }
        }
    }

  // Now list1 contains the degree 1 primes and list2 the degree 2
  // primes, each ordered by norm; we merge these lists so that they
  // are still sorted by norm:

  alpha=list1.begin(); beta=list2.begin();
  while ((alpha!=list1.end()) && (beta!=list2.end()))
    if (quadnorm(*alpha)<quadnorm(*beta))
      quadprimes.push_back(*alpha++);
    else
      quadprimes.push_back(*beta++);

  while (alpha!=list1.end()) quadprimes.push_back(*alpha++);
  while ( beta!=list2.end()) quadprimes.push_back(*beta++);

  nquadprimes = quadprimes.size();
}

Quad primdiv(const Quad& a)
{
  long na=quadnorm(a);
  if (na<2) return 0;   // must return something!
  vector<Quad>::const_iterator pr;
  for (pr=quadprimes.begin(); pr!=quadprimes.end(); pr++)
    {
      Quad p=*pr;
      if (div(p,a)) return p;
      long np=quadnorm(p);
      if (np*np>na) return makepos(a);
    }
  cout<<"No prime divisor found for "<<a<<" so assuming prime!\n";
  return makepos(a);
}

vector<Quad> pdivs(const Quad& aa)
{ Quad a=aa; long norma=quadnorm(a);
  vector<Quad> plist; // will hold prime factors
  if (norma<2) return plist;
  vector<Quad>::const_iterator pr;
  for (pr=quadprimes.begin(); pr!=quadprimes.end(); pr++)
     { Quad p = *pr;
       if (div(p,a)) 
	 {
	   plist.push_back(p);
	   while (div(p,a)) a/=p; 
	   norma=quadnorm(a);	 
	   if (norma==1) return plist;
	 }
       else 
	 {
	   long normp=quadnorm(p);
	   if (normp*normp>norma) 
	     {
	       plist.push_back(makepos(a));
	       return plist;
	     }
	 }
     }
  //In case of p-factors outside range, assume the cofactor is prime:
  if (quadnorm(a)>1)
    plist.push_back(makepos(a));  
  return plist;
}

vector<Quad> posdivs(const Quad& a)   // all "positive" divisors (up to units)
{
  vector<Quad> plist=pdivs(a); Quad p; 
  int e, nu = 1; int nd=nu;
  vector<int> elist;
  vector<Quad>::iterator pr;
  for (pr=plist.begin(); pr!=plist.end(); pr++)
    {
      e=val(*pr,a); 
      elist.push_back(e);
      nd*=(1+e);
    }
  vector<Quad> dlist(nd);
  dlist[0]=1;
  nd=nu;
  vector<int>::iterator ei;
  for (pr=plist.begin(), ei=elist.begin(); pr!=plist.end(); pr++, ei++)
   {
     p = *pr; 
     e = *ei;
     for (int j=0; j<e; j++)
       for (int k=0; k<nd; k++)
         dlist[nd*(j+1)+k] = makepos(p*dlist[nd*j+k]);
     nd*=(e+1);
   } 
 return dlist;
}

vector<Quad> alldivs(const Quad& a)       // all divisors
{
  vector<Quad> plist=pdivs(a); Quad p; 
  int e, nu = Quad::nunits; int nd=nu;
  vector<int> elist;
  vector<Quad>::iterator pr;
  for (pr=plist.begin(); pr!=plist.end(); pr++)
   {
     e=val(*pr,a); 
     elist.push_back(e);
     nd*=(1+e);
   }
  vector<Quad> dlist(nd);
  for(int i=0; i<nu; i++) dlist[i]=quadunits[i];
  nd=nu;
  vector<int>::iterator ei;
  for (pr=plist.begin(), ei=elist.begin(); pr!=plist.end(); pr++, ei++)
   {
     p = *pr; 
     e = *ei;
     for (int j=0; j<e; j++)
       for (int k=0; k<nd; k++)
         dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
     nd*=(e+1);
   } 
 return dlist;
}

vector<Quad> sqdivs(const Quad& a) // all divisors whose square divides a, up to +/-
{
  vector<Quad> plist=pdivs(a); Quad p;
  int e, nu = Quad::nunits/2; int nd=nu;
  vector<int> elist;
  vector<Quad>::iterator pr;
  for (pr=plist.begin(); pr!=plist.end(); pr++)
   {
     e=val(*pr,a)/2; 
     elist.push_back(e);
     nd*=(1+e);
   }
  vector<Quad> dlist(nd);
  for(int i=0; i<nu; i++) dlist[i]=quadunits[i];
  nd=nu;
  vector<int>::iterator ei;
  for (pr=plist.begin(), ei=elist.begin(); pr!=plist.end(); pr++, ei++)
   {
     p = *pr; 
     e = *ei;
     for (int j=0; j<e; j++)
       for (int k=0; k<nd; k++)
         dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
     nd*=(e+1);
   } 
 return dlist;
}

vector<Quad> sqfreedivs(const Quad& a)       // all square-free divisors
{
  vector<Quad> plist=pdivs(a); Quad p;
  int e, nu = 2; int nd=nu;
  vector<int> elist;
  vector<Quad>::iterator pr;
  for (pr=plist.begin(); pr!=plist.end(); pr++)
   {
     e=1; 
     elist.push_back(e);
     nd*=(1+e);
   }
  vector<Quad> dlist(nd);
  for(int i=0; i<nu; i++) dlist[i]=quadunits[i];
  nd=nu;
  vector<int>::iterator ei;
  for (pr=plist.begin(), ei=elist.begin(); pr!=plist.end(); pr++, ei++)
   {
     p = *pr; 
     e = *ei;
     for (int j=0; j<e; j++)
       for (int k=0; k<nd; k++)
         dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
     nd*=(e+1);
   } 
 return dlist;
}

Quad quadgcd1(const Quad& aa, const Quad& bb)  //Only works for Euclidean fields!
{Quad a=aa,b=bb,c,q;
 while (b.nm != 0)
   {
     //     cout<<"a="<<a<<" (norm "<<a.nm<<"), b="<<b<<" (norm "<<b.nm<<")\n";
     q = a/b; c = a-q*b; a = b; b = c;
     //     cout<<"quotient "<<q<<", remainder "<<c<<" (norm "<<c.nm<<")"<<endl;
     //     cout<<"a="<<a<<" (norm "<<a.nm<<"), b="<<b<<" (norm "<<b.nm<<")\n";
     if(b.nm >= a.nm){cout<<"error--norm not reduced!\n";exit(1);}
   }
 while (!pos(a)) a*=fundunit;
 return a;
}

Quad quadbezout1(const Quad& alpha, const Quad& beta, Quad& coeff1, Quad& coeff2)
{Quad a=alpha,b=beta,c,x=0,oldx=1,newx,y=1,oldy=0,newy,q,g;
 while (b!=0)
 { q = a/b; 
   c    = a    - q*b; a    = b; b = c;
   newx = oldx - q*x; oldx = x; x = newx;
   newy = oldy - q*y; oldy = y; y = newy;
  }
 coeff1=oldx; coeff2=oldy; g=a;
//Next two lines try to get coeff1,coeff2 as small as possible
 if(beta!=0)
 {
  coeff1 = coeff1%beta;            //reduced
  coeff2 = (g-coeff1*alpha)/beta;     //should be exact
 }
 while (!pos(g))
   { g*=fundunit; coeff1*=fundunit; coeff2*=fundunit;
   }
#ifdef testbezout
//CHECK:
  if (div(g,alpha) && div(g,beta) && (g==coeff1*alpha+coeff2*beta)) {;}  //OK
  else 
    {cerr<<"Error in quadbezout1!"<<endl;
     cerr<<"alpha = "<<alpha<<endl;
     cerr<<"beta  = "<<beta<<endl;
     cerr<<"coeff1 = "<<coeff1<<endl;
     cerr<<"coeff2 = "<<coeff2<<endl;
     cerr<<"g   = "<<g<<endl;
     }
#endif
 return g;
}

Quad invmod(const Quad& a, const Quad& p)
{Quad x,y;
 Quad g=quadbezout(a,p,x,y);
 if (g==1) return x;
 else {cerr<<"invmod called with "<<a<<" and "<<p<<" -- not coprime!"<<endl;
       return 0;
      }
}

int coprime(const Quad& a, const Quad& b) 
{
  Quad g = quadgcd(a,b); 
  return g==1;
}

int invertible(const Quad& a, const Quad& b, Quad& inverse)
{ Quad y; Quad g = quadbezout(a,b,inverse,y);
  return g==1;
}

long nearest_long_to_Quad_quotient ( const Quad& alpha, const Quad& beta)
// returns round(true_real_part_of(alpha/beta))
{
  return roundover((alpha*quadconj(beta)).re(), beta.norm());
}

// there follow 4 "flavours" of findminQuad returning different amounts of data

vector<long> findminquadcoeffs(const Quad&al, const Quad&be, Quad& beta, Quad& alpha)
{ alpha=al;
  beta=be;
  long n, v;
  Quad temp;
  vector<long> c = {1, 0};
  vector<long> d = {0, 1};

  while (
	 n = nearest_long_to_Quad_quotient(alpha,beta),
         //= round(true_real_part_of(alpha/beta))
         alpha -= n*beta,
//       d     -= n*c,
         d[0]     -= n*c[0],
         d[1]     -= n*c[1],
         quadnorm(beta) > quadnorm(alpha)
         )
    {
      temp = alpha; alpha = -beta; beta = temp;
      v    = d[0];     d[0]     = -c[0];    c[0]    = v;
      v    = d[1];     d[1]     = -c[1];    c[1]    = v;
    }
// beta is now the non-zero Quad of least norm in the lattice [al,be]
// alpha is another Quad such that [al,be]=[alpha,beta]
#ifdef testbezout
  if (beta != c[1]*al + c[0]*be)
    {
      cerr << "Error in findminquadscoeffs!" << endl;
      cerr << "[al,be] = ["<<al<<","<<be<<"]"<<endl;
      cerr << "c1,c0 = "<<c[1]<<","<<c[0]<<endl;
      cerr << "beta = "<<beta<< " not equal to "<<c[1]*al + c[0]*be<<endl;
    }
#endif
  return c;
}

vector<long> findminquadcoeffs(const Quad& alpha, const Quad& beta, Quad& gen0)
{ Quad gen1;
  return findminquadcoeffs(alpha,beta,gen0,gen1);
}

void findminquad(const Quad&al, const Quad&be, Quad& beta, Quad& alpha)
// same as findminquadscoeffs but don't need coeffs
{ alpha=al;
  beta=be;
  long n;
  Quad temp;
  while (
	 n = nearest_long_to_Quad_quotient(alpha,beta),
         //= round(true_real_part_of(alpha/beta))
         alpha -= n*beta,
	 quadnorm(beta) > quadnorm(alpha)
	 )
    {
      temp = alpha; alpha = -beta; beta = temp;
    }
// beta is now the non-zero Quad of least norm in the lattice [al,be]
// alpha is another Quad such that [al,be]=[alpha,beta]
}
 
void findminquad(const Quad& alpha, const Quad& beta, Quad& gen0)
// same as findminquadcoeffs but don't need coeffs
{ Quad gen1;
  findminquad(alpha,beta,gen0,gen1);
}

Quad quadbezout2(const Quad& alpha, const Quad& beta, Quad& coeff1, Quad& coeff2)
{
  Quad g;
  if (div(beta, alpha)) { g=beta; coeff1=0; coeff2=1;}
  else if (div(alpha, beta)) { g=alpha; coeff1=1; coeff2=0;}
  else
    {
      long n = Quad::n; long t = Quad::t;
      vector<long> rv(4), iv(4), basis(3), x(4), y(4), z(4);
      rv[0] = alpha.re(); iv[0] = alpha.im();
      rv[1] = beta.re();  iv[1] = beta.im();
      rv[2] = -n*iv[0];    iv[2] = rv[0] + t*iv[0];
      rv[3] = -n*iv[1];    iv[3] = rv[1] + t*iv[1];
      findzbasiscoeffs(rv,iv,basis,x,y);
      Quad al=Quad(basis[0],basis[1]), be=basis[2];
      vector<long> coeff = findminquadcoeffs(al,be,g);
      for(int i=0; i<4; i++) z[i] = coeff[0]*y[i] + coeff[1]*x[i];
      coeff1 = Quad(z[0],z[2]);    coeff2 = Quad(z[1],z[3]);
//Next two lines try to get coeff1,2 as small as possible
      if(beta!=0)
      {
       coeff1 = coeff1%beta;            //reduced
       coeff2 = (g-coeff1*alpha)/beta;  //should be exact
      }
    }
  while (!pos(g)) { g*=fundunit; coeff1*=fundunit; coeff2*=fundunit;} 
#ifdef testbezout
//CHECK:
  if (div(g,alpha) && div(g,beta) && (g==coeff1*alpha+coeff2*beta)) {;}  //OK
  else 
    {cerr<<"Error in quadbezout2!"<<endl;
     cerr<<"alpha = "<<alpha<<endl;
     cerr<<"beta  = "<<beta<<endl;
     cerr<<"coeff1 = "<<coeff1<<endl;
     cerr<<"coeff2 = "<<coeff2<<endl;
     cerr<<"g   = "<<g<<endl;
     }
#endif
  return g;
}

Quad quadgcd2(const Quad& alpha, const Quad& beta)
//Same as quadbezout2 except don't need coeff1, coeff2
{
  if (div(beta, alpha)) return beta;
  if (div(alpha, beta)) return alpha;
  long n = Quad::n, t=Quad::t;
  vector<long> rv(4), iv(4), basis(3);
  rv[0] = alpha.r; iv[0] = alpha.i;
  rv[1] = beta.r;  iv[1] = beta.i;
  rv[2] = -n*iv[0];    iv[2] = rv[0] + t*iv[0];
  rv[3] = -n*iv[1];    iv[3] = rv[1] + t*iv[1];
//cout<<"About to call findzbasis with rv="<<rv<<", iv="<<iv<<", basis="<<basis<<endl;
  findzbasis(rv,iv,basis);
  Quad g, al=Quad(basis[0],basis[1]), be=basis[2];
  vector<long> v = findminquadcoeffs(al,be,g);
  while (!pos(g)) g*=fundunit;
#ifdef testbezout
//CHECK:
  if (div(g,alpha) && div(g,beta)) {;}  //OK
  else
    {cerr<<"Error in quadgcd2!"<<endl;
     cerr<<"alpha = "<<alpha<<endl;
     cerr<<"beta  = "<<beta<<endl;
     cerr<<"g   = "<<g<<endl;
   }
#endif
  return g;
}
 

// HNF of ideal (alpha) as a triple [a c d] where [a,c+d*w] is a Z-basis with
//
// a,d>0; c>=0
// N=a*d = Norm(alpha)
// d|a and d|b
// 0 <=c < a

vector<long> HNF(const Quad& alpha)
{
  long N = alpha.norm(), xa = alpha.re(), ya = alpha.im(), u, v;
  long g = bezout(xa,ya,u,v);  // g=u*xa+v*ya=gcd(xa,ya)
  long x = xa/g, y = ya/g;
  // Now the HNF is g*[a, b+w] for some b mod a=N/g
  long t = Quad::t, n = Quad::n;
  long a = N/(g*g);
  long b = posmod(((v-t*u)*x - n*u*y), a);
  vector<long> ans;
  ans.push_back(a*g);
  ans.push_back(b*g);
  ans.push_back(g);
  return ans;
}

// Ideal label: formed from the Norm and HNF of the ideal (alpha)
// (subject to change!)

string old_ideal_label(const Quad& alpha)  // returns label of ideal (alpha)
{
  vector<long>H = HNF(alpha);
  stringstream s;
  s << "[" << alpha.norm() << "," << H[1] << "," << H[2] << "]";
  return s.str();
}

string ideal_label(const Quad& alpha)  // returns label of ideal (alpha)
{
  vector<long>H = HNF(alpha);
  stringstream s;
  s << alpha.norm() << "." << H[1] << "." << H[2];
  return s.str();
}

string field_label() // returns field label, e.g. '2.0.4.1'
{
  stringstream s;
  s << "2.0." << Quad::disc << ".1";
  return s.str();
}

int are_associate(const Quad& a, const Quad& b)
{
  if(a==0) return (b==0);
  if(a.norm() != b.norm()) return 0;
  vector<Quad>::const_iterator eps;
  for(eps=quadunits.begin(); eps!=quadunits.end(); eps++)
    if(a*(*eps)==b) return 1;
  return 0;
}

int is_ideal_Galois_stable(const Quad& a)
{
  Quad b(quadconj(a));
  return are_associate(a, b);
}

matop::matop(const Quad& p, const Quad& n)
{
 if (p==n)
   {
     mats.resize(1, mat22(0,-1,n,0));
   }
 else
 if (div(p,n))   // W involution, 1 term
   {
      Quad u,v,a,b;
      for (u=1, v=n; div(p,v); v/=p, u*=p) ;
      quadbezout(u,v,a,b);
      mats.resize(1, mat22(u*a,-b,n,u));
   }
else                 // Hecke operator, p+1 terms
  {
    vector<Quad> resmodp = residues(p);
    vector<Quad>::const_iterator r=resmodp.begin();
    while(r!=resmodp.end())
      mats.push_back(mat22(1,*r++,0,p));
    mats.push_back(mat22(p,0,0,1));
  }
}

