// FILE QUADS.CC

#define testbezout    // define this to turn on self-verification of bezout

#include <iostream>

//#include "vector.h"  // needed for zbasis functions for non-Euclidean fields.
#include <eclib/arith.h>
#include "quads.h"

//Declare static data members of class Quad:

int Quad::d;
int Quad::disc;
int Quad::t;
int Quad::n;
char Quad::name;
int Quad::maxnorm;
int Quad::nunits;

//Primes
vector<Quad> quadprimes;  //Initialised by initquadprimes, see below
long nquadprimes;         //The number of them.
vector<Quad> quadunits, squareunits;
Quad fundunit;


// declaration of "extern" functions declared in quads.h:
long (*quadnorm)(const Quad& a);
Quad (*mult)(const Quad& a, const Quad& b);
Quad (*qdivi)(const Quad& a, long c);
int (*pos)(const Quad& a);
Quad (*quadgcd)(const Quad& aa, const Quad& bb);
Quad (*quadbezout)(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy);
Quad (*quadconj)(const Quad& a);

void Quad::field(int dd, int max)
{d=dd; 
 if ((d+1)%4) {t=0; disc=4*d; n=d;       
               quadconj=&quadconj0; quadnorm=&quadnorm0; 
               mult=&mult0; qdivi=&qdivi0;
              }
 else         {t=1; disc=d;   n=(d+1)/4; 
               quadconj=&quadconj1; quadnorm=&quadnorm1; 
               mult=&mult1; qdivi=&qdivi1;
              }
 switch (d) {
 case 1:  pos=&pos13; name='i'; nunits=4; fundunit=Quad(0,1); break;
 case 2:  pos=&pos2;  name='t'; nunits=2; fundunit=Quad(-1);  break;
 case 3:  pos=&pos13; name='w'; nunits=6; fundunit=Quad(0,1); break;
 default: pos=&pos2;  name='a'; nunits=2; fundunit=Quad(-1); 
 }
 switch (d) {
 case 1: case 2: case 3: case 7: case 11:       // Euclidean
   quadgcd=&quadgcd1; quadbezout=&quadbezout1;  
   break;
 case 19: case 43: case 67: case 163:           // Non-Euclidean
   quadgcd=&quadgcd2; quadbezout=&quadbezout2;
   break;
 default: 
   cerr << "Field does not have class number one!" << endl;
   exit(1);    //abort program instantly!
 }
 int i;
 quadunits.push_back(1);
 quadunits.push_back(fundunit);
 for(i=2; i<nunits; i++)
   quadunits.push_back(fundunit*quadunits[i-1]);
 for(i=0; 2*i<nunits; i++)
   squareunits.push_back(quadunits[2*i]);
 maxnorm=max; initquadprimes();
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

Quad::Quad(const bigcomplex& z)
{bigfloat x=real(z), y=imag(z);
 if(d>1) y/=sqrt(to_bigfloat(d));
 if(d>2) {x-=y; y*=2.0;}
 longify(x, r); longify(y, i);    //Rounded
}

Quad::operator bigcomplex() const
{bigfloat x=to_bigfloat(r), y=to_bigfloat(i);
 if(d>2) {y/=2.0; x+=y;}
 if(d>1) y*=sqrt(to_bigfloat(d));
 return bigcomplex(x,y);
}

int div(const Quad& a, const Quad& b)
{
 if (a==0) return (b==0);
 if (b==0) return 1;
 long nb = quadnorm(a);
 Quad c = b*quadconj(a);
 return ((real(c)%nb)==0) && ((imag(c)%nb)==0);
}

int ndiv(const Quad& a, const Quad& b)
{
 return !div(a,b);
}

int val(const Quad& factor, const Quad& number)
{ 
  if ((number==0) || (quadnorm(factor)<=1)) 
    {
      cout << "Error in val(): factor = "<<factor<< " should have norm>1"<<endl;
      exit(1);
    }
  int e = 0; Quad n=number, f=factor, nf;
  while (nf=n/f, f*nf==n) {e++; n=nf;}
  return e;
}

vector<Quad> residues(const Quad& a)
{
  long norma = quadnorm(a), m = gcd(real(a),imag(a));
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
    pi =  (d==1 ? Quad(1,1) :
           d==2 ? Quad(0,1) :
           d==3 ? Quad(1,1) :
           Quad(-1,2));
    list.push_back(pi);
    break;
  case -1: // inert
    pi = Quad(p,0);
    list.push_back(pi);
    break;
  case 1: // split
    if(t==0) factorp0(p,a,b,d); else factorp1(p,a,b,d);
    pi = makepos(Quad(a,b));
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
  vector<Quad>::const_iterator pi, alpha, beta;
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
  plist.push_back(makepos(a));  
  return plist;
}

vector<Quad> posdivs(const Quad& a)   // all "positive" divisors (up to units)
{
  vector<Quad> plist=pdivs(a); Quad p; 
  int e, nu = 1; int nd=nu;
  vector<int> elist;
  vector<Quad>::const_iterator pr;
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
  vector<Quad>::const_iterator pr;
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
  vector<Quad>::const_iterator pr;
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
  vector<Quad>::const_iterator pr;
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
 while (b!=0)
   {
     q = a/b; c = a-q*b; a = b; b = c;
     if(quadnorm(b)>=quadnorm(a)){cout<<"error--norm not reduced!\n";break;}
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
  Quad g=quadgcd(a,b); 
  return g==1;
}

int invertible(const Quad& a, const Quad& b, Quad& inverse)
{ Quad y; Quad g = quadbezout(a,b,inverse,y);
  return g==1;
}

//functions needed for non-euclidean fields to compute bezout/quadgcd

long vecbezout(int n, long* a, long* c) 
//returns g = content(a) = a.c
{
  long x=0,ci=0,g=0;       //This does not initialise them properly
                           //for the call to bezout.  Don't know why.
  x++; ci++;               //This does the trick!
//cout<<"In vecbezout with a="<<a<<", c="<<c<<endl;
  for(int i=0; i<n; i++)
    { 
//cout<<"...calling bezout with g="<<g<<", a[i]="<<a[i]<<", x="<<x<<", ci="<<ci<<endl;
      g=bezout(g,a[i],x,ci);
//cout<<"...returns x="<<x<<", ci="<<ci<<endl;
      c[i]=ci;
      for(int j=0; j<i; j++) c[j]*=x;
    }
  return g;
}
 
long vecgcd(int n, long* a)
//returns g = content(a)
{
  long g=0;
  for(int i=0; (i<n)&(g!=1); i++) g=gcd(g,a[i]);
  return g;
}
 
long dot(int n, long* a, long* c) 
//returns g = a.c
{
  long g=0;
  for(int i=0; i<n; i++) g+=a[i]*c[i];
  return g;
}
 
void findzbasiscoeffs(int n, long* first, long* second, 
                      long* basis, long* x, long* y)
//Given: a (2xn) matrix a with rows "first", "second"
//Returns: "x","y": cols of coeffs (nx2), and
//         "basis" [e1,e2,f1] such that the cols of
//         (e1 f1)  are a Z-basis for the cols of a.
//         (e2  0)  
{
  long* u=new long[n];
  long* newfirst=new long[n];                    //temps
  long e2=vecbezout(n,second,x);
  long e1=dot(n,first,x);  //dot product
//Now [e1,e2] is the x-combination of the data, with e2=gcd(second)
//newfirst = first-e1*(second/e2);
  int i;
  for(i=0; i<n; i++) newfirst[i]=first[i]-e1*second[i]/e2;
  long f1 = vecbezout(n,newfirst,u);
  basis[0] = e1;   basis[1] = e2;   basis[2] = f1;
//  y = u - ((u*second)/e2)*x;
  long t = dot(n,u,second);
  for(i=0; i<n; i++) y[i]=u[i]-(t*x[i])/e2;
#ifdef testbezout
//Check:
  if( ! (  (e1==dot(n,first,x))   &&  
           (f1==dot(n,first,y))   &&
           (e2==dot(n,second,x))  &&
           (0==dot(n,second,y)) ))
  {cerr<<"Error in findzbasis!"  <<endl; }
#endif                            
 delete[] u; delete[] newfirst;
}
 
void findzbasis(int n, long* first, long* second, long* basis)
//Same as findzbasiscoeffs except don't need x,y
{
  long* x=new long[n];   
  long* newfirst=new long[n];     //temps
//cout<<"In findzbasis with first="<<first<<", second="<<second<<", basis="<<basis<<endl;
//cout<<"About to call vecbezout with second and x="<<x<<endl;
  long e2=vecbezout(n,second,x);
  long e1=dot(n,first,x);  //dot product
//Now [e1,e2] is the x-combination of the data, with e2=gcd(second)
//newfirst = first-e1*(second/e2);
  for(int i=0; i<n; i++) newfirst[i]=first[i]-e1*second[i]/e2;
  long f1 = vecgcd(n,newfirst);
  basis[0] = e1;   basis[1] = e2;   basis[2] = f1;
  delete[] x; delete[] newfirst;
}
 
long* findminquad(Quad alpha, Quad beta, Quad& gen)
{
  long n,normalpha,normbeta=quadnorm(beta); Quad temp;
  long* c = new long[2];
  long* d = new long[2];  
  long v;
  c[0] = 1; c[1] = 0; d[0] = 0; d[1] = 1;
  while (
	 longify(real(bigcomplex(alpha)/bigcomplex(beta)), n),
         alpha -= n*beta,
//       d     -= n*c,
         d[0]     -= n*c[0],
         d[1]     -= n*c[1],
         normalpha = quadnorm(alpha),
         (normbeta > normalpha)
         )
    {
      temp = alpha; alpha = -beta; beta = temp;
      normbeta = normalpha;
      v    = d[0];     d[0]     = -c[0];    c[0]    = v;
      v    = d[1];     d[1]     = -c[1];    c[1]    = v;
    }
  gen = beta; 
  delete[] d;
  return c;
}
 
Quad quadbezout2(const Quad& alpha, const Quad& beta, Quad& coeff1, Quad& coeff2)
{
  Quad g;  
  if (div(beta, alpha)) { g=beta; coeff1=0; coeff2=1;}
  else if (div(alpha, beta)) { g=alpha; coeff1=1; coeff2=0;}
  else
    {
      long n = Quad::n; long t = Quad::t;
      long* rv = new long[4];
      long* iv = new long[4];
      long* basis = new long[3];
      long* x = new long[4];
      long* y = new long[4];   
      long* z = new long[4];
      rv[0] = real(alpha); iv[0] = imag(alpha);
      rv[1] = real(beta);  iv[1] = imag(beta);
      rv[2] = -n*iv[0];    iv[2] = rv[0] + t*iv[0];
      rv[3] = -n*iv[1];    iv[3] = rv[1] + t*iv[1];
      findzbasiscoeffs(4,rv,iv,basis,x,y);
      Quad al=Quad(basis[0],basis[1]), be=basis[2];
      long* coeff = findminquad(al,be,g);
      for(int i=0; i<4; i++) z[i] = coeff[0]*y[i] + coeff[1]*x[i];
      coeff1 = Quad(z[0],z[2]);    coeff2 = Quad(z[1],z[3]);
//Next two lines try to get coeff1,2 as small as possible
      if(beta!=0)
      {
       coeff1 = coeff1%beta;            //reduced
       coeff2 = (g-coeff1*alpha)/beta;  //should be exact
      }
      delete[] rv; delete[] iv; delete[] basis; 
      delete[] x; delete[] y; delete[] z; delete[] coeff;
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
  long* rv = new long[4];
  long* iv = new long[4];
  long* basis = new long[3];
  rv[0] = real(alpha); iv[0] = imag(alpha);   
  rv[1] = real(beta);  iv[1] = imag(beta);
  rv[2] = -n*iv[0];    iv[2] = rv[0] + t*iv[0];
  rv[3] = -n*iv[1];    iv[3] = rv[1] + t*iv[1];
//cout<<"About to call findzbasis with rv="<<rv<<", iv="<<iv<<", basis="<<basis<<endl;
  findzbasis(4,rv,iv,basis);
  Quad g, al=Quad(basis[0],basis[1]), be=basis[2];
  long* v = findminquad(al,be,g);
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
  delete[] rv; delete[] iv; delete[] basis; delete[] v;
  return g;
}
 

//-------------------------------------------------------------------------
// N.B.  The point of the following ghastly function is that the
// built-in gcc division truncates towards 0, while we need
// rounding, with a consistent behaviour for halves (they go up here). 
//-------------------------------------------------------------------------
long roundover(long aa, long bb)
{ //cout<<"roundover("<<aa<<","<<bb<<") = ";
  long ans;
  if(aa%bb==0) ans = aa/bb;   //exact case always OK
  else 
    {
      long a=aa,b=bb;
      if(b<0){a=-a; b=-b;}
      long c = a+a+b, b2=b+b; 
      ans = c/b2;
      if((c<0)&&(c%b2))ans--;
    }
  //cout<<ans<<endl;
  return ans;
}

// HNF of ideal (alpha) as a triple [a c d] where [a,c+d*w] is a Z-basis with
//
// a,d>0; c>=0
// N=a*d = Norm(alpha)
// d|a and d|b
// 0 <=c < a

vector<long> HNF(const Quad& alpha)
{
  long N = quadnorm(alpha);
  long xa = real(alpha), ya = imag(alpha), u, v;
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
  s << "[" << quadnorm(alpha) << "," << H[1] << "," << H[2] << "]";
  return s.str();
}

string ideal_label(const Quad& alpha)  // returns label of ideal (alpha)
{
  vector<long>H = HNF(alpha);
  stringstream s;
  s << quadnorm(alpha) << "." << H[1] << "." << H[2];
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
  if(quadnorm(a)!=quadnorm(b)) return 0;
  vector<Quad>::const_iterator eps;
  for(eps=quadunits.begin(); eps!=quadunits.end(); eps++)
    if(a*(*eps)==b) return 1;
  return 0;
}

int is_ideal_Galois_stable(const Quad& a)
{
  return are_associate(a, quadconj(a));
}
