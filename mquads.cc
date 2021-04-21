// FILE MQUADS.CC

#define testbezout    // define this to turn on self-verification of bezout

#include <sstream>

#include "eclib/marith.h"
#include "mquads.h"

// (pointer) functions declared extern in mquads.h

mQuad (*mquadconj)(const mQuad& a);
bigint (*mquadnorm)(const mQuad& a);
mQuad (*mult)(const mQuad& a, const mQuad& b);
mQuad (*qdivi)(const mQuad& a, const bigint& c);
int (*pos)(const mQuad& a);
mQuad (*mquadgcd)(const mQuad& aa, const mQuad& bb);
mQuad (*mquadbezout)(const mQuad& aa, const mQuad& bb, mQuad& xx, mQuad& yy);

//Declare static data members of class mQuad:

int mQuad::d;
int mQuad::disc;
bigint mQuad::t;
bigint mQuad::n;
char mQuad::name;
bigint mQuad::maxnorm;
int mQuad::nunits;

//Primes
mQuadlist mquadprimes;  //Initialised by initmquadprimes, see below
long nmquadprimes;     //The number of them.
mQuadlist mquadunits;
mQuad fundunit;

void mQuad::field(int dd, const bigint& max)
{d=dd; 
 if ((d+1)%4) {t=0; disc=4*d; n=d;       
               mquadconj=&mquadconj0; mquadnorm=&mquadnorm0; 
               mult=&mult0; qdivi=&qdivi0;
              }
 else         {t=1; disc=d;   n=(d+1)/4; 
               mquadconj=&mquadconj1; mquadnorm=&mquadnorm1; 
               mult=&mult1; qdivi=&qdivi1;
              }
 switch (d) {
 case 1:  pos=&pos13; name='i'; nunits=4; fundunit=mQuad(0,1); break;
 case 2:  pos=&pos2;  name='t'; nunits=2; fundunit=mQuad(-1);  break;
 case 3:  pos=&pos13; name='w'; nunits=6; fundunit=mQuad(0,1); break;
 default: pos=&pos2;  name='a'; nunits=2; fundunit=mQuad(-1); 
 }
 switch (d) {
 case 1: case 2: case 3: case 7: case 11:       // Euclidean
   mquadgcd=&mquadgcd1; mquadbezout=&mquadbezout1;  
   break;
 case 19: case 43: case 67: case 163:           // Non-Euclidean
   mquadgcd=&mquadgcd2; mquadbezout=&mquadbezout2;
   break;
 default: 
   cerr << "Field does not have class number one!" << endl;
   exit(1);    //abort program instantly!
 }
 mquadunits=mQuadlist(nunits); mquadunits[0]=1;
 for(int i=1; i<nunits; i++) mquadunits[i]=fundunit*mquadunits[i-1];
 maxnorm=max; initmquadprimes();
}

void mQuad::displayfield(ostream& s)
{s<<"Q(sqrt("<<-d<<"))\tdiscriminant = "<<disc;
 s<<"\tmin poly("<<name<<") = "<<name<<"^2"; if(t!=0) s<<"-"<<name; 
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
 s<<nmquadprimes<<" primes initialised, max norm = " << maxnorm << endl;
}

mQuad::mQuad(const complex& z)
{double x=real(z), y=imag(z);
 if(d>1) y/=sqrt((double)d);
 if(d>2) {x-=y; y*=2.0;}
 r=(bigint)(long)(x+0.5); i=(bigint)(long)(y+0.5);    //Rounded
}

mQuad::operator complex()
{double x=I2double(r), y=I2double(i);
 if(d>2) {y/=2.0; x+=y;}
 if(d>1) y*=sqrt((double)d);
 return complex(x,y);
}

int div(const mQuad& a, const mQuad& b)
{
 if (a==0) return (b==0);
 if (b==0) return 1;
 bigint nb = mquadnorm(a);
 mQuad c = b*mquadconj(a);
 return ((real(c)%nb)==0) && ((imag(c)%nb)==0);
}

int ndiv(const mQuad& a, const mQuad& b)
{
 return !div(a,b);
}

int val(const mQuad& factor, const mQuad& number)
{ if ((number==0) || (mquadnorm(factor)<=1)) return 9999;
  int e = 0; mQuad n=number, f=factor, nf;
  while (nf=n/f, f*nf==n) {e++; n=nf;}
  return e;
}

mQuadlist residues(const mQuad& a)
{
  bigint norma = mquadnorm(a), m = gcd(real(a),imag(a));
  bigint rednorma = (norma/m)/m;
  mQuadlist ans(I2long(norma)); int count=0;
  for(int j=0; j<m*rednorma; j++)
    for(int k=0; k<m; k++) 
      {mQuad res = mQuad(j,k)%a; //cout<<"res "<<count<<" = "<<res;
       ans[count++]= res ;
      }
  return ans;
}

ostream& operator<<(ostream& s, const mQuad& a)
{
bigint  i = a.i, r = a.r;
if (i==0) s<<r;
else {if (r==0) 
       {
        if (i==1) ;
        else if (i==-1) s << "-";
             else s << i;
       }
       else 
       {s<<r; 
        if(i>0) s<<"+"; else {s<<"-"; i.negate();}
        if(i>1) s<<i;
       }
      s<<(mQuad::name);
     }
return s;
}

char* to_string(const mQuad& a)
{
  char* ans=new char[25];    //allows 10 digits per part plus signs etc.
  ostringstream s(ans);      //lets us output to "ans" via s, as for cout etc.
  s << a; s.put(char(0));
  return ans;
}

//Functions for computing mQuad-primes, initializing the mquadlist
//mquadprimes.  NB all primes are "pos" i.e. normalized w.r.t. units

void factorp0(long p, long& a, long& b, int d)
{ int found=0;
  for (b=1; !found; b++)
  { long a2 = p - d*b*b;
    a = (long)round(sqrt((double)a2));
    found = (a*a == a2);
  }
  b--;  
}
 
void factorp1(long p, long& a, long& b, int d)
{ int found=0; long fourp = 4*p;
  for (b=1; !found; b++)
  { long a2 = fourp -d*b*b;
    a = (long)round(sqrt((double)a2));
    found = (a*a == a2);
  }
  b--;
  a=(a-b)/2;
}

void mQuad::initmquadprimes()
{
  int d=mQuad::d, disc=-mQuad::disc, t=I2long(mQuad::t);
  long imaxnorm = I2long(maxnorm);
  mQuadlist list1(imaxnorm),list2(imaxnorm);
  int i, count1=0, count2=0;
  long p ; long a,b;

  for (primevar pr; pr.ok()&&pr<imaxnorm; pr++)
    { p=pr;
      int sig = kronecker(disc,p);
      switch (sig) {
      case  0: 
               if((d==1)||(d==3)) list1[count1++]=mQuad(1,1); break;
               if (d==2)          list1[count1++]=mQuad(0,1); break;
               if (d>3)           list1[count1++]=mQuad(-1,2);
               break;
      case -1: if(p*p<=maxnorm) list2[count2++]=mQuad(p,0);
               break;
      case 1:
        if(t==0) factorp0(p,a,b,d); else factorp1(p,a,b,d);
        list1[count1++] = makepos(mQuad(a,b));
        list1[count1++] = makepos(mquadconj(mQuad(a,b)));
      }
    }
  
  long npr = count1+count2;  //Total number of mQuad-primes found
  mquadprimes=mQuadlist(npr);
  long ipr=0, point1=0, point2=0;
  mQuad alpha=list1[0], beta=list2[0];
  while ((point2<count2) && (point1<count1))
    {  if (mquadnorm(alpha)<mquadnorm(beta))
         {mquadprimes[ipr++]=alpha;  alpha=list1[++point1];}
    else 
      {mquadprimes[ipr++]=beta;   beta=list2[++point2];}
     }
  
//Only one of the following (almost certainly the first) will run:

  while (point1<count1) mquadprimes[ipr++]=list1[point1++];
  while (point2<count2) mquadprimes[ipr++]=list2[point2++];

  nmquadprimes = npr;
}



mQuad primdiv(const mQuad& a)
{
 mQuad p=0,q;
 if (mquadnorm(a)<2) return 0;   //Well, it has to return something!
 for (mQuadvar pr(mquadprimes); pr.ok() && p==0; ++pr)
   {q=pr; bigint nq;
    if (div(q,a)) p = q;
    else if (nq=mquadnorm(q),nq*nq>mquadnorm(a)) p=makepos(a);
   }
 if (p==0) {p=makepos(a); 
            cerr<<"No prime divisor found for "<<a<<" so assuming prime!\n";
           }
 return p;
}


mQuadlist pdivs(const mQuad& aa)
{ mQuad a=aa; bigint normp;
  int np = 0; mQuadlist plist(10);              //Hope this is enough
  for (mQuadvar pr(mquadprimes); (mquadnorm(a)>1) && pr.ok(); ++pr)
     { mQuad p = pr;
       if (div(p,a)) {plist[np++] = p;
                      while (div(p,a)) a/=p; 
                    }
       else if (normp=mquadnorm(p),normp*normp>mquadnorm(a)) 
         {plist[np++] = makepos(a);   a=1;
        }
     }
 //In case of p-factors outside range, assume the cofactor is prime:
 if (mquadnorm(a)>1) plist[np++]=makepos(a);  
 plist.truncate(np);
 return plist;
}


mQuadlist posdivs(const mQuad& a)       // all "positive" divisors (up to units)
{
 mQuadlist plist=pdivs(a); 
 int np = plist.length;
 int e, nu = 1; int nd=nu;
 int* elist = new int[np];
 mQuadvar pr(plist); 
 for(; pr.ok(); ++pr) 
   {elist[pr.index]=e=val(pr,a); nd*=(1+e);}
 mQuadlist dlist(nd);
 dlist[0]=1;
 mQuad p; nd=nu;
 for(pr.init(plist); pr.ok(); ++pr)
   {
     p = pr; int e = elist[pr.index];
     for (int j=0; j<e; j++)
       for (int k=0; k<nd; k++)
         dlist[nd*(j+1)+k] = makepos(p*dlist[nd*j+k]);
     nd*=(e+1);
   } 
 return dlist;
}

mQuadlist alldivs(const mQuad& a)       // all divisors
{
 mQuadlist plist=pdivs(a); 
 int np = plist.length;
 int e, nu = mQuad::nunits; int nd=nu;
 int* elist = new int[np];
 mQuadvar pr(plist);
 for(; pr.ok(); ++pr) 
   {elist[pr.index]=e=val(pr,a); nd*=(1+e);}
 mQuadlist dlist(nd);
 for(int i=0; i<nu; i++) dlist[i]=mquadunits[i];
 mQuad p; nd=nu;
 for(pr.init(plist); pr.ok(); ++pr)
   {
     p = pr; int e = elist[pr.index];
     for (int j=0; j<e; j++)
       for (int k=0; k<nd; k++)
         dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
     nd*=(e+1);
   } 
 return dlist;
}

mQuadlist sqdivs(const mQuad& a)     // all divisors whose square divides a, up to +/-
{
 mQuadlist plist=pdivs(a); 
 int np = plist.length;
 int e, nu = mQuad::nunits/2; int nd=nu;
 int* elist = new int[np];
 mQuadvar pr(plist); 
 for(; pr.ok(); ++pr) 
   {elist[pr.index]=e=val(pr,a)/2; nd*=(1+e);}
 mQuadlist dlist(nd);
 for(int i=0; i<nu; i++) dlist[i]=mquadunits[i];
 mQuad p; nd=nu;
 for(pr.init(plist); pr.ok(); ++pr)
   {
     p = pr; int e = elist[pr.index];
     for (int j=0; j<e; j++)
       for (int k=0; k<nd; k++)
         dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
     nd*=(e+1);
   } 
 return dlist;
}

mQuadlist sqfreedivs(const mQuad& a)       // all square-free divisors
{
 mQuadlist plist=pdivs(a); 
 int np = plist.length;
 int e, nu = 2; int nd=nu;
 int* elist = new int[np];
 mQuadvar pr(plist);
 for(; pr.ok(); ++pr) 
   {elist[pr.index]=e=1; nd*=(1+e);}
 mQuadlist dlist(nd);
 for(int i=0; i<nu; i++) dlist[i]=mquadunits[i];
 mQuad p; nd=nu;
 for(pr.init(plist); pr.ok(); ++pr)
   {
     p = pr; int e = elist[pr.index];
     for (int j=0; j<e; j++)
       for (int k=0; k<nd; k++)
         dlist[nd*(j+1)+k] = p*dlist[nd*j+k];
     nd*=(e+1);
   } 
 return dlist;
}

//#define TRACE

mQuad mquadgcd1(const mQuad& aa, const mQuad& bb) 
  //Only works for Euclidean fields!
{mQuad a=aa,b=bb,c,q;
 while (b!=0)
   {
#ifdef TRACE
cout<<"a = "<<a<<", b = "<<b<<endl;
#endif
     q = a/b; c = a-q*b; a = b; b = c;
#ifdef TRACE
cout<<"q = "<<q<<", a = "<<a<<", b = "<<b<<endl;
#endif
     if(mquadnorm(b)>=mquadnorm(a)){cout<<"error--norm not reduced!\n";break;}
   }
 while (!pos(a)) a*=fundunit;
 return a;
}

mQuad mquadbezout1(const mQuad& alpha, const mQuad& beta, mQuad& coeff1, mQuad& coeff2)
{mQuad a=alpha,b=beta,c,x=0,oldx=1,newx,y=1,oldy=0,newy,q,g;
 while (b!=0)
 { q = a/b; 
#ifdef TRACE
cout<<"a = "<<a<<", b = "<<b<<endl;
#endif
   c    = a    - q*b; a    = b; b = c;
#ifdef TRACE
cout<<"q = "<<q<<", a = "<<a<<", b = "<<b<<endl;
#endif
   newx = oldx - q*x; oldx = x; x = newx;
   newy = oldy - q*y; oldy = y; y = newy;
   if(mquadnorm(b)>=mquadnorm(a)){cout<<"error--norm not reduced!\n";break;}
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
    {cerr<<"Error in mquadbezout1!"<<endl;
     cerr<<"alpha = "<<alpha<<endl;
     cerr<<"beta  = "<<beta<<endl;
     cerr<<"coeff1 = "<<coeff1<<endl;
     cerr<<"coeff2 = "<<coeff2<<endl;
     cerr<<"g   = "<<g<<endl;
     }
#endif
 return g;
}
 
mQuad invmod(const mQuad& a, const mQuad& p)
{mQuad x,y;
 mQuad g=mquadbezout(a,p,x,y);
 if (g==1) return x;
 else {cerr<<"invmod called with "<<a<<" and "<<p<<" -- not coprime!"<<endl;
       return 0;
      }
}

int coprime(const mQuad& a, const mQuad& b) 
{
  mQuad g=mquadgcd(a,b); 
  return g==1;
}

int invertible(const mQuad& a, const mQuad& b, mQuad& inverse)
{ mQuad y; mQuad g = mquadbezout(a,b,inverse,y);
  return g==1;
}

//functions needed for non-euclidean fields to compute bezout/mquadgcd

bigint vecbezout(int n, bigint* a, bigint* c) 
//returns g = content(a) = a.c
{
  bigint x=0,ai=0,ci=0,g=0;  //This does not initialise them properly
                              //for the call to bezout.  Don't know why.
  x++; ci++;                  //This does the trick!
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
 
bigint vecgcd(int n, bigint* a)
//returns g = content(a)
{
  bigint g=0;
  for(int i=0; (i<n)&(g!=1); i++) g=gcd(g,a[i]);
  return g;
}
 
bigint dot(int n, bigint* a, bigint* c) 
//returns g = a.c
{
  bigint g=0;
  for(int i=0; i<n; i++) g+=a[i]*c[i];
  return g;
}
 

void findzbasiscoeffs(int n, bigint* first, bigint* second, 
                      bigint* basis, bigint* x, bigint* y)
//Given: a (2xn) matrix a with rows "first", "second"
//Returns: "x","y": cols of coeffs (nx2), and
//         "basis" [e1,e2,f1] such that the cols of
//         (e1 f1)  are a Z-basis for the cols of a.
//         (e2  0)  
{
  bigint* u=new bigint[n];
  bigint* newfirst=new bigint[n];                    //temps
  bigint e2=vecbezout(n,second,x);
  bigint e1=dot(n,first,x);  //dot product
//Now [e1,e2] is the x-combination of the data, with e2=gcd(second)
//newfirst = first-e1*(second/e2);
  int i;
  for(i=0; i<n; i++) newfirst[i]=first[i]-e1*second[i]/e2;
  bigint f1 = vecbezout(n,newfirst,u);
  basis[0] = e1;   basis[1] = e2;   basis[2] = f1;
//  y = u - ((u*second)/e2)*x;
  bigint t = dot(n,u,second);
  for(i=0; i<n; i++) y[i]=u[i]-(t*x[i])/e2;
#ifdef testbezout
//Check:
  if( ! (  (e1==dot(n,first,x))   &&  
           (e2==dot(n,first,y))   &&
           (f1==dot(n,second,x))  &&
           (0==dot(n,second,y)) ))
  {cerr<<"Error in findzbasis!"  <<endl; }
#endif                            
 delete[] u; delete[] newfirst;
}
 
void findzbasis(int n, bigint* first, bigint* second, bigint* basis)
//Same as findzbasiscoeffs except don't need x,y
{
  bigint* x=new bigint[n];   
  bigint* newfirst=new bigint[n];     //temps
//cout<<"In findzbasis with first="<<first<<", second="<<second<<", basis="<<basis<<endl;
//cout<<"About to call vecbezout with second and x="<<x<<endl;
  bigint e2=vecbezout(n,second,x);
  bigint e1=dot(n,first,x);  //dot product
//Now [e1,e2] is the x-combination of the data, with e2=gcd(second)
//newfirst = first-e1*(second/e2);
  for(int i=0; i<n; i++) newfirst[i]=first[i]-e1*second[i]/e2;
  bigint f1 = vecgcd(n,newfirst);
  basis[0] = e1;   basis[1] = e2;   basis[2] = f1;
  delete[] x; delete[] newfirst;
}
 
//(redundant old version superceded by findminmquad)
/*
bigint* findminvec(const bigint* basis)
{
  bigint normw=mQuad::n;
  bigint e1,e2,f1,f2,n,norm1,norm2,temp;
  bigint* c(2),d(2),v;
  c[1] = 1; c[2] = 0; d[1] = 0; d[2] = 1;
  e1 = basis[1]; f1 = basis[3]; 
  e2 = basis[2]; f2 = 0;
  while (
         n = roundover(((2*e1 + e2)*f1 + (e1 + 2*normw*e2)*f2), 
                        (2*(f1*(f1+f2) + normw*f2*f2))),
         e1 -=  n*f1, e2 -=  n*f2,
         d  -=  n*c,
         norm1 = e1*(e1 + e2) + normw*e2*e2,
         norm2 = f1*(f1 + f2) + normw*f2*f2,
         (norm2 > norm1)
         )
    {
      temp = e1; e1 = -f1; f1 = temp; 
      temp = e2; e2 = -f2; f2 = temp;
      v = d; d = -c; c = v;
    }
  return c;
}
*/
 
bigint* findminmquad(mQuad alpha, mQuad beta, mQuad& gen)
{
  bigint n,normalpha,normbeta=mquadnorm(beta); mQuad temp;
  bigint* c = new bigint[2];
  bigint* d = new bigint[2];  
  bigint v;
  c[0] = 1; c[1] = 0; d[0] = 0; d[1] = 1;
  while (
         n = roundover(real(alpha*mquadconj(beta)),normbeta),
         //= round(real(alpha/beta))
         alpha -= n*beta,
//       d     -= n*c,
         d[0]     -= n*c[0],
         d[1]     -= n*c[1],
         normalpha = mquadnorm(alpha),
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
 
mQuad mquadbezout2(const mQuad& alpha, const mQuad& beta, mQuad& coeff1, mQuad& coeff2)
{
  mQuad g;  
  if (div(beta, alpha)) { g=beta; coeff1=0; coeff2=1;}
  else if (div(alpha, beta)) { g=alpha; coeff1=1; coeff2=0;}
  else
    {
      bigint n = mQuad::n; bigint t = mQuad::t;
      bigint* rv = new bigint[4];
      bigint* iv = new bigint[4];
      bigint* basis = new bigint[3];
      bigint* x = new bigint[4];
      bigint* y = new bigint[4];   
      bigint* z = new bigint[4];
      rv[0] = real(alpha); iv[0] = imag(alpha);
      rv[1] = real(beta);  iv[1] = imag(beta);
      rv[2] = -n*iv[0];    iv[2] = rv[0] + t*iv[0];
      rv[3] = -n*iv[1];    iv[3] = rv[1] + t*iv[1];
      findzbasiscoeffs(4,rv,iv,basis,x,y);
      mQuad al=mQuad(basis[0],basis[1]), be=basis[2];
      bigint* coeff = findminmquad(al,be,g);
      for(int i=0; i<4; i++) z[i] = coeff[0]*y[i] + coeff[1]*x[i];
      coeff1 = mQuad(z[0],z[2]);    coeff2 = mQuad(z[1],z[3]);
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
    {cerr<<"Error in mquadbezout2!"<<endl;
     cerr<<"alpha = "<<alpha<<endl;
     cerr<<"beta  = "<<beta<<endl;
     cerr<<"coeff1 = "<<coeff1<<endl;
     cerr<<"coeff2 = "<<coeff2<<endl;
     cerr<<"g   = "<<g<<endl;
     }
#endif
  return g;
}

mQuad mquadgcd2(const mQuad& alpha, const mQuad& beta)
//Same as mquadbezout2 except don't need coeff1, coeff2
{
  if (div(beta, alpha)) return beta;
  if (div(alpha, beta)) return alpha;
  bigint n = mQuad::n, t=mQuad::t;
  bigint* rv = new bigint[4];
  bigint* iv = new bigint[4];
  bigint* basis = new bigint[3];
  rv[0] = real(alpha); iv[0] = imag(alpha);   
  rv[1] = real(beta);  iv[1] = imag(beta);
  rv[2] = -n*iv[0];    iv[2] = rv[0] + t*iv[0];
  rv[3] = -n*iv[1];    iv[3] = rv[1] + t*iv[1];
//cout<<"About to call findzbasis with rv="<<rv<<", iv="<<iv<<", basis="<<basis<<endl;
  findzbasis(4,rv,iv,basis);
  mQuad g, al=mQuad(basis[0],basis[1]), be=basis[2];
  bigint* v = findminmquad(al,be,g);
  while (!pos(g)) g*=fundunit;
#ifdef testbezout
//CHECK:
  if (div(g,alpha) && div(g,beta)) {;}  //OK
  else 
    {cerr<<"Error in mquadgcd2!"<<endl;
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
#if(0)
bigint roundover(const bigint& aa, const bigint& bb)
{ //cout<<"roundover("<<aa<<","<<bb<<") = ";
  bigint ans;
  if(aa%bb==0) ans = aa/bb;   //exact case always OK
  else 
    {
      bigint a=aa,b=bb;
      if(b<0){a=-a; b=-b;}
      bigint c = a+a+b, b2=b+b; 
      ans = c/b2;
      if((c<0)&&((c%b2)!=0))ans--;
    }
//  cout<<ans<<endl;
  return ans;
}
#endif
