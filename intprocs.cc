// FILE INTPROCS.CC

#include <eclib/arith.h>
#include <assert.h>
#include <numeric>
#include "intprocs.h"

// The point of the following function is that the built-in gcc
// division truncates towards 0, while we need rounding, with a
// consistent behaviour for halves (they go up here).
//
// For b>0, roundover(a,b) = q such that a/b = q + r/b with -1/2 <= r/b < 1/2

long roundover_old(long aa, long bb)
{
  long a=aa, b=bb, q, r;
  assert (b>0); // the following code requires b>0
  r = (a<0? b-(-a)%b : a%b);
  if (2*r>=b) r-=b;
  // Now    -b   <= 2*r    < b
  assert ((-b<=2*r) && (2*r<b));
  q = (a-r)/b;
  return q;
}

long roundover(long a, long b)
{
  std::ldiv_t qr = ldiv(a, b);
  long r = qr.rem, q = qr.quot;
  long r2 = r<<1;
  return (r2<-b? q-1: (r2>=b? q+1: q));
}

// eclib only has the following for bigints
long sqrt_mod_p(long a, long p) // p odd prime, a quadratic residue
{
  bigint rr;
  sqrt_mod_p(rr, BIGINT(posmod(a,p)), BIGINT(p));
  return I2long(rr);
}

long vecgcd(const vector<long>& a)
{
  long g=0;
  for(vector<long>::const_iterator ai=a.begin(); ai!=a.end() && (g!=1); ai++)
    g = gcd(g, *ai);
  return g;
}

// returns content and sets c such that content(a) = a.c
long vecbezout(const vector<long>& a, vector<long>& c)
{
  long x = 1, g = vecgcd(a);
  vector<long> a0=a;
  if (g>1)
    for(vector<long>::iterator ai=a0.begin(); ai!=a0.end(); ai++)
      (*ai) /= g;
  // Now a0 is primitive: we do this to make numbers smaller in what follows
  int n=(int)a.size();
  c = vector<long>(n, 0);
  long g1=0;
  for(int i=0; i<n &&g1!=1; i++)
    {
      g1=bezout(g1,a0[i],x,c[i]);
      for(int j=0; j<i; j++) c[j]*=x;
    }
  return g;
}

// returns g = content(a) and c such that g = a.c, with c reduced mod s
long xmodvecbezout(long s, const vector<long>& a, vector<long>& c)
{
  int n=(int)a.size();
  long i, j, x=1, ci=1, g=0;
  for(i=0; i<n; i++)
    {
      g = bezout(g,a[i],x,ci);
      c[i] = xmod(ci,s);
      for(j=0; j<i; j++)
        c[j] = xmodmul(c[j],x,s);
    }
  return g;
}

long dot(const vector<long>& a, const vector<long>& c)
//returns g = a.c
{
  return std::inner_product(a.begin(), a.end(), c.begin(), 0);
}

// returns a.c computed as a double
double ddot(const vector<long>& a, const vector<long>& c)
{
  double g;
  vector<long>::const_iterator ai,ci;
  for(g=0, ai=a.begin(), ci=c.begin(); ai!=a.end(); ++ai, ++ci)
    g += (double(*ai) * (double)(*ci));
  return g;
}

// inner product mod s
long xmoddot(long s, const vector<long>& a, const vector<long>& c)
{
  long g;
  vector<long>::const_iterator ai,ci;
  for(g=0, ai=a.begin(), ci=c.begin(); ai!=a.end(); ++ai, ++ci)
    g = xmod( g + xmodmul(*ai, *ci,s), s);
  return g;
}

//#define testbezout

// If (b,d)!=(0,0), return an HNF basis [(aa,bb), (cc,0)] for the
// Z-module [[a,b],[c,d]]:
vector<long> hnf22(long a, long b, long c, long d)
{
#ifdef testbezout
  cout<<"  - hnf22("<<a<<", "<<b<<", "<<c<<", "<<d<<") = "<<flush;
#endif
  long x,y;
  long bb = bezout(b,d, x,y);
  long cc = abs(a*d-b*c)/bb;
  long aa = a*x+c*y;
  if (cc) aa%=cc;
  vector<long> v = {aa, bb, cc};
#ifdef testbezout
  cout<<v<<endl;
#endif
  return v;
}

//Sets basis={e1,e2,f1} such that [[e1,f1], [e2,0]] is a Z-basis for
//the Z-module spanned by [first[i], second[i]]

void findzbasis(const vector<long>& first, const vector<long>& second, vector<long>& basis)
{
#ifdef testbezout
  cout<<"findzbasis("<<first<<", "<<second<<"): "<<endl;
#endif
  long a, b, c, d;
  vector<long>::const_iterator ai, bi, aj, bj;
  // Find a nonsingular 2x2 block:
  int t=1;
  for (ai=first.begin(), bi=second.begin(); t && ai!=first.end(); ai++, bi++)
    {
      a = *ai;
      b = *bi;
      for (aj=ai+1, bj=bi+1; t && aj!=first.end(); aj++, bj++)
        {
          c = *aj;
          d = *bj;
          if (a*d-b*c)
            {
              basis = hnf22(a,b, c,d);
              t = 0;
            }
        }
    }
  a = basis[0];
  b = basis[1];
  c = basis[2];
  assert (!t);
  assert (basis[1]*basis[2]!=0);
#ifdef testbezout
  cout<<" - after step 0, {a,b,c} = "<<basis<<endl;
#endif
  // process all the rest
  for (ai=first.begin(), bi=second.begin(); ai!=first.end(); ai++, bi++)
    {
      basis = hnf22(a,b, *ai, *bi);
      c = gcd(c, basis[2]);
      a = basis[0]%c;
      b = basis[1];
      basis = {a, b, c};
#ifdef testbezout
      cout<<" - after one step using ("<<(*ai)<<","<<(*bi)<<"), {a,b,c} = "<<basis<<endl;
#endif
    }
}

// No longer used:

// sets basis={e1,e2,f1} such that [[e1,e2], [f1,0]] is a Z-basis
// for the Z-module spanned by [first[i], second[i]], and also sets x, y to be vectors such that
//
// [e1, e2] = [first.x, second.x]
// [f1,  0] = [first.y, second.y]

void findzbasiscoeffs(const vector<long>& first, const vector<long>& second,
                      vector<long>& basis, vector<long>& x, vector<long>& y)
{
  int i, n=(int)first.size();
  vector<long> u(n), newfirst(n);
  basis.resize(3);
  basis[1] = vecbezout(second,x);
  basis[0] = dot(first,x);  //dot product
//Now [basis[0], basis[1]] is the x-combination of the data, with basis[1]=gcd(second)
//newfirst = first-e1*(second/e2);
  for(i=0; i<n; i++)
    newfirst[i] = first[i] - basis[0]*(second[i]/basis[1]);
  basis[2] = vecbezout(newfirst,u);
  //  y = u - ((u*second)/basis[1])*x;
  long t = dot(u,second);
  for(i=0; i<n; i++)
    y[i]=u[i]-(t*x[i])/basis[1];
  // reduce e1 mod f1
  std::ldiv_t qr = ldiv(basis[0], basis[2]);
  long q = qr.quot;
  if (q!=0)
    {
#ifdef testbezout
      cout<<"findzbasis("<<first<<","<<second<<") --> basis="<<basis<<"i.e. [e1,e2]=["<<basis[0]<<","<<basis[1]<<"], [f1,0]=["<<basis[2]<<",0]"<<endl;
      cout<<" reducing e1="<<basis[0]<< " mod f1="<<basis[2]<<": new e1="<<qr.rem<<endl;
#endif
      basis[0] = qr.rem;
      for(i=0; i<n; i++)
        x[i] -= q*y[i];
    }
#ifdef testbezout
  cout<<"findzbasis("<<first<<","<<second<<") --> basis="<<basis<<"i.e. [e1,e2]=["<<basis[0]<<","<<basis[1]<<"], [f1,0]=["<<basis[2]<<",0]"<<endl;
  cout<<"coefficient vectors x="<<x<<", y="<<y<<endl;
//Check:
  if( ! (  (basis[0]==dot(first,x))   &&
           (basis[1]==dot(second,x))  &&
           (basis[2]==dot(first,y))   &&
           (       0==dot(second,y)) ))
  {
    cerr<<"Error in findzbasis!"  <<endl;
    exit(1);
  }
#endif
}

