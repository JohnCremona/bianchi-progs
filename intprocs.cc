// FILE INTPROCS.CC

#include <assert.h>
#include <numeric>
#include "intprocs.h"

#ifdef INT_IS_ZZ

INT rounded_division(INT a, INT b, int round_down)
{
  INT q, r;
  ::divides(a,b, q,r);
  assert (a==b*q+r);
  INT r2 = 2*r;
  if (round_down)
    {
      // We want -b <= r2 < +b
      if (r2<-b)
        q-=1;
      else
        if (r2>=b)
          q+=1;
    }
  else
    {
      // We want -b < r2 <= +b
      if (r2<=-b)
        q-=1;
      else
        if (r2>b)
          q+=1;
    }
  return q;
}

int divides(const INT& aa, const INT& bb)
{
  return div(aa,bb);
}

int divrem(const bigint& a, const bigint& b, bigint& q, bigint& r)
{
  return divides(a, b, q, r);
}

INT isqrt(const INT& a)
{
  INT r; isqrt(a,r); return r;
}

// This is the heart of vecbzout3; if we are not using multiprecision
// integers we do things differently using bigfloats
void get_lambda_mu(const INT& x, const INT& y, const INT& z, const INT& w,
                   const INT& a1, const INT& b1, const INT& c1, const INT& h1,
                   INT& lambda, INT& mu)
{
  INT x2y2 = x*x+y*y;
  lambda = rounded_division(x2y2*c1*z-h1*w, x2y2*c1*c1+h1*h1);
  mu = rounded_division((b1*x-a1*y)*z, a1*a1+b1*b1);
}

#endif

#ifdef FLINT
void get_lambda_mu(INT x, INT y, INT z, INT w,
                   INT a1, INT b1, INT c1, INT h1,
                   INT& lambda, INT& mu)
{
  INT c1x2y2 = c1 * fmma(x,x,y,y); // c1 * (x*x+y*y)
  INT p = fmms(c1x2y2, z, h1, w);
  INT q = fmma(c1x2y2, c1, h1, h1);
  lambda = rounded_division(p, q);
  mu = rounded_division(fmms(b1,x,a1,y)*z, fmma(a1,a1,b1,b1));
}
#endif

#ifdef INT_IS_long

// This is the heart of vecbzout3; if we are not using multiprecision
// integers we do things differently using bigfloats
void get_lambda_mu(const INT& x, const INT& y, const INT& z, const INT& w,
                   const INT& a1, const INT& b1, const INT& c1, const INT& h1,
                   INT& lambda, INT& mu)
{
  bigfloat rx2y2 = pow(to_bigfloat(x),2) + pow(to_bigfloat(y),2);
  bigfloat rlambda = (rx2y2*c1*z-h1*w) / (rx2y2*pow(to_bigfloat(c1),2)+pow(to_bigfloat(h1),2));
  bigfloat rmu = to_bigfloat(b1*x-a1*y)*z /  (pow(to_bigfloat(a1),2)+pow(to_bigfloat(b1),2));
  longify(rlambda, lambda);
  longify(rmu, mu);
}
#endif

INT vecgcd(const vector<INT>& a)
{
  INT g(0);
  for(vector<INT>::const_iterator ai=a.begin(); ai!=a.end() && (g!=1); ++ai)
    g = gcd(g, *ai);
  return g;
}

// returns content and sets c such that content(a) = a.c
INT vecbezout2(vector<INT>& a, vector<INT>& c)
{
  int n=(int)a.size();
  if (n!=2) return vecbezout(a, c);
  INT x, y;
  INT g = bezout(a[0], a[1], x, y);
  c = {-y, x};
  assert (a[0]*c[0]+a[1]*c[1]==g);
  return g;
}

//#define testbezout3

// returns content and sets c such that content(a) = a.c
INT vecbezout3(vector<INT>& a, vector<INT>& c)
{
  int n=(int)a.size();
  if (n!=3) return vecbezout(a, c);

  INT aa=a[0], bb=a[1], cc=a[2];
#ifdef testbezout3
  cout<<"Computing vecbezout3("<<a<<")"<<endl;
#endif
  if ((is_zero(aa))&&(is_zero(bb)))
    {
      if (cc<0)
        {
          c = {ZERO,ZERO,MONE};
          return -cc;
        }
      else
        {
          c = {ZERO,ZERO,ONE};
          return cc;
        }
    }
  INT x, y, z, w;
  INT h = bezout(aa,bb, x, y);
  INT a1 = aa/h, b1 = bb/h;
  INT g = bezout(h,cc, z, w);
  INT h1=h/g, c1=cc/g;
  c = {x*z, y*z, w};
  assert (a[0]*c[0]+a[1]*c[1]+a[2]*c[2]==g);
#ifdef testbezout3
  cout<<"   g="<<g<<", first solution is "<<c<<"; "<<flush;
  vector<INT> perp1 = {-c1*x, -c1*y, h1};
  vector<INT> perp2 = {-b1, a1, 0};
  cout<<"   primitive basis of perp: "<<perp1<<", "<<perp2<<endl;
  assert (a[0]*perp1[0]+a[1]*perp1[1]+a[2]*perp1[2]==0);
  assert (a[0]*perp2[0]+a[1]*perp2[1]+a[2]*perp2[2]==0);
#endif
  // now minimize
  INT lambda, mu;
  get_lambda_mu(x,y,z,w,a1,b1,c1,h1, lambda, mu);
#ifdef testbezout3
  cout << " (lambda,mu)=("<<lambda<<","<<mu<<")"<<endl;
#endif

  c[0] -= lambda*c1*x+mu*b1;
  c[1] -= lambda*c1*y-mu*a1;
  c[2] += lambda*h1;
  assert (a[0]*c[0]+a[1]*c[1]+a[2]*c[2]==g);
#ifdef testbezout3
  cout<<" vecbezout3("<<a<<") = "<<c<<endl;
#endif
  return g;
}

//#define testbezout

// returns content and sets c such that content(a) = a.c
INT vecbezout(vector<INT>& a, vector<INT>& c)
{
#ifdef testbezout
  cout<<"Computing vecbezout("<<a<<")"<<endl;
#endif
  int n=(int)a.size();
  if (n==2) return vecbezout2(a, c);
  if (n==3) return vecbezout3(a, c);
  INT x(1), g = vecgcd(a);
  vector<INT> a0=a;
  if (g>1)
    for(vector<INT>::iterator ai=a0.begin(); ai!=a0.end(); ++ai)
      (*ai) /= g;
  // Now a0 is primitive: we do this to make numbers smaller in what follows
  c = vector<INT>(n, INT(0));
  INT g1(0);
  for(int i=0; i<n &&g1!=1; i++)
    {
      g1=bezout(g1,a0[i],x,c[i]);
      for(int j=0; j<i; j++) c[j]*=x;
    }
#ifdef testbezout
  cout<<"vecbezout("<<a<<") returns "<<g<<", with coefficients "<<c<<endl;
#endif
  return g;
}

INT dot(vector<INT>& a, vector<INT>& c)
//returns g = a.c
{
  return std::inner_product(a.begin(), a.end(), c.begin(), INT(0));
}

//#define testbezout

// If (b,d)!=(0,0), return an HNF basis [(aa,bb), (cc,0)] for the
// Z-module [[a,b],[c,d]]:
vector<INT> hnf22(INT a, INT b, INT c, INT d)
{
#ifdef testbezout
  cout<<"  - hnf22("<<a<<", "<<b<<", "<<c<<", "<<d<<") = "<<flush;
#endif
  INT x,y;
  INT bb = bezout(b,d, x,y);
  INT cc = abs(a*(d/bb)-(b/bb)*c);
  INT aa = a*x+c*y;
  if (!is_zero(cc)) aa = posmod(aa,cc);
  vector<INT> v = {aa, bb, cc};
#ifdef testbezout
  cout<<v<<endl;
#endif
  return v;
}

//#define testbezout

//Sets basis={e1,e2,f1} such that [[e1,f1], [e2,0]] is a Z-basis for
//the Z-module spanned by [first[i], second[i]]

void findzbasis(const vector<INT>& first, const vector<INT>& second, vector<INT>& basis)
{
#ifdef testbezout
  cout<<"findzbasis("<<first<<", "<<second<<"): "<<endl;
#endif
  INT a, b, c, d;
  vector<INT>::const_iterator ai, bi, aj, bj;
  // Find a nonsingular 2x2 block:
  int t=1;
  for (ai=first.begin(), bi=second.begin(); t && ai!=first.end(); ++ai, ++bi)
    {
      a = *ai;
      b = *bi;
      for (aj=ai+1, bj=bi+1; t && aj!=first.end(); ++aj, ++bj)
        {
          c = *aj;
          d = *bj;
          if (a*d != b*c)
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
  assert (b*c!=0);
#ifdef testbezout
  cout<<" - after step 0, {a,b,c} = "<<basis<<endl;
#endif
  // process all the rest
  for (ai=first.begin(), bi=second.begin(); ai!=first.end(); ++ai, ++bi)
    {
      INT e = *ai %c, f = *bi;
      basis = hnf22(a,b, e, f);
#ifdef testbezout
      cout<<" - sub-basis using ("<<(*ai)<<","<<(*bi)<<"), {a,b,c} = "<<basis<<endl;
#endif
      c = gcd(c, basis[2]);
      a = basis[0]%c;
      b = basis[1];
      basis = {a, b, c};
#ifdef testbezout
      cout<<" - after one step using ("<<(*ai)<<","<<(*bi)<<"), {a,b,c} = "<<basis<<endl;
#endif
    }
}

int div_disc(INT D1, INT D)
{
  if (!divides(D1,D))
    return 0;
  INT d = posmod(D/D1, INT(4));
  return (d==0 || d==1);
}
