// FILE INTPROCS.CC

#include <eclib/arith.h>
#include <assert.h>
#include <numeric>
#include "intprocs.h"

int is_nonnegative(QUINT a)
{
  return is_positive(a) || is_zero(a);
}

long squarefree_part(long d)
{
  if (d==0) return d;
  vector<long> sd = sqdivs(d);
  long maxd = sd[sd.size()-1];
  long ans = d/(maxd*maxd);
  //cout << "d has max square divisor "<<maxd<<"^2"<<" and squarefree part "<<ans<<endl;
  return ans;
}

// The point of the following function is that the built-in gcc
// division truncates towards 0, while we need rounding, with a
// consistent behaviour for halves (they go up here).
//
// For b>0, rounded_division(a,b) = q such that a/b = q + r/b with -1/2 <= r/b < 1/2
//
// bigint version of a similar function roundover(a,b) is in marith.h implmented as
// {bigint a0=(a%b); bigint c = (a-a0)/b; if(2*a0>b) c+=1; return c;}
// which is not quite the same

long rounded_division(long a, long b)
{
  std::ldiv_t qr = ldiv(a, b);
  long r = qr.rem, q = qr.quot;
  long r2 = r<<1;
  return (r2<-b? q-1: (r2>=b? q+1: q));
}

QUINT rounded_division(QUINT a, QUINT b)
{
  QUINT q, r;
  ::divides(a,b, q,r);
  assert (a==b*q+r);
  QUINT r2 = 2*r;
  // We want -b <= r2 < +b
  if (r2<-b)
    q-=1;
  else
    if (r2>=b)
      q+=1;
  return q;
}

QUINT vecgcd(const vector<QUINT>& a)
{
  QUINT g(0);
  for(vector<QUINT>::const_iterator ai=a.begin(); ai!=a.end() && (g!=1); ++ai)
    g = gcd(g, *ai);
  return g;
}

// returns content and sets c such that content(a) = a.c
QUINT vecbezout2(const vector<QUINT>& a, vector<QUINT>& c)
{
  int n=(int)a.size();
  if (n!=2) return vecbezout(a, c);
  QUINT x, y;
  QUINT g = bezout(a[0], a[1], x, y);
  c = {-y, x};
  assert (a[0]*c[0]+a[1]*c[1]==g);
  return g;
}

//#define testbezout3

// returns content and sets c such that content(a) = a.c
QUINT vecbezout3(const vector<QUINT>& a, vector<QUINT>& c)
{
  int n=(int)a.size();
  if (n!=3) return vecbezout(a, c);

  QUINT aa=a[0], bb=a[1], cc=a[2];
  QUINT ZERO(0), ONE(1);
#ifdef testbezout3
  cout<<"Computing vecbezout3("<<a<<")"<<endl;
#endif
  if ((is_zero(aa))&&(is_zero(bb)))
    {
      if (is_negative(cc))
        {
          c = {ZERO,ZERO,-ONE};
          return -cc;
        }
      else
        {
          c = {ZERO,ZERO,ONE};
          return cc;
        }
    }
  QUINT x, y, z, w;
  QUINT h = bezout(aa,bb, x, y);
  QUINT a1 = aa/h, b1 = bb/h;
  QUINT g = bezout(h,cc, z, w);
  QUINT h1=h/g, c1=cc/g;
  c = {x*z, y*z, w};
  assert (a[0]*c[0]+a[1]*c[1]+a[2]*c[2]==g);
#ifdef testbezout3
  cout<<"   g="<<g<<", first solution is "<<c<<"; "<<flush;
  vector<QUINT> perp1 = {-c1*x, -c1*y, h1};
  vector<QUINT> perp2 = {-b1, a1, 0};
  cout<<"   primitive basis of perp: "<<perp1<<", "<<perp2<<endl;
  assert (a[0]*perp1[0]+a[1]*perp1[1]+a[2]*perp1[2]==0);
  assert (a[0]*perp2[0]+a[1]*perp2[1]+a[2]*perp2[2]==0);
#endif
  // now minimize
  QUINT lambda, mu;
#ifdef testbezout3
  QUINT x2y2 = x*x+y*y;
  lambda = rounded_division(x2y2*c1*z-h1*w, x2y2*c1*c1+h1*h1);
  mu = rounded_division((b1*x-a1*y)*z, a1*a1+b1*b1);
  cout << " (lambda,mu)=("<<lambda<<","<<mu<<")"<<endl;
#endif
  // bigfloat rx2y2 = pow(to_bigfloat(x),2) + pow(to_bigfloat(y),2);
  // bigfloat rlambda = (rx2y2*c1*z-h1*w) / (rx2y2*pow(to_bigfloat(c1),2)+pow(to_bigfloat(h1),2));
  // bigfloat rmu = to_bigfloat(b1*x-a1*y)*z /  (pow(to_bigfloat(a1),2)+pow(to_bigfloat(b1),2));
  // longify(rlambda, lambda);
  // longify(rmu, mu);
  //#ifdef testbezout3
  //  cout << "--using bigfloats, (lambda,mu)=("<<lambda<<","<<mu<<")"<<endl;
  //#endif
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
QUINT vecbezout(const vector<QUINT>& a, vector<QUINT>& c)
{
#ifdef testbezout
  cout<<"Computing vecbezout("<<a<<")"<<endl;
#endif
  int n=(int)a.size();
  if (n==2) return vecbezout2(a, c);
  if (n==3) return vecbezout3(a, c);
  QUINT x(1), g = vecgcd(a);
  vector<QUINT> a0=a;
  if (g>1)
    for(vector<QUINT>::iterator ai=a0.begin(); ai!=a0.end(); ++ai)
      (*ai) /= g;
  // Now a0 is primitive: we do this to make numbers smaller in what follows
  c = vector<QUINT>(n, BIGINT(0));
  QUINT g1(0);
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

// returns g = content(a) and c such that g = a.c, with c reduced mod s
QUINT xmodvecbezout(QUINT s, const vector<QUINT>& a, vector<QUINT>& c)
{
  int n=(int)a.size();
  int i, j;
  QUINT x(1), ci(1), g(0);
  for(i=0; i<n; i++)
    {
      g = bezout(g,a[i],x,ci);
      c[i] = mod(ci,s);
      for(j=0; j<i; j++)
        c[j] = mod(c[j]*x,s);
    }
  return g;
}

QUINT dot(const vector<QUINT>& a, const vector<QUINT>& c)
//returns g = a.c
{
  return std::inner_product(a.begin(), a.end(), c.begin(), BIGINT(0));
}

//#define testbezout

// If (b,d)!=(0,0), return an HNF basis [(aa,bb), (cc,0)] for the
// Z-module [[a,b],[c,d]]:
vector<QUINT> hnf22(QUINT a, QUINT b, QUINT c, QUINT d)
{
#ifdef testbezout
  cout<<"  - hnf22("<<a<<", "<<b<<", "<<c<<", "<<d<<") = "<<flush;
#endif
  QUINT x,y;
  QUINT bb = bezout(b,d, x,y);
  QUINT cc = abs(a*(d/bb)-(b/bb)*c);
  QUINT aa = a*x+c*y;
  if (!is_zero(cc)) aa = posmod(aa,cc);
  vector<QUINT> v = {aa, bb, cc};
#ifdef testbezout
  cout<<v<<endl;
#endif
  return v;
}

//#define testbezout

//Sets basis={e1,e2,f1} such that [[e1,f1], [e2,0]] is a Z-basis for
//the Z-module spanned by [first[i], second[i]]

void findzbasis(const vector<QUINT>& first, const vector<QUINT>& second, vector<QUINT>& basis)
{
#ifdef testbezout
  cout<<"findzbasis("<<first<<", "<<second<<"): "<<endl;
#endif
  QUINT a, b, c, d;
  vector<QUINT>::const_iterator ai, bi, aj, bj;
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
      QUINT e = *ai %c, f = *bi;
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
