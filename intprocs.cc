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

//functions needed for non-euclidean fields to compute bezout/quadgcd

// returns content and sets c such that content(a) = a.c
long vecbezout(const vector<long>& a, vector<long>& c)
{
  long x=1,g=0;
  int n=(int)a.size();
  c.resize(n);
  for(int i=0; i<n; i++)
    {
      g=bezout(g,a[i],x,c[i]);
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

// content
long vecgcd(const vector<long>& a)
{
  //  return std::accumulate(a.begin(), a.end(), 0, gcd);
  long g=0;
  int n=(int)a.size();
  for(int i=0; (i<n)&(g!=1); i++) g=gcd(g,a[i]);
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
  for(g=0, ai=a.begin(), ci=c.begin(); ai!=a.end(); ai++, ci++)
    g += (double(*ai) * (double)(*ci));
  return g;
}

// inner product mod s
long xmoddot(long s, const vector<long>& a, const vector<long>& c)
{
  long g;
  vector<long>::const_iterator ai,ci;
  for(g=0, ai=a.begin(), ci=c.begin(); ai!=a.end(); ai++, ci++)
    g = xmod( g + xmodmul(*ai, *ci,s), s);
  return g;
}

// sets basis={e1,e2,f1} such that [[e1,f1], [e2,0]] is a Z-basis
// for the Z-module spanned by [first[i], second[i]], and also sets x, y to be vectors such that
//
// first.x = e1
// first.y = e2
// second.x = f1
// second.y = 0

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
  for(i=0; i<n; i++) newfirst[i] = first[i] - basis[0]*second[i]/basis[1];
  basis[2] = vecbezout(newfirst,u);
//  y = u - ((u*second)/basis[1])*x;
  long t = dot(u,second);
  for(i=0; i<n; i++) y[i]=u[i]-(t*x[i])/basis[1];
#ifdef testbezout
//Check:
  if( ! (  (basis[0]==dot(first,x))   &&
           (basis[2]==dot(first,y))   &&
           (basis[1]==dot(second,x))  &&
           (0==dot(second,y)) ))
  {cerr<<"Error in findzbasis!"  <<endl; }
#endif
}

//Same as findzbasiscoeffs except don't need x,y: sets
//basis={e1,e2,f1} such that [[e1,f1], [e2,0]] is a Z-basis for the
//Z-module spanned by [first[i], second[i]]

void findzbasis(const vector<long>& first, const vector<long>& second, vector<long>& basis)
{
  vector<long> x(first.size()), y(first.size());
  findzbasiscoeffs(first, second, basis, x, y);
}

// returns a s.t.  s*a + (f dot x) equals d
long special_solve_for_last_coeff(long d, vector<long>& f, vector<long>& x, long s)
{
  double a = d - ddot(f, x);
  if ( fmod(a, s) != 0 )
    {
      cerr << "Error: inexact division in special_solve_for_last_coeff!"<<endl;
    }
  a/=s;
  if (fabs(a)>MAXLONG)
    {
      cerr << "Error: overflow in special_solve_for_last_coeff!"<< endl;
    }
  return (long)a;
}

void specialfindzbasiscoeffsmod(const vector<long>& first, const vector<long>& second,
                                vector<long>& basis, vector<long>& x, vector<long>& y)
{
  int n=(int)first.size();
  if ((n<2)|| (second[n-2] != 0) || (first[n-1]!=0) ||
      (first[n-2]==0) || (first[n-2] != second[n-1]))
    { cerr << "Invalid call to specialfindzbasiscoeffsmod!" << endl;
      exit(1);
    }
  long s=first[n-2];
  vector<long> u(n), newfirst(n), newsecond(n);

  long e2=xmodvecbezout(s,second,x);
//Now second.x=e2 ...
// ... except that x[n-2] is arbitrary (currently zero) and x[n-1] is wrong !

  long e1=xmoddot(s,first,x);                  //dot product
//Now  first.x=e1,  for suitable choice of x[n-2] (depending on correct
// value of x[n-1]

//That is, [e1,e2] is the x-combination of the data, with e2=gcd(second) ...
// ... except that the currently held values of x[n-2], x[n-1] are incorrect !

//newfirst = first-e1*(second/e2);
// (It was to stop this step overflowing that we chose a small e1.)
  for(long i=0; i<n; i++)
    {
      newsecond[i]=second[i]/e2;
      if (i==n-2) newfirst[i]=s;  // prevent reduction of s to zero !
      else newfirst[i]=xmod(first[i]-xmodmul(e1,newsecond[i],s),s);
    }
// Since newfirst[n-2]=s, f1 comes out OK even though rest of row is mod s

  long f1 = xmodvecbezout(s,newfirst,u);
//Now newfirst.u=f1, except that u[n-2] is wrong!

  basis = {e1, e2, f1};

  long t = xmoddot(s,u,newsecond);
// correct despite u[n-2] being wrong - because newsecond[n-2]=0 !

//  y = u - ((second/e2)*u)*x;
  for(long i=0; i<n-2; i++) y[i]=xmod(u[i]-xmodmul(t,x[i],s),s);

// So far, so good.
// Now try computing x[n-2], x[n-1], y[n-2], y[n-1] !!

// The trouble with the naive method is that it overflows too easily.
// x[n-2] = (e1 - dot(n-2, first,x))/s;
// x[n-1] = (e2 - dot(n-2,second,x))/s;

// y[n-2] = (f1 - dot(n-2, first,y))/s;
// y[n-1] = ( 0 - dot(n-2,second,y))/s;

  vector<long> short_x(x.begin(), x.begin()+n-2);
  vector<long> short_y(y.begin(), y.begin()+n-2);
  vector<long> short_first(first.begin(), first.begin()+n-2);
  vector<long> short_second(second.begin(), second.begin()+n-2);
  x[n-2] = special_solve_for_last_coeff(e1,  short_first,  short_x, s);
  x[n-1] = special_solve_for_last_coeff(e2,  short_second, short_x, s);
  y[n-2] = special_solve_for_last_coeff(f1,  short_first,  short_y, s);
  y[n-1] = special_solve_for_last_coeff( 0,  short_second, short_y, s);

if (0)
    { cerr <<"Reporting from specialfindzbasiscoeffsmod!"  <<endl;
      cerr <<"modulus s = "<<s<<endl;

      cerr <<"first = " << first <<endl;
      cerr << "second= " << second << endl;

      cerr << "e1 = " << e1 << endl;
      cerr << "e2 = " << e2 << endl;
      cerr << "f1 = " << f1 << endl;

      cerr <<" x = " << x << endl;

      cerr <<" y = " << y << endl;

      cerr <<"Warning: x[n-2], x[n-1], y[n-2], y[n-1] are less robust!"<<endl;
    }

}

