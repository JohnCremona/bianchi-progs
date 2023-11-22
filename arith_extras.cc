// FILE ARITH_EXTRAS.CC: miscellaneous integer (int/long) functions which could be in eclib/arith.cc

#include <assert.h>
#include <numeric>
#include "arith_extras.h"

int isqrt(long a, long& root)
{
  if (a<0) {return 0;}
  root = round(sqrt(a));
  return a==root*root;
}

long isqrt(const long a) {long r; isqrt(a,r); return r;}

int divrem(long a, long b, long& q, long& r)
{
  std::ldiv_t qr = ldiv(a, b);
  r = qr.rem;
  q = qr.quot;
  return (r==0);
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

void sqrt_mod_p(long & x, long a, long p)
{
  bigint rr;
  sqrt_mod_p(rr, BIGINT(posmod(a,p)), BIGINT(p));
  x = I2long(rr);
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

// return list of integers from first to last inclusive
vector<long> range(long first, long last)
{
  long n = last-first+1;
  vector<long> ans(n);
  std::iota(ans.begin(), ans.end(), first);
  return ans;
}

// return 1 with r=sqrt(a) if a is square, else return 0:
long is_square(long a, long& r)
{
  if (a<0) return 0;
  r = (long)(sqrt((double)a)+0.001);
  return (a == r*r);
}

// return the dot product (0/1) of a and b bitwise, with 0<=a,b<2^r
int dotbits(long a, long b, int r)
{
  int x=0;
  for (int i=0; i<r; i++)
    x ^= (bit(a,i)*bit(b,i));
  return x;
}

// return list of bits of a
vector<int> bits(long a, int r)
{
  vector<int> ans(r);
  for(int i=0; i<r; i++)
    ans[i] = bit(a,i);
  // assert (a==from_bits(ans,r));
  return ans;
}

// recover a from its bit vector of length r
long from_bits(vector<int> aa, int r)
{
  long a=0;
  for(int i=0; i<r; i++)
    if (aa[i])
      a |= (1<<i); // sets the i'th bit to 1
  // assert (aa==bits(a,r));
  return a;
}

// return a basis for the orthogonal complement of a<2^r (viewed as a bit vector of length r)
vector<long> dotperp(long a, int r)
{
  if (a==0) // trivial special case
    {
      vector<long> ans(r);
      for (int i=0; i<r; i++)
        ans[i] = 1<<i;
      return ans;
    }
  else
    {
      mat_i m(1,r);
      for (int j=1; j<=r; j++)
        m.set(1,j,bit(a,j-1));
      subspace_i ker = pkernel(m,2); // right kernel mod 2
      assert (dim(ker)==r-1);
      mat_i bas = basis(ker);
      vector<long> ans(r-1, 0);
      for (int i=0; i<r-1; i++)
        {
          vec_i coli = bas.col(i+1);
          for (int j=0; j<r; j++)
            if (coli[j+1])
              ans[i] |= 1<<j;
        }
      return ans;
    }
}

// return a basis for the orthogonal complement of the span of a in alist (viewed as bit vectors of length r)
vector<long> dotperp(vector<long> alist, int r)
{
  int s = alist.size();
  mat_i m(s,r);
  for (int i=1; i<=s; i++)
    for (int j=1; j<=r; j++)
      m.set(i,j,bit(alist[i-1],j-1));
  subspace_i ker = pkernel(m,2); // right kernel mod 2
  assert (dim(ker)==r-s);
  mat_i bas = basis(ker);
  vector<long> ans(r-s, 0);
  for (int i=0; i<r-s; i++)
    {
      vec_i coli = bas.col(i+1);
      for (int j=0; j<r; j++)
        if (coli[j+1])
          ans[i] |= 1<<j;
    }
  return ans;
}

vec reduce_modp(const vec& v, const scalar& p)
{
  if (p==0) return v;
  long i, d=dim(v);
  vec ans(d);
  for(i=1; i<=d; i++)
    ans[i] = mod(v[i], p);
  return ans;
}

mat reduce_modp(const mat& m, const scalar& p)
{
  if (p==0) return m;
  long i, j, nr=m.nrows(), nc=m.ncols();
  mat ans(nr,nc);
  for(i=1; i<=nr; i++)
    for(j=1; j<=nc; j++)
      ans(i,j) = mod(m(i,j),p);
  return ans;
}
