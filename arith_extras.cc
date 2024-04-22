// FILE ARITH_EXTRAS.CC: miscellaneous integer (int/long) functions

#include <assert.h>
#include <numeric>
#include "arith_extras.h"

void sqrt_mod_p(long & x, long a, long p)
{
  bigint rr;
  sqrt_mod_p(rr, BIGINT(posmod(a,p)), BIGINT(p));
  x = I2long(rr);
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
