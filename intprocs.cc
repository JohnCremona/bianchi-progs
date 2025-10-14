// FILE INTPROCS.CC

#include <assert.h>
#include <numeric>
#include "intprocs.h"
#include "flint_snf.h"

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

INT content(const vector<INT>& a)
{
  INT g(0);
  for( const auto& ai : a)
    {
      g = gcd(g, ai);
      if (g.is_one())
        return g;
    }
  return g;
}

// returns content and sets c such that content(a) = a.c
INT vecbezout2(const vector<INT>& a, vector<INT>& c)
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
INT vecbezout3(const vector<INT>& a, vector<INT>& c)
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
INT vecbezout(const vector<INT>& a, vector<INT>& c)
{
#ifdef testbezout
  cout<<"Computing vecbezout("<<a<<")"<<endl;
#endif
  int n=(int)a.size();
  if (n==2) return vecbezout2(a, c);
  if (n==3) return vecbezout3(a, c);
  INT x(1), g = content(a);
  vector<INT> a0=a;
  if (g>1)
    std::for_each(a0.begin(), a0.end(), [g](INT& ai) {ai/=g;});
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

// If (b,d)!=(0,0), return {aa,bb,cc} where [(aa,bb), (cc,0)] is an HNF basis
// for the Z-module [[a,b],[c,d]]:
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

// Test if the Z-module spanned by [coords[0][i], coords[1][i]] is all of Z^2 (for coprimality testing)
int span_Z2(const pair< vector<INT>, vector<INT>> & coords)
{
  if (coords.first.size()!=coords.second.size()) return 0;
  if (coords.first.size()<2) return 0;
  INT g(0);
  for (auto ai=coords.first.begin(), bi=coords.second.begin(); g!=1 && ai!=coords.first.end(); ++ai, ++bi)
    for (auto aj=ai+1, bj=bi+1; g!=1 && aj!=coords.first.end(); ++aj, ++bj)
      g = gcd(g, ((*ai)*(*bj)-(*aj)*(*bi)));
  return g==1;
}

// Returns {a,b,c} where [(a,b), (c,0)] is an HNF basis for the
// Z-module spanned by [coords.first[i], coords.second[i]]

vector<INT> Zbasis(const pair<vector<INT>, vector<INT>>& coords)
{
  vector<vector<INT>> M = {coords.second, coords.first};
  M = HNF(M);
  return {M[1][0], M[0][0], M[1][1]};
}

vector<INT> oldZbasis(const pair<vector<INT>, vector<INT>>& coords)
{
#ifdef testbezout
  cout<<"findzbasis("<<coords.first<<", "<<coords.second<<"): "<<endl;
#endif
  INT a, b, c, d;
  vector<INT> basis;
  // Find a nonsingular 2x2 block:
  int t=1;
  for (auto ai=coords.first.begin(), bi=coords.second.begin(); t && ai!=coords.first.end(); ++ai, ++bi)
    {
      a = *ai;
      b = *bi;
      for (auto aj=ai+1, bj=bi+1; t && aj!=coords.first.end(); ++aj, ++bj)
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
  for (auto ai=coords.first.begin(), bi=coords.second.begin(); ai!=coords.first.end(); ++ai, ++bi)
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
  return basis;
}

int div_disc(INT D1, INT D)
{
  if (!divides(D1,D))
    return 0;
  INT d = posmod(D/D1, INT(4));
  return (d==0 || d==1);
}

// Return a list of all vectors of length dim which are primitive,
// with all entries <= bound (in absolute value), modulo
// multiplication by -1 (the first nonzero entry in each will be
// positive).
vector<vector<int>> all_linear_combinations(int dim, int bound)
{
  if (dim<1) return {};
  if (dim==1) return {{1}};
  vector<int> v(dim,0);
  v[0] = 1;
  vector<vector<int>> ans = {v};
  // recurse
  vector<vector<int>> ans1 = all_linear_combinations(dim-1, bound);
  for (auto v1: ans1)
    {
      int maxv1 = * std::max_element(v1.begin(), v1.end(), [](int a, int b){return abs(a)<abs(b);});
      int max_scale = bound / maxv1; // round down
      v = v1;
      v.insert(v.begin(), 0);
      ans.push_back(v);
      //
      // for (int i=0; i<dim; i++)
      //   cout << " ";
      // cout << v << endl;
      //
      for (int b=-max_scale; b<=max_scale; b++)
        {
          if (b==0)
            continue;
          vector<int> vb = v;
          std::transform(v.begin(), v.end(), vb.begin(), [b](int c){return b*c;});
          for (int a=1; a<=bound; a++)
            {
              if (gcd(a,b)==1)
                {
                  vb[0] = a;
                  ans.push_back(vb);
                  //
                  // for (int i=0; i<dim; i++)
                  //   cout << " ";
                  // cout << vb << endl;
                  //
                }
            }
        }
    }
  return ans;
}
