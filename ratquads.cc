// ratquads.cc

#include "ratquads.h"
#include "qideal.h"
#include "mat22.h"
#include "geometry.h"

// reduce to lowest terms: when non-principal this method only divides
// n and d by the gcd of the content.  Returns 1 iff principal.

int RatQuad::reduce()
{
  // first divide out by content
  long c = gcd(n.content(), d.content());
  if (c>1) {n/=c; d/=c;}

  // find gcd(n,d) when ideal (n,d) is principal (will return 0 if ideal not principal):
  Quad g = quadgcd(n,d);
  long ng = quadnorm(g);
  if (ng>1) {n/=g; d/=g;}

  // final adjustment by units:
  while (!pos(d)) {n*=fundunit; d*=fundunit;}
  return (ng>0);
}

void RatQuad::reduce(const Qideal& N)
{
  if (reduce()) return; // if principal then after reduce() the ideal is (1), so nothing more to do
  Qideal I = ideal();
  Quad a,b;
  I.equivalent_coprime_to(N, a, b); // (a/b)*(n,d) is coprime to N
  n *= a; n /= b;
  d *= a; d /= b;
  while (!pos(d)) {n*=fundunit; d*=fundunit;}
}

Qideal RatQuad::ideal() const
{
  return Qideal({n,d});
}

Qideal RatQuad::denominator_ideal() const
{
  return Qideal(d)/Qideal({n,d});
}

int RatQuad::is_principal() const
{
  return ideal().is_principal();
}

//#define DEBUG_CUSP_EQ

// The following, where the modulus is a Quad, only works for two
// principal cusps: we assume that c1 and c2 are reduced (coprime
// numerator and denominator).

int cuspeq(const RatQuad& c1, const RatQuad& c2, const Quad& N, int plusflag)
{
#ifdef DEBUG_CUSP_EQ
  cout<<"Testing equivalence of cusps "<<c1<<" and "<<c2;
  cout<<" (N="<<N<<")"<<endl;
#endif
  Quad q1 = c1.d, q2 = c2.d, q3, s1,r1,s2,r2;
  quadbezout(c1.n,q1,s1,r1);  s1*=q2;
  quadbezout(c2.n,q2,s2,r2);  s2*=q1;
  q3 = quadgcd(q1*q2,N);
#ifdef DEBUG_CUSP_EQ
  cout<<"s1 =  "<<s1<<", s2 = " << s2 << ", q3 = "<<q3<<endl;
#endif
  int equiv=0;
  vector<Quad> units = (plusflag? quadunits: squareunits);
  for (vector<Quad>::const_iterator u = units.begin(); u!= units.end() && !equiv; u++)
    equiv = div(q3,(s1-(*u)*s2));
#ifdef DEBUG_CUSP_EQ
  cout<<"Returning "<<equiv<<endl;
#endif
  return equiv;
}

// General cusp equivalence modulo Gamma_0(N) where N is an ideal:

int cuspeq(const RatQuad& c1, const RatQuad& c2, const Qideal& N, int plusflag)
{
#ifdef DEBUG_CUSP_EQ
  cout<<"Testing equivalence of cusps "<<c1<<" and "<<c2;
  cout<<" (N="<<N<<")"<<endl;
#endif
  // test ideals are in the same class:
  Qideal I1 = c1.ideal(), I2 = c2.ideal();
  if (!I1.is_equivalent(I2)) return 0;

  // denominator test:
  Qideal D = c1.denominator_ideal()+N;
  if (D != c2.denominator_ideal()+N) return 0;

  // adjust representations so that ideals are coprime to N:
  RatQuad cc1(c1), cc2(c2);
  cc1.reduce(N);
  cc2.reduce(N);

  // test whether there is a unit u such that
  // (1) d2 = u*d1 (mod N)
  // (2) n1 = u*n2 (mod D)

  vector<Quad> units = (plusflag? quadunits: squareunits);
  for (vector<Quad>::const_iterator ui = units.begin(); ui!= units.end(); ui++)
    {
      Quad u = *ui;
      if (N.divides(cc2.d-u*cc1.d) && D.divides(cc1.n-u*cc2.n))
        return 1;
    }
  return 0;
}

// if type = t>=0 , return U{alpha[i],oo}
// if type = -s<0 , return U{sigmas[s],oo}
modsym::modsym(const mat22& U, int type) // U in SL2
{
  RatQuad alpha;
  if(type==0) // always true for Euclidean fields: apply to {0,oo}
    alpha = RatQuad(0);
  else
   if (type>0) // apply to {alpha,oo} where alpha = alphas[type]
     {
       mat22 M = M_alphas[type];
       alpha = RatQuad(-M.d, M.c);
     }
   else // type=t<0 means apply to {sigma,oo} where sigma = sigmas[-t]
     {
       alpha = sigmas[-type];
     }
  a = U(alpha);
  b = U(RatQuad(1,0)); // = U(oo); no need to reduce as n,d are coprime
}
