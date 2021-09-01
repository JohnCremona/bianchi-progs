// ratquads.cc

#include "ratquads.h"
#include "qideal.h"
#include "mat22.h"
#include "geometry.h"

// Constants
RatQuad RatQuad::oo(1,0,1, 0);
RatQuad RatQuad::one(1,1,1, 0);
RatQuad RatQuad::zero(0,1,1, 0);

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
  if (d==0)
    return Qideal(0);
  return Qideal(d) / Qideal({n,d});
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
  if (c1==c2) return 1;
#ifdef DEBUG_CUSP_EQ
  cout<<"Testing equivalence of cusps "<<c1<<" and "<<c2;
  cout<<" (N="<<N<<")"<<endl;
#endif
  Quad q1 = c1.d, q2 = c2.d, q3, s1,r1,s2,r2;
  quadbezout(c1.n,q1,s1,r1);  s1*=q2;
  quadbezout(c2.n,q2,s2,r2);  s2*=q1;
  q3 = quadgcd(q1*q2,N);
#ifdef DEBUG_CUSP_EQ
  cout<<" - s1 =  "<<s1<<", s2 = " << s2 << ", q3 = "<<q3<<endl;
#endif
  int equiv=0;
  vector<Quad>& units = (plusflag? quadunits: squareunits);
  for (vector<Quad>::const_iterator u = units.begin(); u!= units.end() && !equiv; ++u)
    {
#ifdef DEBUG_CUSP_EQ
      cout<<" - testing unit "<<(*u)<<endl;
#endif
      equiv = div(q3,(s1-(*u)*s2));
    }
#ifdef DEBUG_CUSP_EQ
  cout<<"Returning "<<equiv<<endl;
#endif
  return equiv;
}

// General cusp equivalence modulo Gamma_0(N) where N is an ideal:

int cuspeq(const RatQuad& c1, const RatQuad& c2, const Qideal& N, int plusflag)
{
  if (c1==c2) return 1;
#ifdef DEBUG_CUSP_EQ
  cout<<"Testing equivalence of cusps "<<c1<<" and "<<c2;
  cout<<" (N="<<N<<")"<<endl;
#endif
  // test ideals are in the same class:
  Qideal I1 = c1.ideal(), I2 = c2.ideal();
  if (!I1.is_equivalent(I2))
    {
#ifdef DEBUG_CUSP_EQ
      cout << " - ideals "<<I1<<", "<<I2<<" are not equivalent"<<endl;
      cout<<" - Returning 0"<<endl;
#endif
      return 0;
    }
#ifdef DEBUG_CUSP_EQ
  cout << " - cusps are in the same ideal class"<<endl;
  cout << " - absolute denominator ideals: " << c1.denominator_ideal() << " and " << c2.denominator_ideal() << endl;
#endif
  // denominator test:
  Qideal D1 = c1.denominator_ideal()+N,
    D2 = c2.denominator_ideal()+N;
#ifdef DEBUG_CUSP_EQ
  cout << " - relative denominator ideals: " << D1 << " and " << D2 << endl;
#endif
  if (D1 != D2)
    {
#ifdef DEBUG_CUSP_EQ
      cout << " - denominator ideals "<<D1<<", "<<D2<<" are not equal"<<endl;
      cout<<" - Returning 0"<<endl;
#endif
      return 0;
    }

#ifdef DEBUG_CUSP_EQ
  cout<<" - denominator test passes, denominator ideal = "<<D1<<endl;
#endif
  // adjust representations so that ideals are coprime to N and equal:
  RatQuad cc1 = c1, cc2 = c2;
  cc1.reduce(N);
  cc2.reduce(N);
  I1 = cc1.ideal();
  I2 = cc2.ideal();
  assert (I1==I2);
  assert (N.is_coprime_to(I1));
#ifdef DEBUG_CUSP_EQ
  cout<<" - adjusted representations: "<<cc1<<", "<<cc2<<" with equal ideals "<<I1<<" coprime to N"<<endl;
#endif

  // Use the criterion of Cor.2 in "Manin symbols over number fields"

  // Form ab-matrices with first columns equal to the two cusp representations

  mat22 M1 = AB_matrix(cc1.n, cc1.d), // [a1,b1;a2,b2]
    M2 = AB_matrix(cc2.n, cc2.d);     // [a1',b1';a2',b2']
  assert (M1.det()==M2.det());
  Quad a2db2 = M2.entry(1,0)*M1.entry(1,1), a2b2d = M1.entry(1,0)*M2.entry(1,1);
#ifdef DEBUG_CUSP_EQ
  cout<<" - A = "<<M1<<", det = "<<M1.det()<<endl;
  cout<<" - B = "<<M2<<", det = "<<M2.det()<<endl;
  cout<<" - a2'*b2 = "<<a2db2<<", a2*b2' = "<<a2b2d<<endl;
#endif

  Qideal M = D1*D1+N; // D1==D2
  M = I1*I2*M;
#ifdef DEBUG_CUSP_EQ
  cout<<" - M = I1*I2*(D^2+N) = "<<M<<endl;
#endif

  // test whether there is a unit u such that
  // (1) a2'*b2 = u*a2*b2' (mod M)

  if (M.divides(a2db2-a2b2d))
    {
#ifdef DEBUG_CUSP_EQ
      cout<<" - Returning 1"<<endl;
#endif
      return 1;
    }
  vector<Quad>& units = (plusflag? quadunits: squareunits);
  for (vector<Quad>::const_iterator ui = units.begin()+1; ui!= units.end(); ++ui)
    {
      Quad u = *ui;
#ifdef DEBUG_CUSP_EQ
      cout<<" - testing unit "<<u<<endl;
#endif
      if (M.divides(a2db2-u*a2b2d))
        {
#ifdef DEBUG_CUSP_EQ
          cout<<" - Returning 1"<<endl;
#endif
          return 1;
        }
    }
#ifdef DEBUG_CUSP_EQ
  cout<<" - Returning 0"<<endl;
#endif
  return 0;
}

// if type = t>=0 , return U{alpha[t],oo}
// if type = -s<0 , return U{sigmas[s],oo}

modsym::modsym(const mat22& U, int type) // U in SL2
{
  RatQuad alpha = (type>=0? alphas[type]: sigmas[-type]);
  a = U(alpha);
  b = U.image_oo();
}
