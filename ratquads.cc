// ratquads.cc

#include "ratquads.h"
#include "qideal.h"
#include "mat22.h"
#include "geometry.h"

Cusp_comparison Cusp_cmp;

// reduce to lowest terms: when non-principal this method only divides
// n and d by the gcd of the content.  Returns 1 iff principal.

int RatQuad::reduce()
{
  // first divide out by content
  INT c = gcd(n.content(), d.content());
  if (c>1) {n/=c; d/=c;}

  // find gcd(n,d) when ideal (n,d) is principal (will return 0 if ideal not principal):
  Quad g = quadgcd(n,d);
  INT ng = quadnorm(g);
  if (ng>1) {n/=g; d/=g;}

  // final adjustment by units:
  while (!pos(d)) {n*=fundunit; d*=fundunit;}
  return (ng>0);
}

void RatQuad::reduce(long n)
{
  reduce(Qideal(Quad(INT(n))));
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
  if (d.is_zero())
    return Qideal(ZERO);
  return Qideal(d) / Qideal({n,d});
}

int RatQuad::is_principal() const
{
  if (Quad::class_number==1)
    return 1;
  else
    // quadgcd(n,d) returns 0 if n,d do not generate a principal ideal
    return coprime(n,d); //quadgcd(n,d) != Quad::zero;
}

void RatQuad::normalise()                              // scale so ideal is a standard class rep
{
  Qideal I({n,d});
  Qideal J = class_representative(I);
  I *= J.conj();
  Quad g;
  I.is_principal(g);
  n *= J.norm(); n /= g;
  d *= J.norm(); d /= g;
  while (!pos(d)) {n*=fundunit; d*=fundunit;}
}

// return rational x,y s.t. this = x+y*sqrt(-d)
vector<RAT> RatQuad::coords(int rectangle) const
{
  Quad a = n*d.conj();
  INT b = d.norm();
  RAT x(a.r, b), y(a.i, b);
  if (rectangle && Quad::t)
    {
      y /= 2;
      x += y;
    }
  return {x,y};
}

// True iff x in (-1/2,1/2] and y in (-1/2,1/2] or (-1/4,1/4]
int RatQuad::in_rectangle() const
{
  vector<RAT> xy = coords(1); // so this = x+y*sqrt(-d)
  RAT x=xy[0], y=xy[1], half(1,2);
  if (Quad::t) y *=2;
  return -half<x && x<=half && -half<y && y<=half;
}

// x in [0,1/2] and y in [0,1/2] or [0,1/4]
int RatQuad::in_quarter_rectangle() const
{
  vector<RAT> xy = coords(1); // so this = x+y*sqrt(-d)
  RAT x=xy[0], y=xy[1], half(1,2);
  if (Quad::t) y *=2;
  return 0<=x && x<=half && 0<=y && y<=half;
}

// subtract Quad to put into rectangle
RatQuad reduce_to_rectangle(const RatQuad& a, Quad& shift)
{
  vector<RAT> xy = a.coords(1);  // so a = x+y*sqrt(-d)
  RAT x = xy[0], y=xy[1];
  INT xshift, yshift;
  if (Quad::t) y *= 2;
  yshift = y.round();
  xshift = (Quad::t? ((2*x-yshift)/2).round() : x.round());
  shift = Quad(xshift, yshift);
  // cout << " a = " << a <<endl;
  // cout << " shift = "<<shift<<" with norm "<<shift.norm()<<endl;
  // cout << " a-shift = " << a-shift << " with rect coords "<<(a-shift).coords(1)<<endl;
  assert ((a-shift).in_rectangle());
  return a-shift;
}

// list of Quad(s) a s.t. N(z-a)<1; at most 1 if just_one, otherwise
// we return all

// NB if a1 and a2 work then N(a1-a2)<4; only the Euclidean fields
// have elements of norm 2 or 3; our main use of this is in the
// non-Euclidean case.
vector<Quad> nearest_quads(const RatQuad& z, int just_one)
{
  vector<Quad> ans;
  Quad r;
  RatQuad z0 = reduce_to_rectangle(z, r);
  assert (z-r==z0);
  if (z0.norm()<1)
    ans.push_back(r);
  else
    return ans; // empty
  if (just_one)
    return ans;
  if (!Quad::is_Euclidean) // only possible shifts are by +-1
    {
      if ((z0-Quad::one).norm()<1)
        {
          ans.push_back(r+Quad::one);
          return ans;
        }
      if ((z0+Quad::one).norm()<1)
        {
          ans.push_back(r-Quad::one);
          return ans;
        }
      return ans;
    }
  vector<Quad> shifts = {Quad::one}; // adjustments up to sign
  if (Quad::d < 7) shifts.push_back(Quad(1,1)); // norms 2, 3, 3 for d=1,2,3
  if (Quad::nunits == 2) shifts.push_back(Quad(0,1)); // norms 2, 2, 3 for d=2,7,11
  for ( const auto& s : shifts)
    {
      for ( const auto& u : quadunits)
        {
          Quad us = u*s;
          if ((z0-us).norm()<1)
            ans.push_back(us+r);
        }
    }
  return ans;
}

vector<Quad> nearest_quads_to_quotient(const Quad& a, const Quad& b, int just_one)
{
  vector<Quad> ans;
  Quad q = a/b; // rounded
  Quad r=a; r.subprod(q,b); // a=q*b+r
  INT bnorm=b.norm();
  if (r.norm()<bnorm)
    ans.push_back(q);
  else
    return ans; // empty
  if (just_one)
    return ans;
  if (!Quad::is_Euclidean) // only possible shifts are by +-1
    {
      if ((r-b).norm()<bnorm)
        {
          q +=1;
          ans.push_back(q);
          return ans;
        }
      if ((r+b).norm()<bnorm)
        {
          q -=1;
          ans.push_back(q);
          return ans;
        }
      return ans;
    }
  vector<Quad> shifts = {Quad::one}; // adjustments up to sign
  if (Quad::d < 7) shifts.push_back(Quad(1,1)); // norms 2, 3, 3 for d=1,2,3
  if (Quad::nunits == 2) shifts.push_back(Quad(0,1)); // norms 2, 2, 3 for d=2,7,11
  for ( const auto& s : shifts)
    {
      for ( const auto& u : quadunits)
        {
          Quad us = u*s;
          if ((r-b*us).norm()<bnorm)
            ans.push_back(us+q);
        }
    }
  return ans;
}

// Finding cusp in list, with or without translation

// Return index of c in clist, or -1 if not in list
int cusp_index(const RatQuad& c, const vector<RatQuad>& clist)
{
  auto ci = std::find(clist.begin(), clist.end(), c);
  if (ci==clist.end())
    return -1;
  return ci-clist.begin();
}

// Return index i of c mod O_K in clist, with a=c-clist[i], or -1 if not in list
int cusp_index_with_translation(const RatQuad& c, const vector<RatQuad>& clist, Quad& t)
{
  int i=0;
  for ( const auto& ci : clist)
    {
      RatQuad diff = c-ci;
      if (diff.is_integral())
        {
          diff.is_integral(t);
          assert (clist[i]+t==c);
          return i;
        }
      i++;
    }
  return -1;
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
  if ((c1-c2).is_integral()) return 1;
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
#ifdef DEBUG_CUSP_EQ
  cout<<"Testing equivalence of cusps "<<c1<<" and "<<c2;
  cout<<" (N="<<N<<")"<<endl;
#endif
  // Quick tests
  if (c1==c2) return 1;
  if ((c1-c2).is_integral()) return 1;

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
  Qideal I = cc1.ideal();
  assert (I==cc2.ideal());
  assert (N.is_coprime_to(I));
#ifdef DEBUG_CUSP_EQ
  cout<<" - adjusted representations: "<<cc1<<", "<<cc2<<" with equal ideals "<<I<<" coprime to N"<<endl;
#endif

  // Use the criterion of Cor.2 in "Manin symbols over number fields"

  // Form ab-matrices with first columns equal to the two cusp representations

  mat22 M1 = AB_matrix(cc1.n, cc1.d); // [a1,b1;a2,b2]
#ifdef DEBUG_CUSP_EQ
  cout<<" - A = "<<M1<<", det = "<<M1.det()<<endl;
#endif
  mat22 M2 = AB_matrix(cc2.n, cc2.d); // [a1',b1';a2',b2']
#ifdef DEBUG_CUSP_EQ
  cout<<" - B = "<<M2<<", det = "<<M2.det()<<endl;
#endif
  assert (M1.det()==M2.det());
  Quad a2db2 = M2.entry(1,0)*M1.entry(1,1), a2b2d = M1.entry(1,0)*M2.entry(1,1);
#ifdef DEBUG_CUSP_EQ
  cout<<" - a2'*b2 = "<<a2db2<<", a2*b2' = "<<a2b2d<<endl;
#endif

  Qideal M = (D1*D1+N)*I.norm(); // recall D1==D2
#ifdef DEBUG_CUSP_EQ
  cout<<" - M = N(I)*(D^2+N) = "<<M<<endl;
#endif

  // test whether there is a unit u such that
  // (1) a2'*b2 = u*a2*b2' (mod M)

  if (M.divides(a2db2-a2b2d))
    {
#ifdef DEBUG_CUSP_EQ
      cout<<" - testing unit 1 - equivalent"<<endl;
#endif
      return 1;
    }
  vector<Quad>& units = (plusflag? quadunits: squareunits);
  for (vector<Quad>::const_iterator ui = units.begin()+1; ui!= units.end(); ++ui)
    {
      Quad u = *ui;
#ifdef DEBUG_CUSP_EQ
      cout<<" - testing unit "<<u<<flush;
#endif
      if (M.divides(a2db2-u*a2b2d))
        {
#ifdef DEBUG_CUSP_EQ
          cout<<" - equivalent"<<endl;
#endif
          return 1;
        }
    }
#ifdef DEBUG_CUSP_EQ
  cout<<" - not equivalent"<<endl;
#endif
  return 0;
}

// test function
int cuspeq_conj(const RatQuad& c1, const RatQuad& c2, const Qideal& N, int plusflag)
{
  int t1 = cuspeq(c1, c2, N, plusflag);
  int t2 = cuspeq(c1.conj(), c2.conj(), N.conj(), plusflag);
  if (t1!=t2)
    {
      cout<<"Problem testing equivalence of cusps "<<c1<<" and "<<c2<<" modulo "<<N<<endl;
      cout<<" - direct test yields "<<t1<<" while conjugate test yields "<<t2<<endl;
      exit(1);
    }
  return t1;
}

// if type = t>=0 , return U{alpha[t],oo}
// if type = -s<0 , return U{sigmas[s],oo}

modsym::modsym(const mat22& U, int type) // U in SL2
  : a(U(type>=0? alphas[type]: sigmas[-type])),
    b(U.image_oo())
{;}
