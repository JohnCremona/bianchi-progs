// FILE MAT22.CC

#include "primes.h"
#include "mat22.h"

// Implementation of Hecke and Atkin-Lehner matrices

vector<mat22> Hecke(const Quad& p)  // P=(p) principal prime
{
  vector<Quad> resmodp = residues(p);
  vector<mat22> mats(1+resmodp.size());
  for (vector<Quad>::const_iterator r=resmodp.begin(); r!=resmodp.end(); ++r)
    mats.push_back(mat22(Quad::one,*r,Quad::zero,p));
  mats.push_back(mat22(p,Quad::zero,Quad::zero,Quad::one));
  return mats;
}

mat22 AtkinLehner(const Quad& p, const Quad& n) // P=(p) principal prime
{
  Quad u,v,a,b;
  for (u=Quad::one, v=n; div(p,v); v/=p, u*=p) ;
  quadbezout(u,v,a,b);
  return mat22(u*a,-b,n,u);
}

mat22 AtkinLehner(const Quad& p, Qideal& N) // P=(p) principal prime dividing N
{
  Qideal P(p);
  return AtkinLehner(P, N);
}


mat22 AtkinLehner(Qideal& M1, Qideal& M2) // assume [M1] square and M1,M2 coprime
{
  if ((M1.is_principal() && M2.is_principal()))
    {
      Quad u = M1.gen(), v = M2.gen(), a,b;
      quadbezout(u,v,a,b);
      mat22 W(u*a,-b,u*v,u);
      assert (W.det()==u);
      return W;
    }

  Qideal I = M1.sqrt_coprime_to(M2);
  if (I.is_zero())
    {
      cerr<<"Cannot compute AtkinLehner("<<ideal_label(M1)<<") as its class is not a square"<<endl;
      return mat22();
    }
  Quad g, a, c, x, y;
  Qideal IM1 = I*M1;
  Qideal IsqM1 = I*IM1;
  int i = IsqM1.is_principal(g);
  assert (i && "I^2*M1 is principal");
  Qideal IN = IM1*M2;
  Qideal C = IN.equivalent_coprime_to(IN, c, x, 1); // CIN=(c)
  Qideal M2C = M2*C;
  Qideal A = IM1.equivalent_coprime_to(M2C, a, x, 1); // AIM1=(a)
  Qideal M1A = M1*A;
  i = M1A.is_coprime_to(M2C, x, y);
  assert(i && "AI is coprime to M2C");
  Quad d = g*x/a;
  assert (a*d==g*x);
  Quad b = -g*y/c;
  assert (b*c==-g*y);
  assert (IN.contains(c));
  assert (IM1.contains(a));
  assert (IM1.contains(d));
  assert (I.contains(b));
  assert (a*d-b*c==g);
  return mat22(a,b,c,d);
}

mat22 AtkinLehnerP(Quadprime& P, const Qideal& N) // =AL(P^e,N/P^e) where P^e||N
{
  Qideal M1(Quad::one), M2(N);
  while (P.divides(M2))
    {
      M1 *= P;
      M2 /= P;
    }
  return AtkinLehner(M1,M2);
}

vector<mat22> Hecke(const Quad& p, Qideal& N)
{
  vector<Quad> resmodp = residues(p);
  vector<Quad>::const_iterator r=resmodp.begin();
  vector<mat22> mats;
  mats.reserve(I2long(quadnorm(p)));
  while(r!=resmodp.end())
    mats.push_back(mat22(Quad::one,*r++,Quad::zero,p));
  mats.push_back(mat22(p,Quad::zero,Quad::zero,Quad::one));
  return mats;
}

//#define DEBUG_HECKE

vector<mat22> Hecke(Quadprime& P, Qideal& N) // assume [P] square
{
#ifdef DEBUG_HECKE
  cout<<"In Hecke("<<P<<"), level "<<N<<endl;
#endif
  if (P.is_principal())
    return Hecke(P.gen(), N);
  vector<mat22> mats;
  Qideal A = P.sqrt_coprime_to(N); // A^2*P is principal, with A coprime to N, or A=0
  if (A.is_zero())
    {
      cerr<<"Cannot compute Hecke("<<ideal_label(P)<<") as its class is not a square"<<endl;
      return mats;
    }
  Qideal AP = A*P;
  int i = ((A*AP).is_principal());
  assert(i && "A^2*P is principal");
  Quad g;
  mat22 M = AP.AB_matrix_of_level(A, N, g);
#ifdef DEBUG_HECKE
  cout<<" A = "<<A<<", A*A*P= "<<A*AP<<endl;
  cout<<" base M = "<<M<<endl;
#endif

  vector<Quad> Ngens = N.gens();
  Quad nu = Ngens[0];
  if (P.divides(nu))
    nu = Ngens[1];
  assert (!P.divides(nu));

  vector<Quad> resmodp = P.residues();
  mats.reserve(I2long(P.norm()));
  for(vector<Quad>::const_iterator r=resmodp.begin(); r!=resmodp.end(); ++r)
    {
      Quad a = *r;
      mats.push_back(M*mat22(Quad::one,a,nu,Quad::one+a*nu));
    }
  mats.push_back(M);
#ifdef DEBUG_HECKE
  cout<<" Hecke matrices are "<<mats<<endl;
#endif
  return mats;
}

RatQuad mat22::operator()(const RatQuad& q)const
{
  Quad r = q.num(), s = q.den();
  apply_left(r, s);
  return RatQuad(r,s, 1);
}

RatQuad mat22::image_oo() const {return RatQuad(a,c);}
RatQuad mat22::preimage_oo() const {return RatQuad(-d,c);}
RatQuad mat22::image_0() const {return RatQuad(b,d);}
RatQuad mat22::preimage_0() const {return RatQuad(b,-a);}
