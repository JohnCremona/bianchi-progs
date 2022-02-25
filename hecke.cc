// FILE HECKE.CC: Implemention of Hecke operators for class homspace

#include "hecke.h"

// Implementation of Hecke and Atkin-Lehner matrices

// These functions return either one mat22 (2x2 matrix) or a vector of
// them.  Usage will normally be via the associated functions HeckeOp
// (etc.)  which return an object of class matop which is just a
// wrapper around a list of mat22's, with an associated string, the
// operator name.

//////////////////////////////////////////////////////////////////////
//
// For use only over fields of class number 1, probably now redundant.
//
//////////////////////////////////////////////////////////////////////

// T(P) for P=(p) principal prime
vector<mat22> Hecke(const Quad& p)  // P=(p) principal prime
{
  vector<Quad> resmodp = residues(p);
  vector<mat22> mats(1+resmodp.size());
  for (vector<Quad>::const_iterator r=resmodp.begin(); r!=resmodp.end(); ++r)
    mats.push_back(mat22(Quad::one,*r,Quad::zero,p));
  mats.push_back(mat22(p,Quad::zero,Quad::zero,Quad::one));
  return mats;
}

// W(P) for P=(p) principal prime dividing n
mat22 AtkinLehner(const Quad& p, const Quad& n)
{
  Quad u,v,a,b;
  for (u=Quad::one, v=n; div(p,v); v/=p, u*=p) ;
  quadbezout(u,v,a,b);
  return mat22(u*a,-b,n,u);
}

// W(P) for P=(p) principal prime dividing N
mat22 AtkinLehner(const Quad& p, Qideal& N)
{
  Qideal P(p);
  return AtkinLehner(P, N);
}

//////////////////////////
//
// Atkin-Lehner operators
//
//////////////////////////

// For [M1] square and M1,M2 coprime; the level is N=M1*M2.
//
// This is a matrix representing the operator T(I,I)*W(M1) where
// I^2*M1 is principal

mat22 AtkinLehner(Qideal& M1, Qideal& M2)
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
  assert (!I.is_zero() && "[M1] should be square in  AtkinLehner(M1,M2)");
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

// W(P^e) at level N where P^e||N

mat22 AtkinLehnerP(Quadprime& P, const Qideal& N)
{
  Qideal M1(Quad::one), M2(N);
  while (P.divides(M2))
    {
      M1 *= P;
      M2 /= P;
    }
  return AtkinLehner(M1,M2);
}

//#define DEBUG_HECKE

// N(P)+1 matrices representing T(A,A)*T(P) where P does not divide N
// and the class [P] is square with A^2P principal, A coprime to N.

vector<mat22> Hecke(Quadprime& P, Qideal& N) // assume [P] square
{
#ifdef DEBUG_HECKE
  cout<<"In Hecke("<<P<<"), level "<<N<<endl;
#endif
  if (P.is_principal())
    return Hecke(P.gen());
  vector<mat22> mats;
  Qideal A = P.sqrt_coprime_to(N); // A^2*P is principal, with A coprime to N, or A=0
  assert (!A.is_zero() && "[P] should be square in Hecke(P)");
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
  long normP = I2long(P.norm());
  assert ((long)mats.size()==1+normP);
  return mats;
}

// N(P)^2+N(P)+1 matrices representing T(A,A)*T(P^2), where P does
// not divide N, with AP principal and A coprime to N.

vector<mat22> HeckeSq(Quadprime& P, Qideal& N)
{
#ifdef DEBUG_HECKE
  cout<<"In Hecke("<<P<<"^2), level "<<N<<endl;
#endif
  vector<mat22> mats;
  Quad g,h;
  Qideal A = P.equivalent_coprime_to(N,g,h,1); // A*P=(g) is principal, with A coprime to N
  Qideal P2 = P*P;
  Qideal AP2 = A*P2;
  assert (Qideal(g)==A*P);
  // Make an (AP^2,A)-matrix of level N:
  mat22 M = AP2.AB_matrix_of_level(A, N, h); //  sets h to a generator of (AP)^2, not needed

  // (1) M * lift(1:a) for a mod P^2                (N(P)^2 matrices)
  // (2) M * lift(a:1) for a mod P^2 non-invertible (N(P) matrices)
  // with the second factor a lift from P^1(P^2) to Gamma_0(N)
  vector<Quad> resmodp2 = P2.residues();
  Quad u, v;
  N.is_coprime_to(P2, u, v); // u+v=1, u in N, v in P2
  for(vector<Quad>::const_iterator r=resmodp2.begin(); r!=resmodp2.end(); ++r)
    {
      Quad a = *r;
      mats.push_back(M*lift_to_Gamma_0(N, P2, Quad::one, a, u, v));
      if (P.contains(a))
        mats.push_back(M*lift_to_Gamma_0(N, P2, a, Quad::one, u, v));
    }

  // (3) an (AP,AP) matrix of level N, i.e. just diag(g,g) (1 matrix, so N(P)^2_N(P)+1 in all)
  mats.push_back(mat22::scalar(g));

#ifdef DEBUG_HECKE
  cout<<" Hecke matrices are "<<mats<<endl;
#endif
  long normP = I2long(P.norm());
  assert ((long)mats.size()==1+normP*(1+normP));
  return mats;
}

// (N(P)+1)(N(Q)+1) matrices representing T(A,A)*T(PQ) where P,Q do
// not divide N and class [PQ] square with A^2PQ principal, A coprime
// to N.  When PQ is principal we take A=1.

vector<mat22> HeckePQ(Quadprime& P, Quadprime& Q, Qideal& N)
{
#ifdef DEBUG_HECKE
  cout<<"In Hecke("<<P<<"*"<<Q<<"), level "<<N<<endl;
#endif
  vector<mat22> mats;
  Quad g;
  Qideal PQ = P*Q;
  Qideal A = PQ.sqrt_coprime_to(N); // APQ is principal, with A coprime to N, or A=0
  assert (!A.is_zero() && "[PQ] should be square in HeckePQ()");
  Qideal AP=A*P, AQ=A*Q, APQ = A*PQ;
  mat22 M = APQ.AB_matrix_of_level(A, N, g); // sets g to a generator of A^2PQ,
                                             //raising an error if not principal

  // (1) M*lift(1:a) for a mod PQ                (N(P)N(Q) matrices)
  // (2) M*lift(a:1) for a mod PQ not invertible (N(P)+N(Q)-1 matrices)
  vector<Quad> resmodpq = PQ.residues();
  Quad u, v;
  N.is_coprime_to(PQ, u, v); // u+v=1, u in N, v in PQ
  for(vector<Quad>::const_iterator r=resmodpq.begin(); r!=resmodpq.end(); ++r)
    {
      Quad a = *r;
      mats.push_back(M*lift_to_Gamma_0(N, PQ, Quad::one, a, u, v));
      if (P.contains(a) or Q.contains(a)) // or both
        mats.push_back(M*lift_to_Gamma_0(N, PQ, a, Quad::one, u, v));
    }

  // (3) (AP,AQ) and (AQ,AP) matrices of level N (2 matrices, so (N(P)+1)(N(Q)+1) in all)
  mats.push_back(AP.AB_matrix_of_level(AQ, N, g));
  mats.push_back(AQ.AB_matrix_of_level(AP, N, g));

#ifdef DEBUG_HECKE
  cout<<" Hecke matrices are "<<mats<<endl;
#endif
  long normP = I2long(P.norm()), normQ = I2long(Q.norm());
  assert ((long)mats.size()==(1+normP)*(1+normQ));
  return mats;
}

// T(P)W(M1) for P*M1 principal, P not dividing N=M1*M2, (M1,M2)=1.

// Later we'll implement a more general version giving T(A,A)T(P)W(M1) when [P*M1] is square

//#define DEBUG_HECKE
vector<mat22> HeckeAL(Quadprime& P, Qideal& M1, Qideal& M2)
{
#ifdef DEBUG_HECKE
  cout<<"In HeckeAL() with P="<<P<<", M1="<<ideal_label(M1)<<", M2="<<ideal_label(M2)<<endl;
#endif
  vector<mat22> mats;
  Qideal PM1 = P*M1;
  Quad g, h, h1, x, u, v;
  int i = PM1.is_principal(g);
  assert (i && "P*M1 is principal");
  i = P.is_coprime_to(M1, u, v); // u+v=1, u in P, v in M1
  assert (i && "P and M1 are coprime");
  Qideal M3 = M2.equivalent_coprime_to(PM1, h, x, 1); // M2*M3=(h)
  assert (Qideal(h) == M2*M3);
  // i = PM1.is_coprime_to(h, h1); // h*h1=1 mod PM1
  // assert (i && "M2*M3 is coprime to P*M1");
  // h *= h1;
  // M3 *= h1;
  Qideal M2M3 = M2*M3;
  assert (Qideal(h) == M2M3);

  Quad a, b, c, d, r, s;
  i = PM1.is_coprime_to(M2M3, r, s); // r+s=1, r in PM1, s in M2M3
  assert (i);

  // First handle (1:0) mod P*M1, finding a lift [a,b;c,d] with c in M2M3 so h|c
  mat22 m = lift_to_Gamma_0(M2M3, PM1, 1, 0, s, r);
#ifdef DEBUG_HECKE
  cout<<" Lift of (1:0) mod "<<PM1<<" to Gamma_0("<<M2M3<<") is "<<m<<endl;
#endif
  a = m.entry(0,0);
  b = m.entry(0,1)*h;
  c = m.entry(1,0)/h;
  d = m.entry(1,1);
  m = mat22(d, -c, -b*g, g*a);
  assert (M1.contains(d));
  assert (M1.contains(g));
  assert ((M1*M2).contains(b*g));
  assert (m.det()==g);
  mats.push_back(m);
  vector<Quad> resmodp = P.residues();
  for(vector<Quad>::const_iterator x=resmodp.begin(); x!=resmodp.end(); ++x)
    {
      c = v*(*x)+u;
      d = v;
      m = lift_to_Gamma_0(M2M3, PM1, c, d, s, r);
      a = m.entry(0,0);
      b = m.entry(0,1)*h;
      c = m.entry(1,0)/h;
      d = m.entry(1,1);
      m = mat22(d, -c, -b*g, g*a);
      assert (M1.contains(d) && M1.contains(g));
      assert ((M1*M2).contains(b*g));
      assert (m.det()==g);
      mats.push_back(m);
    }
#ifdef DEBUG_HECKE
  cout<<" Hecke matrices are "<<mats<<endl;
#endif
  return mats;
}

// T(P)W(Q^e) for P*Q^e principal, at level N where Q^e||N and P does not divide N

// Later we'll implement a more general version giving T(A,A)T(P)W(Q^e) when [P*Q^e] is square

vector<mat22> HeckeALP(Quadprime& P, Quadprime& Q, Qideal& N)
{
#ifdef DEBUG_HECKE
  cout<<"In HeckeALP() with P="<<P<<", Q="<<Q<<", N="<<ideal_label(N)<<endl;
#endif
  Qideal M1(Quad::one), M2(N);
  while (Q.divides(M2))
    {
      M1 *= Q;
      M2 /= Q;
    }
  return HeckeAL(P,M1,M2);
}


// Matrix inducing T(A,A) at level N, when A^2 is principal and A+N=1

mat22 Char(Qideal& A, const Qideal& N)
{
  Quad g;
  return A.AB_matrix_of_level(A, N, g);
}

