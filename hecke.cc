// FILE HECKE.CC: Implemention of Hecke operators for class homspace

#include "hecke.h"
#include "P1N.h"

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
vector<mat22> HeckeP(const Quad& p)  // P=(p) principal prime
{
  vector<Quad> resmodp = residues(p);
  vector<mat22> mats; mats.reserve(1+resmodp.size());
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

// Level N=M1*M2 and M1,M2 coprime and M1 principal, a matrix
// representing W(M1,M2).

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

  // we require M1 to be principal

  Quad g, a, c, x, y;
  int i = M1.is_principal(g);
  assert (i && "M1 should be principal in  AtkinLehner(M1,M2)");
  Qideal N = M1*M2;
  Qideal C = N.equivalent_coprime_to(N, c, x, 1); // CN=(c)
  Qideal M2C = M2*C;
  Qideal A = M1.equivalent_coprime_to(M2C, a, x, 1); // AM1=(a)
  Qideal M1A = M1*A;
  i = M1A.is_coprime_to(M2C, x, y);
  assert(i && "AI is coprime to M2C");
  Quad d = g*x/a;
  assert (a*d==g*x);
  Quad b = -g*y/c;
  assert (b*c==-g*y);
  assert (M1.contains(a));
  assert (N.contains(c));
  assert (M1.contains(d));
  assert (a*d-b*c==g);
  return mat22(a,b,c,d);
}

// Level N=M1*M2 and M1,M2 coprime and [M1] square with A^2*M1
// principal: a matrix representing T(A,A)W(M1)

mat22 AtkinLehner_Chi(Qideal& M1, Qideal& M2, Qideal& A)
{
  Quad g, a, c, x, y;
  Qideal AM1 = A*M1;
  Qideal AsqM1 = A*AM1;
  int i = AsqM1.is_principal(g);
  assert (i && "A^2*M1 must be principal in AtkinLehner_Chi(M1,M2,A)");
  Qideal AN = AM1*M2;
  Qideal C = AN.equivalent_coprime_to(AN, c, x, 1); // CAN=(c)
  Qideal M2C = M2*C;
  Qideal B = AM1.equivalent_coprime_to(M2C, a, x, 1); // BAM1=(a)
  Qideal M1B = M1*B;
  i = M1B.is_coprime_to(M2C, x, y);
  assert(i && "BA is coprime to M2C");
  Quad d = g*x/a;
  assert (a*d==g*x);
  Quad b = -g*y/c;
  assert (b*c==-g*y);
  assert (AN.contains(c));
  assert (AM1.contains(a));
  assert (AM1.contains(d));
  assert (A.contains(b));
  assert (a*d-b*c==g);
  return mat22(a,b,c,d);
}

// Level N, Q|N prime, Q^e||N where Q^e must be principal: then this
// is a matrix representing W(Q^e).

mat22 AtkinLehnerQ(Quadprime& Q, const Qideal& N)
{
  Qideal M1(Quad::one), M2(N);
  while (Q.divides(M2))
    {
      M1 *= Q;
      M2 /= Q;
    }
  return AtkinLehner(M1,M2);
}

// Level N, Q|N prime, Q^e||N where [Q^e] is square with A^2*Q^e
// principal, A coprime to N: this is a matrix representing
// T(A,A)W(Q^e).

mat22 AtkinLehnerQ_Chi(Quadprime& Q, Qideal& A, const Qideal& N)
{
  Qideal M1(Quad::one), M2(N);
  while (Q.divides(M2))
    {
      M1 *= Q;
      M2 /= Q;
    }
  return AtkinLehner_Chi(M1,M2, A);
}

//#define DEBUG_HECKE

// Level N, P prime not dividing N.

// N(P)+1 matrices representing T(P) when P is principal. These do not depend on N.

vector<mat22> HeckeP(Quadprime& P)
{
#ifdef DEBUG_HECKE
  cout<<"In HeckeP("<<P<<"), level "<<N<<endl;
#endif
  Quad g;
  int i = P.is_principal(g);
  assert (i && "HeckeP(P) called with non-principal P");
  return HeckeP(g);
}

// N(P)+1 matrices representing T(A,A)*T(P) when the class [P] is
// square with A^2*P principal, A coprime to N.

vector<mat22> HeckeP_Chi(Quadprime& P, Qideal& A, Qideal& N)
{
#ifdef DEBUG_HECKE
  cout<<"In HeckeP_Chi("<<P<<","<<A<<"), level "<<N<<endl;
#endif
  vector<mat22> mats;
  Qideal AP = A*P;
  int i = ((A*AP).is_principal());
  assert(i && "A^2*P must be principal in HeckeP_chiA(A,P)");
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
  long normP = I2long(P.norm());
  mats.reserve(1+normP);
  for(vector<Quad>::const_iterator r=resmodp.begin(); r!=resmodp.end(); ++r)
    {
      Quad a = *r;
      mats.push_back(M*mat22(Quad::one,a,nu,Quad::one+a*nu));
    }
  mats.push_back(M);
#ifdef DEBUG_HECKE
  cout<<" Hecke matrices are "<<mats<<endl;
#endif
  assert ((long)mats.size()==1+normP);
  return mats;
}

// Level N, P prime not dividing N.

// Returns N(P)^2+N(P)+1 matrices representing T(P^2), when P^2 is
// principal.

vector<mat22> HeckeP2(Quadprime& P, Qideal& N)
{
#ifdef DEBUG_HECKE
  cout<<"In HeckePSq(P,N) with P="<<P<<", N="<<N<<endl;
#endif
  Qideal P2 = P*P;
  Quad g, u,v,a;
  int i = P2.is_principal(g);
  assert (i && "P^2 must be principal in HeckePSq(P,N)");
  N.is_coprime_to(P2, u, v); // u+v=1, u in N, v in P2

  mat22 M1 = mat22(g,0,0,1);                // a (P^2,O)-matrix of level N
  mat22 M2 = P.AB_matrix_of_level(P, N, g); // a (P,P)-matrix of level N

  vector<mat22> mats;
  vector<Quad> resmodp2 = P2.residues();

  // (1) M1 * lift(1:a) for a mod P^2                (N(P)^2 matrices)
  // (2) M1 * lift(a:1) for a mod P^2 non-invertible (N(P) matrices)
  // with the second factor a lift from P^1(O/P^2) to Gamma_0(N)
  for(vector<Quad>::const_iterator r=resmodp2.begin(); r!=resmodp2.end(); ++r)
    {
      a = *r;
      mats.push_back(M1*lift_to_Gamma_0(N, P2, Quad::one, a, u, v));
      if (P.contains(a))
        mats.push_back(M1*lift_to_Gamma_0(N, P2, a, Quad::one, u, v));
    }

  // (3) M2, a (P,P) matrix of level N  (1 matrix, so N(P)^2_N(P)+1 in all):
  mats.push_back(M2);
#ifdef DEBUG_HECKE
  cout<<" Hecke matrices are "<<mats<<endl;
#endif
  long normP = I2long(P.norm());
  assert ((long)mats.size()==1+normP*(1+normP));
  return mats;
}

// Returns N(P)^2+N(P)+1 matrices representing T(A,A)T(P^2) with
// (AP)^2 principal and A coprime to N.

vector<mat22> HeckeP2_Chi(Quadprime& P, Qideal& A, Qideal& N)
{
#ifdef DEBUG_HECKE
  cout<<"In HeckePSq_Chi(P,N) with P="<<P<<", A="<<A<<", N="<<N<<endl;
#endif
  Qideal P2 = P*P;
  Qideal AP = A*P, AP2 = A*P2;
  Qideal A2P2 = A*AP2;
  Quad g, u, v, a;
  int i = A2P2.is_principal(g);
  assert (i && "(AP)^2 must be principal in HeckePSq_Chi(P,A,N)");
  N.is_coprime_to(P2, u, v); // u+v=1, u in N, v in P2

  vector<mat22> mats;
  vector<Quad> resmodp2 = P2.residues();
  mat22 M1 = AP2.AB_matrix_of_level(A, N, g); // An (AP^2,A)-matrix of level N
  mat22 M2 = AP.AB_matrix_of_level(AP, N, g); // An (AP,AP)-matrix of level N

  // (1) M1 * lift(1:a) for a mod P^2                (N(P)^2 matrices)
  // (2) M1 * lift(a:1) for a mod P^2 non-invertible (N(P) matrices)
  // with the second factor a lift from P^1(P^2) to Gamma_0(N)
  for(vector<Quad>::const_iterator r=resmodp2.begin(); r!=resmodp2.end(); ++r)
    {
      a = *r;
      mats.push_back(M1*lift_to_Gamma_0(N, P2, Quad::one, a, u, v));
      if (P.contains(a))
        mats.push_back(M1*lift_to_Gamma_0(N, P2, a, Quad::one, u, v));
    }

  // (3) M2, an (AP,AP) matrix of level N  (1 matrix, so N(P)^2_N(P)+1 in all):
  mats.push_back(M2);

#ifdef DEBUG_HECKE
  cout<<" Hecke matrices are "<<mats<<endl;
#endif
  long normP = I2long(P.norm());
  assert ((long)mats.size()==1+normP*(1+normP));
  return mats;
}

// Level N, P,Q distinct primes not dividing N.

// Returns (N(P)+1)(N(Q)+1) matrices representing T(PQ) when P*Q is
// principal.

vector<mat22> HeckePQ(Quadprime& P, Quadprime& Q, Qideal& N)
{
#ifdef DEBUG_HECKE
  cout<<"In HeckePQ(P,Q,N) with P="<<P<<", Q="<<Q<<", N="<<N<<endl;
#endif
  Quad g, u, v, a;
  Qideal PQ = P*Q;
  int i = PQ.is_principal(g);
  assert (i && "PQ must be principal in HeckePQ(P,Q,N)");
  vector<Quad> resmodpq = PQ.residues();
  vector<mat22> mats;
  mat22 M = mat22(g,0,0,1);
  N.is_coprime_to(PQ, u, v); // u+v=1, u in N, v in PQ

  // (1) M*lift(1:a) for a mod PQ                (N(P)N(Q) matrices)
  // (2) M*lift(a:1) for a mod PQ not invertible (N(P)+N(Q)-1 matrices)
  for(vector<Quad>::const_iterator r=resmodpq.begin(); r!=resmodpq.end(); ++r)
    {
      a = *r;
      mats.push_back(M*lift_to_Gamma_0(N, PQ, Quad::one, a, u, v));
      if (P.contains(a) or Q.contains(a)) // or both
        mats.push_back(M*lift_to_Gamma_0(N, PQ, a, Quad::one, u, v));
    }

  // (3) (P,Q) and (Q,P) matrices of level N (2 matrices, so (N(P)+1)(N(Q)+1) in all)
  mats.push_back(P.AB_matrix_of_level(Q, N, g));
  mats.push_back(Q.AB_matrix_of_level(P, N, g));

#ifdef DEBUG_HECKE
  cout<<" Hecke matrices are "<<mats<<endl;
#endif
  long normP = I2long(P.norm()), normQ = I2long(Q.norm());
  assert ((long)mats.size()==(1+normP)*(1+normQ));
  return mats;
}

// Level N, B square-free, principal and coprime to N

// Returns \psi(B) = \prod_{P|B}(N(P)+1) matrices representing T(B).

vector<mat22> HeckeB(Qideal& B, Qideal& N)
{
#ifdef DEBUG_HECKE
  cout<<"In HeckeB(B,N) with B="<<B<<" = "<<B.factorization()<", N="<<N<<endl;
#endif
  Quad g, u, v, c, d;
  int i = B.is_principal(g);
  assert (i && "B must be principal in HeckePQ(B,N)");
  mat22 M = mat22(g,0,0,1);
  N.is_coprime_to(B, u, v); // u+v=1, u in N, v in B

  // The matrices are
  // M*lift(c:d) for all (c:d) in P1(O/B), lifted to Gamma_0(N)    (psi(B) matrices)
  P1N P1B(B);
  vector<mat22> mats;
  mats.reserve(P1B.size());
  for(i=0; i<P1B.size(); i++)
    {
      P1B.make_symb(i,c,d);
      mats.push_back(M*lift_to_Gamma_0(N, B, c, d, u, v));
    }

#ifdef DEBUG_HECKE
  cout<<" Hecke matrices are "<<mats<<endl;
#endif
  return mats;
}

// Returns (N(P)+1)(N(Q)+1) matrices representing T(A,A)T(PQ) when
// [P*Q] is square, with A^2PQ principal, A coprime to N.

vector<mat22> HeckePQ_Chi(Quadprime& P, Quadprime& Q, Qideal&A, Qideal& N)
{
#ifdef DEBUG_HECKE
  cout<<"In HeckePQ_Chi(P,Q,A,N) with P="<<P<<", Q="<<Q<<", A="<<A<<", N="<<N<<endl;
#endif
  Quad g, u, v, a;
  Qideal PQ = P*Q, AP=A*P, AQ=A*Q;
  Qideal APQ = A*PQ, A2PQ=AP*AQ;
  int i = A2PQ.is_principal(g);
  assert (i && "A^2PQ must be principal in HeckePQ(P,Q,A,N)");
  vector<mat22> mats;
  mat22 M = APQ.AB_matrix_of_level(A, N, g); // sets g to a generator of A^2PQ,
                                             //raising an error if not principal
  vector<Quad> resmodpq = PQ.residues();
  N.is_coprime_to(PQ, u, v); // u+v=1, u in N, v in PQ

  // (1) M*lift(1:a) for a mod PQ                (N(P)N(Q) matrices)
  // (2) M*lift(a:1) for a mod PQ not invertible (N(P)+N(Q)-1 matrices)
  for(vector<Quad>::const_iterator r=resmodpq.begin(); r!=resmodpq.end(); ++r)
    {
      a = *r;
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

// Returns psi(B) matrices representing T(A,A)T(B) when
// [B] is square, with A^2*B principal, A coprime to N.

vector<mat22> HeckeB_Chi(Qideal& B, Qideal&A, Qideal& N)
{
#ifdef DEBUG_HECKE
  cout<<"In HeckeB_Chi(B,A,N) with B="<<B<<" = "<<B.factorization()<<", N="<<N<<endl;
#endif
  Quad g, u, v, c,d;
  Qideal AB = A*B;
  Qideal A2B = A*AB;
  int i = A2B.is_principal(g);
  assert (i && "A^2B must be principal in HeckePQ(B,A,N)");
  N.is_coprime_to(B, u, v); // u+v=1, u in N, v in B
  mat22 M = AB.AB_matrix_of_level(A, N, g); // sets g to a generator of A^2B,
                                             //raising an error if not principal

  // The matrices are
  // M*lift(c:d) for all (c:d) in P1(O/B), lifted to Gamma_0(N)    (psi(B) matrices)
  P1N P1B(B);
  vector<mat22> mats;
  mats.reserve(P1B.size());
  for(i=0; i<P1B.size(); i++)
    {
      P1B.make_symb(i,c,d);
      mats.push_back(M*lift_to_Gamma_0(N, B, c, d, u, v));
    }

#ifdef DEBUG_HECKE
  cout<<" Hecke matrices are "<<mats<<endl;
#endif
  assert ((long)mats.size()==P1B.size());
  return mats;
}

// Level N=M1*M2, P prime not dividing N, M1,M2 coprime

// Returns N(P)+1 matrices for T(P)W(M1) if P*M1 is principal

// Later we'll implement a more general version giving
// T(A,A)T(P)W(M1) when [P*M1] is square

//#define DEBUG_HECKE
vector<mat22> HeckePAL(Quadprime& P, Qideal& M1, Qideal& M2)
{
#ifdef DEBUG_HECKE
  cout<<"In HeckePAL(P,M1,M1) with P="<<P<<", M1="<<ideal_label(M1)<<", M2="<<ideal_label(M2)<<endl;
#endif
  vector<mat22> mats;
  Qideal PM1 = P*M1;
  Quad g, h, h1, x, u, v;
  int i = PM1.is_principal(g);
  assert (i && "P*M1 must be principal in HeckePAL(P,M1,M2");
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
#ifdef DEBUG_HECKE
      cout<<"lifting (c:d)=("<<c<<":"<<d<<") from "<<PM1<<" to Gamma_0("<<M2M3<<")"<<endl;
#endif
      m = lift_to_Gamma_0(M2M3, PM1, c, d, s, r);
#ifdef DEBUG_HECKE
      cout<<" --> "<<m<<endl;
#endif
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

// Level N, P prime not dividing N, Q^e||N

// Returns matrices for T(P)W(Q^e) for P*Q^e principal

// Later we'll implement a more general version giving T(A,A)T(P)W(Q^e) when [P*Q^e] is square

vector<mat22> HeckePALQ(Quadprime& P, Quadprime& Q, Qideal& N)
{
#ifdef DEBUG_HECKE
  cout<<"In HeckePALQ() with P="<<P<<", Q="<<Q<<", N="<<ideal_label(N)<<endl;
#endif
  Qideal M1(Quad::one), M2(N);
  while (Q.divides(M2))
    {
      M1 *= Q;
      M2 /= Q;
    }
  return HeckePAL(P,M1,M2);
}


// Matrix inducing T(A,A) at level N, when A^2 is principal and A+N=1

mat22 Char(Qideal& A, const Qideal& N)
{
  Quad g;
  return A.AB_matrix_of_level(A, N, g);
}

