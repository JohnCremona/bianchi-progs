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
  vector<mat22> mats(resmodp.size()); // will add 1 at end
  Quad zero(0), one(1);
  std::transform(resmodp.begin(), resmodp.end(), mats.begin(),
                 [p, zero, one] (const Quad& r) {return mat22(one,r,zero,p);});
  mats.push_back(mat22(p,zero,zero,one));
  return mats;
}

// W(P) for P=(p) principal prime dividing n
mat22 AtkinLehner(const Quad& p, const Quad& n)
{
  Quad u,v,a,b;
  for (u=Quad::one, v=n; div(p,v,v); u*=p) ;
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

mat22 AtkinLehner_Chi(const Qideal& M1, const Qideal& M2, const Qideal& A)
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

mat22 AtkinLehnerQ(const Quadprime& Q, const Qideal& N)
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

mat22 AtkinLehnerQ_Chi(const Quadprime& Q, const Qideal& A, const Qideal& N)
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
  cout<<"In HeckeP("<<P<<")"<<endl;
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
  Quad nu = Ngens[0], one(1);
  if (P.divides(nu))
    nu = Ngens[1];
  assert (!P.divides(nu));

  vector<Quad> resmodp = P.residues();
  long normP = I2long(P.norm());
  vector<mat22> mats(normP); // will add 1 at end
  std::transform(resmodp.begin(), resmodp.end(), mats.begin(),
                 [nu, M, one] ( const Quad& a) {return M*mat22(one,a,nu,one+a*nu);});
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
  Qideal NP2 = N*P2;
  Quad g, u,v, one(1);
  int i = P2.is_principal(g);
  assert (i && "P^2 must be principal in HeckePSq(P,N)");
  N.is_coprime_to(P2, u, v); // u+v=1, u in N, v in P2

  mat22 M1 = mat22::diag(g,Quad(1));                // a (P^2,O)-matrix of level N
  mat22 M2 = P.AB_matrix_of_level(P, N, g); // a (P,P)-matrix of level N

  vector<mat22> mats;
  vector<Quad> resmodp2 = P2.residues();

  // (1) M1 * lift(1:a) for a mod P^2                (N(P)^2 matrices)
  // (2) M1 * lift(a:1) for a mod P^2 non-invertible (N(P) matrices)
  // with the second factor a lift from P^1(O/P^2) to Gamma_0(N)
  for( const auto& a : resmodp2)
    {
      mats.push_back(M1*lift_to_Gamma_0(NP2, one, a, u, v));
      if (P.contains(a))
        mats.push_back(M1*lift_to_Gamma_0(NP2, a, one, u, v));
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
  Qideal NP2 = N*P2;
  Qideal AP = A*P, AP2 = A*P2;
  Qideal A2P2 = A*AP2;
  Quad g, u, v, one(1);
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
  for( const auto& a : resmodp2)
    {
      mats.push_back(M1*lift_to_Gamma_0(NP2, one, a, u, v));
      if (P.contains(a))
        mats.push_back(M1*lift_to_Gamma_0(NP2, a, one, u, v));
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
  Quad g, u, v, one(1);
  Qideal PQ = P*Q;
  Qideal NPQ = N*PQ;
  int i = PQ.is_principal(g);
  assert (i && "PQ must be principal in HeckePQ(P,Q,N)");
  vector<Quad> resmodpq = PQ.residues();
  vector<mat22> mats;
  mat22 M = mat22::diag(g,Quad(1));
  N.is_coprime_to(PQ, u, v); // u+v=1, u in N, v in PQ

  // (1) M*lift(1:a) for a mod PQ                (N(P)N(Q) matrices)
  // (2) M*lift(a:1) for a mod PQ not invertible (N(P)+N(Q)-1 matrices)
  for( const auto& a : resmodpq)
    {
      mats.push_back(M*lift_to_Gamma_0(NPQ, one, a, u, v));
      if (P.contains(a) or Q.contains(a)) // or both
        mats.push_back(M*lift_to_Gamma_0(NPQ, a, one, u, v));
    }

  // (3) (P,Q) and (Q,P) matrices of level N (2 matrices, so (N(P)+1)(N(Q)+1) in all)
  mats.push_back(P.AB_matrix_of_level(Q, N, g));
  mats.push_back(Q.AB_matrix_of_level(P, N, g));

#ifdef DEBUG_HECKE
  cout<< mats.size() << " Hecke matrices" << endl;
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
  cout<<"In HeckeB(B,N) with B="<<B<<" = "<<B.factorization()<<", N="<<N<<endl;
#endif
  Quad g, u, v, c, d;
  int i = B.is_principal(g);
  assert (i && "B must be principal in HeckePQ(B,N)");
  mat22 M = mat22::diag(g,Quad(1));
  N.is_coprime_to(B, u, v); // u+v=1, u in N, v in B

  // The matrices are
  // M*lift(c:d) for all (c:d) in P1(O/B), lifted to Gamma_0(N)    (psi(B) matrices)
  Qideal NB = N*B;
  P1N P1B(B);
  vector<mat22> mats;
  mats.reserve(P1B.size());
  for(i=0; i<P1B.size(); i++)
    {
      P1B.make_symb(i,c,d);
      mats.push_back(M*lift_to_Gamma_0(NB, c, d, u, v));
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
  Quad g, u, v, one(1);
  Qideal PQ = P*Q, AP=A*P, AQ=A*Q;
  Qideal APQ = A*PQ, A2PQ=AP*AQ;
  Qideal NPQ = N*PQ;
  int i = A2PQ.is_principal(g);
  assert (i && "A^2PQ must be principal in HeckePQ(P,Q,A,N)");
  vector<mat22> mats;
  mat22 M = APQ.AB_matrix_of_level(A, N, g); // sets g to a generator of A^2PQ,
                                             //raising an error if not principal
  vector<Quad> resmodpq = PQ.residues();
  N.is_coprime_to(PQ, u, v); // u+v=1, u in N, v in PQ

  // (1) M*lift(1:a) for a mod PQ                (N(P)N(Q) matrices)
  // (2) M*lift(a:1) for a mod PQ not invertible (N(P)+N(Q)-1 matrices)
  for( const auto& a : resmodpq)
    {
      mats.push_back(M*lift_to_Gamma_0(NPQ, one, a, u, v));
      if (P.contains(a) or Q.contains(a)) // or both
        mats.push_back(M*lift_to_Gamma_0(NPQ, a, one, u, v));
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
// [B] is squarefree, with A^2*B principal, A coprime to N.

vector<mat22> HeckeB_Chi(Qideal& B, Qideal&A, Qideal& N)
{
#ifdef DEBUG_HECKE
  cout<<"In HeckeB_Chi(B,A,N) with B="<<B<<" = "<<B.factorization()<<", N="<<N<<endl;
#endif
  Quad g, u, v, c,d;
  Qideal AB = A*B;
  Qideal NB = N*B;
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
      mats.push_back(M*lift_to_Gamma_0(NB, c, d, u, v));
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
  Quad g, h, h1, u, v, t;
  int i = PM1.is_principal(g);
  assert (i && "P*M1 must be principal in HeckePAL(P,M1,M2");
  i = P.is_coprime_to(M1, u, v); // u+v=1, u in P, v in M1
  assert (i && "P and M1 are coprime");
  Qideal M3 = M2.equivalent_coprime_to(PM1, h, t, 1); // M2*M3=(h)
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
  Qideal PM1M2M3 = PM1*M2M3;

  // First handle (1:0) mod P*M1, finding a lift [a,b;c,d] with c in M2M3 so h|c
  mat22 m = lift_to_Gamma_0(PM1M2M3, Quad(1), Quad(0), s, r);
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
  for( const auto& x : resmodp)
    {
      c = v*x+u;
      d = v;
#ifdef DEBUG_HECKE
      cout<<"lifting (c:d)=("<<c<<":"<<d<<") from "<<PM1<<" to Gamma_0("<<M2M3<<")"<<endl;
#endif
      m = lift_to_Gamma_0(PM1M2M3, c, d, s, r);
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

vector<mat22> HeckePALQ(Quadprime& P, const Quadprime& Q, Qideal& N)
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

void matop::set_long_name()
{
  long_name = short_name + (char_name.size()?
                            "*chi(" + char_name + ")" :
                            "");
}

string gmatop::name() const
{
  ostringstream s;
  auto ci=coeffs.begin();
  int first = 1;
  for (auto T: ops)
    {
      scalar c = *ci++;
      if (c!=0)
        {
          if (!first)
            {
              s << "+";
            }
          first = 0;
          if (c!=1)
            s << "[" << c << "]";
          s<< T.name();
        }
    }
  return s.str();
}

// constructor for one of the following, where P is prime not dividing N:
// (1) T_P if P principal, else
// (2) T_P*T_{A,A} if P*A^2 principal, else
// (3) T_{P^2} if P^2 principal, else
// (4) T(A,A)*T(P^2); all at level N
matop AutoHeckeOp(Quadprime& P, Qideal& N)
{
  if (P.is_principal())
    {
      return HeckePOp(P, N);
    }
  if (P.has_square_class())
    {
      Qideal A = P.sqrt_coprime_to(N);
      return HeckePChiOp(P, A, N);
    }
  if ((P*P).is_principal())
    {
      return HeckeP2Op(P, N);
    }
  Qideal A = P.equivalent_mod_2_coprime_to(N,1);
  return HeckeP2ChiOp(P, A, N);
}

// Utilities for contructing names

string opname(const Quad& p, const Quad& n)
{
  ostringstream ans;
  ans << (div(p,n) ? "W" : "T") << "(" << p << ")";
  return ans.str();
}

string opname(const string& Tname, const string& Pname)
{
  return Tname + "(" + Pname + ")";
}

string opnameW(const Quadprime& P)
{
  return opname("W", prime_label(P));
}

string opnameT(const Quadprime& P)
{
  return opname("T", prime_label(P));
}

string opname(const Quadprime& P, const Qideal& N)
{
  return opname((P.divides(N) ? "W" : "T"), prime_label(P));
}

string opname(Qideal& N)
{
  return opname("W", ideal_label(N));
}

string opnameAA(Qideal& A)
{
  return opname("chi", ideal_label(A));
}

// Constructors for various matops

matop AtkinLehnerOp(const Quad& p, const Quad& n)
{
  return matop(AtkinLehner(p,n), opname(p,n));
}

// For M1 principal and M1,M2 coprime:

// operator W(M1,M1) at level  N=M1*M2

matop AtkinLehnerOp(Qideal& M1, Qideal& M2)
{
  return matop(AtkinLehner(M1,M2), opname("W", ideal_label(M1)));
}

// For [M1] square with A^2*M1 principal and M1,M2 coprime, A,N coprime:

// operator T(A,A)*W(M1,M2) at level N=M1*M2

matop AtkinLehner_ChiOp(Qideal& M1, const Qideal& M2, Qideal& A)
{
  return matop(AtkinLehner_Chi(M1,M2, A), opname("W", ideal_label(M1)), ideal_label(A));
}

// For Q prime, Q^e||N, Q^e principal:

// operator W(Q^e) at level N

matop AtkinLehnerQOp(const Quadprime& Q, const Qideal& N)
{
  return matop(AtkinLehnerQ(Q,N), opnameW(Q));
}

// For Q prime, Q^e||N, [Q^e] square with A^2*Q^e principal, A coprime
// to N:

// operator T(A,A)W(Q^e) at level N

matop AtkinLehnerQChiOp(const Quadprime& Q, Qideal& A, const Qideal& N)
{
  return matop(AtkinLehnerQ_Chi(Q,A,N), opnameW(Q), ideal_label(A));
}

// For P prime not dividing N, P principal:

// The operator T(P) at level N

matop HeckePOp(Quadprime& P, const Qideal& N)
{
  return matop(HeckeP(P), opname(P,N), "", P);
}

// For P prime not dividing N, [P] square with A^2*P principalm A
// coprime to N:

// The operator T(A,A)T(P) at level N

matop HeckePChiOp(Quadprime& P, Qideal& A, Qideal& N)
{
  ostringstream s;
  s << opname(P,N) << " * " + opnameAA(A);
  return matop(HeckeP_Chi(P,A,N), opnameT(P), ideal_label(A), P, P, A);
}

// For P prime not dividing N with P^2 principal:

// The operator T(P^2) at level N

matop HeckeP2Op(Quadprime& P, Qideal& N)
{
  ostringstream s;
  s << "T(" << P << "^2)";
  return matop(HeckeP2(P,N), opname("T", prime_label(P) + "^2"), "", P);
}

// For P prime not dividing N and A coprime to N with (AP)^2 principal:

// The operator T(A,A)*T(P^2) at level N

matop HeckeP2ChiOp(Quadprime& P, Qideal& A, Qideal& N)
{
  return matop(HeckeP2_Chi(P,A,N), opname("T", prime_label(P) + "^2"), ideal_label(A), P, P, A);
}

// For P,Q distinct primes not dividing N, with P*Q principal:

// The operator T(PQ) at level N

matop HeckePQOp(Quadprime& P, Quadprime& Q, Qideal& N)
{
  return matop(HeckePQ(P,Q,N), opname("T", prime_label(P)+"*"+prime_label(Q)), "", P, Q);
}

// For P,Q distinct primes not dividing N, with [P*Q] square, A^2*P*Q
// principal with A coprime to N:

// The operator T(A,A) T(PQ) at level N

matop HeckePQChiOp(Quadprime& P, Quadprime& Q, Qideal& A, Qideal& N)
{
  return matop(HeckePQ_Chi(P,Q,A,N), opname("T", prime_label(P)+"*"+prime_label(Q)), ideal_label(A), P, Q, A);
}

// For B squarefree principal coprime to N:

// The operator T(B) at level N

matop HeckeBOp(Qideal& B, Qideal& N)
{
  return matop(HeckeB(B,N), opname("T", ideal_label(B)));
}

// For B squarefree coprime to N, with [B] square, A^2*B
// principal with A coprime to N:

// The operator T(A,A) T(B) at level N

matop HeckeBChiOp(Qideal& B, Qideal& A, Qideal& N)
{
  return matop(HeckeB_Chi(B,A,N), opname("T", ideal_label(B)), ideal_label(A));
}

// The operator T(P)W(Q) where P does not divide N, Q^e||N,
// P*Q^e principal.

// We could implement a more general version giving T(A,A)T(P)W(Q^e)
// when [P*Q^e] is square

matop HeckePALQOp(Quadprime& P, const Quadprime& Q, Qideal& N)
{
  return matop(HeckePALQ(P,Q,N), opnameT(P)+"*"+opnameW(Q));
}

// The operator T(P)W(M1) where P does not divide N=M1*M2, M1,M2 coprime and P*M1 principal

// Later we'll implement a more general version giving T(A,A)T(P)W(M1) when [P*M1] is square

matop HeckePALOp(Quadprime& P, Qideal& M1, Qideal& M2)
{
  return matop(HeckePAL(P,M1,M2), opnameT(P)+"*"+opname("W",ideal_label(M1)));
}

matop FrickeOp(Qideal& N)
{
  return matop(Fricke(N), opname("W",ideal_label(N)));
}

matop CharOp(Qideal& A, const Qideal& N)
{
  return matop(Char(A,N), opname("chi", ideal_label(A)));
}


// constructor for one of the following, where Q is prime dividing N
// with Q^e||N:
// (1) W_Q^e if Q^e principal, ( and set t=0);
// (2) W_Q^e*T_{A,A} if Q^e*A^2 principal (and set t=1 and A)
// (3) W_Q*T_P for P  good with Q^e*P principal (and set t=2 and P)
matop AutoALOp(Quadprime& Q, Qideal& N, int& t, Qideal& A, Quadprime& P)
{
  t = 0;
  int e = val(Q,N);
  if (e==0) return matop(); // identity
  Qideal Qe = Q;
  while (--e) Qe*=Q;
  if (Qe.is_principal())
    return AtkinLehnerQOp(Q,N);
  if (Qe.has_square_class())
    {
      t = 1;
      A = Qe.sqrt_coprime_to(N);
      return AtkinLehnerQChiOp(Q,A,N);
    }
  for (auto Pi : Quadprimes::list)
    if (!Pi.divides(N) && (Pi*Qe).is_principal())
      {
        t = 2;
        P = Pi;
        return HeckePALQOp(P, Q, N);
      }
  return matop(); // fall-back to identity (not useful)
}
