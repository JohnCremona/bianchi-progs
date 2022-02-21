// FILE MAT22.CC

#include "primes.h"
#include "mat22.h"

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

// Same as above with M1=P^e, M2=N/P^e where P^e||N

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
    return Hecke(P.gen(), N);
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

// Matrix inducing T(A,A) at level N, when A^2 is principal and A+N=1

mat22 Char(Qideal& A, const Qideal& N)
{
  Quad g;
  return A.AB_matrix_of_level(A, N, g);
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

//#define DEBUG_LIFT

// return a matrix [a, b; c, d] with det=1 and (c:d)=(cc:dd) in P^1(N)
mat22 lift_to_SL2(Qideal& N, const Quad& cc, const Quad& dd)
{
  Quad a, b, c(cc), d(dd), inv, x, y, z, h;
#ifdef DEBUG_LIFT
  cout<<"Lifting symbol (c:d)=("<<c<<":"<<d<<") mod "<<N<<" to SL2"<<endl;
#endif
  // Special cases (1): (c:1), (1:d) need no work:
  Quad one=Quad::one, zero=Quad::zero;
  if (d==one) return mat22(one,zero,c,one);
  if (c==one) return mat22(zero,-one,one,d);

#ifdef DEBUG_LIFT
  cout<<"Neither c nor d is invertible modulo "<<N<<": testing whether ideal (c,d) is principal"<<endl;
#endif
  // General case: neither c nor d is invertible.

  // Test if (c,d)=(h), principal:
  h = quadbezout(c,d, x, y);
  if (!h.is_zero()) // then is principal with c*x+d*y=h, and h=1 since reduced
    {
#ifdef DEBUG_LIFT
      cout<<"ideal (c,d)=("<<h<<"), success"<<endl;
#endif
      a = y;
      b = -x;
      if (h.norm()>1) // should not happen as (c:d) was reduced
        {
          c /= h;
          d /= h;
        }
    }
  else
    {  // Now we must work harder.
#ifdef DEBUG_LIFT
      cout<<" (c,d) not principal, working harder..."<<endl;
#endif
      int t = N.is_coprime_to(c, d, x, y, 1);   // c*x+d*y = 1 mod N with y invertible
      assert (t==1);
#ifdef DEBUG_LIFT
      cout<<" c*x+d*y=1 mod N with x = "<<x<<" and y = "<<y<<endl;
#endif
      t = N.is_coprime_to(y, z);            // y*z = 1 mod N
      assert (t==1);
#ifdef DEBUG_LIFT
      cout<<" inverse of y mod N is z = "<<z<<" with y*z="<<y*z<<endl;
#endif
      a = Quad::one;
      b = N.reduce(-x*z);
      c = N.reduce(c*y); // so b*c == -x*c = d*y-1
      d = a + b*c; // = d*y mod N
    }
#ifdef DEBUG_LIFT
  cout<<" replacing c by "<<c<<" and d by "<<d<<", which are coprime"<<endl;
#endif
  assert (a*d-b*c==one);
  mat22 M(a,b,c,d);
#ifdef DEBUG_LIFT
  cout<<" returning  "<< M <<endl;
#endif
  return M;
}

// return a matrix [a, b; c, d] with det=1 and c in M and (c:d)=(cc:dd) in P^1(N)
// If (u,v)!=(0,0) they should satisfy u+v=1 with u in N, v in M,
// otherwise such u,v will be computed and returned.
mat22 lift_to_Gamma_0(Qideal& M, Qideal& N, const Quad& cc, const Quad& dd, const Quad& u, const Quad& v)
{
  // CRT: lift (cc:dd) in P^1(N) to (c:d) in P^1(MN) which also lifts (0:1) in P^1(M), then lift that
  Quad uu(u), vv(v);
  if (uu.is_zero() && vv.is_zero())
    M.is_coprime_to(N, uu, vv); // u+v=1, u in M, v in N
  Qideal MN = M*N;
  return lift_to_SL2(MN, cc*uu, dd*uu+vv);
}
