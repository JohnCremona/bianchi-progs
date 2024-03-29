// FILE HECKE.H: class matop for Hecke and other operators

#if     !defined(_HECKE_H)
#define _HECKE_H      1       //flags that this file has been included

#include <assert.h>
#include "mat22.h"

// Atkin-Lehner and Hecke operators are implemented in the class matop
// which consists of a list of 2x2 matrices and a string holding the
// operator's name.

class matop;  // details below

// We first define various functions returning the matrix lists.
// These are where all the hard work is done: the later functions
// returning matop objects are simple wrappers which include the
// operator name.

// The general principle is that we can only implement *principal*
// operators on a homspace: for the general case we would have to
// implement a more general homspace with h (=class number)
// components, each the quotient of H_3 by a twisted form of
// Gamma_0(N).  These principal operators (at level N) are generated
// by these three kinds:

// Unramified character operators:

// - T(A,A, N) where A^2=(g) is principal and A is coprime to N: one
//   matrix of determinant g (an (A,A)-matrix of level N);

mat22 Char(Qideal& A, const Qideal& N);

// - T(B, N) where B=(g) is principal and coprime to N, where the
//   matrices all have determinant g;

// Pure Hecke operators:

// T(P) =T(P,N) where P is a principal prime (the N(P)+1 matrices here
// do not depend on N)
vector<mat22> HeckeP(Quadprime& P);

// T(PQ, N) for PQ principal, P,Q distinct primes not dividing N
vector<mat22> HeckePQ(Quadprime& P, Quadprime& Q, Qideal& N);

// more general T(B, N) with B square-free
vector<mat22> HeckeB(Qideal& B, Qideal& N);

// T(P^2), when P^2 is principal (and P not), P not dividing N,
vector<mat22> HeckeP2(Quadprime& P, Qideal& N);

// Pure Atkin-Lehner operators

// - W(M1,M2) where M1=(g) is principal and coprime to M2.  Here,
//   N=M1*M2: one matrix of determinant g. As a special case we have
//   W(Q,N) where Q is a prime with Q^e||N, which is W(Q^e,N/Q^e).

// W(M1) at level N=M1*M2, where M1 is principal and M1,M2 coprime
mat22 AtkinLehner(Qideal& M1, Qideal& M2);

// W(Q^e) at level N where Q^e||N and Q^e is principal
mat22 AtkinLehnerQ(Quadprime& Q, const Qideal& N);

// We also need certain products of these.  All the operators above
// have version where the relevant ideal is not principal but has
// square ideal class; these can be made into principal operators by
// composing with T(A,A) for suitable A, for use in situations where
// the central character values (eigenvalues of T(A,A)) are known, in
// particular when it is known that all the T(A,A) act trivially.

// Adjusted Hecke operators:

// T(A,A)T(P,N) where [P] is square, A^2*Pprincipal, P prime and A
// coprime to N
vector<mat22> HeckeP_Chi(Quadprime& P, Qideal&A, Qideal& N);

// T(A,A)T(PQ, N) for P,Q distinct primes not dividing N with [PQ]
// square, A coprime to N with A^2PQ principal
vector<mat22> HeckePQ_Chi(Quadprime& P, Quadprime& Q, Qideal&A, Qideal& N);

// T(A,A)T(B, N) with B square-free, coprime to N, and A^2B principal

vector<mat22> HeckeB_Chi(Qideal& B, Qideal& A, Qideal& N);

// T(A,A)T(P^2), for any prime P not dividing N, (AP)^2 principal
vector<mat22> HeckeP2_Chi(Quadprime& P, Qideal& A, Qideal& N);

// Adjusted Atkin-Lehner operators

// - T(A,A)W(M1,M2) where M1 is coprime to M2, where A^2*M1 is
//   principal.  Here, N=M1*M2: one matrix of determinant g. As a
//   special case we have T(A,A)W(Q,N) where Q is a prime with Q^e||N,
//   where A^2*Q^e is principal, which is W(Q^e,N/Q^e).

// W(M1) at level N=M1*M2, where M1 is principal and M1,M2 coprime
mat22 AtkinLehner_Chi(Qideal& M1, Qideal& M2, Qideal& A);

// W(Q^e) at level N where Q^e||N and Q^e is principal
mat22 AtkinLehnerQ_Chi(Quadprime& Q, Qideal& A, const Qideal& N);

// Products of Hecke and Atkin-Lehner operators

// T(P)W(M1) at level N for P*M1 principal, P not dividing N=M1*M2
vector<mat22> HeckePAL(Quadprime& P, Qideal& M1, Qideal& M2);
// T(P)W(Q^e) at level N for P*Q^e principal, P not dividing N, Q^e||N
vector<mat22> HeckePALQ(Quadprime& P, Quadprime& Q, Qideal& N);

// NB We will also need adjusted versions of these: not yet implemented.

// Utilities for contructing names

inline string opname(const Quad& p, const Quad& n)
{
  ostringstream ans;
  ans << (div(p,n) ? "W" : "T") << "(" << p << ")";
  return ans.str();
}

inline string opname(const Quadprime& P, const Qideal& N)
{
  ostringstream ans;
  ans << (P.divides(N) ? "W" : "T") << "(" << P << ")";
  return ans.str();
}

inline string opname(Qideal& N)
{
  ostringstream ans;
  ans << "W(" << ideal_label(N) << ")";
  return ans.str();
}

inline string opnameAA(Qideal& A)
{
  ostringstream ans;
  string s = ideal_label(A);
  ans << "T(" << s << "," << s << ")";
  return ans.str();
}

// For use only over fields of class number 1, probably now redundant.

// T(P) for P=(p) principal prime
vector<mat22> HeckeP(const Quad& p);
// W(P) for P=(p) principal prime dividing n
mat22 AtkinLehner(const Quad& p, const Quad& n);
// W(P) for P=(p) principal prime dividing N
mat22 AtkinLehner(const Quad& p, Qideal& N);

// For use over general fields

inline mat22 Fricke(const Quad& n)
{
  return mat22(Quad::zero,-Quad::one, n,Quad::zero);
}

inline mat22 Fricke(Qideal& N) // assumes N principal
{
  Qideal One(Quad::one);
  return AtkinLehner(N, One);
}

class matop {  // formal sum of 2x2 matrices
 public:
  vector<mat22> mats;
  string the_name;
  matop() {;}
  explicit matop(const mat22& m, const string& n="") :mats({m}), the_name(n) {;}
  explicit matop(const vector<mat22>& mlist, const string& n="") :mats(mlist), the_name(n) {;}
  mat22 operator[](int i) const {return mats[i];}
  int length() const {return mats.size();}
  string name() const {return the_name;}
};

// Constructors for various matops

inline matop AtkinLehnerOp(const Quad& p, const Quad& n)
{
  return matop(AtkinLehner(p,n), opname(p,n));
}

// For M1 principal and M1,M2 coprime:

// operator W(M1,M1) at level  N=M1*M2

inline matop AtkinLehnerOp(Qideal& M1, Qideal& M2)
{
  ostringstream s;
  s << "W(" << ideal_label(M1) << ")";
  return matop(AtkinLehner(M1,M2), s.str());
}

// For [M1] square with A^2*M1 principal and M1,M2 coprime, A,N coprime:

// operator T(A,A)*W(M1,M2) at level N=M1*M2

inline matop AtkinLehner_ChiOp(Qideal& M1, Qideal& M2, Qideal& A)
{
  ostringstream s;
  s << "W(" << ideal_label(M1) << ") * " + opnameAA(A);
  return matop(AtkinLehner_Chi(M1,M2, A), s.str());
}

// For Q prime, Q^e||N, Q^e principal:

// operator W(Q^e) at level N

inline matop AtkinLehnerQOp(Quadprime& Q, const Qideal& N)
{
  return matop(AtkinLehnerQ(Q,N), opname(Q,N));
}

// For Q prime, Q^e||N, [Q^e] square with A^2*Q^e principal, A coprime
// to N:

// operator T(A,A)W(Q^e) at level N

inline matop AtkinLehnerQChiOp(Quadprime& Q, Qideal& A, const Qideal& N)
{
  ostringstream s;
  s << "W(" << Q << ") * " + opnameAA(A);
  return matop(AtkinLehnerQ_Chi(Q,A,N), opname(Q,N));
}

// For P prime not dividing N, P principal:

// The operator T(P) at level N

inline matop HeckePOp(Quadprime& P, const Qideal& N)
{
  return matop(HeckeP(P), opname(P,N));
}

// For P prime not dividing N, [P] square with A^2*P principalm A
// coprime to N:

// The operator T(A,A)T(P) at level N

inline matop HeckePChiOp(Quadprime& P, Qideal& A, Qideal& N)
{
  ostringstream s;
  s << "T(" << P << ") * " + opnameAA(A);
  return matop(HeckeP_Chi(P,A,N), opname(P,N));
}

// For P prime not dividing N with P^2 principal:

// The operator T(P^2) at level N

inline matop HeckeP2Op(Quadprime& P, Qideal& N)
{
  ostringstream s;
  s << "T(" << P << "^2)";
  return matop(HeckeP2(P,N), s.str());
}

// For P prime not dividing N and A coprime to N with (AP)^2 principal:

// The operator T(A,A)*T(P^2) at level N

inline matop HeckeP2ChiOp(Quadprime& P, Qideal& A, Qideal& N)
{
  ostringstream s;
  s << "T(" << P << "^2) * " + opnameAA(A);
  return matop(HeckeP2_Chi(P,A,N), s.str());
}

// For P,Q distinct primes not dividing N, with P*Q principal:

// The operator T(PQ) at level N

inline matop HeckePQOp(Quadprime& P, Quadprime& Q, Qideal& N)
{
  ostringstream s;
  s << "T(" << P << "*" << Q << ")";
  return matop(HeckePQ(P,Q,N), s.str());
}

// For P,Q distinct primes not dividing N, with [P*Q] square, A^2*P*Q
// principal with A coprime to N:

// The operator T(A,A) T(PQ) at level N

inline matop HeckePQChiOp(Quadprime& P, Quadprime& Q, Qideal& A, Qideal& N)
{
  ostringstream s;
  s << "T(" << P << "*" << Q << ") * " + opnameAA(A);
  return matop(HeckePQ_Chi(P,Q,A,N), s.str());
}

// For B squarefree principal coprime to N:

// The operator T(B) at level N

inline matop HeckeBOp(Qideal& B, Qideal& N)
{
  ostringstream s;
  s << "T(" << B << ")";
  return matop(HeckeB(B,N), s.str());
}

// For B squarefree coprime to N, with [B] square, A^2*B
// principal with A coprime to N:

// The operator T(A,A) T(B) at level N

inline matop HeckeBChiOp(Qideal& B, Qideal& A, Qideal& N)
{
  ostringstream s;
  s << "T(" << B << ") * " + opnameAA(A);
  return matop(HeckeB_Chi(B,A,N), s.str());
}

// The operator T(P)W(Q) where P does not divide N, Q^e||N,
// P*Q^e principal.

// Later we'll implement a more general version giving T(A,A)T(P)W(Q^e) when [P*Q^e] is square

inline matop HeckePALQOp(Quadprime& P, Quadprime& Q, Qideal& N)
{
  ostringstream s;
  s << "T(" << P << ")*W(" << Q << ")";
  return matop(HeckePALQ(P,Q,N), s.str());
}

// The operator T(P)W(M1) where P does not divide N=M1*M2, M1,M2 coprime and P*M1 principal

// Later we'll implement a more general version giving T(A,A)T(P)W(M1) when [P*M1] is square

inline matop HeckePALOp(Quadprime& P, Qideal& M1, Qideal& M2)
{
  ostringstream s;
  s << "T(" << P << ")*W(" << ideal_label(M1) << ")";
  return matop(HeckePAL(P,M1,M2), s.str());
}

inline matop FrickeOp(Qideal& N)
{
  return matop(Fricke(N), opname(N));
}

inline matop CharOp(Qideal& A, const Qideal& N)
{
  ostringstream s;
  s << "nu";
  return matop(Char(A,N), s.str());
}

#endif
