// FILE HECKE.H: class matop for Hecke and other operators

#if     !defined(_HECKE_H)
#define _HECKE_H      1       //flags that this file has been included

#include <assert.h>
#include "mat22.h"

// Atkin-Lehner and Hecke operators are implemented in the class matop
// which consists of a list of 2x2 matrices and a string holding the
// operator's name.

// We first define various functions returning the matrix lists.
// These are where all the hard work is done: the later functions
// returning matop objects are simple wrappers which include the
// operator name.

// The general principle is that we can only implement *principal*
// operators on a homspace: for the general case we would have to
// implement a more general homspace with h (=class number)
// components, each the quotient of H_3 by a twisted form of
// Gamma_0(N).  These principal operators (at level N) are of three kinds:

// - T(A) where A=(g) is principal and coprime to N, where the matrices
//   all have determinant g;

// - T(A,A) where A^2=(g) is principal and A is coprime to M: one matrix
//   of determinant g (an (A,A)-matrix of level N);

// - W(M) where M=(g) is principal, M|N and M is coprime to N/M; one
//   matrix of determinant g.

// The most general principal operator then has the form
// T(A1)*T(A2,A2)*W(M). We have not implemented the general case, only
// the ones we have (so far) needed.

// When we are operating on the subspace of the homspace at level N
// with trivial central (unramified nebentypus) character \chi, the
// operators T(A,A) act trivially (by definition), in which case we
// can compute T(A) and W(M) for a wider set of A and M:

//  - If A is coprime to N and the ideal class [A] is square, say
//    A*B^2=(g) with B coprime to AN, then we can compute the action of
//    the principal operator T(A)T(B,B) as a proxy for the action of
//    T(A).

// - Similarly, if M|N is soprime to N/M and has square class, say
//   M*B^2=(g) with B coprime to N, then we can compute the action of
//   the principal operator W(M)T(B,B) as a proxy for the action of
//   W(M).

// In particular, when the class number is odd we can always use these
// proxy operators, since every eigenform is a twist of one with
// trivial character. In this case the action of the principal
// operators T(A)T(B,B) and W(M)T(B,B) (as above) on the principal H_3
// quotient will give the eigenvalues of T(A) and W(M) on this
// trivial-character twist.

// The case of even class number is more complicated, and not fully
// implemented.

// Operators implemented so far:

// HeckePOp(P, N) is T(P) at level N when P is principal
// HeckePSqOp(P, N) is T(P^2) at level N when P^2 is principal (P coprime to N)
// HeckePQOp(P, Q, N) is T(PQ) at level N where P*Q is principal (P,Q coprime to N and distinct)

// AtkinLehnerOp(M1, M2) is W(M1) at level N=M1*M2 when M1 is principal
// AtkinLehnerQOp(Q, N) is W(Q^e) at level N where Q is prime, Q^e||N and Q^e is principal
// FrickeOp(N) is W(N) at level N where N is principal

// HeckePALQOp(P, Q, N) is T(P)W(Q^e) at level N where P*Q^e is principal (P,Q as above)
// HeckePALOp(P, M1, M2) is T(P)W(M1) at level N=M1*M2 where P*M1 is principal (P,M1,M2 as above)

// CharOp(A, N) is T(A,A) at level N where A^2 is principal, A coprime to N

// For all but the last, there is an extended version where the
// relevant ideal is only assumed to be of square class, not
// necessarily principal, and the operator returned is T(B,B) times
// the advertised operator for suitable B coprime to N.  It is up to
// the user whether this is sensible!  It always will be in odd class
// number, where the extended version is the default.  It is also
// relevant when the class group has no 4-torsion, i.e. is
// (odd)*(elementary abelian 2-group).

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

inline string opname(const Qideal& N)
{
  ostringstream ans;
  ans << "W(" << N << ")";
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

// Atkin-Lehner operators

// W(M1) at level N=M1*M2, where M1 is principal (or [M1] square if extend=1) and M1,M2 coprime
mat22 AtkinLehner(Qideal& M1, Qideal& M2, int extend=Quad::class_number%2);
// W(Q^e) at level N where Q^e||N and Q^e is principal (or [Q^e] square if extend=1)
mat22 AtkinLehnerQ(Quadprime& Q, const Qideal& N, int extend=Quad::class_number%2);

// Pure Hecke operators

// T(P) where P does not divide N and P principal (or [P] square if extend=1)
vector<mat22> HeckeP(Quadprime& P, Qideal& N, int extend=Quad::class_number%2);
// T(P^2), when P^2 is principal (and P not), P not dividing N,
// or T(A,A)T(P^2) where A*P is principal if extend=1
vector<mat22> HeckePSq(Quadprime& P, Qideal& N, int extend=Quad::class_number%2);
// T(PQ) for [PQ] square, P,Q not dividing N
vector<mat22> HeckePQ(Quadprime& P, Quadprime& Q, Qideal& N, int extend=Quad::class_number%2);

// Products of Hecke and Atkin-Lehner operators

// T(P)W(M1) at level N for P*M1 principal, P not dividing N=M1*M2
vector<mat22> HeckePAL(Quadprime& P, Qideal& M1, Qideal& M2, int extend=Quad::class_number%2);
// T(P)W(Q^e) at level N for P*Q^e principal, P not dividing N, Q^e||N
vector<mat22> HeckePALQ(Quadprime& P, Quadprime& Q, Qideal& N, int extend=Quad::class_number%2);

inline mat22 Fricke(const Quad& n)
{
  return mat22(Quad::zero,-Quad::one, n,Quad::zero);
}

inline mat22 Fricke(Qideal& N, int extend=Quad::class_number%2) // assumes [N] square
{
  Qideal One(Quad::one);
  return AtkinLehner(N, One);
}

// Matrix inducing T(A,A) at level N, when A^2 is principal and A,N coprime

mat22 Char(Qideal& A, const Qideal& N);

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

inline matop AtkinLehnerOp(const Quad& p, const Quad& n)
{
  return matop(AtkinLehner(p,n), opname(p,n));
}

//  operator W(Q^e) or T(A,A)*W(Q^e) where Q^e||N, either Q^e
//  principal of extend=1 and class [Q^e] square with A^2*P^e
//  principal, A coprime to N

inline matop AtkinLehnerQOp(Quadprime& Q, const Qideal& N, int extend=Quad::class_number%2)
{
  return matop(AtkinLehnerQ(Q,N, extend), opname(Q,N));
}

// assume [M1] square and M1,M2 coprime
//  operator T(A,A)*W(M1) where N=M1*M2, class [M1] square with
//  A^2*M1 principal, A coprime to N, M1,M2 coprime

inline matop AtkinLehnerOp(Qideal& M1, Qideal& M2, int extend=Quad::class_number%2)
{
  ostringstream s;
  s << "W(" << ideal_label(M1) << ")";
  return matop(AtkinLehner(M1,M2, extend), s.str());
}

// The operator T(A,A)*T(P) where P does not divide N, class [P]
// square with P*A^2 principal, A coprime to N

inline matop HeckePOp(Quadprime& P, Qideal& N, int extend=Quad::class_number%2)
{
  return matop(HeckeP(P,N, extend), opname(P,N));
}

// The operator T(A,A)*T(P^2), where P does not divide N, with AP
// principal and A coprime to N

inline matop HeckePSqOp(Quadprime& P, Qideal& N, int extend=Quad::class_number%2)
{
  ostringstream s;
  s << "T(" << P << "^2)";
  return matop(HeckePSq(P,N, extend), s.str());
}

// The operator T(A,A)*T(PQ) where P,Q do not divide N, class [PQ]
// square with A^2*P*Q principal, A coprime to N

inline matop HeckePQOp(Quadprime& P, Quadprime& Q, Qideal& N, int extend=Quad::class_number%2)
{
  ostringstream s;
  s << "T(" << P << "*" << Q << ")";
  return matop(HeckePQ(P,Q,N, extend), s.str());
}

// The operator T(P)W(Q) where P does not divide N, Q^e||N,
// P*Q principal.

// Later we'll implement a more general version giving T(A,A)T(P)W(Q^e) when [P*Q^e] is square

inline matop HeckePALQOp(Quadprime& P, Quadprime& Q, Qideal& N, int extend=Quad::class_number%2)
{
  ostringstream s;
  s << "T(" << P << ")*W(" << Q << ")";
  return matop(HeckePALQ(P,Q,N, extend), s.str());
}

// The operator T(P)W(M1) where P does not divide N=M1*M2, M1,M2 coprime and P*M1 principal

// Later we'll implement a more general version giving T(A,A)T(P)W(M1) when [P*M1] is square

inline matop HeckePALOp(Quadprime& P, Qideal& M1, Qideal& M2, int extend=Quad::class_number%2)
{
  ostringstream s;
  s << "T(" << P << ")*W(" << ideal_label(M1) << ")";
  return matop(HeckePAL(P,M1,M2, extend), s.str());
}

inline matop FrickeOp(Qideal& N, int extend=Quad::class_number%2)
{
  return matop(Fricke(N, extend), opname(N));
}

inline matop CharOp(Qideal& A, const Qideal& N)
{
  ostringstream s;
  s << "nu";
  return matop(Char(A,N), s.str());
}

#endif
