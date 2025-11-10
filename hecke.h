// FILE HECKE.H: class matop for Hecke and other operators

#if     !defined(_HECKE_H)
#define _HECKE_H      1       //flags that this file has been included

#include <assert.h>
#include "mat22.h"

// Atkin-Lehner and Hecke operators are implemented in the class matop
// which consists of a list of 2x2 matrices and a string holding the
// operator's name. The class gmatop is a formal Z-linear combination
// of matops.

class matop;  // details below
class gmatop;  // details below

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
mat22 AtkinLehnerQ(const Quadprime& Q, const Qideal& N);

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
mat22 AtkinLehner_Chi(const Qideal& M1, const Qideal& M2, const Qideal& A);

// W(Q^e) at level N where Q^e||N and Q^e is principal
mat22 AtkinLehnerQ_Chi(const Quadprime& Q, const Qideal& A, const Qideal& N);

// Products of Hecke and Atkin-Lehner operators

// T(P)W(M1) at level N for P*M1 principal, P not dividing N=M1*M2
vector<mat22> HeckePAL(Quadprime& P, Qideal& M1, Qideal& M2);
// T(P)W(Q^e) at level N for P*Q^e principal, P not dividing N, Q^e||N
vector<mat22> HeckePALQ(Quadprime& P, const Quadprime& Q, const Qideal& N);

// Utilities for contructing names

string opname(const Quad& p, const Quad& n);          // T(p) or W(p)
string opname(const Quadprime& P, const Qideal& N);   // T(P) or W(P)
string opname(Qideal& N);                             // W(N)
string opnameAA(Qideal& A);                           // T(A,A)

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

// a matop is a formal sum of 2x2 matrices representing a principal
// Hecke operator
class matop {
 public:
  vector<mat22> mats;
  // For an operator like T(A) the long and short names are "T(A)" and
  // the char name is "".  For an operator like T(A)T(B,B) the short
  // name is "T(A)", the char name is the label of B and the full name
  // is the concatenation of these with " * "
  string long_name;
  string short_name;
  string char_name;
  vector<int> genus_char; // list of chi(B)
  matop() // default is identity
    :mats({mat22::identity}), short_name("I"), char_name("") {set_long_name();}
  explicit matop(const mat22& m, const string& n="", const string& c="")
    :mats({m}), short_name(n), char_name(c)  {set_long_name();}
  explicit matop(const vector<mat22>& mlist, const string& n="", const string& c="")
    :mats(mlist), short_name(n), char_name(c) {set_long_name();}
  void set_long_name();
  mat22 operator[](int i) const {return mats[i];}
  int length() const {return mats.size();}
  string sname() const {return short_name;}
  string cname() const {return char_name;}
  string name() const  {return long_name;}
};

// a gmatop is a linear combination of matops
class gmatop {
 public:
  vector<matop> ops;
  vector<scalar> coeffs;
  gmatop() {;}
  explicit gmatop(const vector<matop>& Tlist, const vector<scalar> clist)
    :ops(Tlist), coeffs(clist)  {;}
  explicit gmatop(const vector<matop>& Tlist)
    :ops(Tlist), coeffs(vector<scalar>(Tlist.size(), scalar(1)))  {;}
  explicit gmatop(const matop& T)
    :ops({T}), coeffs({scalar(1)})  {;}
  void set_coeff(int i, const scalar& c) {coeffs[i] = c;}
  string name() const;
};

// Constructors for various matops

matop AtkinLehnerOp(const Quad& p, const Quad& n);

// For M1 principal and M1,M2 coprime:
// operator W(M1,M1) at level  N=M1*M2

matop AtkinLehnerOp(Qideal& M1, Qideal& M2);

// For [M1] square with A^2*M1 principal and M1,M2 coprime, A,N coprime:
// operator T(A,A)*W(M1,M2) at level N=M1*M2

matop AtkinLehner_ChiOp(Qideal& M1, const Qideal& M2, Qideal& A);

// For Q prime, Q^e||N, Q^e principal:
// operator W(Q^e) at level N

matop AtkinLehnerQOp(const Quadprime& Q, const Qideal& N);

// For Q prime, Q^e||N, [Q^e] square with A^2*Q^e principal, A coprime
// to N: operator T(A,A)W(Q^e) at level N

matop AtkinLehnerQChiOp(const Quadprime& Q, Qideal& A, const Qideal& N);

// For principal P prime not dividing N: the operator T(P) at level N

matop HeckePOp(Quadprime& P, const Qideal& N);

// For P prime not dividing N, [P] square with A^2*P principal, A
// coprime to N: the operator T(A,A)T(P) at level N

matop HeckePChiOp(Quadprime& P, Qideal& A, Qideal& N);

// For P prime not dividing N with P^2 principal: the operator T(P^2)
// at level N

matop HeckeP2Op(Quadprime& P, Qideal& N);

// For P prime not dividing N and A coprime to N with (AP)^2 principal:
// The operator T(A,A)*T(P^2) at level N

matop HeckeP2ChiOp(Quadprime& P, Qideal& A, Qideal& N);

// For P,Q distinct primes not dividing N, with P*Q principal: the
// operator T(PQ) at level N

matop HeckePQOp(Quadprime& P, Quadprime& Q, Qideal& N);

// For P,Q distinct primes not dividing N, with [P*Q] square, A^2*P*Q
// principal with A coprime to N: the operator T(A,A) T(PQ) at level N

matop HeckePQChiOp(Quadprime& P, Quadprime& Q, Qideal& A, Qideal& N);

// For squarefree principal B coprime to N: the operator T(B) at level N

matop HeckeBOp(Qideal& B, Qideal& N);

// For B squarefree coprime to N, with [B] square, A^2*B principal
// with A coprime to N: the operator T(A,A) T(B) at level N

matop HeckeBChiOp(Qideal& B, Qideal& A, Qideal& N);

// The operator T(P)W(Q) where P does not divide N, Q^e||N, P*Q^e
// principal.

// (We have not yet implemented a more general version giving
// T(A,A)T(P)W(Q^e) when [P*Q^e] is square)

matop HeckePALQOp(Quadprime& P, const Quadprime& Q, const Qideal& N);

// The operator T(P)W(M1) where P does not divide N=M1*M2, M1,M2
// coprime and P*M1 principal

// (We have not yet implemented a more general version giving
// T(A,A)T(P)W(M1) when [P*M1] is square)

matop HeckePALOp(Quadprime& P, Qideal& M1, Qideal& M2);

matop FrickeOp(Qideal& N);

matop CharOp(Qideal& A, const Qideal& N);

// constructor for one of the following, where P is prime not dividing N:
// (1) T_P if P principal, else
// (2) T_P*T_{A,A} if P*A^2 principal, else
// (3) T_{P^2} if P^2 principal, else
// (4) T(A,A)*T(P^2); all at level N
matop AutoHeckeOp(Quadprime& P, Qideal& N);

// constructor for one of the following, where Q is prime dividing N
// with Q^e||N:
// (1) W_Q^e if Q^e principal, ( and set t=0);
// (2) W_Q^e*T_{A,A} if Q^e*A^2 principal (and set t=1 and A)
// (3) W_Q*T_P for P  good with Q^e*P principal (and set t=2 and P)
matop AutoALOp(Quadprime& Q, Qideal& N, int& t, Qideal& A, Quadprime& P);

#endif
