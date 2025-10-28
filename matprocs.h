// FILE MATPROCS.H:  interface with NTL mat_ZZ and related classes

#if     !defined(_MATPROCS_H)
#define _MATPROCS_H      1       //flags that this file has been included

// include stuff from eclib
#include "intprocs.h"

#include "qidloop.h"

// convert an eclib mat to an NTL mat_ZZ:
template<class T>
mat_ZZ mat_to_mat_ZZ(Zmat<T> A);

// convert an eclib mat to an NTL mat_ZZ_p:
mat_ZZ_p mat_to_mat_ZZ_p(mat A);

// compute char poly of A:
ZZX charpoly(const mat_ZZ& A);

// compute char poly of A:
inline ZZX charpoly(const mat_m& A) {return charpoly(mat_to_mat_ZZ(A));}

// compute char poly of A/den mod m:
ZZX scaled_charpoly(const mat_ZZ& A, const ZZ& den = to_ZZ(1), const scalar& m = scalar(0));

// return A mod m (or just A if m==0)
mat_ZZ reduce_mat(const mat_ZZ& A, const ZZ& m);

// evaluate f(A) (assumes f monic)
mat_ZZ evaluate(const ZZX& f, const mat_ZZ& A);
// evaluate f(A) (assumes f monic)
template<class T>
mat_m evaluate(const ZZX& f, const Zmat<T>& A);

// p should be monic:
mat_ZZ CompanionMatrix(const ZZX& p);

// check that a matrix is a scaled involution (modulo m, if m!=0):
int check_involution(const mat_ZZ& A, scalar den, const scalar& m, int verbose=0);

// check that two matrices commute (modulo m, if m!=0):
int commute(const mat_ZZ& A, const mat_ZZ& B, const scalar& m);

// check that a matrix commutes (modulo m, if m!=0) with all those in a list:
int check_commute(const mat_ZZ& A, const vector<mat_ZZ>& Blist, const scalar& m);
// rank of an NTL matrix:
long rank(mat_ZZ A);

// nullity of an NTL matrix:
long nullity(mat_ZZ A);

// Linear combinarion of n>0 matrices, all dxd
mat_m lin_comb_mats(const vec_m& co, const vector<mat_m>& mats);
// Linear combinarion of n>0 matrices, all dxd
mat_m lin_comb_mats(const vector<ZZ>& co, const vector<mat_m>& mats);

// same as m.output(cout) except no newlines between rows
template<class T>
void output_flat_matrix(const Zmat<T>& m, ostream&s = cout);

#endif
