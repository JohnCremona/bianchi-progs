// HECKETEST.CC  -- Test for Hecke operators

#include <eclib/mmatrix.h>
#include <NTL/LLL.h>
#include <NTL/mat_poly_ZZ.h>
#include <NTL/ZZXFactoring.h>
#include "qidloop.h"
#include "newforms.h"
//#define LOOPER

#define MAXPRIME 10000

// convert an eclib mat to an NTL mat_ZZ:
mat_ZZ mat_to_mat_ZZ(mat A);
// compute char poly of A/den:
ZZX scaled_charpoly(const mat_ZZ& A, const ZZ& den);
// check that a matrix is a scaled involution:
int check_involution(const mat_ZZ& A, long den=1, int verbose=0);
// check that a matrix commutes with all those in a list:
int check_commute(const mat_ZZ& A, const vector<mat_ZZ>& Blist);
// display factors of a polynomial:
void display_factors(const ZZX& f);
// rank of an NTL matrix:
long rank(mat_ZZ A);
// nullity of an NTL matrix:
long nullity(mat_ZZ A);

// function to sort a factorization vector, first by degree of factor
// then exponent of factor then lexicographically

struct factor_comparison {
  bool operator()(pair_ZZX_long& fac1, pair_ZZX_long& fac2)
  {
    // first sort by degree of the factor
    int s = deg(fac1.a) - deg(fac2.a);
    if(s) return (s<0); // true if fac1 has smaller degree

    // then sort by exponent of the factor
    s = fac1.b - fac2.b;
    if(s) return (s<0); // true if fac1 is to a lower exponent

    // finally lexicographically compare the coefficient lists
    return std::lexicographical_compare(fac1.a.rep.begin(), fac1.a.rep.end(), fac2.a.rep.begin(), fac2.a.rep.end());
  }
}
    fact_cmp;

int main(void)
{
  long d, max(MAXPRIME);
  int np, ntp;
 Quad n; int show_mats, show_pols, show_factors, plusflag, cuspidal=1;
 cerr << "Enter field (one of "<<valid_fields<<"): " << flush;  cin >> d;
 if (!check_field(d))
   {
     cerr<<"field must be one of: "<<valid_fields<<endl;
     exit(1);
   }
 Quad::field(d,max);
 Quad::displayfield(cout);
 int n2r = Quad::class_group_2_rank;
 cerr << "Plus space (0/1)? "; cin>>plusflag;
 cerr << "Cuspidal subspace (0/1)? "; cin>>cuspidal;
 cerr << "See the hecke matrices (0/1)? "; cin >> show_mats;
 cerr << "See the char polys (0/1)? "; cin >> show_pols;
 cerr << "Factor the char polys (0/1)? "; cin >> show_factors;
 Qideal N;
#ifdef LOOPER
 QUINT firstn, lastn;
 cerr<<"Enter first and last norm for Quad loop: ";
 cin >> firstn >> lastn;
 cerr << "How many Hecke matrices T(P)? ";
 cin >> np;
 Qidealooper loop(firstn, lastn, 1, 1); // sorted within norm
 while( loop.not_finished() )
   {
     N = loop.next();
#else
 while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
   {
#endif
  QUINT normn = N.norm();
  cout << ">>>> Level " << ideal_label(N) <<" = "<<gens_string(N)<<", norm = "<<normn<<" <<<<" << endl;
  homspace h(N,plusflag,cuspidal,0);  //level, plusflag, cuspidal, verbose
  int dim = h.h1dim();
  int den = h.h1denom();
  int cden = h.h1cdenom();
  if (cuspidal)
    {
      cout << "Cuspidal dimension = " << dim << endl;
      den=cden;
    }
  else
    {
      cout << "Dimension = " << dim << endl;
    }
  if(den!=1) cout << " denominator = " << den << endl;
  long hmod = h.h1hmod();
  if(hmod)
    {
      cout << "Failed to lift basis from Z/"<<hmod<<" to Z!" << endl;
      cout << "Hence characteristic polynomials are only correct modulo "<<hmod
           <<" and their factorizations are not useful."<<endl;
    }

  vector<Quadprime>::const_iterator pr;
  if (dim>0)
    {
      vector<mat_ZZ> tplist, tpqlist, wqlist, nulist, tpwqlist;

      if (n2r)
        {
          // Compute unramified quadratic characters and check that they are involutions and commute

          vector<Qideal> t2ideals = make_nulist(N);
          mat_ZZ I = ident_mat_ZZ(dim);

          for (vector<Qideal>::iterator t2i=t2ideals.begin(); t2i!=t2ideals.end(); ++t2i)
            {
              Qideal A = *t2i;
              cout << "Computing nu_"<< A <<"..." << flush;
              mat_ZZ nu = mat_to_mat_ZZ(h.calcop(CharOp(A, N), 0, 0));
              cout << "done. " << endl;
              if (show_mats)
                cout << "Matrix is \n" << nu << endl;

              ZZX charpol = scaled_charpoly(nu, to_ZZ(den));
              if (show_pols)
                {
                  cout << "Coefficients of characteristic polynomial are " << charpol << endl;
                }
              if(show_factors)
                {
                  display_factors(charpol);
                }
              cout << endl;
              if (!check_involution(nu,den, 1))
                {
                  exit(1);
                }
              if (!check_commute(nu, nulist))
                {
                  cout << "********* unramified character matrices do not commute with each other ***********" << endl;
                  exit(1);
                }
              nulist.push_back(nu);
            }

          for (int i=0; i<(1<<n2r); i++)
            {
              mat_ZZ A = I;
              string sgs = "";
              for (int j=0; j<n2r; j++)
                {
                  int s = bit(i,j);
                  A = A * (nulist[j] + (s? -den: den)*I);
                  sgs += (s? "-": "+");
                }
              cout << "nu-eigenspace " << sgs << ": " << flush;
              int n = rank(A);
              cout << n <<endl;
            }
        }

      // Compute Atkin-Lehner operators W(Q) and check that they are
      // involutions and commute with each other and with the
      // character matrices.  Restricted to Q such that the power of Q
      // dividing N has square ideal class

      vector<Quadprime> badprimes = N.factorization().sorted_primes();
      vector<Quadprime> squarebadprimes = make_squarebadprimes(N, badprimes);
      vector<Qideal> badprimepowers;
      for (pr=badprimes.begin(); pr!=badprimes.end(); ++pr)
        {
          Quadprime Q = *pr;
          Qideal Qe = Q;
          int e = val(Q, N);
          while(--e) Qe*=Q;
          badprimepowers.push_back(Qe);
        }

      for (pr=squarebadprimes.begin(); pr!=squarebadprimes.end(); ++pr)
        {
          Quadprime Q = *pr;
          int e = val(Q,N);
          Qideal Qe = Q;
          matop op;
          while (--e) Qe *= Q;
          if (Qe.is_principal()) // we compute W_Q directly
            {
              cout << "Computing W("<<Q<<")..." << flush;
              op = AtkinLehnerQOp(Q, N);
            }
          else
            if (Qe.has_square_class()) // we compute T(A,A)*W_Q
              {
                Qideal A = Qe.sqrt_coprime_to(N);
                cout << "Computing W("<<Q<<")*T(A,A) for A="<<A<<"..." << flush;
                op = AtkinLehnerQChiOp(Q, A, N);
              }
            else  // we have an odd power of an ideal with non-square ideal class and compute nothing
              continue;
          mat_ZZ wq =  mat_to_mat_ZZ(h.calcop(op,0,0));
	  cout << "done. " << endl;
          if (show_mats)
            cout << "Matrix is \n" << wq << endl;

          ZZX charpol = scaled_charpoly(wq, to_ZZ(den));
          if (show_pols)
            {
              cout << "Coefficients of characteristic polynomial are " << charpol << endl;
            }
          if(show_factors)
            {
              display_factors(charpol);
            }
          cout << endl;
          if (!check_involution(wq,den, 1))
            {
              exit(1);
            }
          if (!check_commute(wq, nulist))
            {
              cout << "********* W(Q) does not commute with unramified character matrices ***********" << endl;
              exit(1);
            }
          if (!check_commute(wq, wqlist))
            {
              cout << "********* W(Q) matrices do not commute with each other ***********" << endl;
              exit(1);
            }
	  wqlist.push_back(wq);
	}

      // Compute Hecke operators T(P) or T(P^2) and check that they
      // commute with each other and with the character matrices and
      // with A-Ls

#ifndef LOOPER
      cerr << "How many Hecke matrices T(P)? ";
      cin >> np;
      cout<<endl;
#endif
      Quadprime P0; // when h=2 this will be the first good nonprincipal prime.
      int P0_set=0;
      for (pr=Quadprimes::list.begin(), ntp=0;
	   pr!=Quadprimes::list.end() && (ntp<np);
	   ++pr)
	{
          Quadprime P = *pr;
          mat_ZZ tp, tpq, tpwq;
	  if (P.divides(N))
            continue;
          int P_class_order=1; // will hold order of [P] mod squares
          matop op;
          if (P.is_principal())
            {
              op = HeckePOp(P, N);
              cout << "Computing T(" << P << ")..."<<flush;
            }
          else
            {
              if (P.has_square_class())
                {
                  Qideal A = P.sqrt_coprime_to(N);
                  op = HeckePChiOp(P, A, N);
                  cout << "Computing T(" << P << ")*T(A,A) for A="<<A<<"..."<<flush;
                }
              else
                {
                  P_class_order = 2;
                  if ((!P0_set) && (Quad::class_number==2))
                    {
                      P0 = P;
                      P0_set = 1;
                      //cout<<"Setting P0 to "<<P0<<endl;
                    }
                  if ((P*P).is_principal())
                    {
                      op = HeckeP2Op(P, N);
                      cout << "Computing T(" << P << "^2)..."<<flush;
                    }
                  else
                    {
                      Qideal A = P.equivalent_coprime_to(N,1);
                      op = HeckeP2ChiOp(P, A, N);
                      cout << "Computing T(" << P << "^2)*T(A,A) for A="<<A<<"..."<<flush;
                    }
                }
            }

          tp = mat_to_mat_ZZ(h.calcop(op, 0, show_mats));
          cout << "done. " << endl;
          if (show_mats)
            cout << "Matrix is \n" << tp <<endl;

          ntp++;
          tplist.push_back(tp);

          ZZX charpol = scaled_charpoly(tp, to_ZZ(den));
          if (show_pols)
            {
              cout << "Coefficients of characteristic polynomial are " << charpol << endl;
            }
          if(show_factors)
            {
              display_factors(charpol);
            }
          cout << endl;

          if (!check_commute(tp, nulist))
            {
              cout << "********* T(P) does not commute with unramified character matrices ***********" << endl;
              exit(1);
            }
          if (!check_commute(tp, wqlist))
            {
              cout << "********* T(P) does not commute with W(Q) matrices ***********" << endl;
              exit(1);
            }
          if (!check_commute(tp, tplist))
            {
              cout << "********* T(P) does not commute with other T(P) matrices ***********" << endl;
              exit(1);
            }

          // Computing T(P)T(P0) when P*P0 is principal (P, P0 not principal), or just of square class

          if (!(P0_set && (P0!=P) && (P_class_order==2)))
            continue;

          Qideal PP0 = P*P0;
          if (PP0.is_principal())
            {
              op = HeckePQOp(P,P0,N);
              cout << "Computing T(" << P << ") T(" << P0 << ")..."<<flush;
            }
          else
            if (PP0.has_square_class())
              {
                Qideal A = PP0.sqrt_coprime_to(N);
                op = HeckePQChiOp(P,P0,A,N);
                cout << "Computing T(" << P << ") T(" << P0 << ")*T(A,A) with A = "<<A<<"..."<<flush;
              }
            else
              continue;

          tpq = mat_to_mat_ZZ(h.calcop(op, 0, 0));
          cout << "done. " << endl;
          if (show_mats)
            cout << "Matrix is \n" << tpq << endl;
          tpqlist.push_back(tpq);

          charpol = scaled_charpoly(tpq, to_ZZ(den));
          if (show_pols)
            {
              cout << "Coefficients of characteristic polynomial are " << charpol << endl;
            }
          if(show_factors)
            {
              display_factors(charpol);
            }
          cout << endl;

          if (!check_commute(tpq, nulist))
            {
              cout << "********* T(PQ) does not commute with unramified character matrices ***********" << endl;
              exit(1);
            }
          if (!check_commute(tpq, wqlist))
            {
              cout << "********* T(PQ) does not commute with all W matrices ***********" << endl;
              exit(1);
            }
          if (!check_commute(tpq, tplist))
            {
              cout << "********* T(PQ) does not commute with all T matrices ***********" << endl;
              exit(1);
            }
          if (!check_commute(tpq, tpqlist))
            {
              cout << "********* T(PQ) does not commute with other T(PQ) matrices ***********" << endl;
              exit(1);
            }

          // Computing T(P)W(Q) for suitable Q^e||N

          vector<Qideal>::iterator qe=badprimepowers.begin();
          vector<Quadprime>::iterator q=badprimes.begin();
          for (; q!=badprimes.end(); ++q,++qe)
            {
              Quadprime Q = *q;
              Qideal Qe = *qe;
              if ((P*Qe).is_principal())
                {
                  cout<<"Computing T("<<P<<")W("<<ideal_label(Qe)<<")..."<<flush;
                  tpwq = mat_to_mat_ZZ(h.calcop(HeckePALQOp(P,Q,N), 0, 0));
                  cout << "done. " << endl;
                  if (show_mats)
                    cout << "Matrix is \n" << tpwq << endl;
                  tpwqlist.push_back(tpwq);
                  charpol = scaled_charpoly(tpwq, to_ZZ(den));
                  if (show_pols)
                    {
                      cout << "Coefficients of characteristic polynomial are " << charpol << endl;
                    }
                  if(show_factors)
                    {
                      display_factors(charpol);
                    }
                  cout << endl;
                  if (!check_commute(tpwq, nulist))
                    {
                      cout << "********* T(P)W(Q) does not commute with unramified character matrices ***********" << endl;
                      exit(1);
                    }
                  if (!check_commute(tpwq, wqlist))
                    {
                      cout << "********* T(P)W(Q) does not commute with all W matrices ***********" << endl;
                      exit(1);
                    }
                  if (!check_commute(tpwq, tplist))
                    {
                      cout << "********* T(P)W(Q) does not commute with all T matrices ***********" << endl;
                      exit(1);
                    }
                  if (!check_commute(tpwq, tpqlist))
                    {
                      cout << "********* T(P)W(Q) does not commute with all T(PQ) matrices ***********" << endl;
                      exit(1);
                    }
                  if (!check_commute(tpwq, tpwqlist))
                    {
                      cout << "********* T(P)W(Q) does not commute with other T(P)W(Q) matrices ***********" << endl;
                      exit(1);
                    }
                }
            }

        }
    }      // end of if(dim>0)

}       // end of while()
exit(0);
}       // end of main()

mat_ZZ mat_to_mat_ZZ(mat A)
{
  int i, j, d = A.nrows();

  // copy into an NTL matrix:
  mat_ZZ ntl_A;
  ntl_A.SetDims(d,d);
  for(i=1; i<=d; i++)
    for(j=1; j<=d; j++)
      ntl_A(i,j)=conv<ZZ>(A(i,j));
  return ntl_A;
}

ZZX scaled_charpoly(const mat_ZZ& A, const ZZ& den)
{
  ZZX charpol;
  CharPoly(charpol, A);
  long d = deg(charpol);
  if (den>1)
    {
      bigint dpow(1);
      for(int i=0; i<=d; i++)
        {
          SetCoeff(charpol, d-i, coeff(charpol, d-i)/dpow);
          dpow *= den;
        }
    }
  return charpol;
}

int check_involution(const mat_ZZ& A, long den, int verbose)
{
  int res = IsDiag(power(A,2), A.NumRows(), to_ZZ(den*den));
  if (verbose)
    cout << (res? "Involution!": "NOT an involution....") << "\n";
  return res;
}

// check that a matrix commutes with all those in a list:
int check_commute(const mat_ZZ& A, const vector<mat_ZZ>& Blist)
{
  for(vector<mat_ZZ>::const_iterator B = Blist.begin(); B!=Blist.end(); B++)
    {
      if ((A*(*B))!=((*B)*A))
        return 0;
    }
  return 1;
}

// display factors of a polynomaial:
void display_factors(const ZZX& f)
{
  ZZ content; vec_pair_ZZX_long factors;
  factor(content, factors, f);
  ::sort(factors.begin(), factors.end(), fact_cmp);
  cout<<"Factors of characteristic polynomial are:"<<endl;
  long nf = factors.length();
  for(int i=0; i<nf; i++)
    {
      cout<<(i+1)<<":\t"<<factors[i].a
          <<"\t(degree "<<deg(factors[i].a)<<")";
      cout<<"\t to power "<<factors[i].b;
      cout<<endl;
    }
}

// rank of an NTL matrix:
long rank(mat_ZZ A)
{
  ZZ d2;
  return image(d2, A);
}

// nullity of an NTL matrix:
long nullity(mat_ZZ A)
{
  return A.NumRows()-rank(A);
}
