// BASECHANGE.CC  -- dimension and splitting of (non-)base-change forms for Galois stable levels

#include "matprocs.h"
#include "qidloop.h"
#include "newforms.h"
//#define LOOPER

#define MAXPRIME 10000

mat galois_conjugate_matrix(homspace& h, int cuspidal=0, int dual=1, int display=0);
mat galois_conjugate_matrix(homspace& h, const ssubspace& s, int dual=1, int display=0);

int main(void)
{
  long d, maxpnorm(MAXPRIME);
  int np, ntp;
  Quad n; int show_mats=0, show_pols=1, show_factors=1, plusflag=1, cuspidal=1;
 cerr << "Enter field: " << flush;  cin >> d;
 if (!check_field(d))
   {
     cerr<<"field must be one of: "<<valid_fields<<endl;
     exit(1);
   }
 Quad::field(d,maxpnorm);
 Quad::displayfield(cout);
 // cerr << "Plus space (0/1)? "; cin>>plusflag;
 // cerr << "Cuspidal subspace (0/1)? "; cin>>cuspidal;
 // cerr << "See the hecke matrices (0/1)? "; cin >> show_mats;
 // cerr << "See the char polys (0/1)? "; cin >> show_pols;
 // cerr << "Factor the char polys (0/1)? "; cin >> show_factors;
 Qideal N;
#ifdef LOOPER
 long firstn, lastn;
 cerr<<"Enter first and last norm for Quad loop: ";
 cin >> firstn >> lastn;
 cerr << "How many Hecke matrices T(P)? ";
 cin >> np;
 Qidealooper loop(firstn, lastn, 1, 1); // sorted within norm
 while( loop.not_finished() )
   {
     N = loop.next();
#else
 while(cerr<<"Enter Galois stable level (ideal label or generator): ", cin>>N, !N.is_zero())
   {
#endif
     if (!N.is_Galois_stable())
       {
         cout << "Skipping level "<< ideal_label(N) << " as not Galois stable" << endl;
         continue;
       }
  cout << ">>>> Level " << ideal_label(N) <<" = "<<gens_string(N)<<", norm = "<<N.norm()<<" <<<<" << endl;
  homspace h(N,plusflag,0);  //level, plusflag, verbose
  int h1dim = (cuspidal? h.h1cuspdim(): h.h1dim());
  scalar den = (cuspidal? h.h1cdenom(): h.h1denom());
  cout << (cuspidal?"Cuspidal dimension = ":"Dimension = ") << h1dim << endl;
  if(den!=1) cout << " denominator = " << den << endl;
  scalar hmod = h.h1hmod();
  if(hmod!=0)
    {
      cout << "Failed to lift basis from Z/"<<hmod<<" to Z!" << endl;
      cout << "Hence characteristic polynomials are only correct modulo "<<hmod
           <<" and their factorizations are not useful."<<endl;
    }

  if (h1dim>0)
    {
      vector<mat_ZZ> tplist, tpqlist, wqlist, tpwqlist;

      // Compute trivial character subspace
      ssubspace V = h.trivial_character_subspace();

      cout << "Dimension of trivial character subspace = " << dim(V) <<endl;
      pair<int,int> dd = h.trivial_character_subspace_dimensions();
      cout << "Dimension of trivial character subspaces: total = " << dd.first <<", cuspidal = " <<dd.second <<endl;

      // Compute matrix of Galois conjugation restricted to the trivial character subspace:
      mat C = galois_conjugate_matrix(h, V, 1, 1);
      mat_ZZ CC =  mat_to_mat_ZZ(C);
      if (!check_involution(CC, den, hmod, 1))
        cout<<"Problem -- not an involution!"<<endl;

      ZZX charpol = scaled_charpoly(CC, to_ZZ(den), hmod);
      if (show_pols)
        {
          cout << "Coefficients of characteristic polynomial are " << charpol << endl;
        }
      if(show_factors)
        {
          display_factors(charpol);
        }
      cout << endl;

      smat SC(C);
      scalar modulus(default_modulus<scalar>());
      int dplus = SC.nullity(scalar(den), modulus);
      int dminus = SC.nullity(scalar(-den), modulus);
      cout << "Conjugation has eigenvalue multiplicities:\n";
      cout << "+1: "<<dplus<<endl;
      cout << "-1: "<<dminus<<endl;

      // Compute Atkin-Lehner operators W(Q) and check that they are
      // involutions and commute with each other and with the
      // character matrices.  Restricted to Q such that the power of Q
      // dividing N has square ideal class

      vector<Quadprime> badprimes = N.factorization().sorted_primes();
      vector<Quadprime> squarebadprimes = make_squarebadprimes(N, badprimes);
      vector<Qideal> badprimepowers;
      for ( const auto& Q : badprimes)
        {
          Qideal Qe = Q;
          int e = val(Q, N);
          while(--e) Qe*=Q;
          badprimepowers.push_back(Qe);
        }

      for ( auto& Q : squarebadprimes)
        {
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
          mat m = h.calcop(op,cuspidal,0,0);
          mat_ZZ wq =  mat_to_mat_ZZ(m);
	  cout << "done. " << endl;
          if (show_mats)
            cout << "Matrix is \n" << m << "\n="<< m << endl;

          ZZX charpol = scaled_charpoly(wq, to_ZZ(den), hmod);
          if (show_pols)
            {
              cout << "Coefficients of characteristic polynomial are " << charpol << endl;
            }
          if(show_factors)
            {
              display_factors(charpol);
            }
          cout << endl;
          if (!check_involution(wq,den, hmod, 1))
            {
              exit(1);
            }
          if (!check_commute(wq, wqlist, hmod))
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
      ntp = 0;
      for ( auto& P : Quadprimes::list)
	{
          if (ntp>=np) break;
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
                      Qideal A = P.equivalent_mod_2_coprime_to(N,1);
                      op = HeckeP2ChiOp(P, A, N);
                      cout << "Computing T(" << P << "^2)*T(A,A) for A="<<A<<"..."<<flush;
                    }
                }
            }

          mat m = h.calcop(op,cuspidal, 0, show_mats);
          tp = mat_to_mat_ZZ(m);
          cout << "done. " << endl;
          if (show_mats)
            cout << "Matrix is \n" << m <<endl;

          ntp++;
          tplist.push_back(tp);

          ZZX charpol = scaled_charpoly(tp, to_ZZ(den), hmod);
          if (show_pols)
            {
              cout << "Coefficients of characteristic polynomial are " << charpol << endl;
            }
          if(show_factors)
            {
              display_factors(charpol);
            }
          cout << endl;

          if (!check_commute(tp, wqlist, hmod))
            {
              cout << "********* T(P) does not commute with W(Q) matrices ***********" << endl;
              exit(1);
            }
          if (!check_commute(tp, tplist, hmod))
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

          m = h.calcop(op,cuspidal, 0, 0);
          tpq = mat_to_mat_ZZ(m);
          cout << "done. " << endl;
          if (show_mats)
            cout << "Matrix is \n" << m << endl;
          tpqlist.push_back(tpq);

          charpol = scaled_charpoly(tpq, to_ZZ(den), hmod);
          if (show_pols)
            {
              cout << "Coefficients of characteristic polynomial are " << charpol << endl;
            }
          if(show_factors)
            {
              display_factors(charpol);
            }
          cout << endl;

          if (!check_commute(tpq, wqlist, hmod))
            {
              cout << "********* T(PQ) does not commute with all W matrices ***********" << endl;
              exit(1);
            }
          if (!check_commute(tpq, tplist, hmod))
            {
              cout << "********* T(PQ) does not commute with all T matrices ***********" << endl;
              exit(1);
            }
          if (!check_commute(tpq, tpqlist, hmod))
            {
              cout << "********* T(PQ) does not commute with other T(PQ) matrices ***********" << endl;
              exit(1);
            }

          // Computing T(P)W(Q) for suitable Q^e||N

          auto qe=badprimepowers.begin();
          auto q=badprimes.begin();
          for (; q!=badprimes.end(); ++q,++qe)
            {
              Quadprime Q = *q;
              Qideal Qe = *qe;
              if ((P*Qe).is_principal())
                {
                  cout<<"Computing T("<<P<<")W("<<ideal_label(Qe)<<")..."<<flush;
                  mat m1 = h.calcop(HeckePALQOp(P,Q,N),cuspidal, 0, 0);
                  tpwq = mat_to_mat_ZZ(m1);
                  cout << "done. " << endl;
                  if (show_mats)
                    cout << "Matrix is \n" << m1 << endl;
                  tpwqlist.push_back(tpwq);
                  charpol = scaled_charpoly(tpwq, to_ZZ(den), hmod);
                  if (show_pols)
                    {
                      cout << "Coefficients of characteristic polynomial are " << charpol << endl;
                    }
                  if(show_factors)
                    {
                      display_factors(charpol);
                    }
                  cout << endl;
                  if (!check_commute(tpwq, wqlist, hmod))
                    {
                      cout << "********* T(P)W(Q) does not commute with all W matrices ***********" << endl;
                      exit(1);
                    }
                  if (!check_commute(tpwq, tplist, hmod))
                    {
                      cout << "********* T(P)W(Q) does not commute with all T matrices ***********" << endl;
                      exit(1);
                    }
                  if (!check_commute(tpwq, tpqlist, hmod))
                    {
                      cout << "********* T(P)W(Q) does not commute with all T(PQ) matrices ***********" << endl;
                      exit(1);
                    }
                  if (!check_commute(tpwq, tpwqlist, hmod))
                    {
                      cout << "********* T(P)W(Q) does not commute with other T(P)W(Q) matrices ***********" << endl;
                      exit(1);
                    }
                }
            }

        }
    }      // end of if(h1dim>0)

}       // end of while()
exit(0);
}       // end of main()

 mat galois_conjugate_matrix(homspace& h, int cuspidal, int dual, int display)
{
  mat m(h.h1dim(), h.h1dim());
  for (long j=0; j<h.h1dim(); j++)
     {
       m.setcol(j+1,h.chain(h.freemods[j].conj()));
     }
  if(cuspidal) m = restrict_mat(smat(m),h.kern).as_mat();
  if(dual) m = transpose(m);
  if (display)
    {
      cout << "Matrix of Galois conjugation = " << m;
      if (h.h1dim()>1)
        cout << endl;
    }
  return m;
}

mat galois_conjugate_matrix(homspace& h, const ssubspace& s, int dual, int display)
{
  long d=dim(s);
  smat sm(d, h.h1dim());
  for (long j=0; j<d; j++)
     {
       long jj = pivots(s)[j+1]-1;
       sm.setrow(j+1,svec(h.chain(h.freemods[jj].conj())));
     }
  mat m = mult_mod_p(sm,basis(s), default_modulus<scalar>()).as_mat();
  if(!dual) m=transpose(m); // as above code computes the transpose

  if (display)
    {
      cout << "Matrix of Galois conjugation (restricted to subspace of dimension "<<d<<") = " << m;
      if (d>1)
        cout << endl;
    }
  return m;
}
