// NEWHECKE.CC  -- Hecke operators: factored characteristic polynomials on the new cuspidal subspace

#include "eclib.h"

#include "matprocs.h"
#include "qidloop.h"
#include "newforms.h"

//#define LOOPER

#define MAXPRIME 10000

int main(void)
{
  scalar modulus = default_modulus<scalar>();
#if (SCALAR_OPTION==3)
  //  NextPrime(modulus, power_ZZ(2,256));
    NextPrime(modulus, power_ZZ(2,512));
#endif
  long d, maxpnorm(MAXPRIME);
  cerr << "Enter field: " << flush;
  cin >> d;
  if (!check_field(d))
   {
     cerr<<"field must be one of: "<<valid_fields<<endl;
     exit(1);
   }
  Quad::field(d,maxpnorm);
  Quad::displayfield(cout);

  int show_pols=1;
  cerr << "See the full char polys (0/1)? ";
  cin >> show_pols;
  int np;
  cerr << "How many Hecke matrices T(P)? ";
  cin >> np;
  Qideal N;
#ifdef LOOPER
  long firstn, lastn;
  cerr<<"Enter first and last norm for Quad loop: ";
  cin >> firstn >> lastn;
  Qidealooper loop(firstn, lastn, 1, 1); // both conjugates, sorted within norm
  while( loop.not_finished() )
    {
      N = loop.next();
#else
      while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
        {
#endif
          cout << ">>>> Level " << ideal_label(N) <<" = "<<gens_string(N)<<", norm = "<<N.norm()<<" <<<<" << endl;
  homspace* h = get_homspace(N, modulus);
  int dim = h->h1cuspdim();
  cout << "Cuspidal dimension = " << dim << endl;
  //scalar den = h->h1cdenom();
  //if(den!=1) cout << " denominator = " << den << endl;

  if (dim>0)
    {
      scalar hmod = h->h1hmod();
      if(hmod!=0)
        {
          cout << "Failed to lift basis from Z/"<<hmod<<" to Z -- "
               << "homology is modulo "<<hmod<<endl;
          ZZ_p::init(to_ZZ(hmod));
        }
      int ntp = 0;
      for ( auto& P : Quadprimes::list)
	{
	  if (P.divides(N)) // || !P.is_principal())
            continue;
          ntp++;
          if (ntp>np) break;
          matop T = AutoHeckeOp(P, N);
          cout << "P = " << P << " is in ideal class " << P.ideal_class()
               << ", using T = "<<T.name()<<endl;
          if (show_pols)
            {
              ZZX charpol = get_full_poly(N, T, modulus);
              cout << "Full characteristic polynomial: "
                   << polynomial_string(charpol)
                //                   << "\tCoefficients: " << charpol
                   << endl;
              cout <<"Factors:"<<endl;
              display_factors(charpol);
              cout << endl;
            }
          ZZX newpol = get_new_poly(N, T, modulus);
          int dimnew = deg(newpol);
          if (dimnew==0)
            {
              cout << "Newspace is trivial" << endl;
              break;
            }
          if (ntp==1)
            {
              cout << "Newspace has dimension " << dimnew << endl;
            }
          cout << "New characteristic polynomial: "
               << polynomial_string(newpol)
            //               << "\tCoefficients: " << newpol
               << endl;
          cout << "Factors:" <<endl;
          display_factors(newpol);
          cout << endl;
        }
    }      // end of if(dim>0)
}       // end of while()
exit(0);
}       // end of main()
