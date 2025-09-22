// NEWHECKE.CC  -- Hecke operators: factored characteristic polynomials on the new cuspidal subspace

// Currently only implemented for class number 1

#include "eclib.h"

#include "matprocs.h"
#include "qidloop.h"
#include "newforms.h"
#include "nfd.h"

//#define LOOPER

#define MAXPRIME 10000

int main(void)
{
  scalar modulus = default_modulus<scalar>();
  //  NextPrime(modulus, power_ZZ(2,256));

  long d, maxpnorm(MAXPRIME);
  int np, ntp;
  Quad n; int show_pols=1;
  cerr << "Enter field (one of "<<class_number_one_fields<<"): " << flush;  cin >> d;
  if (!is_class_number_one(d))
    {
      cerr<<"field must be one of: "<<class_number_one_fields<<endl;
      exit(1);
    }
  Quad::field(d,maxpnorm);
  Quad::displayfield(cout);
  cerr << "See the full char polys (0/1)? "; cin >> show_pols;
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
#ifndef LOOPER
      cerr << "How many Hecke matrices T(P)? ";
      cin >> np;
      cout<<endl;
#endif
      ntp = 0;
      for ( auto& P : Quadprimes::list)
	{
	  if (P.divides(N))
            continue;
          ntp++;
          if (ntp>np) break;
          if (show_pols)
            {
              cout << "Characteristic polynomial of T(" << P << ")"<<endl;
              ZZX charpol = get_full_poly(N,P,modulus);
              cout << "Coefficients: " << charpol << endl;
              cout<<"Factors:"<<endl;
              display_factors(charpol);
              cout << endl;
            }
          cout << "New characteristic polynomial of T(" << P << ")"<<endl;
          ZZX newpol = get_new_poly(N,P,modulus);
          cout << "Coefficients: " << newpol << endl;
          cout<<"Factors:"<<endl;
          display_factors(newpol);
          cout << endl;
        }
    }      // end of if(dim>0)
}       // end of while()
exit(0);
}       // end of main()
