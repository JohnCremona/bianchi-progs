// FILE TNFD.CC:  test program for nfd (d-dimensional newform) class
//////////////////////////////////////////////////////////////////////////
//
// Adapted from the similar program (over Q) in eclib
//
//////////////////////////////////////////////////////////////////////////

#include "nfd.h"

//#define DEBUG
//#define COMPARE_OLD
//#define LOOPER
#define MAXPRIME 10000

int main()
{
  // init_time();
  cout << "Program tnfd." << endl;
  scalar modulus = default_modulus<scalar>();
  NextPrime(modulus, power_ZZ(2,256));

  long d, maxpnorm(MAXPRIME);
  cerr << "Enter field (one of "<<class_number_one_fields<<"): " << flush;  cin >> d;
  if (!is_class_number_one(d))
    {
      cerr<<"field must be one of: "<<class_number_one_fields<<endl;
      exit(1);
    }
  Quad::field(d,maxpnorm);
  Quad::displayfield(cout);

  int verbose=1;
  cerr << "Verbose output? (0/1) "; cin >> verbose;

  Quad n;
  Qideal N;
#ifdef LOOPER
  long firstn, lastn;
  cerr<<"Enter first and last norm for Quad loop: ";
  cin >> firstn >> lastn;
  Qidealooper loop(firstn, lastn, 1, 1); // sorted within norm
  while( loop.not_finished() )
    {
      N = loop.next();
#else
  while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
    {
#endif
     cout << ">>>> Level " << ideal_label(N) <<" = "<<gens_string(N)<<", norm = "<<N.norm()<<" <<<<" << endl;
     homspace hplus(N, modulus, 1);
     int dimH = hplus.h1cuspdim();
     cout << "dimension = " << dimH << endl;

     nfd form = nfd(&hplus, verbose);
     int dimS = form.dimS;
     cout<<"Finished constructing nfd object, with an irreducible subspace of dimension "<<dimS
         <<" and defining polynomial "<<form.f<<endl;
     if(dimS==0)
       continue;
     bigint den=form.dHS;
     int i, ip, nap=5;
     cout<<"Number of ap? ";  cin>>nap;
     //     start_time();
      ip = 0;
      for ( auto& P : Quadprimes::list)
	{
	  if (P.divides(N))
            continue;
          ip++;
          if (ip>nap)
            break;
          if(0) // verbose)
	   {
	     mat_m tp = form.heckeop(P);
	     if(den>1) cout<<den<<"*";
	     cout<<"Matrix of ";
	     cout<<"T(" <<ideal_label(P)<<") = ";
	     tp.output_pari(cout);
             cout << endl;
	     vector<bigint> cptp = tp.charpoly();
	     for(i=0; i<dimS; i++)
	       {
		 bigint temp = cptp[i];
		 divide_exact(temp,form.Hscales[dimS-i],temp);
		 divide_exact(temp,form.Sscales[dimS-i],temp);
		 cptp[i]=temp;
	       }
	     cout<<"char poly = "<<cptp<<endl;
	   }

          vec_m apvec = form.ap(P);
          cout<<"a_"<<ideal_label(P)<<" = ";
          if (dimS==1)
            cout << apvec[1];
          else
            cout << apvec;
          cout<<endl;
        } // end of prime loop
    }     // end of level loop
  //     stop_time();
  //     show_time();
  cout<<endl;
  exit(0);
}      // end of main()

