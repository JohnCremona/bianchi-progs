// FILE TNFD.CC:  test program for Newforms (d-dimensional newform) class
//////////////////////////////////////////////////////////////////////////
//
// Adapted from the similar program (over Q) in eclib
//
//////////////////////////////////////////////////////////////////////////

//#define LOOPER
#ifdef LOOPER
#include "qidloop.h"
#endif
#include "nfd.h"
#include "heckefield.h"

#define MAXPRIME 10000

int main()
{
  cout << "Program tnfd: constructing Bianchi newforms of arbitrary dimension." << endl;
  scalar modulus = default_modulus<scalar>();
#if (SCALAR_OPTION==3)
  //  NextPrime(modulus, power_ZZ(2,256));
  NextPrime(modulus, power_ZZ(2,512));
#endif
  long d, maxpnorm(MAXPRIME);
  cerr << "Enter field: " << flush;
  cin >> d;
  Quad::field(d,maxpnorm);
  Quad::displayfield(cout);
  //int n2r = Quad::class_group_2_rank;

  int verbose=1;
  cerr << "Verbose output? (0/1) ";
  cin >> verbose;

  int triv_char_only = 0;
  if (Quad::class_group_2_rank > 0)
    {
      cout << "Trivial character only? (0/1) ";
      cin >> triv_char_only;
    }

  int nap=5;
  cerr<<"Number of ap? ";
  cin>>nap;

  Quad n;
  Qideal N;
#ifdef LOOPER
  long firstn, lastn;
  cerr<<"Enter first and last norm for Quad loop: ";
  cin >> firstn >> lastn;
  Qidealooper loop(firstn, lastn, 0, 1); // 0 = not both conjugates; 1 = sorted within norm
  while( loop.not_finished() )
    {
      N = loop.next();
#else
  while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
    {
#endif
      cout << endl;
      cout << ">>>> Level " << ideal_label(N) <<" = "<<gens_string(N)<<", norm = "<<N.norm()<<" <<<<" << endl;
     homspace* H1 = get_homspace(N, modulus);
     if (verbose)
       cout << "Constructed homspace of dimension " << H1->h1dim()
            << ", cuspidal dimension " << H1->h1cuspdim()
            << ", denominator " << H1->h1cdenom() << endl;

     int maxnp = 7, maxc = 2;
     Newforms forms(H1, maxnp, maxc, verbose);
     if (!forms.ok())
       {
         cout << "Failed to find a splitting operator using lnear combnations of " << maxnp
              << " operators with coefficients up to" << maxc
              << endl;
         continue; // to next level
       }
     cout << "Splitting using " << forms.splitopname() << endl;
     int nnf = forms.nforms();
     cout << "Found " << nnf << " newforms";
     if (nnf)
       cout << " with dimensions " << forms.dimensions();
     cout << endl;
     if (!nnf)
       continue;
     forms.display_newforms(triv_char_only);
     int nnf_triv_char = std::count_if(forms.newforms.begin(), forms.newforms.end(),
                                       [](Newform F){return F.trivial_char()==1;});
     if (triv_char_only&& nnf_triv_char==0)
       {
         cout << "No newforms have trivial character"<<endl;
         continue;
       }
     cout << "Hecke eigenvalues:" << endl;
     int ip = 0;
     for ( auto& P : Quadprimes::list)
       {
         if (P.divides(N))
           continue;
         ip++;
         if (ip>nap)
           break;
         matop T = AutoHeckeOp(P,N);
         vector<HeckeFieldElement> apvec = forms.eig(T);
         cout<<T.name() << ":\t";
         auto F = forms.newforms.begin();
         for (auto ap: apvec)
           {
             if ((!triv_char_only) || F->trivial_char())
               {
                 cout << ap << "\t";
               }
             ++F;
           }
         cout<<endl;
       } // end of prime loop
     //  cout<<endl;
    }     // end of level loop
  exit(0);
}   // end of main()

