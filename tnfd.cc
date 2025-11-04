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
#include "field.h"

#define MAXPRIME 10000

int main()
{
  cout << "Program tnfd: constructing Bianchi newforms of arbitrary dimension." << endl;
  scalar modulus = default_modulus<scalar>();
#if (SCALAR_OPTION==3)
  //  NextPrime(modulus, power_ZZ(to_ZZ(2),256));
  NextPrime(modulus, pow(to_ZZ(2),512));
#endif
  long d, maxpnorm(MAXPRIME);
  cerr << "Enter field: " << flush;
  cin >> d;
  if (d<1)
    {
      cout << "Field parameter d must be positive for Q(sqrt(-d))" << endl;
      exit(0);
    }
  Quad::field(d,maxpnorm);
  Quad::displayfield(cout);
  int n2r = Quad::class_group_2_rank;

  int verbose=1;
  cerr << "Verbose output? (0/1) ";
  cin >> verbose;

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
     if (verbose)
       cout << "Splitting using " << forms.splitopname() << endl;
     int nnf = forms.nforms();
     cout << "Found " << nnf << " homological newform";
     if (nnf!=1) cout << "s";
     if (nnf)
       {
         cout << " of dimension";
         if (nnf>1)
           cout << "s " << forms.dimensions();
         else
           cout << " " << forms.dimensions()[0];
       }
     cout << endl;
     if (!nnf)
       continue;
     int nnf_triv_char = std::count_if(forms.newforms.begin(), forms.newforms.end(),
                                       [](Newform F){return F.trivial_char()==1;});
     if (n2r)
       cout << "Of these, " << nnf_triv_char
            << (nnf_triv_char==1? " has": " have")
            << " trivial character"<<endl<<endl;

     cout << "Newform data" << endl;
     forms.display_newforms();

     cout << "Hecke eigenvalues:" << endl << endl;

     // For homological forms with trivial character we find the full
     // eigensystem with trivial character:

     if (nnf_triv_char > 0)
       {
         if (n2r>0)
           cout << "Full eigensystems for forms with trivial character" << endl;
         int inf=1;
         for (auto F: forms.newforms)
           {
             if (F.trivial_char())
               {
                 cout << endl;
                 if (verbose)
                   cout << "Computing eigenvalues for newform #" << inf <<endl;
                 map<Quadprime, Eigenvalue> eigs = F.eigs(nap, verbose);
                 F.display(1);
                 cout << endl;
                 cout << "Eigenvalues for first " << nap << " good primes:" << endl;
                 for (auto x: eigs)
                   cout << x.first << ":\t" << x.second << endl;
                 if (F.is_self_twist()==+1)
                   cout << "*** form appears to have self-twist ***" << endl;
               }
             inf++;
           }
       }

     cout<<endl;

     if (nnf > nnf_triv_char)
       {
         cout << "Principal eigenvalues for forms with non-trivial character" << endl;

         for (auto F: forms.newforms)
           {
             if (F.trivial_char())
               continue;
             cout << endl;
             if (verbose)
               cout << "Computing principal eigenvalues for newform #" << F.get_index() <<endl;
             map<string, FieldElement> eigs = F.principal_eigs(nap, verbose);
             F.display(1);
             cout << endl;
             cout << "Principal eigenvalues for first " << nap << " good primes:" << endl;
             for (auto x: eigs)
               cout << x.first << ":\t" << x.second << endl;
             if (F.is_self_twist()==+1)
               cout << "*** form appears to have self-twist ***" << endl;
           }
       }
    }     // end of level loop
  exit(0);
}   // end of main()

