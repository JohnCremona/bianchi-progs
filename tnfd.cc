// FILE TNFD.CC:  test program for Newspace (d-dimensional newform) class
//////////////////////////////////////////////////////////////////////////

//#define LOOPER
#ifdef LOOPER
#include "qidloop.h"
#endif
#include "nfd.h"
#include "field.h"
#include <eclib/pari_init.h>

#define MAXPRIME 10000

int main()
{
  cout << "Program tnfd: constructing Bianchi newforms of arbitrary dimension." << endl;
  eclib_pari_init();

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
  int C4 = is_C4();

  int verbose=1;
  cerr << "Verbose output? (0/1) ";
  cin >> verbose;

  int nap;
  cerr<<"Number of ap? ";
  cin>>nap;
  // Make sure the first nap primes are closed under conjugation:
  Quadprime P=Quadprimes::list[nap-1]; // this is the last prime if we use nap
  if (P.norm() == Quadprimes::list[nap].norm()) // then the next one has the same norm so should also be included
    {
      cerr<<"Increasing nap from "<<nap<<" to "<<nap+1<<" to include the conjugate of "<<P<<endl;
      nap+=1;
    }

  Quad n;
  Qideal N;
#ifdef LOOPER
  long firstn, lastn;
  cerr<<"Enter first and last norm for levels: ";
  cin >> firstn >> lastn;
 // This loop only covers one of each conjugate pair.  For levels not
 // Galois stable, if both_conj is true and output is true, we'll
 // output the conjugate data too.
  Qidealooper loop(firstn, lastn, 0, 1); // 0 = not both conjugates; 1 = sorted within norm
  while( loop.not_finished() )
    {
      N = loop.next();
#else
  while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
    {
#endif
      cout << endl;
      Qideal Nbar = N.conj();
      int conjugate_level_equal = (N==Nbar);

      string Nlabel = ideal_label(N);
      cout << ">>>> Level " << Nlabel <<" = "<<gens_string(N)<<", norm = "<<N.norm()<<" <<<<" << endl;
     homspace* H1 = get_homspace(N, modulus);
     if (verbose)
       cout << "Constructed homspace of dimension " << H1->h1dim()
            << ", cuspidal dimension " << H1->h1cuspdim()
            << ", denominator " << H1->h1cdenom() << endl;
     // if (H1->h1cuspdim()==0)
     //   {
     //     continue;
     //   }
     int maxnp = 7, maxc = 2;
     Newspace forms(H1, maxnp, maxc, verbose);
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
       {
         //     cout << "Outputting newform data to files" << endl;
         forms.output_to_file();
         if (!conjugate_level_equal)
           forms.output_to_file(1);
         continue;
       }
     int nnf_triv_char = std::count_if(forms.newforms.begin(), forms.newforms.end(),
                                       [](Newform F){return F.is_char_trivial()==1;});
     if (n2r)
       cout << "Of these, " << nnf_triv_char
            << (nnf_triv_char==1? " has": " have")
            << " trivial character"<<endl<<endl;

     // cout << "Newform data" << endl;
     // forms.display_newforms();

     cout << "Hecke eigenvalues:" << endl << endl;

     // For homological forms with trivial character we find the full
     // eigensystem with trivial character:

     if (C4 || (nnf_triv_char > 0))
       {
         if (n2r>0)
           {
             if (C4)
               cout << "Full eigensystems for forms with character chi_0 (trivial)"
                << " and character chi_1" << endl;
           }
         else
           {
             cout << "Full eigensystems for forms with trivial character" << endl;
           }
         int inf=1;
         for (auto& F: forms.newforms)
           {
             int triv_char = F.is_char_trivial();
             if (C4 || triv_char)
               {
                 cout << endl;
                 if (verbose)
                   cout << "Computing eigenvalues for newform #" << inf <<endl;
                 map<Quadprime, Eigenvalue> aP = F.aPeigs(nap, verbose);
                 map<Quadprime, int> ALeigs = F.ALeigs(verbose);
                 // display newform data including AL and aP (and all principal eigs if verbose):
                 F.display(1, 1, verbose);
                 cout << endl;
               } // end of triv_char or C4 test
             inf++;
           }
       }

     cout<<endl;
     if (!C4 && (nnf > nnf_triv_char))
       {
         cout << "Principal eigenvalues for forms with non-trivial character" << endl;

         for (auto& F: forms.newforms)
           {
             if (F.is_char_trivial())
               continue;
             cout << endl;
             if (verbose)
               cout << "Computing principal eigenvalues for newform #" << F.get_index() <<endl;
             map<pair<int,string>, Eigenvalue> eigs = F.principal_eigs(nap, verbose);
             F.display(0, 0, 1);
             cout << endl;
             // cout << "Principal eigenvalues involving first " << nap << " good primes:" << endl;
             // for (auto x: eigs)
             //   cout
             //     // << "(" << x.first.first << ") "
             //     << x.first.second << ":\t" << x.second << endl;
           }
       }
     if (verbose)
       cout << "Outputting Newspace data" << endl;
     forms.output_to_file();
     if (!conjugate_level_equal)
       {
         if (verbose)
           cout << "Outputting conjugate Newspace data" << endl;
         forms.output_to_file(1);
       }
     if (verbose)
       cout << "Finished outputting Newspace data"<< endl;
    }     // end of level loop
  cout << endl;
  exit(0);
}   // end of main()

