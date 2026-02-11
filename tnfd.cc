// FILE TNFD.CC:  test program for Newspace (d-dimensional newform) class
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

  // We only compute full newform data for forms with trivial
  // character (which cover all newforms over fields with odd class
  // number) except (at present) for C4 fields such as Q(sqrt(-17)).
  int triv_char_only = n2r && !C4;
  if (triv_char_only)
    cout << "Only computing newforms with trivial character" <<endl;

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
       cout << "Constructed level " << Nlabel << " homspace of dimension " << H1->h1dim()
            << ", cuspidal dimension " << H1->h1cuspdim()
            << ", denominator " << H1->h1cdenom() << endl;

     // The Newspace constructor looks for a suitable splitting
     // operator which is a linear combination of up to maxnp prime
     // Hecke operators with coefficients up to maxc:
     int maxnp = 7, maxc = 2;
     Newspace NS(H1, maxnp, maxc, triv_char_only, verbose);
     if (!NS.ok())
       {
         cout << "Failed to find a splitting operator using linear combnations of " << maxnp
              << " operators with coefficients up to" << maxc
              << endl;
         continue; // to next level
       }
     if (verbose)
       cout << "Splitting using " << NS.splitopname() << endl;
     int nnf = NS.nforms();
     cout << "Found " << nnf << " homological newform";
     if (nnf!=1) cout << "s";
     cout << " at level " << Nlabel;
     if (nnf)
       {
         cout << " of dimension";
         if (nnf>1)
           cout << "s " << NS.dimensions();
         else
           cout << " " << NS.dimensions()[0];
       }
     cout << endl;

     if (!nnf)
       {
         //     cout << "Outputting newform data to files" << endl;
         NS.output_to_file();
         if (!conjugate_level_equal)
           NS.output_to_file(1);
         continue;
       }
     int nnf_triv_char = std::count_if(NS.newforms.begin(), NS.newforms.end(),
                                       [](Newform F){return F.is_char_trivial()==1;});
     if (n2r)
       cout << "Of these, " << nnf_triv_char
            << (nnf_triv_char==1? " has": " have")
            << " trivial character"<<endl<<endl;

     // We compute eigenvalues, coefficients and traces and then resort

     // For homological forms with trivial character we find the full
     // eigensystem with trivial character:

     int inf=1;
     for (auto& F: NS.newforms)
       {
         if (C4 || F.is_char_trivial())
           {
             if (verbose)
               cout << "Computing eigenvalues for newform #" << inf <<endl;
             F.compute_eigs(nap, verbose);
             if (verbose)
               cout << endl;
           }
         else
           {
             if (verbose)
               cout << "NOT computing eigenvalues for newform #" << inf << " (non-trivial character)" << endl;
           }
         inf++;
       }

     // After computing eigenvalues, resort the newforms: the forms
     // with trivial character will now be sorted by their sequence of
     // traces.

     if (nnf_triv_char > 1)
       {
         if (verbose)
           cout << "Resorting the newforms with trivial character using traces" << endl;
         NS.sort_newforms();
       }

     if (C4 || (nnf_triv_char > 0))
       {
         if (C4)
           cout << "Full eigensystems for forms with character chi_0 (trivial)"
                << " and character chi_1 (up to unramified quadratic twist)"
                << endl;
         else
           cout << "Full eigensystems for forms with trivial character (up to unramified quadratic twist)"
                << endl;
         for (auto& F: NS.newforms)
           {
             F.display(1, 1, verbose, 1);
             cout << endl;
           }

         if (verbose)
           cout << "Outputting Newspace data" << endl;
         NS.output_to_file();
         if (!conjugate_level_equal)
           {
             if (verbose)
               cout << "Outputting conjugate Newspace data" << endl;
             NS.output_to_file(1);
           }
         if (verbose)
           cout << "Finished outputting Newspace data"<< endl;
       }

     if (n2r)
       {
         cout << "Adding quadratic twists..." << flush;
         NS.add_unram_quadratic_twists();
         cout << "done" << endl << endl;

         if (C4 || (nnf_triv_char > 0))
           {
             if (C4)
               cout << "Full eigensystems for forms with character chi_0 (trivial)"
                    << " and character chi_1 (including unramified quadratic twist)" << endl;
             else
               cout << "Full eigensystems for forms with trivial character"
                    << " (including unramified quadratic twist)" << endl;
             for (auto& F: NS.newforms)
               {
                 F.display(1, 1, verbose, 1);
                 cout << endl;
               }
           }
       }
    }     // end of level loop
  cout << endl;
  exit(0);
}   // end of main()

