// FILE modularity.cc: checks individual ap against input values for given level(s)

#include <fstream>
#include "newforms.h"
#include "eclib/curvesort.h" // for letter codes

long prime_index(const Quad& p)
{
  return find(quadprimes.begin(), quadprimes.end(), makepos(p)) - quadprimes.begin();
}

int main(void)
{
  cerr << "Program modularity: for given field and level, reads newforms from\n";
  cerr << "file (which should exist already), and computes Hecke eigenvalues for\n";
  cerr << "a list of input primes, checking the values against input values, modulo a prime p.\n";
  cerr << "----------------------------------------------------------------\n\n";
  cerr << "Input format:\n";
  cerr << "<field> <level> <modulus> <nforms> <nprimes>\n";
  cerr << "followed by nprimes lines each containing a prime's coefficients (2 integers) and nforms integers:\n";
  cerr << "<prime> <ap_1> <ap_2> ... <ap_nforms>\n";
  cerr << "----------------------------------------------------------------\n\n";
  // max is the maximum norm of precomputed primes.  It should be
  // large enough to include all prime factors of levels computed,
  int d,max=250000;
  cerr<<"Enter field: \n";
  cin >> d;
  Quad::field(d,max);
  Qideal N;
  int verbose=0, showforms=0;
  //cerr << "Verbose? "; cin>>verbose;
  //cerr << "Display newforms (1/0)? "; cin>>showforms;

  cerr<<"Enter level (ideal label or generator): \n";
  cin>>N;
  long normn = N.norm();
  if (verbose)
    cout << ">>>> Level " << ideal_label(N) <<" = "<<gens_string(N)<<", norm = "<<normn<<" <<<<" << endl;

  cerr<<"Enter prime p: \n";
  int p;
  cin>>p;
  if (verbose)
    cout << "(working modulo " << p << ")" << endl;

  int nforms, nprimes;
  cerr<<"Enter number of newforms and number of primes to check: \n";
  cin >> nforms >> nprimes;

  newforms nf(N,verbose>1);
  nf.createfromdata();
  if (verbose && showforms)
    nf.display();
  int nnf = nf.n1ds;
  int nap = nf.nap;
  if((nnf==0)||(nap==0))
    {
      if (verbose)
        cout<<"No newforms."<<endl;
      else
        cout<<"?"<<endl;
      exit(0);
    }
  if (verbose)
    cerr << "There are " << nnf << " newforms on file, with a_P for the first " << nap << " primes P (with index 0.." << (nap-1) << ")" << endl;
  // primes_needed will be a list of the (quad) primes for which
  // values of a_p will be input, and prime_indexes will be a list of
  // the index (starting at 0) of these in the standard list
  // quadprimes of all primes

  Quad P;
  vector<Quad> primes_needed(nprimes);
  vector<int> prime_indexes(nprimes);
  long maxnormp=0, maxip=0;
  long ip, np, ap, nform, kform;
  int computation_needed = 0;
  vector< vector<long> > apvecs_in(nforms);
  vector< vector<long> > apvecs_comp(nnf);
  for (nform=0; nform<nforms; nform++)
    {
      apvecs_in[nform].resize(nprimes);
    }
  for (nform=0; nform<nnf; nform++)
    {
      apvecs_comp[nform].resize(nprimes);
    }

  // Read in primes and ap:

  for(np=0; np<nprimes; np++)
    {
      // if(verbose)
      //   cerr << "Enter a prime p followed by "<<nforms<<" ap: "<<endl;
      cin >> P;
      for (nform=0; nform<nforms; nform++)
        {
          cin >> ap;
          apvecs_in[nform][np] = ap%p;
        }
      primes_needed[np] = P;
      prime_indexes[np] = ip = prime_index(P);
      if (ip>maxip)
        maxip = ip;
      long normp = quadnorm(P);
      if (normp>maxnormp)
        maxnormp = normp;
      // if(verbose)
      //   cerr << "p=" << P <<" (index "<<ip<<", norm "<<normp<<"): a_P = "<<ap<<endl;
        }       // end of prime loop

  // See whether we need to compute more ap:

  if(verbose)
    cerr << "Largest prime index (based at 0) for which we need ap is " << maxip <<"."<<endl;
  if (maxip>nap-1)
    {
      computation_needed = 1;
      if(verbose)
        {
          cout << "No stored ap for p = " << p << " which is has index " << ip << "(starting from 0): only " << nap << " a_P are on file." << endl;
          cout << "We'll have to compute the modular symbol space and eigenspaces in order to compute a-P" << endl;
        }
    }
  else
    {
      if(verbose)
        {
          cout << "All required a_P are on file (the last is for the " << maxip << "th prime, and we have " << nap << ")" << endl;
        }
    }

  if (computation_needed)   // Compute ap for these primes
    {

      if(verbose>1) cout << "Making homspace and bases..."<<flush;
      nf.makebases();
      if(verbose>1)
        cout << "done."<<endl;

      for(np=0; np<nprimes; np++)
        {
          P = primes_needed[np];
          Quadprime PP = Qideal(P).factorization().prime(0);
          vector<long> apv = nf.apvec(PP);
          if(verbose)
            cerr << "List of a_P for P="<<P<<": "<<apv<<endl;
          for (nform=0; nform<nnf; nform++)
            apvecs_comp[nform][np] = apv[nform]%p;
        }       // end of prime loop
    }
  else   // Extract the ap for these primes from the newform data
    {
      for(np=0; np<nprimes; np++)
        for (nform=0; nform<nnf; nform++)
          apvecs_comp[nform][np] = nf.nflist[nform].aplist[prime_indexes[np]]%p;
    }

  // Find each input ap list in the newforms aplists:

  for (nform=0; nform<nforms; nform++)
    {
      int n_matches = 0;
      if (verbose)
        cout << "Input Hecke eigenvalue data a_P (mod "<<p<<") = " << apvecs_in[nform] << endl;
      for (kform=0; kform<nnf; kform++)
        {
          if (verbose)
            cout << "Comparing with computed form "<<codeletter(kform)<<" with a_P (mod "<<p<<") = " << apvecs_comp[kform] << endl;
          if (apvecs_comp[kform] == apvecs_in[nform])
            {
              if (verbose)
                cout << " MATCHES form # "<<(kform+1)<<endl;
              else
                {
                  if (n_matches>0) cout << " ";
                  cout << codeletter(kform);
                }
              n_matches++;
            }
        }
      if (!verbose) cout <<endl;
      if (n_matches==0)
        {
          if (verbose)
            cout << " HAS NO MATCH" << endl;
          else
            cout << "?"<<endl;
        }
    }
}       // end of main()
