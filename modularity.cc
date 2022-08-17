// FILE modularity.cc: checks individual ap against input values for given level(s)

#include <fstream>
#include "newforms.h"

long prime_index(const Quadprime& P)
{
  return find(Quadprimes::list.begin(), Quadprimes::list.end(), P) - Quadprimes::list.begin();
}

int main(void)
{
  cerr << "Program modularity: for given field and level, reads newforms from\n";
  cerr << "file (which should exist already), and computes Hecke eigenvalues for\n";
  cerr << "a list of input primes, checking the values agains input values.\n";
  cerr << "----------------------------------------------------------------\n\n";
  cerr << "Input format:\n";
  cerr << "<field> <level> <nforms> <nprimes>\n";
  cerr << "followed by nprimes lines each containing a prime (2 integers if principal, or a label) and nforms integers:\n";
  cerr << "<prime> <ap_1> <ap_2> ... <ap_nforms>\n";
  cerr << "----------------------------------------------------------------\n\n";
  // max is the maximum norm of precomputed primes.  It should be
  // large enough to include all prime factors of levels computed,
  long d, max(250000);
  cerr<<"Enter field: \n";
  cin >> d;
  Quad::field(d,max);
  Qideal N;
  int verbose=0, showforms=0;
  // cerr << "Verbose? "; cin>>verbose;
  // cerr << "Display newforms (1/0)? "; cin>>showforms;

  cerr<<"Enter level (ideal label or generator): \n";
  cin>>N;
  QUINT normn = N.norm();
  if (verbose)
    cout << ">>>> Level " << ideal_label(N) <<" = "<<gens_string(N)<<", norm = "<<normn<<" <<<<" << endl;

  int nforms, nprimes;
  cerr<<"Enter number of newforms and number of primes to check: \n";
  cin >> nforms >> nprimes;

  newforms nf(N,verbose>1);
  nf.read_from_file_or_find();
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
    {
      cerr << "There are " << nnf << " newforms on file";
      if (Quad::class_group_2_rank)
        cerr << " (up to unramified quadratic twist)";
      cerr << ", with ap for the first " << nap << " primes (with index 0.." << (nap-1) << ")" << endl;
    }
  Quadprime P;
  // primes_needed will be a list of the prime ideals P for which
  // values of a_P will be input, and prime_indexes will be a list of
  // the index (starting at 0) of these in the standard list of primes

  vector<Quadprime> primes_needed(nprimes);
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
          apvecs_in[nform][np] = ap;
        }
      primes_needed[np] = P;
      prime_indexes[np] = ip = prime_index(P);
      if (ip>maxip)
        maxip = ip;
      long normp = I2long(P.norm());
      if (normp>maxnormp)
        maxnormp = normp;
      // if(verbose)
      //   cerr << "P = " << P <<" (index "<<ip<<", norm "<<normp<<"): ap = "<<ap<<endl;
        }       // end of prime loop

  // See whether we need to compute more ap:

  if(verbose)
    cout << "Largest prime index (based at 0) for which we need ap is " << maxip <<"."<<endl;
  if (maxip>nap-1)
    {
      computation_needed = 1;
      if(verbose)
        {
          cout << "No stored ap for P = " << P << " which is has index " << ip << "(starting from 0): only " << nap << " ap are on file." << endl;
          cout << "We'll have to compute the modular symbol space and eigenspaces in order to compute ap" << endl;
        }
    }
  else
    {
      if(verbose)
        {
          cout << "All required ap are on file (the last is for the " << maxip << "th prime, and we have " << nap << ")" << endl;
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
          vector<long> apv = nf.apvec(P);
          if(verbose)
            cerr << "List of a_P for P = "<<P<<": "<<apv<<endl;
          for (nform=0; nform<nnf; nform++)
            apvecs_comp[nform][np] = apv[nform];
        }       // end of prime loop
    }
  else   // Extract the ap for these primes from the newform data
    {
      for(np=0; np<nprimes; np++)
        for (nform=0; nform<nnf; nform++)
          apvecs_comp[nform][np] = nf.nflist[nform].aplist[prime_indexes[np]];
    }

  // Find each input ap list in the newforms aplists:

  for (nform=0; nform<nforms; nform++)
    {
      int not_found = 1;
      vector<long> apvec_in = apvecs_in[nform];
      if (verbose)
        cout << "Input Hecke eigenvalue data ap = " << apvec_in << endl;
      if (Quad::class_group_2_rank==0)
        {
          for (kform=0; (kform<nnf) &&not_found; kform++)
            {
              string code = codeletter(kform);
              if (verbose)
                cout << "Comparing with computed form "<<code<<", ap = " << apvecs_comp[kform] << endl;
              if (apvecs_comp[kform] == apvec_in)
                {
                  if (verbose)
                    cout << "  - input data MATCHES newform "<< code <<endl;
                  else
                    cout << code << endl;
                  not_found = 0;
                }
            }
        }
      else // even class number must look at unramified quadratic twists
        {
          for (kform=0; (kform<nnf) &&not_found; kform++)
            {
              QUINT D = nf.nflist[kform].CMD;
              vector<QUINT> twists = disc_factors_mod_D((D==ZERO?ONE:D));
              int ntwists = twists.size();
              vector<long> apvec = apvecs_comp[kform];
              for (int jtwist = 0; jtwist<ntwists; jtwist++)
                {
                  int nform = kform*ntwists+jtwist;
                  string code = codeletter(nform);
                  vector<long> apvec_twist = apvec;
                  // Twist the ap:
                  auto Pi = primes_needed.begin();
                  auto aPi = apvec_twist.begin();
                  for (; aPi!=apvec_twist.end(); ++Pi, ++aPi)
                    (*aPi) *= Pi->genus_character(twists[jtwist]);

                  if (verbose)
                    cout << "  - comparing with computed form "<<code<<", ap = " << apvec_twist << endl;
                  if (apvec_twist == apvec_in)
                    {
                      if (verbose)
                        cout << "  - input data MATCHES newform "<< code <<endl;
                      else
                        cout << code << endl;
                      not_found = 0;
                    }
                }
            }
        }
      if (not_found)
        {
          if (verbose)
            cout << " - has NO matching form" << endl;
          else
            cout << "?"<<endl;
        }
    }
}       // end of main()
