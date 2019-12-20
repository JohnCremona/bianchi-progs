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
  cerr << "a list of input primes, checking the values agains input values.\n";
  cerr << "----------------------------------------------------------------\n\n";
  cerr << "Input format:\n";
  cerr << "<field> <level> <nforms> <nprimes>\n";
  cerr << "followed by nprimes lines each containing a prime's coefficients (2 integers) and nforms integers:\n";
  cerr << "<prime> <ap_1> <ap_2> ... <ap_nforms>\n";
  cerr << "----------------------------------------------------------------\n\n";
  // max is the maximum norm of precomputed primes.  It should be
  // large enough to include all prime factors of levels computed,
  int d,max=250000;
  cerr<<"Enter field (1, 2, 3, 7 or 11): \n";
  cin >> d;
  Quad::field(d,max);
  Quad n;
  int verbose=0, showforms=0;
  //cerr << "Verbose? "; cin>>verbose;
  //cerr << "Display newforms (1/0)? "; cin>>showforms;

  cerr<<"Enter level: \n";
  cin>>n;
  n = makepos(n);
  long normn = quadnorm(n);
  string efilename = eigfile(n);
  if (verbose)
    cout << ">>>> Level " << ideal_label(n) <<" = ("<<n<<"), norm = "<<normn<<" <<<<" << endl;

  int nforms, nprimes;
  cerr<<"Enter number of newforms and number of primes to check: \n";
  cin >> nforms >> nprimes;

  newforms nf(n,verbose>1);
  nf.createfromdata();
  if (verbose && showforms)
    nf.display();
  int nnf = nf.n1ds;
  int nap = nf.nap;
  if((nnf==0)||(nprimes==0))
    {
      if (verbose)
        cout<<"No newforms."<<endl;
      else
        cout<<"?"<<endl;
      exit(0);
    }

  Quad p;
  // primes_needed will be a list of the (quad) primes for which
  // values of a_p will be input, and prime_indexes will be a list of
  // the index (starting at 0) of these in the standard list
  // quadprimes of all primes

  vector<Quad> primes_needed(nprimes);
  vector<int> prime_indexes(nprimes);
  long normp, maxnormp=0, maxip=0;
  long ip, np, nform, kform;
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
      if(verbose)
	cerr << "Enter a prime p followed by "<<nforms<<" ap: "<<endl;
      cin >> p;
      for (nform=0; nform<nforms; nform++)
	cin >> apvecs_in[nform][np];
      primes_needed[np] = p;
      prime_indexes[np] = ip = prime_index(p);
      if (ip>maxip)
        maxip = ip;
      normp = quadnorm(p);
      if (normp>maxnormp)
        maxnormp = normp;
    }       // end of prime loop

  // See whether we need to compute more ap:

  if (maxip>nap)
    {
      computation_needed = 1;
      if(verbose)
        {
          cout << "No stored ap for p = " << p << " which is the " << ip << "'th prime (only " << nap << " ap are on file)" << endl;
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
          p = primes_needed[np];
          vector<long> apv = nf.apvec(p);
          if(verbose)
            cerr << "List of a_p for p="<<p<<": "<<apv<<endl;
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
      if (verbose)
        cout << "Hecke eigenvalue data ap = " << apvecs_in[nform];
      for (kform=0; (kform<nnf) &&not_found; kform++)
        {
          if (apvecs_comp[kform] == apvecs_in[nform])
            {
              if (verbose)
                cout << " MATCHES form # "<<(kform+1)<<endl;
              else
                cout << codeletter(kform) << endl;
              not_found = 0;
            }
        }
      if (not_found)
        {
          if (verbose)
            cout << "?"<<endl;
          else
            cout << " HAS NO MATCH" << endl;
        }
    }
}       // end of main()
