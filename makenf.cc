#include "qidloop.h"
#include "newforms.h"   // which includes quads.h & moddata.h & etc.

//#define LOOPER // this is set in Makefile to build makenf_loop,
                 // which loops over a range of norms and skips any level for which
                 // data already exists.

int main ()
{
  long d, max(200000);
  cerr << "Enter field (one of "<<valid_fields<<"): " << flush;  cin >> d;
  if (!check_field(d))
   {
     cerr<<"field must be one of: "<<valid_fields<<endl;
     exit(1);
   }
 Quad::field(d,max);
 Quad::displayfield(cout);
 Quad n; int verbose=0;
 int nap;
 cerr << "Verbose? "; cin>>verbose;
  cerr << "How many primes for Hecke eigenvalues? ";
  cin >> nap; cerr << endl;
  int nQP = Quadprimes::list.size();
  if (nap>nQP)
    {
      cerr<<"Reducing from "<<nap<<" to "<<nQP<<", the number of Quadprimes initialized"<<endl;
      nap=nQP;
    }
  int output=1;
  cerr << "Output Hecke eigenvalues? (0/1) ";  cin >> output;
#ifdef LOOPER
 long firstn, lastn;
 int both_conj;
 cerr<<"Both conjugates? (0/1) "; cin >> both_conj;
 cerr<<"Enter first and last norm for Quads: ";
 cin >> firstn >> lastn;

 // // Open file for output of dimension table
 // stringstream dimtabfilename;
 // dimtabfilename << getenv("NF_DIR");
 // if (dimtabfilename.str().empty())
 //   {
 //     dimtabfilename.clear();
 //     dimtabfilename<<"./newforms";
 //   }
 // dimtabfilename << "/dimtabeis."<<d<<"."<<firstn<<"-"<<lastn;
 // ofstream dimtab(dimtabfilename.str().c_str());

 Qidealooper loop(firstn, lastn, both_conj, 1); // sorted within norm
 while( loop.not_finished() )
   {
     Qideal N = loop.next();
#else
     Qideal N;
     while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
       {
#endif
     cout<<endl;
     QUINT normN = N.norm();
     string efilename = eigfile(N);
     string label = ideal_label(N);
     cout << ">>>> Level " << label <<" = "<<gens_string(N)<<", norm = "<<normN<<" <<<<" << endl;
     newforms nf(N,verbose);
#ifdef LOOPER // skip this level if we already have a newforms file
     if(nf.read_from_file())
       {
         if (verbose)
           cout<<"data file "<<efilename<<" already exists, skipping this level"<<endl;
         continue;
       }
#endif
     nf.find();

// #ifdef LOOPER
//      int dimcusp, dimeis, dimall;
//      // output lines as in dimtabeis:
//      dimtab << d << "\t2\t";           // field and weight
//      dimtab << label<<"\t\t"; // level and norm
//      dimcusp = nf.h1->h1cuspdim();
//      dimall = nf.h1->h1dim();
//      dimeis = dimall-dimcusp;
//      dimtab << dimall << "\t\t"
//             << dimcusp << "\t\t"
//             << dimeis << endl;
// #endif

     // So far the newforms may include some "fake rationals" so don't
     //display yet
     // nf.display();

     // Compute the required number of ap. This also deletes any false
     // rationals and sets various data.

     nf.getap(1, nap, verbose);
     nf.sort_lmfdb();
     nf.display(verbose);

     if(output)
       {
	 cout << "Writing data to file "<<efilename<<"..."<<flush;
	 nf.output_to_file(efilename);
	 cout << "done." << endl;
       }
     cout<<"==========================================="<<endl;
   }
 cout<<endl;
// #ifdef LOOPER
//  dimtab.close();
// #endif
}
