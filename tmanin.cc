#include "qidloop.h"
#include "newforms.h"   // which includes quads.h & moddata.h & etc.
//#define LOOPER

// List of fields for which this has been implemented so far:
vector<long> fields = {1,2,3,7,11,19,43,67,163,23,31};

int main ()
{
  long d;
  QUINT max=200000;
  cerr << "Enter field (one of "<<fields<<"): " << flush;  cin >> d;
  if (!check_field(d, fields))
   {
     cerr<<"field must be one of: "<<fields<<endl;
     exit(1);
   }
 Quad::field(d,max);
 Quad::displayfield(cout);
 Quad n; int verbose=0;
 int startp, stopp;
 cerr << "Verbose? "; cin>>verbose;
  cerr << "Which primes for Hecke eigenvalues (first#, last#)? ";
  cin >> startp >> stopp; cerr << endl;
  int nQP = Quadprimes::list.size();
  if (stopp>nQP)
    {
      cerr<<"Reducing last# to "<<nQP<<", the number of Quadprimes initialized"<<endl;
      stopp=nQP;
    }
  int output=1;
  cerr << "Output Hecke eigenvalues? (0/1) ";  cin >> output;
#ifdef LOOPER
 QUINT firstn, lastn;
 int both_conj;
 cerr<<"Both conjugates? (0/1) "; cin >> both_conj;
 cerr<<"Enter first and last norm for Quads: ";
 cin >> firstn >> lastn;
 stringstream dimtabfilename;
 dimtabfilename << "dimtabeis."<<d<<"."<<firstn<<"-"<<lastn;
 ofstream dimtab(dimtabfilename.str().c_str());

 Qidealooper loop(firstn, lastn, both_conj, 1); // sorted within norm
 while( loop.not_finished() )
   {
     Qideal N = loop.next();
#else
     Qideal N;
     while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
       {
#endif
     cerr<<endl;
     QUINT normN = N.norm();
     string efilename = eigfile(N);
     string label = ideal_label(N);
     cout << ">>>> Level " << label <<" = "<<gens_string(N)<<", norm = "<<normN<<" <<<<" << endl;
     newforms nf(N,verbose);
     nf.createfromscratch();
#ifdef LOOPER
     int dimcusp, dimeis, dimall;
     // output lines as in dimtabeis:
     dimtab << d << "\t2\t";           // field and weight
     dimtab << label<<"\t\t"; // level and norm
     dimcusp = nf.h1->h1cuspdim();
     dimall = nf.h1->h1dim();
     dimeis = dimall-dimcusp;
     dimtab << dimall << "\t\t"
            << dimcusp << "\t\t"
            << dimeis << endl;
#endif
     //nf.display();
     nf.getap(startp,stopp,verbose);
     //cout << "After sort_lmfdb():\n";
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
#ifdef LOOPER
 dimtab.close();
#endif
}
