#include "qidloop.h"
#include "newforms.h"   // which includes quads.h & moddata.h & etc.
#define LOOPER

int main ()
{
  long d, max(200000);
  cerr << "Enter field (one of "<<valid_fields<<"): " << flush;  cin >> d;
  if (!check_field(d))
   {
     cerr<<"field must be one of: "<<valid_fields<<endl;
     exit(1);
   }
  long ch=0;
  cerr << "Enter characteristic p (prime): " << flush;  cin >> ch;
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
 long firstn, lastn;
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
     string efilename = eigfile(N, ch);
     string label = ideal_label(N);
     cout << ">>>> Level " << label <<" = "<<gens_string(N)<<", norm = "<<normN<<" <<<<" << endl;
     newforms nf(N,verbose, ch);
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
