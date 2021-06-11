#include <fstream>
#include "newforms.h"   // which includes quads.h & moddata.h & etc.
//#define LOOPER
#ifdef LOOPER
#include "looper.h"
#endif

// List of fields for which this has been implemented so far:
vector<int> fields = {1,2,3,7,11,19,43,67};

int main ()
{
 int d,max=200000;
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
 cout << "Verbose? "; cin>>verbose;
  cout << "Which primes for Hecke eigenvalues (first#, last#)? ";
  cin >> startp >> stopp; cout << endl;
  if (stopp>nquadprimes)
    {
      cout<<"Reducing last# to "<<nquadprimes<<endl;
      stopp=nquadprimes;
    }
  int output=1;
  cout << "Output Hecke eigenvalues? (0/1) ";  cin >> output;
 int dimcusp, dimeis, dimall;
#ifdef LOOPER
 long firstn, lastn;
 int both_conj;
 cout<<"Both conjugates? (0/1) "; cin >> both_conj;
 cout<<"Enter first and last norm for Quads: ";
 cin >> firstn >> lastn;
 for(Quadlooper alpha(d,firstn,lastn,both_conj); alpha.ok(); ++alpha)
#else
 Quad alpha;
 while(cout<<"Enter level: ", cin>>alpha, alpha!=0)
#endif
   {
     cout<<endl;
     n = makepos((Quad)alpha);  // makepos normalizes w.r.t. units
     long normn = quadnorm(n);
     string efilename = eigfile(n);
     cout << ">>>> Level " << ideal_label(n) <<" = ("<<n<<"), norm = "<<normn<<" <<<<" << endl;
     newforms nf(n,verbose);
     nf.createfromscratch();
     // output lines as in dimtabeis:
     cout << "DIMTAB " << d << "\t2\t";           // field and weight
     cout << ideal_label(n)<<"\t\t"; // level and norm
     dimcusp = nf.h1->h1cuspdim();
     dimall = nf.h1->h1dim();
     dimeis = dimall-dimcusp;
     cout << dimall << "\t\t" << dimcusp << "\t\t" << dimeis << endl;

     //nf.display();
     nf.getap(startp,stopp,verbose);
     //cout << "After sort_lmfdb():\n";
     nf.sort_lmfdb();
     nf.display();
     if(output)
       {
	 cout << "Writing data to file "<<efilename<<"..."<<flush;
	 nf.output_to_file(efilename);
	 cout << "done." << endl;
       }
     cout<<"==========================================="<<endl;
   }
 cout<<endl;
}
