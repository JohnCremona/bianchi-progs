#include <fstream>
#include "manin.h"   // which includes quads.h & moddata.h & etc.
//#define LOOPER
#ifdef LOOPER
#include "looper.h"
#endif

int main ()
{
 int d,max=1000;
 cout << "Enter field: " << flush;  cin >> d;
 if(!((d==1)||(d==2)||(d==3)||(d==7)||(d==11)))
   {
     cout<<"field must be one of: 1, 2, 3, 7, 11!\n";
     exit(1);
   }
 Quad::field(d,max);
 Quad::displayfield(cout);
 Quad n; int verbose=0;
 int startp, stopp;
 cout << "Verbose? "; cin>>verbose;
  cout << "Which primes for Hecke eigenvalues (first#, last#)? ";
  cin >> startp >> stopp; cout << endl;
  int output=1;
  cout << "Output Hecke eigenvalues? (0/1) ";  cin >> output;
#ifdef LOOPER
 long firstn, lastn;
 int both_conj;
 cout<<"Both conjugates? (0/1) "; cin >> both_conj;
 cout<<"Enter first and last norm for Quads: ";
 cin >> firstn >> lastn;
 if(firstn<2) firstn=2;
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
     cout << ">>>> Level ("<<n<<"), norm = "<<normn<<" <<<<" << endl;
     manin machine(n,0,verbose);   
               // (level, use_old, verbose)
     machine.display();
     machine.getap(startp,stopp,output,efilename,1);
     cout<<"==========================================="<<endl;
   }
 cout<<endl;
}
