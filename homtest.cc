#include <eclib/subspace.h>
#include "homspace.h"   // which includes quads.h & moddata.h
#define LOOPER
#ifdef LOOPER
#include "looper.h"
#endif

int main ()
{
 int d,max=1000;
 cout << "Enter field (1,2,3,7,11): " << flush;  cin >> d;
 if(!((d==1)||(d==2)||(d==3)||(d==7)||(d==11)))
   {
     cout<<"field must be one of: 1, 2, 3, 7, 11!\n";
     exit(1);
   }
 Quad::field(d,max);
 long firstn, lastn; Quad n; int verbose, plusflag=1;
 cout << "Verbose? "; cin>>verbose;
 cout << "Plus space? "; cin>>plusflag;
#ifdef LOOPER
 int both_conj;
 cout<<"Both conjugates? (0/1) "; cin >> both_conj;
 cout<<"Enter first and last norm for Quad loop: ";
 cin >> firstn >> lastn;
 cout<<endl;
 if(firstn<2) firstn=2;
 for(Quadlooper alphaloop(d,firstn,lastn,both_conj); alphaloop.ok(); ++alphaloop)
   {
     Quad alpha = (Quad)alphaloop;
#else
 Quad alpha;
 while(cout<<"Enter level: ", cin>>alpha, alpha!=0)
   {
#endif
     n = makepos((Quad)alpha);  // makepos normalizes w.r.t. units
     long normn = quadnorm(n);
     cout << ">>>> Level ("<<n<<"), norm = "<<normn<<" <<<<";
     if(verbose) cout<<endl;
     else cout << "\t";
     homspace h(n,plusflag,verbose);  //level, plusflag, verbose
     cout << "Dimension = " << h.h1cuspdim() << endl;
   }

}
