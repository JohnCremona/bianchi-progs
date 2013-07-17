#include <fstream>
#include "manin.h"   // which includes quads.h & moddata.h & etc.
#include "lf1.h"
//#define LOOPER
#ifdef LOOPER
#include "looper.h"
#endif

int main ()
{
 cout.precision(10);
 int d,max=1000;
 cout << "Enter field: " << flush;  cin >> d;
 Quad::field(d,max);
 long firstn, lastn; Quad n; int verbose=0, plusflag=1;
 int startp, stopp;
 cout << "Verbose? "; cin>>verbose;
// cout << "Plus space? "; cin>>plusflag;
#ifdef LOOPER
 cout<<"Enter first and last norm for Quads: ";
 cin >> firstn >> lastn;
 int both_conj=0;
 for(Quadlooper alpha(d,firstn,lastn,both_conj); alpha.ok(); alpha++)
#else
 Quad alpha;
 while(cout<<"Enter level: ", cin>>alpha, alpha!=0)
#endif
   {
     n = makepos((Quad)alpha);  long normn = quadnorm(n);
     cout << ">>>> Level " << ideal_label(n) <<" = ("<<n<<"), norm = "<<normn<<" <<<<" << endl;
     newforms nf(n,1,verbose);
     nf.display();
     int denom = nf.h1->h1denom();
     if(denom!=1) cout << "Denom = " << denom << endl;
     for(int i=0; i<nf.n1ds; i++)
       {
         cout << "\nForm number " << i+1 << ": " << endl;
         cout<<"Finding periods -- via L(f_chi)."<<endl;
         Quad lambda = nf.nflist[i].lambda;
         if(lambda!=Quad(1))
           cout<<"Using twisting prime lambda = "<<lambda<<endl;
         period_via_lf1chi per(&(nf.nflist[i]),verbose);
         double period = per.getperiod();
         double lf1chivalue = per.getlf1chivalue();
         rational ratio = per.getratio();
         cout << "Period = " << period << endl;
         cout << "L(f_chi,1) = " << lf1chivalue << endl;
         cout << "ratio = " << ratio << endl;
//       cout << "ratio direct from newform = " << nf.nflist[i].loverp << endl;
       }
   }
}
