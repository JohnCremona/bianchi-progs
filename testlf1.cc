#include <fstream>
#include "newforms.h"
#include "lf1.h"
//#define LOOPER
#ifdef LOOPER
#include "looper.h"
#endif

int main ()
{
 cout.precision(10);
 int d,maxpnorm=10000;
 cout << "Enter field: " << flush;  cin >> d;
 Quad::field(d,maxpnorm);
 Quad n; int verbose=0;
 cout << "Verbose? "; cin>>verbose;
#ifdef LOOPER
 long firstn, lastn;
 cout<<"Enter first and last norm for Quads: ";
 cin >> firstn >> lastn;
 int both_conj=0;
 for(Quadlooper alpha(firstn,lastn,both_conj); alpha.ok(); ++alpha)
#else
   Qideal N;
 while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
#endif
   {
     cout << ">>>> Level " << ideal_label(N) <<" = "<<gens_string(N)<<", norm = "<<N.norm()<<" <<<<" << endl;
     newforms nf(N,verbose);
     nf.read_from_file();
     nf.display();
     // int denom = nf.h1->h1denom();
     // if(denom!=1) cout << "Denom = " << denom << endl;
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
