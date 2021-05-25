#include <fstream>
#include "newforms.h"   // which includes quads.h & moddata.h & etc.
//#define LOOPER
#ifdef LOOPER
#include "looper.h"
#endif

int main ()
{
 int d,max=1000000;
 cout << "Enter field: " << flush;  cin >> d;
 Quad::field(d,max);
 Quad n; int verbose=0;
 cout << "Verbose? "; cin>>verbose;
#ifdef LOOPER
 long firstn, lastn;
 cout<<"Enter first and last norm for Quads: ";
 cin >> firstn >> lastn;
 int both_conj=1;
 for(Quadlooper alpha(d,firstn,lastn,both_conj); alpha.ok(); ++alpha)
#else
 Quad alpha;
 while(cout<<"Enter level: ", cin>>alpha, alpha!=0)
#endif
   {
     n = makepos((Quad)alpha);  long normn = quadnorm(n);
     newforms nf(n,verbose);
     nf.createfromdata();
     int nnf = nf.n1ds;
     if(verbose||nnf)
       {
         cout << "\n>>>> Level " << ideal_label(n) <<" = ("<<n<<"), norm = "<<normn<<" <<<<" << endl;
         nf.display();
       }
   }
}

