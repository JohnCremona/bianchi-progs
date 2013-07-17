#include <fstream>
#include "manin.h"   // which includes quads.h & moddata.h & etc.
#define LOOPER
#ifdef LOOPER
#include "looper.h"
#endif

int main ()
{
 int d,max=1000;
 cout << "Enter field: " << flush;  cin >> d;
 Quad::field(d,max);
 long firstn, lastn; Quad n; int verbose=0;
 cout << "Verbose? "; cin>>verbose;
#ifdef LOOPER
 cout<<"Enter first and last norm for Quads: ";
 cin >> firstn >> lastn;
 if(firstn<2) firstn=2;
 int both_conj=1;
 for(Quadlooper alpha(d,firstn,lastn,both_conj); alpha.ok(); ++alpha)
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
       }
}

