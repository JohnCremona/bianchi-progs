#include "symb.h"   // which includes quads.h & moddata.h
#define LOOPER
#ifdef LOOPER
#include "looper.h"
#endif

int main ()
{
 int d,max=1000;
 cout << "Enter field: " << flush;  cin >> d;
 Quad::field(d,max);
 long firstn, lastn; Quad n; int verbose;
 cout << "Verbose? "; cin>>verbose;
#ifdef LOOPER
 cout<<"Enter first and last norm for Quad loop: ";
 cin >> firstn >> lastn;
 if(firstn<2) firstn=2;
 for(Quadlooper alpha(d,firstn,lastn); alpha.ok(); ++alpha)
#else
 Quad alpha;
 while(cout<<"Enter level: ", cin>>alpha, alpha!=0)
#endif
   {
     n = makepos((Quad)alpha);  // makepos normalizes w.r.t. units
     long normn = quadnorm(n);
     cout << ">>>> Level ("<<n<<"), norm = "<<normn<<" <<<<\t";
     symbdata sd(n);
     if(verbose) cout<<"\n",sd.display();
     int ok = sd.check();
     if(ok)cout<<" OK\n"; else cout<<" !!!!!!CHECK FAILS!!!!!\n";
   }

}
