#include "moddata.h"   // which includes quads.h
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
 cout << "Verbose? "; cin >> verbose;
#ifdef LOOPER
 cout<<"Enter first and last norm for Quad loop: ";
 cin >> firstn >> lastn;
 if(firstn<2) firstn=2;
 int both_conj=0;
 for(Quadlooper alpha(d,firstn,lastn,both_conj); alpha.ok(); ++alpha)
#else
 Quad alpha;
 while(cout<<"Enter level: ", cin>>alpha, alpha!=0)
#endif
   {
     n = makepos((Quad)alpha);  // makepos normalizes w.r.t. units
     long normn = quadnorm(n);
     cout << ">>>> Level " << ideal_label(n) <<" = ("<<n<<"), norm = "<<normn<<" <<<<\t";
     moddata md(n);
     if(verbose) cout<<"\n",md.display();
     int ok = md.check();
     if(ok)cout<<" OK\n"; else cout<<" !!!!!!CHECK FAILS!!!!!\n";
   }

}
