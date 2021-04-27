#include <eclib/subspace.h>
#include "homspace.h"   // which includes quads.h & moddata.h
#define LOOPER
#ifdef LOOPER
#include "looper.h"
#endif

vector<int> valid_fields = {1,2,3,7,11,19,43};

int main ()
{
 int d,max=1000;
 cout << "Enter field (one of "<<valid_fields<<"): " << flush;  cin >> d;
 if (std::find(valid_fields.begin(), valid_fields.end(), d) == valid_fields.end())
   {
     cout<<"field must be one of: "<<valid_fields<<endl;
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
     cout << ">>>> Level " << ideal_label(n) <<" = ("<<n<<"), norm = "<<normn<<" <<<<";
     if(verbose) cout<<endl;
     else cout << "\t";
     homspace h(n,plusflag, 0, verbose);  //level, plusflag, cuspidal, verbose
     cout << "Dimension = " << h.h1cuspdim() << endl;
   }

}
