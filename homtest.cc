#include "homspace.h"
#include "ratquads.h"
#define LOOPER
#ifdef LOOPER
#include "looper.h"
#endif

// List of fields for which this has been implemented so far:
vector<int> fields = {1,2,3,7,11,19,43,67,163, 5};

int main ()
{
 int d,max=1000;
 cerr << "Enter field (one of "<<fields<<"): " << flush;  cin >> d;
 if (!check_field(d, fields))
   {
     cerr<<"field must be one of: "<<fields<<endl;
     exit(1);
   }
 Quad::field(d,max);
 long firstn, lastn; Quad n; int verbose, plusflag=1;
 cerr << "Verbose? "; cin>>verbose;
 cerr << "Plus space? "; cin>>plusflag;
#ifdef LOOPER
 int both_conj;
 cerr<<"Both conjugates? (0/1) "; cin >> both_conj;
 cerr<<"Enter first and last norm for Quad loop: ";
 cin >> firstn >> lastn;
 cerr<<endl;
 for(Quadlooper alphaloop(firstn,lastn,both_conj); alphaloop.ok(); ++alphaloop)
   {
     Quad alpha = (Quad)alphaloop;
#else
 Quad alpha;
 while(cerr<<"Enter level: ", cin>>alpha, alpha!=0)
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
