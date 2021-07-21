#include <fstream>
#include "newforms.h"   // which includes quads.h & moddata.h & etc.
//#define LOOPER
#ifdef LOOPER
#include "looper.h"
#endif

int main ()
{
 int d,max=1000000;
 cerr << "Enter field: " << flush;  cin >> d;
 Quad::field(d,max);
 Quad n;
 long nap;
 cerr << "How many coefficients ap? "<<flush; cin>>nap;
#ifndef LOOPER
 vector<Quad>::const_iterator pr=quadprimes.begin();
 long np=0;
 cout << "Primes: "<<endl;
 while(np<nap)
   {
     Quad p = *pr++;
     np ++;
     cout << p << ", ";
   }
 cout << "..." << endl;
#endif
#ifdef LOOPER
 long firstn, lastn;
 cerr<<"Enter first and last norm for Quads: ";
 cin >> firstn >> lastn;
 cerr<<endl;
 int both_conj=1;
 for(Quadlooper alpha(firstn,lastn,both_conj); alpha.ok(); ++alpha)
#else
 Quad alpha;
 while(cerr<<"Enter level: ", cin>>alpha, cerr<<endl, alpha!=0)
#endif
   {
     n = makepos((Quad)alpha);
     ifstream data(eigfile(n).c_str());
     if(!data)
       {
         cout<<"No data for level " << ideal_label(n) << " = ("<<n<<"), norm = "<< quadnorm(n)<<endl;
       }
     else
       {
         newforms nf(n,0);
         nf.createfromdata();
         nf.list(nap);
       }
   }
}

