#include "qidloop.h"
#include "newforms.h"   // which includes quads.h & moddata.h & etc.
//#define LOOPER

int main ()
{
  long d, maxnorm(10000);
 cerr << "Enter field: " << flush;  cin >> d;
 Quad n; int verbose=0;
 cerr << "Verbose? "; cin>>verbose;
#ifdef LOOPER
 long firstn, lastn;
 cerr<<"Enter first and last norm for Quads: ";
 cin >> firstn >> lastn;
 int both_conj=1;
 if (lastn>maxnorm)  maxnorm = lastn;
#endif

 Quad::field(d,maxnorm);

#ifdef LOOPER
 Qidealooper loop(firstn, lastn, both_conj, 1); // sorted within norm
 while( loop.not_finished() )
   {
     Qideal N = loop.next();
#else
 Qideal N;
 while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
   {
#endif
     newforms nf(N,verbose);
     nf.createfromdata();
     int nnf = nf.n1ds;
     if(verbose||nnf)
       {
         cout << "\n>>>> Level " << ideal_label(N) <<" = "<<gens_string(N)<<", norm = "<<N.norm()<<" <<<<" << endl;
         nf.display();
       }
   }
}

