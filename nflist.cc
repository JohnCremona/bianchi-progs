#include "qidloop.h"
#include "newforms.h"   // which includes quads.h & moddata.h & etc.
//#define LOOPER

int main ()
{
  long d, maxnorm(10000); // must be at least as many as #eigs on file
 cerr << "Enter field: " << flush;  cin >> d;
 Quad n;
 long nap;
 cerr << "How many coefficients ap? "<<flush; cin>>nap;
#ifdef LOOPER
 long firstn, lastn;
 cerr<<"Enter first and last norm for Quads: ";
 cin >> firstn >> lastn;
 cerr<<endl;
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
 if (Quad::class_number==1)
   {
     vector<Quadprime>::const_iterator pr=Quadprimes::list.begin();
     long np=0;
     cout << "Primes: "<<endl;
     while(np<nap)
       {
         Quadprime p = *pr++;
         np ++;
         cout << p.gen() << ", ";
       }
     cout << "..." << endl;
   }
 Qideal N;
 while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
   {
#endif
     string datafilename = eigfile(N);
     ifstream data(datafilename.c_str());
     //cout<<"Opening data file "<<datafilename<<endl;
     if(!data)
       {
         cout<<"No data for level " << ideal_label(N) << " = "<<gens_string(N)<<", norm = "<< N.norm()<<endl;
       }
     else
       {
         //cout<<"Level "<<ideal_label(N)<<" = "<<N<<": "<<flush;
         newforms nf(N,0);
         nf.createfromdata();
         nf.list(nap);
       }
   }
}

