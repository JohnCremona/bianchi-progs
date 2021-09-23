#include "qidloop.h"
#include "newforms.h"   // which includes quads.h & moddata.h & etc.
//#define LOOPER

int main ()
{
 int d,max=1000000;
 cerr << "Enter field: " << flush;  cin >> d;
 Quad::field(d,max);
 Quad n;
 long nap;
 cerr << "How many coefficients ap? "<<flush; cin>>nap;
#ifndef LOOPER
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
#endif
#ifdef LOOPER
 long firstn, lastn;
 cerr<<"Enter first and last norm for Quads: ";
 cin >> firstn >> lastn;
 cerr<<endl;
 int both_conj=1;

 Qidealooper loop(firstn, lastn, both_conj, 1); // sorted within norm
 while( loop.not_finished() )
   {
     Qideal N = loop.next();
#else
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

