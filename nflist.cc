#include "qidloop.h"
#include "newforms.h"   // which includes quads.h & moddata.h & etc.
//#define LOOPER

scalar modulus = default_modulus<scalar>();

int main ()
{
  long d, maxnorm(10000); // must be at least as many as #eigs on file
 cerr << "Enter field: " << flush;  cin >> d;
 Quad n;
 long nap;
 cerr << "How many coefficients ap? "<<flush; cin>>nap;
#ifdef LOOPER
 long firstn, lastn;
 cerr<<"Enter first and last norm for levels: ";
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
     auto pr=Quadprimes::list.begin();
     long np=0;
     cout << "Primes: "<<endl;
     while(np<(nap<0?25:nap))
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
     //cout<<"Level "<<label(N)<<" = "<<N<<": "<<flush;
     newforms nf(N, modulus, 0);
     if(nf.read_from_file())
       {
         nf.list(nap);
         // Uncomment the next two lines if you want to rewrite the newforms files
         // string efilename = eigfile(N);
         // nf.output_to_file(efilename);
       }
     else
       cout << "No data file for level " << label(N) << endl;
   }
}

