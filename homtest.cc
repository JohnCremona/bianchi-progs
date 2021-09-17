#include "qidloop.h"
#include "homspace.h"
#define LOOPER

// List of fields for which this has been implemented so far:
vector<int> fields = {1,2,3,7,11,19,43,67,163, 5, 23, 31};

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
 int verbose, plusflag=1;
 cerr << "Verbose? "; cin>>verbose;
 cerr << "Plus space? "; cin>>plusflag;
#ifdef LOOPER
 long firstn, lastn;
 int both_conj;
 cerr<<"Both conjugates? (0/1) "; cin >> both_conj;
 cerr<<"Enter first and last norm for level: ";
 cin >> firstn >> lastn;
 cerr<<endl;
 Qidealooper loop(firstn, lastn, both_conj, 1); // sorted within norm
 while( loop.not_finished() )
   {
     Qideal N = loop.next();
#else
     Qideal N;
     while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
       {
#endif
         long normn = N.norm();
         cout << ">>>> Level " << ideal_label(N)
              << " = " << gens_string(N)
              <<", norm = "<<normn<<" <<<<";
         if(verbose)
           cout<<endl;
         else
           cout << "\t";
         homspace h(N,plusflag, 0, verbose);  //level, plusflag, cuspidal, verbose
         cout << "Dimension = " << h.h1cuspdim() << endl;
#ifdef CHECK_CONJUGATE
         if (!h.check_conjugate())
           {
             Qideal Nconj = N.conj();
             cout<<"**************************************************"<<endl;
             if (N.is_Galois_stable())
               cout<<"Conjugation map not an isomorphism at level "<<ideal_label(N)<<"="<<N<<endl;
             else
               cout<<"Inconsistency between level "<<ideal_label(N)<<"="<<N<<" and its conjugate "<<ideal_label(Nconj)<<"="<<Nconj<<endl;
             cout<<"**************************************************"<<endl;
           }
#endif
       }
}
