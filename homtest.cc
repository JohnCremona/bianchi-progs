#include "qidloop.h"
#include "homspace.h"
#define LOOPER
//#define CHECK_CONJUGATE
//#define MODP

int main ()
{
  long d, max(1000);
  cerr << "Enter field: " << flush;  cin >> d;
 if (!check_field(d))
   {
     cerr<<"field must be one of: "<<valid_fields<<endl;
     exit(1);
   }
  long ch=0;
#ifdef MODP
  cerr << "Enter characteristic (0 or prime): " << flush;  cin >> ch;
#endif
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
         QUINT normn = N.norm();
         cout << ">>>> Level " << ideal_label(N)
              << " = " << gens_string(N)
              <<", norm = "<<normn<<" <<<<";
         if(verbose)
           cout<<endl;
         else
           cout << "\t";
         homspace h(N,plusflag, verbose, ch);  //level, plusflag, cuspidal, verbose
         cout << "Dimension";
#ifdef MODP
         if (ch) cout << " (mod "<<ch<<")";
#endif
         cout << " = " << h.h1cuspdim() << endl;
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
