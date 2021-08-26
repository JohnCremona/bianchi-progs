#include "qidloop.h"
#include "homspace.h"
//#include "ratquads.h"
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
 long firstn, lastn; Quad n; int verbose, plusflag=1;
 cerr << "Verbose? "; cin>>verbose;
 cerr << "Plus space? "; cin>>plusflag;
#ifdef LOOPER
 int both_conj;
 cerr<<"Both conjugates? (0/1) "; cin >> both_conj;
 cerr<<"Enter first and last norm for level: ";
 cin >> firstn >> lastn;
 cerr<<endl;
 Qidealooper loop(firstn, lastn, both_conj); // not sorted within norm
  while( loop.not_finished() )
   {
     Qideal N = loop.next();
#else
     string Nlabel;
     while(cerr<<"Enter level ideal label: ", cin>>Nlabel, Nlabel[0]!='0')
       {
         Qideal N(Nlabel);
#endif
         long normn = N.norm();
         cout << ">>>> Level " << ideal_label(N)
              << " = " << gens_string(N)
              <<", norm = "<<normn<<" <<<<";
         if(verbose) cout<<endl;
         else cout << "\t";
         homspace h(N,plusflag, 0, verbose);  //level, plusflag, cuspidal, verbose
         cout << "Dimension = " << h.h1cuspdim() << endl;
       }
}
