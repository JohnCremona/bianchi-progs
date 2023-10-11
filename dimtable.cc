#include "qidloop.h"
#include "homspace.h"
//#define MODP

int main ()
{
  long d, max(2000);
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
 long firstn, lastn; Quad n;
 int both_conj, plusflag;
 int dimall, dimcusp;

 cerr<<"Both conjugates? (0/1) "; cin >> both_conj;
 cerr<<"Plus space only? (0/1) "; cin >> plusflag;
 cerr<<"Enter first and last norm for Quad loop: ";
 cin >> firstn >> lastn;
 cerr<<endl;
 Quad::field(d,max);
 int n2r = Quad::class_group_2_rank>0;

 cout << "Table of dimensions of ";
 if (ch) cout<<"mod "<<ch<<" ";
 cout<<"level N homology over Q(sqrt(-"<<d<<"))" << endl;
 if (n2r)
   {
     cout<<"\n";
     cout<<"   (The GL2 dimensions include all unramified quadratic character subspaces, while"<<endl;
     cout<<"   the NGL2 dimensions are of just the trivial character subspace.)"<<endl;
   }
 cout << "\nField\tWeight\tLevel\tNorm\t";
 if (!plusflag)
   cout << "SL2" << "\t\t";
 cout << "GL2";
 if (n2r)
   cout << "\t\tNGL2";
 cout << "\n\t\t\t\t";
 if (!plusflag)
   cout<<"all cuspidal\t";
 cout<<"all cuspidal";
 if (n2r)
   cout<<"\tall cuspidal";
 cout << endl;

 Qidealooper loop(firstn, lastn, both_conj, 1); // sorted within norm
 while( loop.not_finished() )
   {
     Qideal N = loop.next();
     QUINT normn = N.norm();
     cout << d << "\t2\t";                  // field and weight
     cout << ideal_label(N)<<"\t "<<normn<<"\t"; // level and norm
     if (!plusflag)
       {
         homspace h1(N,0, 0, ch);  //level, plusflag, verbose, characteristic
         dimall = h1.h1dim();
         dimcusp = h1.h1cuspdim();
         cout << dimall << "   " << dimcusp << "\t\t";
       }
     homspace h1plus(N,1, 0, ch);  //level, plusflag, verbose, characteristic
     dimall = h1plus.h1dim();
     dimcusp = h1plus.h1cuspdim();
     cout << dimall << "   " << dimcusp;
     if (n2r)
       {
         dimall = h1plus.trivial_character_subspace_dimension(0);
         dimcusp = h1plus.trivial_character_subspace_dimension(1);
         cout << "\t\t" << dimall << "   " << dimcusp;
       }
     cout << endl;
   }

}
