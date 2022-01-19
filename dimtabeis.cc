#include "qidloop.h"
#include "homspace.h"
//#define MODP

int main ()
{
  long d, max(1000);
  cerr << "Enter field (one of "<<valid_fields<<"): " << flush;  cin >> d;
  if (!check_field(d))
    {
      cerr<<"field must be one of: "<<valid_fields<<endl;
      exit(1);
    }
  long ch=0;
#ifdef MODP
  cerr << "Enter characteristic (0 or prime): " << flush;  cin >> ch;
#endif
  cout << "# Table of dimensions of ";
  if (ch) cout<<"mod "<<ch<<" ";
  cout<<"weight 2 Bianchi cuspidal and Eisenstein forms for GL2 over Q(sqrt(-"<<d<<"))" << endl;
  Quad::field(d,max);
 long firstn, lastn; Quad n;
 int both_conj;
 int verbose=0;
 cerr<<"Both conjugates? (0/1) "; cin >> both_conj;
 int plusflag=1;
 // cout<<"Plus space only? (0/1) "; cin >> plusflag;
 cerr<<"Enter first and last norm for Quad loop: ";
 cin >> firstn >> lastn;
 cerr<<endl;
 cout << "# Field\tWeight\tLevel\t";
 cout << "dim(all)\tdim(cuspidal)\tdim(eisenstein)" << endl;

 Qidealooper loop(firstn, lastn, both_conj, 1); // sorted within norm
 while( loop.not_finished() )
   {
     Qideal N = loop.next();
     cout << "\t"<< d << "\t2\t";                  // field and weight
     cout << ideal_label(N)<<"\t"; // level
     homspace hplus(N, plusflag, 0, verbose, ch);  //level, plusflag, cuspidal, verbose
     int dimcusp = hplus.h1cuspdim();
     int dimall = hplus.h1dim();
     int dimeis = dimall-dimcusp;
     cout << dimall << "\t\t" << dimcusp << "\t\t" << dimeis << endl;
   }

}
