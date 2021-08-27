#include "qidloop.h"
#include "homspace.h"

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
 cout << "Table of dimensions of weight 2 Bianchi cusp forms for Q(sqrt(-"<<d<<"))" << endl;
 Quad::field(d,max);
 long firstn, lastn; Quad n;
 int both_conj, plusflag, verbose=0;
 cerr<<"Both conjugates? (0/1) "; cin >> both_conj;
 cerr<<"Plus space only? (0/1) "; cin >> plusflag;
 cerr<<"Enter first and last norm for Quad loop: ";
 cin >> firstn >> lastn;
 cerr<<endl;
 cout << "Field\tWeight\tLevel\tNorm\t";
 if (plusflag)
   cout << "dim(plus)" << endl;
 else
   cout << "dim(all)\tdim(plus)\tdim(minus)" << endl;

 Qidealooper loop(firstn, lastn, both_conj, 1); // sorted within norm
 while( loop.not_finished() )
   {
     Qideal N = loop.next();
     long normn = N.norm();
     cout << d << "\t2\t";                  // field and weight
     cout << ideal_label(N)<<"\t "<<normn<<"\t\t"; // level and norm
     homspace hplus(N,1,0, verbose);  //level, plusflag, cuspidal, verbose
     int dimplus = hplus.h1cuspdim();
     if (!plusflag)
       {
         homspace hall(N,0,0, verbose);  //level, plusflag, cuspidal, verbose
         int dimall = hall.h1cuspdim();
         int dimminus = dimall-dimplus;
         cout << dimall << "\t\t" << dimplus << "\t\t" << dimminus << endl;
       }
     else
       {
         cout << dimplus << endl;
       }
   }

}
