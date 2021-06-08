#include <eclib/subspace.h>
#include "homspace.h"   // which includes quads.h & moddata.h
#include "looper.h"

// List of fields for which this has been implemented so far:
vector<int> fields = {1,2,3,7,11,19,43,67};

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
 int both_conj, plusflag;
 int dimplus, dimminus, dimall;
 cerr<<"Both conjugates? (0/1) "; cin >> both_conj;
 cerr<<"Plus space only? (0/1) "; cin >> plusflag;
 cerr<<"Enter first and last norm for Quad loop: ";
 cin >> firstn >> lastn;
 cerr<<endl;
 cout << "Field\t Weight\t Level\t Norm\t ";
 if (plusflag)
   cout << "dim(plus)" << endl;
 else
   cout << "dim(all)\t dim(plus)\t dim(minus)" << endl;

 for(Quadlooper alphaloop(d,firstn,lastn,both_conj); alphaloop.ok(); ++alphaloop)
   {
     Quad alpha = (Quad)alphaloop;
     n = makepos((Quad)alpha);  // makepos normalizes w.r.t. units
     long normn = quadnorm(n);
     cout << d << "\t2\t";                  // field and weight
     cout << "("<<n<<")\t "<<normn<<"\t\t"; // level and norm
     homspace hplus(n,1,0,0);  //level, plusflag, cuspidal, verbose
     dimplus = hplus.h1cuspdim();
     if (!plusflag)
       {
         homspace hall(n,0,0,0);  //level, plusflag, cuspidal, verbose
         dimall = hall.h1cuspdim();
         dimminus = dimall-dimplus;
         cout << dimall << "\t\t" << dimplus << "\t\t" << dimminus << endl;
       }
     else
       {
         cout << dimplus << endl;
       }
   }

}
