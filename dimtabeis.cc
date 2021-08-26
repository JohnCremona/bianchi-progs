#include <eclib/subspace.h>
#include "homspace.h"   // which includes quads.h & moddata.h
#include "looper.h"

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
 cout << "# Table of dimensions of weight 2 Bianchi cuspidal and Eisenstein forms for GL2 over Q(sqrt(-"<<d<<"))" << endl;
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
 cout << "# Field\t Weight\t Level label\t ";
 cout << "dim(all)\t dim(cuspidal)\t dim(eisenstein)" << endl;

 for(Quadlooper alphaloop(firstn,lastn,both_conj); alphaloop.ok(); ++alphaloop)
   {
     Quad alpha = (Quad)alphaloop;
     n = makepos((Quad)alpha);  // makepos normalizes w.r.t. units
     cout << d << "\t2\t";                  // field and weight
     cout << ideal_label(n)<<"\t\t"; // level and norm
     homspace hplus(Qideal(n), plusflag, 0, verbose);  //level, plusflag, cuspidal, verbose
     int dimcusp = hplus.h1cuspdim();
     int dimall = hplus.h1dim();
     int dimeis = dimall-dimcusp;
     cout << dimall << "\t\t" << dimcusp << "\t\t" << dimeis << endl;
   }

}
