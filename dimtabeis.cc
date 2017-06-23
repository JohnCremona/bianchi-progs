#include <eclib/subspace.h>
#include "homspace.h"   // which includes quads.h & moddata.h
#include "looper.h"

int main ()
{
 int d,max=1000;
 cerr << "Enter field (1,2,3,7,11): " << flush;  cin >> d;
 if(!((d==1)||(d==2)||(d==3)||(d==7)||(d==11)))
   {
     cerr<<"field must be one of: 1, 2, 3, 7, 11!\n";
     exit(1);
   }
 cout << "# Table of dimensions of weight 2 Bianchi cuspidal and Eisenstein forms for GL2 over Q(sqrt(-"<<d<<"))" << endl;
 Quad::field(d,max);
 long firstn, lastn; Quad n;
 int both_conj;
 int dimcusp, dimeis, dimall;
 cerr<<"Both conjugates? (0/1) "; cin >> both_conj;
 // int plusflag;
 // cout<<"Plus space only? (0/1) "; cin >> plusflag;
 cerr<<"Enter first and last norm for Quad loop: ";
 cin >> firstn >> lastn;
 cout<<endl;
 cout << "# Field\t Weight\t Level label\t ";
 cout << "dim(all)\t dim(cuspidal)\t dim(eisenstein)" << endl;

 if(firstn<2) firstn=2;
 for(Quadlooper alphaloop(d,firstn,lastn,both_conj); alphaloop.ok(); ++alphaloop)
   {
     Quad alpha = (Quad)alphaloop;
     n = makepos((Quad)alpha);  // makepos normalizes w.r.t. units
     cout << d << "\t2\t";                  // field and weight
     cout << ideal_label(n)<<"\t\t"; // level and norm
     homspace hplus(n,1,0);  //level, plusflag, verbose
     dimcusp = hplus.h1cuspdim();
     dimall = hplus.h1dim();
     dimeis = dimall-dimcusp;
     cout << dimall << "\t\t" << dimcusp << "\t\t" << dimeis << endl;
   }

}
