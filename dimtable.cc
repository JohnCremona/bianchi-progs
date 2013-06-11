#include <eclib/subspace.h>
#include "homspace.h"   // which includes quads.h & moddata.h
#include "looper.h"

int main ()
{
 int d,max=1000;
 cout << "Enter field (1,2,3,7,11): " << flush;  cin >> d;
 if(!((d==1)||(d==2)||(d==3)||(d==7)||(d==11)))
   {
     cout<<"field must be one of: 1, 2, 3, 7, 11!\n";
     exit(1);
   }
 cout << "Table of dimensions of weight 2 Bianchi cusp forms for Q(sqrt(-"<<d<<"))" << endl;
 Quad::field(d,max);
 long firstn, lastn; Quad n;
 int both_conj;
 int dimplus, dimminus, dimall;
 cout<<"Both conjugates? (0/1) "; cin >> both_conj;
 cout<<"Enter first and last norm for Quad loop: ";
 cin >> firstn >> lastn;
 cout<<endl;
 cout << "Field\t Weight\t Level\t Norm\t dim(all)\t dim(plus)\t dim(minus)" << endl;
 if(firstn<2) firstn=2;
 for(Quadlooper alphaloop(d,firstn,lastn,both_conj); alphaloop.ok(); ++alphaloop)
   {
     Quad alpha = (Quad)alphaloop;
     n = makepos((Quad)alpha);  // makepos normalizes w.r.t. units
     long normn = quadnorm(n);
     cout << d << "\t2\t";                  // field and weight
     cout << "("<<n<<")\t "<<normn<<"\t\t"; // level and norm
     homspace hall(n,0,0);  //level, plusflag, verbose
     dimall = hall.h1cuspdim();
     homspace hplus(n,1,0);  //level, plusflag, verbose
     dimplus = hplus.h1cuspdim();
     dimminus = dimall-dimplus;
     cout << dimall << "\t\t" << dimplus << "\t\t" << dimminus << endl;
   }

}
