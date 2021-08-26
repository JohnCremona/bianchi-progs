// FILE moreap1.cc: computes individual ap for given level(s)

#include <fstream>
#include "newforms.h"

int main(void)
{
 cout << "Program moreap1: for given field and level, assumes that the newforms file exists, and computes more individual Hecke eigenvalues.\n";
 cout << "---------------\n\n";
 int d,max=150000;
 cout << "Enter field: " << flush;  cin >> d;
 Quad::field(d,max);
 Quad n; int verbose=0, showforms=1;
 cout << "Verbose? "; cin>>verbose;
 //cout << "Display newforms (1/0)? "; cin>>showforms;

 while (cout<<"Enter level: \n", cin>>n, n!=0) {
     n = makepos(n);
     long normn = quadnorm(n);
     string efilename = eigfile(n);
     cout << ">>>> Level " << ideal_label(n) <<" = ("<<n<<"), norm = "<<normn<<" <<<<" << endl;
     newforms nf(Qideal(n),verbose);
     nf.createfromdata();
     if (showforms) nf.display();

     int nnf = nf.n1ds;
     if(nnf==0)
       {
         cout<<"No newforms."<<endl;
         continue;
       }
     cout << "Making homspace and bases..."<<flush;
     nf.makebases();
     cout << "done."<<endl;
     long p=1;
     while(p)
       {
         cout << "Enter a rational prime p (0 to finish): "<<endl;
         cin >> p;
         if(p<2) break;
         vector<Quadprime> Plist = Quadprimes_above(p);
         vector<Quadprime>::const_iterator Pi;
         for(Pi=Plist.begin(); Pi!=Plist.end(); ++Pi)
           {
             Quadprime P = *Pi;
             vector<long> apv = nf.apvec(P);
             cout << "ap for "<<P.gen()<<": "<<apv<<endl;
           }
       }       // end of prime loop
 }          // end of level loop
}       // end of main()
