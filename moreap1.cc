// FILE moreap1.cc: computes individual ap for given level(s)

#include <fstream>
#include "newforms.h"

int main(void)
{
 cout << "Program moreap1: for given field and level, assumes that the newforms file exists, and computes more individual Hecke eigenvalues.\n";
 cout << "---------------\n\n";
 int d,max=10000;
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
     newforms nf(n,verbose);
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
     long p=1; int sig;
     while(p)
       {
         cout << "Enter a rational prime p (0 to finish): "<<endl;
         cin >> p;
         if(p<2) break;
         vector<Quad> pilist = Quad::primes_above(p, sig);
         vector<Quad>::const_iterator pij;
         for(pij=pilist.begin(); pij!=pilist.end(); pij++)
           {
             Quad pi = *pij;
             vector<long> apv = nf.apvec(pi);
             cout << "ap for "<<pi<<": "<<apv<<endl;
           }
       }       // end of prime loop
 }          // end of level loop
}       // end of main()
