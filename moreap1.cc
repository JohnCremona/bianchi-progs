// FILE moreap1.cc: computes individual ap for given level(s)

#include <fstream>
#include "newforms.h"

scalar modulus = default_modulus<scalar>();

int main(void)
{
 cout << "Program moreap1: for given field and level, assumes that the newforms file exists, and computes more individual Hecke eigenvalues.\n";
 cout << "---------------\n\n";
 long d, maxpnorm(150000);
 cerr << "Enter field: " << flush;  cin >> d;
 Quad::field(d,maxpnorm);
 Qideal N;
 int verbose=0, showforms=1;
 cerr << "Verbose? "; cin>>verbose;
 //cerr << "Display newforms (1/0)? "; cin>>showforms;

 while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
   {
     string efilename = eigfile(N);
     cout << ">>>> Level " << label(N) <<" = "<<gens_string(N)<<", norm = "<<N.norm()<<" <<<<" << endl;
     newforms nf(N, modulus, verbose);
     nf.read_from_file();
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
         for ( auto& P : Plist)
           {
             vector<long> apv = nf.apvec(P);
             cout << "ap for "<<P<<": "<<apv<<endl;
           }
       }       // end of prime loop
 }          // end of level loop
}       // end of main()
