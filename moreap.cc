// FILE moreap.cc: computes more ap for given level(s)

#include <fstream>
#include "newforms.h"

int main(void)
{
 cout << "Program moreap: for each level, assumes that the newforms file exists, and computes more Hecke eigenvalues.\n";
 cout << "---------------\n\n";
 int d,max=10000;
 cout << "Enter field: " << flush;  cin >> d;
 Quad::field(d,max);
 Quad n; int verbose=0, output, showeigs, showforms, lastp;
 cout << "Verbose? "; cin>>verbose;
 cout << "Output new eigs to file (1/0)? ";  cin>>output;
 cout << "Output new eigs to screen (1/0)? "; cin>>showeigs;
 cout << "Display newforms (1/0)? "; cin>>showforms;

 while (cout<<"Enter level: ", cin>>n, n!=0) {
     n = makepos(n);
     long normn = quadnorm(n);
     string efilename = eigfile(n);
     cout << "How many primes for Hecke eigenvalues? ";
     cin  >> lastp; cout << endl;
     cout << ">>>> Level " << ideal_label(n) <<" = ("<<n<<"), norm = "<<normn<<" <<<<" << endl;
     newforms nf(n,1,verbose);
     if (showforms) nf.display();

     int nnf = nf.n1ds;
     if(nnf==0)
       {
         cout<<"No newforms."<<endl;
         continue;
       }
     int nap = nf.nap;
     if(nap>=lastp)
       {
         cout<<"Already have "<<nap<<" eigenvalues on file, no need to compute more."<<endl;
       }
     else
       {
         cout << "Making homspace and bases..."<<endl;
         nf.makebases();
         cout << "About to start computing ap..."<<endl;
         nf.getap(nap+1,lastp,showeigs);
         cout << "...done."<<endl;
         if(output)
           {
             cout << "Writing data to file "<<efilename<<"..."<<flush;
             nf.output_to_file(efilename);
             cout << "...done." << endl;
           }
       }
 }       // end of while()
}       // end of main()
