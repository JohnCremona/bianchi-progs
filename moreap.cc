// FILE moreap.cc: computes more ap for given level(s)

#include <fstream>
#include "newforms.h"

int main(void)
{
 cout << "Program moreap: for each level, assumes that the newforms file exists, and computes more Hecke eigenvalues.\n";
 cout << "---------------\n\n";
 long d, maxpnorm(150000);
 cerr << "Enter field: " << flush;  cin >> d;
 Quad::field(d,maxpnorm);
 Qideal N;
 int verbose=0, output, showeigs, showforms, lastp;
 cerr << "Verbose? "; cin>>verbose;
 cerr << "Output new eigs to file (1/0)? ";  cin>>output;
 cerr << "Output new eigs to screen (1/0)? "; cin>>showeigs;
 cerr << "Display newforms (1/0)? "; cin>>showforms;

 while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
   {
     string efilename = eigfile(N);
     cerr << "How many primes for Hecke eigenvalues? ";
     cin  >> lastp; cout << endl;
     cout << ">>>> Level " << ideal_label(N) <<" = "<<gens_string(N)<<", norm = "<<N.norm()<<" <<<<" << endl;
     newforms nf(N,verbose);
     int ok = nf.read_from_file();
     if (!ok)
       {
         cout<<"No newform data available for level "<<ideal_label(N)<<endl;
         continue;
       }
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
