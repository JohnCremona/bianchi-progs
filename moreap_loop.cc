// FILE moreap.cc: computes more ap for given level(s)

#include <fstream>
#include "newforms.h"
#include "looper.h"

int main(void)
{
 cout << "Program moreap: for each level, assumes that the newforms file exists, and computes more Hecke eigenvalues.\n";
 cout << "---------------\n\n";
 int d,max=150000;
 cout << "Enter field: " << flush;  cin >> d;
 Quad::field(d,max);
 Quad n; int verbose=0, output=0, showeigs=0, showforms=0, lastp;
 cout << "Verbose? "; cin>>verbose;
 cout << "Output new eigs to file (1/0)? ";  cin>>output;
 cout << "Output new eigs to screen (1/0)? "; cin>>showeigs;
 cout << "Display newforms (1/0)? "; cin>>showforms;
 long firstn, lastn;
 int both_conj;
 cout<<"Both conjugates? (0/1) "; cin >> both_conj;
 cout << "How many primes for Hecke eigenvalues? ";
 cin  >> lastp; cout << endl;
 cout<<"Enter first and last norm for Quads: ";
 cin >> firstn >> lastn;
 if(firstn<2) firstn=2;
 for(Quadlooper alpha(firstn,lastn,both_conj); alpha.ok(); ++alpha)
   {
     n = makepos((Quad)alpha);
     long normn = n.norm();
     string efilename = eigfile(n);
     cout << ">>>> Level " << ideal_label(n) <<" = ("<<n<<"), norm = "<<normn<<" <<<<" << endl;
     newforms nf(n,verbose);
     nf.read_from_file();
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
