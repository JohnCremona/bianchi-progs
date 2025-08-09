// Read in newform data, recompute L/P, cuspidalfactor and matdot, and rewrite

#include <fstream>
#include "newforms.h"
#include "lf1.h"
#define LOOPER
#ifdef LOOPER
#include "looper.h"
#include "qidloop.h"
#endif

int main ()
{
 cout.precision(10);
 int d,maxpnorm=10000;
 cout << "Enter field: " << flush;
 cin >> d;
 Quad::field(d,maxpnorm);
 Quad n; int verbose=0;
 cout << "Verbose? ";
 cin>>verbose;
#ifdef LOOPER
 long firstn, lastn;
 cout<<"Enter first and last norm for Quads: ";
 cin >> firstn >> lastn;
 cerr<<endl;
 Qidealooper loop(firstn, lastn, 0, 1); // sorted within norm
 while( loop.not_finished() )
   {
     Qideal N = loop.next();
#else
     Qideal N;
     while(cerr<<"Enter level (ideal label or generator): ", cin>>N, !N.is_zero())
       {
#endif
     cout<<endl;
     string efilename = eigfile(N);
     string label = ideal_label(N);
     cout << ">>>> Level " << label <<" = "<<gens_string(N)<<", norm = "<<N.norm()<<" <<<<" << endl;
     newforms nf(N, verbose);
     nf.read_from_file();
     nf.makebases();
     if (verbose)
       nf.display();
     cout << "Writing data to file "<<efilename<<"..."<<flush;
     nf.output_to_file(efilename);
     cout << "done." << endl;
     Qideal Nbar = N.conj();
     if (N==Nbar)
       continue;
     string conj_efilename = eigfile(Nbar);
     string conj_label = ideal_label(Nbar);
     cout << "Conjugating data for level "<<label
          <<" into data for conjugate level "<<conj_label
          <<" and resorting"<<endl;
     newforms nfbar(nf);
     nfbar.conjugate();
     nfbar.sort_lmfdb();
     cout << "Writing conjugate data to file "
          <<conj_efilename<<"..."<<flush;
     nfbar.output_to_file(conj_efilename);
     cout << "done." << endl;
   }
 cout << endl;
}
