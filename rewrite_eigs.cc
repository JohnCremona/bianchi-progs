// FILE rewrite_eigs.cc: rewrite newforms/<field_label>/* files for a range of levels

// Only to be used when the schema for these files changes, so the
// data is read under the old schema and written under the new.

#include <fstream>
#include "qidloop.h"
#include "newforms.h"

int main(void)
{
 // cout << "Program rewrite_eigs: for each level, reads an existing newforms file, and rewrites it.\n";
 // cout << "--------------------\n\n";
 long d, max(150000);
 cerr << "Enter field: " << flush;  cin >> d;
 Quad::field(d,max);
 Qideal N;
 int verbose=0, showforms=0;
 if (showforms)
   {
     cerr << "Display newforms (1/0)? ";
     cin>>showforms;
   }
 if (verbose)
   {
     cerr << "Verbose (1/0)? ";
     cin>>verbose;
   }

 long firstn, lastn;
 cerr<<"Enter first and last norm for levels: ";
 cin >> firstn >> lastn;
 cerr<<endl;
 if(verbose)
   cout<<"Rewriting newforms files for field "<<field_label()<<", level norms "<<firstn<<"..."<<lastn<<endl;

 Qidealooper loop(firstn, lastn, 1, 1); // both conjugates, sorted within norm
 while( loop.not_finished() )
   {
     Qideal N = loop.next();
     string efilename = eigfile(N);
     if (verbose)
       cout << ">>>> Level " << ideal_label(N) <<" = "<<gens_string(N)<<", norm = "<<N.norm()<<" <<<<" << endl;
     newforms nf(N,verbose);
     int ok = nf.read_from_file();
     if (!ok)
       {
         cout<<"No newform data available for level "<<ideal_label(N)<<endl;
         continue;
       }
     int nnf = nf.n1ds;
     if(nnf==0)
       {
         if (showforms)
           cout<<"No newforms."<<endl;
         continue;
       }
     nf.list(); // this triggers computation of bc,cm
     if (verbose)
       cout << "Writing new data to file "<<efilename<<"..."<<flush;
     nf.output_to_file(efilename);
     if (verbose)
       cout << "...done." << endl;
   }       // end of while() loop over levels
}       // end of main()
