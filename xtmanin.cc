// XTMANIN.CC -- version of TMANIN.CC which processes all old levels automatically.

#include <fstream>
#include "manin.h"   // which includes quads.h & moddata.h & etc.

int main ()
{
 int d,max=1000;
 cout << "Enter field: " << flush;  cin >> d;
 Quad::field(d,max);
 long firstn, lastn; Quad n; int verbose=0, plusflag=1;
 int startp, stopp;
 cout << "Verbose? "; cin>>verbose;
// cout << "Plus space? "; cin>>plusflag;
  cout << "Which primes for Hecke eigenvalues (first#, last#)? ";
  cin >> startp >> stopp; cout << endl;
  int output=1; 
  cout << "Output Hecke eigenvalues? (0/1) ";  cin >> output;
 Quad alpha;

 while(cout<<"Enter level: ", cin>>alpha, alpha!=0)
   {
     alpha = makepos((Quad)alpha);  // makepos normalizes w.r.t. units
     long normalpha = quadnorm(alpha);
     cout << ">>>> Level ("<<alpha<<"), norm = "<<normalpha<<" <<<<" << endl;

     if(verbose) cout << "Checking that old levels have been processed..." << endl;
     Quadlist alphadivs = posdivs(alpha);
     for(long j=1; j<alphadivs.length; j++)
       {
	 Quad beta = makepos(alphadivs[j]);
	 long normbeta = quadnorm(beta);
	 if((normbeta==1)||(normbeta==normalpha)) continue;

// See if the oldform file for level beta exists; if not, run that level!
	 char* eigfilename = eigfile(beta);
	 ifstream data(eigfilename); 
	 if (!data)
	   {
	     if(verbose)
	       cout << "No data file for beta = " << beta 
		 << "  so constructing newforms at that level..." << endl;
	     manin machine(beta,0,0);   
	     machine.getap(startp,stopp,output,eigfilename,0);
	   }
	 else 
	   {
	     if(verbose) 
	       cout << "Oldform data for level " << beta << " exists."<<endl;
	     data.close();
	   }
	 delete [] eigfilename;
       }
     if(verbose) 
       cout<<"\nFinished processing old levels.  Constructing newforms at level " << alpha << endl <<endl;
     manin machine(alpha,0,verbose);   
     char* eigfilename = eigfile(alpha);
     machine.getap(startp,stopp,output,eigfilename,1);
     delete [] eigfilename;
   }
}
