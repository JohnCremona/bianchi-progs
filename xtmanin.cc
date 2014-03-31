// XTMANIN.CC -- version of TMANIN.CC which processes all old levels automatically.

#include <fstream>
#include "newforms.h"   // which includes quads.h & moddata.h & etc.

int main ()
{
 int d,max=1000;
 cout << "Enter field: " << flush;  cin >> d;
 Quad::field(d,max);
 Quad n; int verbose=0;
 int startp, stopp;
 cout << "Verbose? "; cin>>verbose;
// int plusflag=1;
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
     vector<Quad> alphadivs = posdivs(alpha);
     vector<Quad>::const_iterator betaj;
     for(betaj=alphadivs.begin(); betaj!=alphadivs.end(); betaj++)
       {
	 Quad beta = makepos(*betaj);
	 long normbeta = quadnorm(beta);
	 if((normbeta==1)||(normbeta==normalpha)) continue;

// See if the oldform file for level beta exists; if not, run that level!
	 string eigfilename = eigfile(beta);
	 ifstream data(eigfilename.c_str()); 
	 if (!data)
	   {
	     if(verbose)
	       cout << "No data file for beta = " << beta 
		 << "  so constructing newforms at that level..." << endl;
	     newforms nf(beta,0);
             nf.createfromscratch();
	     nf.getap(startp,stopp,output,eigfilename,0);
	   }
	 else 
	   {
	     if(verbose) 
	       cout << "Oldform data for level " << beta << " exists."<<endl;
	     data.close();
	   }
       }
     if(verbose) 
       cout<<"\nFinished processing old levels.  Constructing newforms at level " << alpha << endl <<endl;
     newforms nf(alpha,0,verbose);   
     string eigfilename = eigfile(alpha);
     nf.getap(startp,stopp,output,eigfilename,1);
   }
}
