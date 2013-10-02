//  manin.cc 

#include <fstream>
#include <iomanip>
#include "manin.h"
#include "looper.h"

manin::manin(const Quad& n, int useolddata, int disp) 
:newforms(n,useolddata, disp)
{
  makeh1plus(); // In case it wasn't done in the newforms constructor
//cout << "In manin constructor, after constructing the newforms." << endl;
//  if(disp)newforms::display();
  allproj();  // Compute homspace::projcoord, so projcycle can be used
  if(n1ds)
    {
      easy=1;
      for (int i=0; i<n1ds; i++) if (nflist[i].pdot==0) easy=0;
      findq(); // Initialize nq, dq, qdotlist, initvec
    }
  //  else if(disp) cout<<"manin: No newforms!" << endl;
}

void manin::findq()
{
  //  Look for a rational q for which {q,infinity} is nontrivial
  int i,stillok,foundq=0, field=Quad::d;
  qdotlist.resize(n1ds);
  if (hmod) qdotlistinv.resize(n1ds);
  for (i=0; i<n1ds; i++)
    {
      pdotlist.push_back(nflist[i].pdot);
      if (hmod)
	pdotlistinv.push_back(invmod(nflist[i].pdot,hmod));
    }
  for (Quadlooper dl(field, 2, 1000, 1); dl.ok()&&!foundq; ++dl)
    { Quad d=(Quad)dl;
      if (coprime(d,modulus))
        { for (Quadlooper nl(field,1,quadnorm(d)/2,1); nl.ok() && !foundq; ++nl)
            { Quad n=(Quad)nl;
	      if (coprime(n,d))
                {   // found a candidate n/d
                  initvec=h1->projcycle(n,d);  //starts at 1
                  for (i=0, stillok=1; (i<n1ds) && stillok; i++)
                    { qdotlist[i]= nflist[i].pdot-nflist[i].dp0*initvec[i+1];
		      if(hmod) 
			{
			  qdotlist[i]=mod(qdotlist[i],h1->h1hmod());
			  if (qdotlist[i])
			    qdotlistinv[i]=invmod(qdotlist[i],hmod);
			}
                      stillok = (qdotlist[i]!=0);
                    }
                  foundq = stillok;
                  if (foundq) { nq=n; dq=d; }
                }
            }
        }
    }
}

void manin::getoneap(const Quad& p, int output, ofstream& out, int verbose)
{
  vec v;
  int i,ap,np; int good = ndiv(p,modulus);
  if(verbose) cout << "p = " << p << "\t"; 
  long maxap=(long)(2*sqrt((double)quadnorm(p))); // for validity check

  if ( easy)  //   Symbol {0,infinity} non-trivial in all cases
    { 
      if (!good)
        { 
          for (i=0; i<n1ds; i++)
            { 
	      int ip = find(plist.begin(),plist.end(),p)-plist.begin();
              ap = nflist[i].aplist[ip];
              if((ap>maxap)||(-ap>maxap))
                {
                  cout<<"Error:  eigenvalue "<<ap<<" for p="<<p
                      <<" for form # "<<(i+1)<<" is outside valid range "
                      <<-maxap<<"..."<<maxap<<endl;
                }
              if(verbose) cout<<setw(5)<<ap<<" ";
              if(output) out<<setw(5)<<ap;
            }
          if(verbose) cout << "   *** bad prime\n";
          if(output) out<<endl;
        }
      else
        {  // good prime
          v = h1->projmaninvector(p);   //starts at 1
          for (i=0; i<n1ds; i++)
            {
	      if (hmod)
		np = mod((v[i+1]*nflist[i].dp0) * pdotlistinv[i], hmod);
	      else
		np = ((v[i+1]*nflist[i].dp0) / pdotlist[i]);
	      ap = 1+quadnorm(p)-np;
              if((ap>maxap)||(-ap>maxap))
                {
                  cout<<"Error:  eigenvalue "<<ap<<" for p="<<p
                      <<" for form # "<<(i+1)<<" is outside valid range "
                      <<-maxap<<"..."<<maxap<<endl;
                }
              if(verbose)cout<<setw(5)<<ap<<" ";
              if(output) out<<setw(5)<<ap;
            }
          if(verbose) cout << endl;
          if(output)out<<endl;
        }
    }        // end of "if easy then..."
  else 
    { 
      if (!good)
        { 
          for (i=0; i<n1ds; i++)
            { 
	      int ip = find(plist.begin(),plist.end(),p)-plist.begin();
              ap = nflist[i].aplist[ip];
              if(verbose)cout<<setw(5)<<ap<<" ";
              if(output) out<<setw(5)<<ap;
            }
          if(verbose)cout << "   *** bad prime\n";
          if(output)out<<endl;
        }
      else
        {
          v = h1->newhecke(p,nq,dq) - (quadnorm(p)+1)*initvec;  //starts at 1
          for (i=0; i<n1ds; i++)
            { 
	      if (hmod)
		np = mod((v[i+1] * nflist[i].dp0) * qdotlistinv[i],hmod);
	      else
		np = ( (v[i+1] * nflist[i].dp0) / qdotlist[i] );
              ap = 1+quadnorm(p)- np;
              if(verbose)cout<<setw(5)<<ap<<" ";
              if(output) out<<setw(5)<<ap;
            }
          if(verbose)cout << endl;
          if(output)out<<endl;
        }
    }      // end of "if easy then ... else ..."
}    

void manin::getap(int first, int last, int output, string eigfile, int verbose)
{
  ofstream out;
  if(output) 
    {
      if(first>1) 
        {out.open(eigfile.c_str(),ios::app);}   //append
      else 
        {
          out.open(eigfile.c_str());
          out<<n1ds<<" "<<n2ds<<endl;
        }
    }
  if(n1ds>0)
    {
      if(output&&(first==1))out<<last<<endl<<endl;
      if (easy&&verbose)     //   Symbol {0,infinity} non-trivial in all cases
        cout << "L(f,1) nonzero for each form"<<endl; 
      vector<Quad>::const_iterator pr;
      for(pr=quadprimes.begin(); pr!=quadprimes.end(); pr++)
	{
	  long index = pr-quadprimes.begin()+1;
	  if ( (index>=first) && (index<=last) )
	    getoneap(*pr,output,out,verbose);
	}
    }
  //  else if(verbose) cout<<"getap: No newforms!"<<endl;
  if(output) out.close();
}
