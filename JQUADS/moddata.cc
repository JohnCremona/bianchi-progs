// FILE moddata.cc: Implementation of member functions for class moddata,
//                  Q(sqrt(-5)) version

#include "moddata.h"
#include "looper.h"

//long level::npdivs, level::ndivs, level::nap, level::normod;
//Quad level::modulus;
//Quadlist level::plist, level::dlist, level::primelist;

set_of_residues::set_of_residues(const Qideal& n)
{
  c=n.cee();
  ac=n.ay()*c;
  bc=n.bee()*c;
  normod=n.norm();
  resl = Quadlist(normod);    // array to hold a rep of minl norm for each res
  int* foundflags = new int[normod];
  long i;
  for (i=0; i<normod; i++) { foundflags[i]=0; }
  long tofind = normod;
  Quadlooper ql(0,0,1,1); // endless, ie "upper limit" 0 (won't ever call ok())
  while (tofind>0)
    {
      Quad q=(Quad)ql;
      i = numres(q);
      if (foundflags[i]==0)
	{ resl[i]=q;
	  foundflags[i]=1;
	  tofind--;
// cout << "Using "<<q<<" of norm "<< q.norm()<<" for ";
// cout << box_resnum(i) << " of norm " << box_resnum(i).norm() << endl;
	}
      ql++;
    }
  delete[] foundflags;
}

long set_of_residues::numres(const Quad& alpha) const
{                                // what number is residue alpha mod modulus?
  long r = posmod(alpha.imag(),c);
  long s = (alpha.imag() - r)/c;
  long t = alpha.real() - s * bc;
  long ans = r*ac + posmod(t,ac);
  return ans;
}

Quad set_of_residues::box_resnum(long i) const 
                                      // which is the i'th residue mod modulus?
{                                     // resnum(0)=0  -- this is important!
  long rdash = i % ac;
  long r = (i-rdash) / ac;
  return Quad(rdash,r);
}

Quad set_of_residues::resnum(long i) const 
                    // which is the smallest rep of i'th residue mod modulus?
{                   // resnum(0)=0  -- this is important!
  return resl(i);
}

int set_of_residues::check(int verbose) const
{                                 // checks whether resnum & numres work OK
                                  // checks both crude and minl norm versions
  int ok=1, ok2=1;
  for(long i=0; i<normod; i++)
    {
      Quad resi = box_resnum(i);
//cout << "Crude (box) residue number " << i << " = " << resi << endl;
      ok&=(i==numres(resi));
      resi = resnum(i);
      ok2&=(i==numres(resi));
    }
  if(verbose)
    {
      if(ok)
	cout << "box residue numbering num ---> res ---> num OK!" << endl;
      else 
	cout << "box residue numbering NOT OK!" << endl;
      if(ok2)
	cout << "minl residue numbering num ---> res ---> num OK!" << endl;
      else 
	cout << "minl residue numbering NOT OK!" << endl;
    }
  return ok&&ok2;
}

void set_of_residues::display(ostream&s, int verbose) const
{
  if (verbose)
    for (long i=0; i<normod; i++)
      { s<< "Res("<<i<<") = "<< resnum(i) << " of norm "<< resnum(i).norm();
	s<< "  to represent standard residue ";
	s<< box_resnum(i) << " of norm " << box_resnum(i).norm() <<endl;
      }
  else
    {
      s << "Residues: "; 
      for(long i=0; i<normod; i++) {if(i) s <<", "; s<<resnum(i);}
    }
  s <<endl;
}

quotient_ring::quotient_ring(const Qideal& n, long phi_n=0): set_of_residues(n)
{
  if (phi_n==0)
    { cerr << "Warning: quotient_ring not implemented for unknown phi";
      exit(1);   // abort program instantly
    }
  invlist = vector<long>(normod);              // codes
  noninvlist = vector<long>(normod-phi_n);     // list of non-units
  Quad resi,x,y;
  long nnoninv = 0;
  for (long i=0; i<normod; i++) {invlist[i]=0;}
  for (long i=0; i<normod; i++)                 //set up codes
    { 
      if (invlist[i] == 0)
	{
	  resi=resnum(i);
	  //cout << "testing residue " << resi;
	  if (comax(resi,n,x,y))
	    {                            // resi is an invertible residue
	      if (i==0)                  // 0 invertible in trivial ring!
		{ invlist[0]=0; }
	      else
		invlist[ invlist[i] = numres(x/resi)] = i;
	      // cout << " --invertible, inverse = " << x/resi << ", number "
	      //  <<invlist[i]<<endl;
	    }
	  else                           // resi is not invertible
	    {
	      //cout << " --not invertible number " << nnoninv << endl;
	      invlist[i]=-nnoninv;
	      noninvlist[nnoninv++]=i;
	    }
	}
    }
}

void quotient_ring::display(ostream& s, int verbose) const
{
  set_of_residues::display(s,verbose);
  s << "invlist: " << invlist << endl;
  s << "noninvlist: " << noninvlist << endl;
}

/*
level::level(const Quad& n, long neigs)
{
  plist=pdivs(n); npdivs=plist.length;
  dlist=posdivs(n); ndivs=dlist.length;
  nap=neigs;
  primelist = Quadlist(nap);
  for(long i=0; i<npdivs; i++) primelist[i]=plist(i);
  for(Quadvar p(quadprimes); p.ok()&&(i<nap); p++) 
    if (ndiv((Quad)p,modulus)) primelist[i++]=(Quad)p;
}
*/

moddata::moddata(const Qideal& n) : modulus(n), nfactn(n)
{
//  modulus = n;
//  nfactn = new prime_factn(n);
  long normod=n.norm();
  phi=psi=normod;
  for(long i=0; i<(nfactn.num_primes()); i++)
    { 
      Qideal p = nfactn.prime(i);
      long np = p.norm();
      phi/=np; phi*=(np-1);
      psi/=np; psi*=(np+1);
    }
}

void moddata::display(ostream& s) const
{
  s << "Level = " << modulus << endl;
  nfactn.display(s);
  s << "phi = " << phi << endl;
  s << "psi = " << psi << endl;
}


// some other general-purpose functions (could be in quads itself)

/*
int squaremod(const Quad& a, const Quad& m, const Quadlist& reslist)
{
  if (div(m,a)) return 0;
  for(Quadvar rvar(reslist); rvar.ok(); rvar++) 
    {
      Quad r=(Quad)rvar; if(div(m,r*r-a)) return +1;
    }
  return -1;
}

int* makechitable(const Quad& lambda, const Quadlist& reslist)
{
  long normlambda = reslist.getlength();
  int* ans = new int[normlambda]; long i=0;
  if(normlambda==1) ans[0]=1;
  else for(Quadvar r(reslist); r.ok(); r++, i++) 
          ans[i]=squaremod((Quad)r,lambda,reslist);
  return ans;
}

double gauss(const Quad& m, const Quadlist& reslist)
{
//cout<<"Computing g(chi) for lambda = " << m << endl;
  double ans1=0; //double ans2=0;
  complex lrd = complex(m)*complex(0,Field::rootdisc);
  for(Quadvar r(reslist); r.ok(); r++)
  {
      double term1 = squaremod((Quad)r,m,reslist)*psif(complex(Quad(r))/lrd);
//    double term2 = squaremod((Quad)r,m,reslist)*psig(complex(Quad(r))/lrd);
//    cout << "term1 = " << term1;
//    cout << ", term2 = " << term2 << endl;
      ans1+=term1;      //    ans2+=term2;
  }
//cout << "Returning g(chi) = " << ans1 << endl;
//cout << "(ans2 = " << ans2 << ")" << endl;
  return ans1;
}
*/

// END OF FILE moddata.cc
