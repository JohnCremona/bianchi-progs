// FILE MODDATA.CC: Implementation of member functions for class moddata,
//                  QUADS version

#include "moddata.h"

level::level(const Quad& n, long neigs)
{
  modulus=makepos(n); normod=quadnorm(n);
  conj_modulus=quadconj(modulus);
  plist=pdivs(n); npdivs=plist.size();
  dlist=posdivs(n); ndivs=dlist.size();
  is_square=1;
  vector<Quad>::const_iterator pr;
  for(pr=plist.begin(); pr!=plist.end() && is_square; pr++)
    if (val(*pr,n)%2) is_square=0;
  is_Galois_stable = are_associate(n, quadconj(n));
  nap=neigs;
  primelist=plist;
  pr = quadprimes.begin();
  while(pr!=quadprimes.end() && (primelist.size()<(unsigned)nap))
    {
      Quad p = *pr++;
      if (ndiv(p,modulus))
	primelist.push_back(p);
    }

  long a = n.r, b = n.i, x, y;
  n0 =  bezout(a,b,x,y);
  n0m0 = normod / n0;
  long a0 = a/n0, b0 = b/n0;
  wmodz = ((Quad::n)*b0 + (Quad::t)*a0)*x - a0*y;
//cout<<"n0, n0m0, wmodz = "<<n0<<", "<<n0m0<<", "<<wmodz<<endl;
}

Quad level::reduce(const Quad& a) const  // return a reduced mod modulus
{
  return a - modulus*qdivi(a*conj_modulus, normod);
}

long level::numres(const Quad& a) const // what number is this residue a mod modulus?
{
  long y = a.i, x = a.r;
  long r = posmod(y,n0);
  long rdash = posmod((x+wmodz*(y-r)) , n0m0);
  return  rdash + r*n0m0;
}

Quad level::resnum(long i) const // which is the i'th residue mod modulus?
{
  std::ldiv_t qr = ldiv(i, n0m0);
  return reduce(Quad(qr.rem, qr.quot));
}

moddata::moddata(const Quad& n) :level(n)
{
 int i,nnoninv; Quad resi,x,y,d; long normd;
 phi=psi=normod;
 vector<Quad>::const_iterator pr = plist.begin();
 while(pr!=plist.end())
   {  
     Quad p = *pr++; 
     long np = quadnorm(p);
     phi*=(np-1); phi/=np;
     psi*=(np+1); psi/=np;
   }
 nsymb = psi;
 nsymb1 = 2*normod-phi;
 nsymb2 = nsymb-nsymb1;
 invlist.resize(normod);          //codes
 noninvlist.resize(normod-phi);   //list of non-units
 noninvdlist.resize(normod-phi);  //list of divisors for each nonunit
//cout<<normod<<" residues, "<<phi<<" invertible and "<<normod-phi<<" non.\n";
//cout<<nsymb<<" symbols, of which "<<nsymb2<<" are special"<<endl;
 nnoninv=0;
 for (i=0; i<normod; i++)            //set up codes
 { 
   resi=resnum(i);
//cout << "testing residue " << resi;
   d = quadbezout(resi,modulus,x,y); normd=quadnorm(d);
   if (normd==1) // this is an invertible residue
     {
       invlist[i] = numres(x); 
//cout << " --invertible, inverse = " << x << ", number "<<invlist[i]<<endl;
     }
   else // it is not invertible
     {
//cout << " --not invertible number "<<nnoninv<<", divisor = "<<d<<endl;
       invlist[i]=-nnoninv;
       noninvlist[nnoninv]=i;
       noninvdlist[nnoninv]=-1;
       if (normd<normod) 
	 noninvdlist[nnoninv]=find(dlist.begin(),dlist.end(),d)-dlist.begin();
       nnoninv++;
     }
 }
 if (ndivs>0) {dstarts.reserve(ndivs);}
}

moddata::~moddata()  { ; }
// No explicit action needed as destructors for members will be called anyway


void moddata::display() const
{
 cout << "Level = " << ideal_label(modulus) <<" = (" << modulus << ")" << endl;
 //if(is_square) cout << "** square level **" << endl;
 cout << "Number of symbols = " << nsymb << endl;
 cout << ndivs << " divisors: " << dlist << endl;
 cout << npdivs << " prime divisors: " << plist << endl;
 cout << "residues: "; 
 for(int i=0; i<normod; i++) {if(i) cout<<","; cout<<resnum(i);}
 cout<<endl;
 cout << "invlist: " << invlist << endl;
 cout << "noninvlist: " << noninvlist << endl;
 cout << "noninvdlist: " << noninvdlist << endl;
}

int moddata::check(int verbose) const //checks whether resnum & numres work OK
{
  int ok=1;
  for(long i=0; i<normod; i++)
    {
      Quad resi = resnum(i);
      if (verbose) cout << "Residue number " << i << " = " << resi << endl;
      ok&=(i==numres(resi));
    }
  if(verbose)
    {
      if(ok)
        cout << "residue numering OK!" << endl;
      else 
        cout << "residue numering NOT OK!" << endl;
    }
  return ok;
}

// brute force test whether a is a square of some element of reslist, mod m

int squaremod(const Quad& a, const Quad& m, const vector<Quad>& reslist)
{
  if (div(m,a)) return 0;
  vector<Quad>::const_iterator r=reslist.begin();
  while(r!=reslist.end())
    {
      Quad res=*r++; 
      if(div(m,res*res-a)) return +1;
    }
  return -1;
}

vector<int> makechitable(const Quad& lambda, const vector<Quad>& reslist)
{
  vector<int> chi;
  if(reslist.size()==1) 
    chi[0]=1;
  else
    {
      vector<Quad>::const_iterator r=reslist.begin();
      while(r!=reslist.end())
	chi.push_back(squaremod(*r++,lambda,reslist));
    }
  return chi;
}

double gauss(const Quad& m, const vector<Quad>& reslist)
{
//cout<<"Computing g(chi) for lambda = " << m << endl;
  double ans1=0; //double ans2=0;
  bigcomplex lrd = bigcomplex(m)*bigcomplex(to_bigfloat(0),sqrt(to_bigfloat(Quad::disc)));
  vector<Quad>::const_iterator r=reslist.begin();
  while(r!=reslist.end())
  {
    Quad res = *r++;
      double term1 = squaremod(res,m,reslist)*psif(bigcomplex(Quad(res))/lrd);
      ans1+=term1;      //    ans2+=term2;
  }
  return ans1;
}

string ideal_code(const Quad& N) // string code for a (principal) ideal N
{
  vector<long>H = HNF(N);
  stringstream s;
  s << quadnorm(N) << "." << H[1] << "." << H[2];
  return s.str();
}

// old versions:

string old_ideal_code(const Quad& N) // string code for a (principal) ideal N
{
  stringstream s;
  long r=N.re(), i=N.im();
  if(r<0)    s << "m";
  s << abs(r);
  s << "i";
  if(i<0)    s << "m";
  s << abs(i);
  s << char(0);
  return s.str();
}

string old_eigfile(const Quad& d)    //returns filename for eigs at level d
{
  stringstream s;
  s << getenv("NF_DIR");
  if (s.str().empty()) {s.clear(); s<<"./newforms";}
  s << "/Qsqrt-" << Quad::d << "/e";
  s << ideal_code(d);
  return s.str();
}

