// FILE SYMB.CC: Implementations for symbols

#include "symb.h"

modsym::modsym(const symb& s, int type) //Constructor for modsym, converting from symb:
{
 Quad c,d,h,x,y;
 c = s.c % ((s.N)->modulus);
 d = s.d % ((s.N)->modulus);
 h = quadbezout(c , d, x, y);
 // matrix [y, -x; c/h, d/h] has det=1 and lifts (c:d)
 c = c/h;
 d = d/h;
 b = RatQuad( y , c);
 if(type==0) // always true for Euclidean fields: apply to {0,oo}
   {
     a = RatQuad(-x , d);
   }
 else // apply to {alpha,oo} where alpha = alphas[type]
   {
     RatQuad alpha = alphas[type];
     a = (y*alpha-x) / (c*alpha+d);
   }
}

//Members of class symblist:

void symblist::display() const
{
  int i;
  vector<symb>::const_iterator s;
  for(i=0, s = symbols.begin(); s!=symbols.end(); i++, s++)
    cout<<i<<":\t"<< *s <<"\n";
}

void symblist::add(symb& s, int start)
{
 if (std::find(symbols.begin()+start, symbols.end(), s) == symbols.end())
   symbols.push_back(s);
}

int symblist::index(symb& s, int start) const
{
 vector<symb>::const_iterator si = std::find(symbols.begin()+start, symbols.end(), s);
 return (si==symbols.end()? -1: si-symbols.begin());
}

symb symblist::item(int n) const
{
  if ((n>=(int)symbols.size())||(n<0)) 
    {
      cerr<<"Error in symblist::item: index out of range!\n";
      exit(1);
    }
  return symbols[n];
}

//Member functions for class symbdata:

vector<RatQuad> alphas; // List of a such that {a,oo} represent edge-orbits
vector<Quad> alpha_denoms; // List of denominators of alphas.
int n_alphas;
vector<mat22> M_alphas;  // List of matrices M_a with det(M_a)=1 such that M_a(a)=oo.

// Given a,b,c,d with ad-bc=2 add alpha=-d/c and M_alpha=[a,b;c,d] to the global lists

void add_alpha(const Quad& a, const Quad& b, const Quad& c, const Quad& d)
{
  static const RatQuad inf(1,0);
  RatQuad alpha(-d,c);
  mat22 M_alpha(a,b,c,d);
  assert (M_alpha.det()==1);
  assert (M_alpha(alpha) == inf);
  alphas.push_back(alpha);
  M_alphas.push_back(M_alpha);
  n_alphas += 1;
}

void define_alphas()
{
  int d = Quad::d;

  // alphas =0 with denominator 1:

  alpha_denoms.push_back(1);
  add_alpha(0,-1,1,0);  // alpha[0] = 0

  if (d<19) return;

  Quad w = Quad::w;

  // alphas with denominator 2:

  alpha_denoms.push_back(2);
  Quad u = (d-3)/8;  // = 2, 5, 8, 20
  add_alpha(w-1,u,2,-w);  // alpha[1] = w/2
  add_alpha(w,u,2,1-w);   // alpha[2] = (w-1)/2

  if (d<43) return;

  alpha_denoms.push_back(3);
  u = -(d+5)/12;  // = -4, -6, -14 so w^2 = w+3*u+1
  add_alpha(1-w,u,3,-w);           // alpha[3] = w/3
  add_alpha(w-1,u,3,w);            // alpha[4] = -w/3
  add_alpha(w,u,3,w-1);            // alpha[5] = (1-w)/3
  add_alpha(-w,u,3,1-w);           // alpha[6] = (w-1)/3
  add_alpha(w+1,-(w+u+1),3,-1-w);  // alpha[7] = (w+1)/3
  add_alpha(-w-1,-(w+u+1),3,1+w);  // alpha[8] = -(w+1)/3

  if (d<67) return;

  if (d==67)
    {
      Quad den(3,-1);
      alpha_denoms.push_back(den);
      alphas.push_back(RatQuad(6+w,den));
      alphas.push_back(RatQuad(-6-w,den));
      alphas.push_back(RatQuad(2+w,den));
      alphas.push_back(RatQuad(-2-w,den));
      den = quadconj(den);
      alpha_denoms.push_back(den);
      alphas.push_back(RatQuad(7-w,den));
      alphas.push_back(RatQuad(w-7,den));
      alphas.push_back(RatQuad(3-w,den));
      alphas.push_back(RatQuad(w-3,den));

      alpha_denoms.push_back(4);
      alphas.push_back(RatQuad(w,4));
      alphas.push_back(RatQuad(-w,4));
      alphas.push_back(RatQuad(w-1,4));
      alphas.push_back(RatQuad(1-w,4));
      alphas.push_back(RatQuad(1+w,4));
      alphas.push_back(RatQuad(-1-w,4));
      alphas.push_back(RatQuad(w-2,4));
      alphas.push_back(RatQuad(2-w,4));
      n_alphas += 16;
      return;
    }
  //  cout << "define_alphas() not yet implemented for field "<<d<<endl;
}


Quad over(const Quad&a, long n)
{
  long x = real(a), y = imag(a);
  long u, v = roundover(y,n);
  if (Quad::t == 0)
    {
      u = roundover(x,n);
    }
  else
    {
      u = roundover(2*x+y-n*v, 2*n);
    }
  return Quad(u, v);
}

Quad over(const Quad&a, const Quad& b)
{
  return over(a*quadconj(b), quadnorm(b));
}

// Find a shift q such that (a/b)-q is close enough to an alpha, and
// set type to the index of that alpha.
Quad translation(const Quad& a, const Quad& b, int& type)
{
  Quad shift(0);
  RatQuad alpha(a,b);
  vector<RatQuad>::iterator t = std::find(alphas.begin(), alphas.end(), alpha);
  if (t!=alphas.end())
    {
      type = std::distance(alphas.begin(), t);
      cout<<"translation("<<a<<","<<b<<") trivially found shift="<<shift<<" and alpha="<<alpha<<" with index "<<type<<endl;
      return shift;
    }
  float maxdist = 0, dist;
  Quad d, u, q, r, best_d=0, best_q;
  long normb = quadnorm(b);
  cout<<"In translaton("<<a<<","<<b<<")"<<endl;
  for (vector<Quad>::iterator di=alpha_denoms.begin(); di!=alpha_denoms.end(); di++)
    {
      d = *di;
      cout<<"Trying d="<<d<<endl;
      u = a*d;
      q = over(u,b);
      r = u-q*b;
      dist = float(1-quadnorm(r)/normb)/quadnorm(d);
      cout<<"dist = "<<dist<<endl;
      if (dist>maxdist)
        {
          maxdist=dist;
          best_d=d;
          best_q=q;
        }
    }
  if (best_d==0)
    {
      cerr<<"translation("<<a<<","<<b<<") found no possible denominators!"<<endl;
    }
  shift = over(best_q,best_d);
  alpha = RatQuad(best_q-best_d*shift, best_d);
  t = std::find(alphas.begin(), alphas.end(), alpha);
  if (t==alphas.end())
    {
      type=-1;
      cerr<<"translation("<<a<<","<<b<<") computed shift="<<shift<<" and alpha="<<alpha<<" which is invalid"<<endl;
      exit(1);
    }
  else
    {
      type = std::distance(alphas.begin(), t);
      cout<<"translation("<<a<<","<<b<<") computed shift="<<shift<<" and alpha="<<alpha<<" with index "<<type<<endl;
    }
  return shift;
}

// index of alpha nearest to a/b, given that a is reduced mod b
// (or -1 if none is near enough: should not happen)
//
// Here alpha = r/s is near enough if M_alpha(a/b) has denominator
// smaller than b, which is iff N(s*a-r*b)<N(b).


int nearest_alpha(const Quad& a, const Quad& b)
{
  long normb = quadnorm(b);
  for (vector<RatQuad>::iterator alpha = alphas.begin(); alpha!=alphas.end(); alpha++)
    {
      Quad r = num(*alpha), s = den(*alpha);
      if (quadnorm(s*a-r*b) < normb)
        return std::distance(alphas.begin(), alpha);
    }
  cerr << "Error in nearest_alpha("<<a<<", "<<b<<"): none of the alphas is near enough"<<endl;
  return -1;
}

void symbdata::init_geometry()
{
  //  cout<<"In init_geometry() with d="<<Quad::d<<endl;
  define_alphas();
  //  cout<<n_alphas<<" alphas: "<<alphas<<endl;
}

symbdata::symbdata(const Quad &n) :moddata(n),specials()
{
  // cout << "In constructor symbdata::symbdata.\n";
  // initialise static data (depending only on the field)
  if (alphas.size()==0)
    init_geometry();
  nsymbx = nsymb*n_alphas;
  // cout << "nsymb2 = " << nsymb2 << "\n";
  dstarts[0]=dstarts[ndivs-1]=0;
//N.B. dlist includes d=1 at 0 and d=mod at end, which we don't want here
 if (nsymb2>0)
 { int ic,id,start; symb s;  Quad c,d;
   for (ic=1; (ic<ndivs-1)&&(specials.count()<nsymb2); ic++)
   { c=dlist[ic];
//cout<<"Looking for specials with c = " << c;
     dstarts[ic]=start=specials.count();
//cout<<" starting with number "<<start<<endl;
     for (id=1; (id<normod-phi)&&(specials.count()<nsymb2); id++)  
     { d = resnum(noninvlist[id]);
       if (coprime(d,c))
         {  s = symb(c,d,this);
//cout<<"s = "<<s<<": "; int oldnum=specials.count();
          specials.add(s,start);     //only adds it if not there already!
//if(oldnum<specials.count())cout<<"new one\n"; else cout<<"old one\n";
       }
     }     // end of d loop
    }      // end of c loop
   if (specials.count()<nsymb2)
   { cerr << "Problem: makesymbols found only " << specials.count() << " symbols ";
     cerr << "out of " << nsymb2 << "\n";
   }
 }
}
 
int symbdata::index2(const Quad& c, const Quad& d) const
{ int kd = code(d);
  if (kd>0)                // d invertible, with inverse res[kd]
    return numres(c*resnum(kd));   // (c:d) = (c*res[kd]:1)
  else
  { int kc = code(c);
    if (kc>0)              // (c:d) = (1:res[kc]*d) if c invertible
      {
	Quad e = resnum(kc);
	e*=d;
	return   normod-code(e);
      }
    else
    {
//cout<<"\nkc="<<kc<<" kd="<<kd;
     int start = dstarts[noninvdlist[-kc]];
     symb s(c,d,this);
//cout<<" start="<<start<<endl;
     return nsymb1+specials.index(s,start);  // should be(?): start);
    }
  }
}

symb symbdata::symbol(int i) const
{ if (i<normod) return symb(resnum(i),1,this);
  else if (i<nsymb1) return symb(1,resnum(noninvlist[i-normod]),this);
  else return specials.item(i-nsymb1);
}

void symbdata::display() const
{ moddata::display();
  cout << "dstarts = " << dstarts << endl;
  cout << "Number of special symbols = " << nsymb2 << endl;
  specials.display();
}

int symbdata::check(int verbose) const
{int moddataok = moddata::check(verbose);
 int i,j,ok=1; symb s;
 for (i=0; i<nsymb; i++)
 {
//  cout<<i<<": "<<flush;
  s = symbol(i);  
//  cout<<s<<": "<<flush;
  j = index(s); 
//  cout<<j<<endl;
  ok&=(i==j);
  if (i!=j) cout << i << "-->" << s << "-->" << j << endl;
 }
 if(verbose)
   {
     if (ok) 
       cout << "symbols check OK!"<<endl;
     else 
       cout << "symbols check found errors!"<<endl;
   }
 return ok&&moddataok;
}

matop::matop(const Quad& p, const Quad& n)
{
 if (p==n)
   {
     mats.resize(1, mat22(0,-1,n,0));
   }
 else
 if (div(p,n))   // W involution, 1 term
   {
      Quad u,v,a,b;
      for (u=1, v=n; div(p,v); v/=p, u*=p) ;
      quadbezout(u,v,a,b);
      mats.resize(1, mat22(u*a,-b,n,u));
   }
else                 // Hecke operator, p+1 terms
  {
    vector<Quad> resmodp = residues(p);
    vector<Quad>::const_iterator r=resmodp.begin();
    while(r!=resmodp.end())
      mats.push_back(mat22(1,*r++,0,p));
    mats.push_back(mat22(p,0,0,1));
  }
}

