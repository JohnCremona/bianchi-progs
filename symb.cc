// FILE SYMB.CC: Implementations for symbols

#include "symb.h"

modsym::modsym(const symb& s) //Constructor for modsym, converting from symb:
{
 Quad c,d,h,x,y;
 c = s.c % ((s.N)->modulus);
 d = s.d % ((s.N)->modulus);
 h = quadbezout(c , d, x, y);
 a=RatQuad(-x , d/h);
 b=RatQuad( y , c/h);
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

vector<RatQuad> symbdata::alphalist; // list of a such that {a,oo} represent edge-orbits
void symbdata::init_geometry()
{
  alphalist.push_back(RatQuad(0));
  int d = Quad::d;
  Quad w = Quad(0,1);

  if (d < 19) return;

  // for d=19 and above, include w/2 and (w-1)/2

  alphalist.push_back(RatQuad(w,2)); // w/2
  alphalist.push_back(RatQuad(w-1,2)); // (w-1)/2
  if (d==19) return;

  // for d=43 and above, include (a+b*w)/3 for a in {-1,0,1}, b in {-1,1}
  alphalist.push_back(RatQuad(w,3));
  alphalist.push_back(RatQuad(w-1,3));
  alphalist.push_back(RatQuad(w+1,3));
  alphalist.push_back(RatQuad(-w,3));
  alphalist.push_back(RatQuad(1-w,3));
  alphalist.push_back(RatQuad(2-w,3));
  if (d==43) return;

  // for d=67, include 16 more
  Quad den(3,-1);
  alphalist.push_back(RatQuad(6+w,den));
  alphalist.push_back(RatQuad(-6-w,den));
  alphalist.push_back(RatQuad(2+w,den));
  alphalist.push_back(RatQuad(-2-w,den));

  den = quadconj(den);
  alphalist.push_back(RatQuad(7-w,den));
  alphalist.push_back(RatQuad(-7+w,den));
  alphalist.push_back(RatQuad(3-w,den));
  alphalist.push_back(RatQuad(-3+w,den));

  alphalist.push_back(RatQuad(w,4));
  alphalist.push_back(RatQuad(-w,4));
  alphalist.push_back(RatQuad(w-1,4));
  alphalist.push_back(RatQuad(1-w,4));
  alphalist.push_back(RatQuad(w+1,4));
  alphalist.push_back(RatQuad(-w-1,4));
  alphalist.push_back(RatQuad(w-2,4));
  alphalist.push_back(RatQuad(2-w,4));
}

symbdata::symbdata(const Quad &n) :moddata(n),specials()
{
  // initialise static data (depending only on the field)
  if (alphalist.size()==0)
    init_geometry();

  // cout << "In constructor symbdata::symbdata.\n";
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

