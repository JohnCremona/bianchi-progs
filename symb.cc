// FILE SYMB.CC: Implementations for symbols

#include "symb.h"
#include "euclid.h"
#include "geometry.h"

// compute a matrix M = [y, -x; c, d] with det=1 lifting (c:d)
mat22 symb::lift_to_SL2() const
{
  // Two special cases: (c:1), (1:d) need no work:
  if (d==1) return mat22(1,0,c,1);
  if (d==-1) return mat22(1,0,-c,1);
  if (c==1) return mat22(0,-1,1,d);
  if (c==-1) return mat22(0,-1,1,-d);
  Quad x, y, sc = c % (N->modulus), sd = d % (N->modulus);
  Quad h = quadbezout(sc , sd, x, y);
  sc /= h;
  sd /= h;
  assert (y*sd+x*sc==1);
  return mat22(y,-x,sc,sd);
}

modsym::modsym(const symb& s, int type) //Constructor for modsym, converting from symb:
{
  mat22 U = s.lift_to_SL2();
  Quad n=1, d=0; // n/d=oo
  U.apply_left(n,d);
  b = RatQuad(n,d); // no need to reduce as n,d are coprime
  if(type==0) // always true for Euclidean fields: apply to {0,oo}
   {
     n=0; d=1; // 0
     U.apply_left(n,d);
     a = RatQuad(n,d); // =M(0) // no need to reduce as n,d are coprime
   }
 else // apply to {alpha,oo} where alpha = alphas[type]
   {
     mat22 M = M_alphas[type];
     n = -M.d; d = M.c; // alpha = r/s
     U.apply_left(n,d);
     a = RatQuad(n,d, 1); //  =M(alpha), reduced
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
  if ((s.N->normod)==1) return 0;
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

symbdata::symbdata(const Quad &n) :moddata(n),specials()
{
  // cout << "In constructor symbdata::symbdata.\n";
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
{
  if (normod==1) return 0;
  int kd = code(d);
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
  // int i;
  // for(i=0; i<nsymb; i++)
  //   cout<<i<<":\t"<< symbol(i) <<"\n";
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
