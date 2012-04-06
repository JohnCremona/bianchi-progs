// FILE SYMB.CC: Implementations for symbols

#include "symb.h"

modsym::modsym(const symb& s) //Constructor for modsym, converting from symb:
{
 Quad c,d,h,x,y;
 c = s.c % level::modulus; 
 d = s.d % level::modulus;
 h = quadbezout(c , d, x, y);
 a=RatQuad(-x , d/h);
 b=RatQuad( y , c/h);
}

//Members of class symblist:
void symblist::add(const symb& s, int start)
{
 if (index(s,start)==-1) 
 {
  if (num<maxnum) list[num++]=s; 
  else cerr << "Error in symblist::add: attempt to add too many symbols to list!\n";
 }
}

int symblist::index(const symb& s, int start) const
{
 int i,ans;
 for (i=start,ans=-1; ((i<num)&&(ans==-1)); i++) if (list[i]==s) ans=i;
 return ans;
}

symb symblist::item(int n) const
{
 if ((n>num)||(n<0)) 
 {cerr<<"Error in symblist::item: index out of range!\n";
  return symb(0,1);
 }
 else return list[n];
}

//Member functions for class symbdata:
symbdata::symbdata(const Quad &n) :moddata(n),specials(nsymb2)
{
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
       {  s = symb(c,d);
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
     symb s(c,d);
//cout<<" start="<<start<<endl;
     return nsymb1+specials.index(s,start);  // should be(?): start);
    }
  }
}

symb symbdata::symbol(int i) const
{ if (i<normod) return symb(resnum(i),1);
  else if (i<nsymb1) return symb(1,resnum(noninvlist[i-normod]));
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
     length=1;  mats=new mat22[1]; mats[0]=mat22(0,-1,n,0);
   }
 else
 if (div(p,n))   // W involution, 1 term
   {
      length=1;
      Quad u,v,a,b;
      for (u=1, v=n; div(p,v); v/=p, u*=p) ;
      quadbezout(u,v,a,b);
      mats = new mat22[1];
      mats[0]=mat22(u*a,-b,n,u);
   }
else                 // Hecke operator, p+1 terms
  {
    vector<Quad> resmodp = residues(p);
    int normp = resmodp.size();
    length=normp+1;
    mats = new mat22[length];
    mat22* m=mats;
    vector<Quad>::const_iterator r=resmodp.begin();
    while(r!=resmodp.end())
      *m++ = mat22(1,*r++,0,p);
    *m++ = mat22(p,0,0,1);
  }
}

