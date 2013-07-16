// FILE SYMB.H: Declarations for M-symbols, modular symbols

#if     !defined(_SYMB_H)
#define _SYMB_H      1       //flags that this file has been included

#include "moddata.h"
#include "ratquads.h"

class symb;
class modsym;
class symblist;
class symbdata;
class symbop;

class symb {
 private:
   Quad c,d;
   const moddata *N;
 public:
   symb() :c(0), d(0)   {;}
   symb(const Quad& cc, const Quad& dd, const moddata * iN) :c(cc), d(dd), N(iN)  {;}
   Quad cee() const        {return c;}
   Quad dee() const        {return d;}
   int operator==(const symb& s) const 
       {Quad sc=d*s.c, sd=c*s.d, m=N->modulus; 
	Quad a = sd-sc; 
	return div(m,a);}
   friend inline ostream& operator<< (ostream&s, const symb&);
   friend class symblist;
   friend class symbdata;
   friend class modsym;
   friend class symbop;
};

class modsym {
 private:
    RatQuad a,b;
 public:
    modsym() :a(0), b(0) {}
    modsym(const RatQuad& ra, const RatQuad& rb) :a(ra),b(rb) {}
    modsym(const symb&);                        //conversion
    RatQuad alpha() const {return a;}
    RatQuad  beta() const {return b;}
    friend ostream& operator<< (ostream& s, const modsym& m); //inline below
};

class symblist {
 private:
    symb *list;
    int num,maxnum;
 public:
    symblist(int n=0) {maxnum=n; num=0; list=new symb[n];}
    ~symblist() {delete[] list;}
    void add(const symb& s, int start=0);
    int index(const symb& s, int start=0) const;
    symb item(int n) const;
    void display() const 
      {for(int i=0; i<num; i++) cout<<i<<":\t"<<list[i]<<"\n";}
    int count() const {return num;}
};

class symbdata :public moddata {
 private:
    symblist specials;         // The list of "special" symbols
 public:
    symbdata(const Quad&);             // The constructor
    int index2(const Quad& c, const Quad& d) const;
    int index(const symb& s) const {return index2(s.c,s.d);}
    symb symbol(int i) const;
    void display() const;
    int check(int verbose=0) const;
};

class mat22 {  //2x2 matrix for linear fractional transformations
private: 
   Quad a,b,c,d;
public: 
   mat22() :a(0),b(0),c(0),d(0) {}
   mat22(const Quad ia, const Quad ib, const Quad ic, const Quad id) 
    :a(ia),b(ib),c(ic),d(id) {}
   RatQuad operator()(const RatQuad& q)const 
    {Quad n=num(q),de=den(q); 
     Quad n1=a*n, n2=b*de;
     Quad d1=c*n, d2=d*de;
     Quad nn=n1+n2, dd=d1+d2;
     return RatQuad(nn,dd);}
   friend class symbop;
};

class matop {  // formal sum of 2x2 matrices
 private: mat22* mats;
 public:  int length;
          matop(const Quad& p, const Quad& n); 
	  // constructor for hecke ops
          ~matop() {delete[] mats;length=0;}
          mat22 operator[](int i) const {return mats[i];}
};

class symbop 
{
private:
  symbdata* sd;
  mat22 m;
public:
  symbop(symbdata* sdi, const mat22& mm) :sd(sdi), m(mm) {}
  symbop(symbdata* sdi, const Quad& a, const Quad& b, const Quad& c, const Quad& d) :sd(sdi), m(a,b,c,d) {}
  int operator()(int i) const 
    {
      symb s = sd->symbol(i);
      Quad c(s.c), d(s.d);
      Quad c1=c*m.a, c2=d*m.c;
      Quad d1=c*m.b, d2=d*m.d;
      int ans = sd->index2(c1+c2,d1+d2);
      return ans;
    }
};
  
inline ostream& operator<< (ostream& s, const symb& sy)
{
   s << "(" << sy.c << ":" << sy.d << ")";
   return s;
}  

inline ostream& operator<< (ostream& s, const modsym& m)
{
   s << "{" << (m.a) << "," << (m.b) << "}";
   return s;
}  
 

#endif
