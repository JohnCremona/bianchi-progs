// FILE SYMB.H: Declarations for M-symbols, modular symbols

#if     !defined(_SYMB_H)
#define _SYMB_H      1       //flags that this file has been included

#include "moddata.h"

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
   friend class mat22;
};

class modsym {
 private:
    RatQuad a,b;
 public:
    modsym() :a(0), b(0) {}
    modsym(const RatQuad& ra, const RatQuad& rb) :a(ra),b(rb) {}
    modsym(const symb&, int type=0);                        //conversion
    RatQuad alpha() const {return a;}
    RatQuad  beta() const {return b;}
    friend ostream& operator<< (ostream& s, const modsym& m); //inline below
};

class symblist {
 private:
  vector<symb> symbols;
 public:
    symblist() {;}
    void add(symb& s, int start=0);
    int index(symb& s, int start=0) const;
    symb item(int n) const;
    void display() const;
  int count() const {return symbols.size();}
};

class symbdata :public moddata {
private:
  symblist specials;         // The list of "special" symbols
protected:
  long nsymbx;               // number of (symb,type) pairs, = nsymb*n_alphas
public:
  symbdata(const Quad&);             // The constructor
  static void init_geometry();       // sets alpha list etc, depending only on the field
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
 // left action on r/s as column vector, changing in place:
  void apply_left(Quad& r, Quad& s) const
  {
    Quad t = a*r+b*s;
    s = c*r+d*s;
    r = t;
  }
  RatQuad operator()(const RatQuad& q)const
  {
    Quad r = num(q), s = den(q);
    apply_left(r, s);
    return RatQuad(r,s);
  }
  // right action on (c:d) symbols as row vectors, changing in place
  void apply_right(Quad& sc, Quad& sd) const
  {
    Quad t = a*sc + c*sd;
    sd = b*sc + d*sd;
    sc = t;
  }
  symb operator()(const symb& s) const
  {
    Quad sc = s.c, sd=s.d;
    apply_right(sc, sd);
    return symb(sc, sd, s.N);
  }
  Quad det() const {return a*d-b*c;}
  friend ostream& operator<< (ostream& s, const mat22& m); // inline below
  friend class symbop;
};

class matop {  // formal sum of 2x2 matrices
private:
  vector<mat22> mats;
 public:
  matop(const Quad& p, const Quad& n);  // constructor for hecke ops
  mat22 operator[](int i) const {return mats[i];}
  int length() const {return mats.size();}
};

class symbop :public mat22
{
private:
  symbdata* sd;
public:
  symbop(symbdata* sdi, const mat22& mm) : mat22(mm), sd(sdi) {}
  symbop(symbdata* sdi, const Quad& a, const Quad& b, const Quad& c, const Quad& d) : mat22(a,b,c,d), sd(sdi)  {}

  // wrapper around the right action of a 2x2 matrix on a
  // (c:d)-symbol, mapping input symbol's index to output symbol's
  // index:

  int operator()(int i) const
  {
    return sd->index(((mat22)*this)(sd->symbol(i)));
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

inline ostream& operator<< (ostream& s, const mat22& m)
{
   s << "[" << (m.a) << "," << (m.b) << "; " << (m.c) << "," << (m.d) << "]";
   return s;
}

extern vector<RatQuad> alphas;  // List of a such that {a,oo} represent edge-orbits.
extern int n_alphas;            // Its length.
extern vector<mat22> M_alphas;  // List of matrices M_a  with det(M_a)=1 such that M_a(a)=oo.
void define_alphas();           // Populate alphas and M_alphas.

int nearest_alpha(const Quad& a, const Quad& b); // index of alpha nearest to a/b, given that a is reduced mod b

#endif
