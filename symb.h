// FILE SYMB.H: Declarations for M-symbols, modular symbols

#if     !defined(_SYMB_H)
#define _SYMB_H      1       //flags that this file has been included

#include "moddata.h"

class symb;
class modsym;
class symblist;
class symbdata;
class symbop;
class edge_relations;

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
       {return div(N->modulus, d*s.c - c*s.d);}
   mat22 lift_to_SL2() const;
   friend inline ostream& operator<< (ostream&s, const symb&);
   friend class symblist;
   friend class symbdata;
   friend class modsym;
   friend class symbop;
   friend class mat22;
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
  int index2(const Quad& c, const Quad& d) const;
  int index(const symb& s) const {return index2(s.c,s.d);}
  symb symbol(int i) const;
  void display() const;
  int check(int verbose=0) const;
  friend class edge_relations;
};

class symbop :public mat22
{
private:
  symbdata* sd;
public:
  symbop(symbdata* sdi, const mat22& mm) : mat22(mm), sd(sdi) {}
  symbop(symbdata* sdi, const Quad& a, const Quad& b, const Quad& c, const Quad& d) : mat22(a,b,c,d), sd(sdi)  {}
  symb operator()(const symb& s) const
  {
    Quad sc = s.c, sd = s.d;
    apply_right(sc, sd);
    return symb(sc, sd, s.N);
  }

  // wrapper around the right action of a 2x2 matrix on a
  // (c:d)-symbol, mapping input symbol's index to output symbol's
  // index:

  int operator()(int i) const
  {
    return sd->index((*this)(sd->symbol(i)));
  }
};

inline ostream& operator<< (ostream& s, const symb& sy)
{
   s << "(" << sy.c << ":" << sy.d << ")";
   return s;
}

#endif
