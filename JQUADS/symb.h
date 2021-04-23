// FILE symb.h   --- Declarations for cusps, M-symbols, modular symbols
// new version based on matrix representation of cusps

#if     !defined(_SYMB_H)
#define _SYMB_H      1       //flags that this file has been included

#include "moddata.h"

class n_mat;    // normaliser matrix, i.e. matrix from the normaliser group
class n_modsym; // new version of modsym

class cdsymb;
class symbdata;
class geometry;
class symbop;

class abstract_mat {
  Quad aa,bb,cc,dd;
  long iclass;       // = 0 for principal 
                     // =-1 if unknown or if not in normaliser group at all

public: 
  abstract_mat() : aa(1), bb(0), cc(0), dd(1), iclass(0) {}

  abstract_mat(const Quad&, const Quad&, const Quad&, const Quad&);

  abstract_mat(const Quad&a,const Quad&b,const Quad&c,const Quad&d,long t) 
    : aa(a), bb(b), cc(c), dd(d), iclass(t) {}

  // member access
  Quad a() const { return aa; }
  Quad b() const { return bb; }
  Quad c() const { return cc; }
  Quad d() const { return dd; }
  long cl() const { return iclass; }

  Quad det() const { return aa*dd - bb*cc ; }
  void operator/=(const long&);                        // divide out by scalar
  abstract_mat conj() const;                           // return complex conj
  abstract_mat operator*(const abstract_mat&) const;   // mx multn
  friend inline ostream& operator<< (ostream& s, const abstract_mat&m);
};

class n_cusp : public abstract_mat {
public: 
  n_cusp() : abstract_mat(1,0,0,1,0) {}   // default is infinity!
  n_cusp(const abstract_mat&m) : abstract_mat(m) {}
//    { cerr << "Given abstract_mat "<<m<< " produced n_cusp "<<*this<<endl; }

  n_cusp(const Quad& ia, const Quad& ib, const Quad& ic, const Quad& id) 
    : abstract_mat(ia,ib,ic,id) {}

  n_cusp(const Quad&ia, const Quad&ib, const Quad&ic, const Quad&id, long t) 
    : abstract_mat(ia,ib,ic,id,t) {}
  friend ostream& operator<< (ostream& s, const n_cusp&z);

  int operator==(const n_cusp&z) const { return ( (c()*z.a()) == (z.c()*a()));}
  int operator!=(const n_cusp&z) const { return !(*this==z); }

  n_cusp infmat() const;

  static n_cusp infty() { return n_cusp(1,0,0,1,0); }
  static n_cusp zero() { return n_cusp(0,-1,1,0,0); }

};

class n_modsym {
  n_cusp a,b;
public:
  n_modsym() {}             // default ctor - calls default ctors for a,b
  n_modsym(const n_cusp& ra, const n_cusp& rb) : a(ra), b(rb) {}
  n_cusp alpha() const { return a; }
  n_cusp  beta() const { return b; }
  n_modsym conj() const;    // return image under complex conjugation
  friend ostream& operator<< (ostream& s, const n_modsym& m); //inline below
};

class n_mat : public abstract_mat {
public: 
  n_mat() : abstract_mat(1,0,0,1,0) {}
  n_mat(const abstract_mat&m) : abstract_mat(m) {}

  n_mat(const Quad& ia, const Quad& ib, const Quad& ic, const Quad& id) 
    : abstract_mat(ia,ib,ic,id) {}

  n_mat(const Quad&ia, const Quad&ib, const Quad&ic, const Quad&id, long t) 
    : abstract_mat(ia,ib,ic,id,t) {}

  int operator==(const n_mat&m) const
    { return (a()==m.a())&&(b()==m.b())&&(c()==m.c())&&(d()==m.d());}
  int operator!=(const n_mat&m) const { return !(*this==m); }

  n_mat inv_mod_scalars() const { return n_mat(d(),-b(),-c(),a(), cl()); }

  n_cusp operator()(const n_cusp&) const;      // action of matrix on cusp

  n_modsym operator()(const n_modsym&s) const  // action of matrix on modsym
    { return n_modsym(operator()(s.alpha()),operator()(s.beta()));}
  // 
  // ASIDE ON C++ SYNTAX:
  // Apparently we don't need to write
  //    this->operator()(s.alpha())
  //

  static n_mat specialmat(const Qideal&);

  static n_mat I() { return n_mat(); }
  static n_mat T(const Quad&q) { return n_mat(1,q,0,1,0); }
};

typedef Tlist<n_mat> matop;

// Obsolete matop class - use general list class
/*
class matop {  // formal sum of 2x2 matrices
private:
  n_mat* mats;
  long lng;
public:
  matop(long n) { lng=n;  mats = new n_mat[n]; }      // c'tor
  ~matop() { delete[] mats; }

  long length() const { return lng; }

  n_mat& operator[] (long i)
    { if ((i>=0)&&(i<lng)) return mats[i];
      else
	{ cerr << "Error: bad index "<<i<< " in matop of length "<<lng<<endl;
	  exit(1);
	}
    }

  n_mat operator() (long i) const
    { if ((i>=0)&&(i<lng)) return mats[i];
      else
	{ cerr << "Error: bad index "<<i<< " in matop of length "<<lng<<endl;
	  exit(1);
	}
    }
};
*/

class cdsymb {
  Quad c,d;
public:
  cdsymb() :c(0), d(0)   {;}
  cdsymb(const Quad& cc, const Quad& dd) :c(cc), d(dd)  {;}
  Quad cee() const {return c;}
  Quad dee() const {return d;}
//
//  There is no longer an equality operator for symbols, since it would
//  require either a global "level" variable, or that each symbol carry
//  a pointer to the level.  Instead, we ask "symbdata" if two cdsymbs
//  are equal.
//
  friend ostream& operator<< (ostream&s, const cdsymb&);
  friend class symbdata;
};

/*
class symblist {
  cdsymb* list;
  long num,maxnum;
public:
  symblist(long n=0);
  ~symblist() {delete[] list;}
  void add(const cdsymb& s, long start=0);
  long index(const cdsymb& s, long start=0) const;
  cdsymb item(long n) const;
  void symblist_display(ostream&s=cout) const 
    {for(long i=0; i<num; i++) s <<i<<":\t"<<list[i]<<"\n";}
  long count() const {return num;}
};
*/

class symbdata : public moddata, public quotient_ring {
  cdsymb* specials;
  long specialsnum;
  long nsymb1,nsymb2;   // with nsymb, moved here from moddata:: by JSB
protected:
  long nsymb;
public:
  symbdata(const Qideal&);             // The constructor
  ~symbdata() {delete[] specials;}
  long numsymb2(const Quad&, const Quad&) const;
  long numsymb(const cdsymb& s) const {return numsymb2(s.c,s.d);}
  cdsymb symbnum(long i) const;
  void display(ostream& s = cout, int verbose=0) const;
  int check(int verbose=0) const;

  int equal(const cdsymb&, const cdsymb&) const;
  n_mat symbtomat(const cdsymb&) const;
  n_mat special_adjustmat(const n_mat&) const;
  n_mat adjustmat(const n_mat&) const;

  long tysofs(long i, long t) const { return nsymb*t + i; }
  long softys(long i) const { return i%nsymb; }
  long toftys(long i) const { return i/nsymb; }
  n_modsym modsymoftys(const long&) const;

  // Hecke ops (as elts of group algebra) for action on modular symbols ...
  matop T_p(const Qideal&) const;                   //... princ prime
  matop U_p(const Qideal&) const;                   //... princ bad prime
  matop T_npnp(const Qideal&) const;                //... sq of non-princ prime
  matop T_npnq(const Qideal&, const Qideal&) const; //... prod of non-pr primes
  matop Chi_p(const Qideal&) const;                 // Normaliser char
};

class pseudo_euclid {

// first the static part, relating to the geometry
// must call init before first constructor call
  static long bubblenum;
  static n_mat* mats;
  static n_mat* imats;
public:
  static void init();
  static void dealloc() { delete[]imats; delete[]mats;
			cerr << "Deallocating pseudo_euclid!"<<endl;}
  static void display(ostream& s = cout, int verbose=0);
  static n_cusp centre(long j);
  static n_cusp co_centre(long j);

// now the "active" part, performing the algorithm for a given cusp
private:
  n_cusp z;
  n_mat m;
  long j;
  void translation_step();       // replaces z by z-q, m by m*T(q), sets j

  n_mat ma;
  long t;
  long eps;
  void reexpress();

public:
  pseudo_euclid(const n_cusp&zz) { z=zz; m=n_mat::I(); }
  int not_infty();               // if not there, takes one step towards infty
//  n_mat mi() const { return m; }
//  n_cusp beta() const { return co_centre(j); }

  n_mat summand_m() const { return ma; }
  long summand_t() const { return t; }
  long summand_eps() const { return eps; }
};

class geometry {
  static long aitch2_nonprinc_det_value;
  static n_cusp* basic_cusps;
  static long num_etypes; // num types of edge, i.e. edge orbits under PGL2(R)
  static n_modsym* basic_edges;
  static long max_length_of_face_rel;
  static long max_relmat_shape_ratio;
  //
  static void init_params();
  static void init_cusp_info();
  static void alloc_edge_info(long n);
  static void init_edge_info();
  //
public:
  static void init();
  static void dealloc() { pseudo_euclid::dealloc(); delete[] basic_edges; }
  static void display(ostream& s = cout, int verbose=0);

  static long aitch2_nonprinc_det_val() { return aitch2_nonprinc_det_value; }
  static n_cusp bc(long i) { return basic_cusps[i]; }
  static long nt() {return num_etypes;}
  static n_modsym be(long i) { return basic_edges[i]; }
  static long max_length_of_rel() {return max_length_of_face_rel;}
  static long relmat_shape() { return max_relmat_shape_ratio; }
};

class symbop {
  symbdata* sd;
  n_mat m;
public:
  symbop(symbdata* sdi, const n_mat& mm) :sd(sdi), m(mm) {}
  symbop(symbdata* sdi, const Quad& a, const Quad& b, const Quad& c, const Quad& d) :sd(sdi), m(a,b,c,d) {}
  long operator()(long i) const 
    // the right action of the symbop on M-symbols
    //
    //    ( ?   ? ) ( a b ) = (       ?                ?       )
    //    (s.c s.d) ( c d )   ( a*s.c + c*s.d    b*s.c + d*s.d )
    //
    {
      cdsymb s = sd->symbnum(i);
      long ans = sd->numsymb2(m.a()* s.cee() + m.c()* s.dee(),
			      m.b()* s.cee() + m.d()* s.dee());
      return ans;
    }
};
  
inline ostream& operator<< (ostream& s, const cdsymb& sy)
{
  s << "(" << sy.c << ":" << sy.d << ")";
  return s;
}  

inline ostream& operator<< (ostream& s, const abstract_mat&m)
{
  s << "((" << m.aa <<","<< m.bb <<"),("<< m.cc <<","<< m.dd <<"),";
  s << " cl=" << m.iclass << ")";
  return s;
}

inline ostream& operator<< (ostream& s, const n_cusp&z)
{
  s << "(" << z.a() << ")/(" << z.c() << ")";
  return s;
}  

inline ostream& operator<< (ostream& s, const n_modsym& m)
{
  s << "{" << (m.a) << "," << (m.b) << "}";
  return s;
}  

// WAS IN SEPARATE FILE cusp.h

class cusplist {
  Qideal modulus;
  int GL_vs_SL_flag;
  n_cusp *list;
  long num,maxnum;
  int eq(const n_cusp&, const n_cusp&) const;
public:
  cusplist(const Qideal&m, int flag, long n=0)
    { modulus=m;
      GL_vs_SL_flag = flag;
      maxnum=n;
      num=0;
      list=new n_cusp[n];
    }
  ~cusplist() {delete[] list;}
  long index(const n_cusp& a);
//  n_cusp item(long n) const {return list[n];}  //should check n really
  void display() const
    {for(long i=0; i<num; i++) cout<<i<<"\t"<<list[i]<<endl;}
  long count() const {return num;}
};   

#endif

// END OF FILE symb.h
