// FILE quadarith.h
#ifndef __QUADARITH_H__
#define __QUADARITH_H__
#include <values.h>
//#include <Complex.h>
#define complex Complex
#include "arith.h"
#include "interface.h"

#include "tlistvar.h"

class Quad;

#include "field.h"

//  various rounding functions

long roundover(long a, long b);
long ldroundover(double a, double b);
long old_roundover_by_JEC(long a, long b);
long roundover_away_from_zero(long a, long b);
bigint floor(const bigint &a, const bigint&b );

//functions assigned by Field::Field constructor

extern Quad (*quadgcd)(const Quad& aa, const Quad& bb);
Quad quadgcd1(const Quad& aa, const Quad& bb); //Euclidean only
Quad quadgcd2(const Quad& aa, const Quad& bb); //General

extern Quad (*quadbezout)(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy);
Quad quadbezout1(const Quad&, const Quad&, Quad&, Quad&); //Euclidean
Quad quadbezout2(const Quad&, const Quad&, Quad&, Quad&); //General

int div(const long& a, const Quad& b);
int ndiv(const long& a, const Quad& b);
int div(const Quad& a, const Quad& b);
int ndiv(const Quad& a, const Quad& b);
long val(const Quad& factor, const Quad& number);

Quad invmod(const Quad& a, const Quad& p);
int coprime(const Quad& a, const Quad& b);
int invertible(const Quad& a, const Quad& b, Quad& inverse);

void specialfindzbasiscoeffsmod(long, long*, long*, long*, long*, long*);
/*
void findzbasiscoeffs(long, long*, long*, long*, long*, long*);
*/
void findzbasis(long, long*, long*, long*);
long* findminquadscoeffs(const Quad&, const Quad&, Quad&, Quad&);
long* findminquadcoeffs(const Quad&, const Quad&, Quad&);
void findminquads(const Quad&, const Quad&, Quad&, Quad&);
void findminquad(const Quad&, const Quad&, Quad&);


class Quad {
private:
    long r,i;
/*    Quad(const complex& z);   //rounds to nearest    */
public:
    Quad(long x=0, long y=0) :r(x),i(y) {;}
    Quad(const Quad& a) :r(a.r), i(a.i) {;}

    long real() const { return r; }
    long imag() const { return i; }

    Quad conj() const 
      { if (Field::t) return Quad(r+i, -i); else return Quad(r, -i); }

    double dnorm() const;     // inline this?
    int norm(long&n) const;   // rtns 1 if oflow, 0 if okay with n=norm
    long norm() const
      {
	long ans;
	int oflow=this->norm(ans);
	if (oflow)
	  { cerr <<"Error: overflow detected in computing norm of Quad!"<<endl;
	    // long d=0,e=0; e/=d;    // for debugging - sets breakpoint
	  }
	return ans;
      }

    int div(const Quad& b) const;
    int small_integer_multiple(long&, Quad&) const;  // rtns 1 if oflow, 0 o/w
                                                     // when this*Quad=long
    void operator=(const Quad& a) {r=a.r; i=a.i;}
    friend long real(const Quad& a);
    friend long imag(const Quad& a);
    
    int operator== (const Quad& b) const {return (r==b.r) && (i==b.i);}
    int operator== (long b) const {return (r==b) && (i==0);}
    int operator!= (const Quad& b) const {return (r!=b.r) || (i!=b.i);}
    int operator!= (long b) const {return (r!=b) || (i!=0);}
//
    Quad operator+ (const Quad& b) const {return Quad(r+b.r,i+b.i);}
    Quad operator+ (long b) const {return Quad(r+b,i);}
    friend Quad operator+(long m, const Quad& a);
    void operator+=(const Quad& b) {r+=b.r; i+=b.i;}
    void operator+=(long b) {r+=b;}
//
    Quad operator- (const Quad& b) const {return Quad(r-b.r,i-b.i);}
    Quad operator- (long b) const {return Quad(r-b,i);}
    friend Quad operator-(long m, const Quad& a);
    void operator-=(const Quad& b) {r-=b.r; i-=b.i;}
    void operator-=(long b) {r-=b;}
    Quad operator- () const {return Quad(-r,-i);}

//  obsolete versions of multn with undiagnosed overflows
//  Quad operator* (const Quad& b) const {return multi(*this,b);}
//  void operator*=(const Quad& b) {*this=multi(*this,b);}
//  Quad operator* (long m) const {return Quad(m*r,m*i);}
//  void operator*=(long m) {r*=m;i*=m;}
//  friend Quad operator*(long m, const Quad& a);

    Quad operator* (long m) const
      {
	double ar=(double)r, ai=(double)i, br=(double)m;
	double prodr = ar * br;
	double prodi = ai * br;
	if ((fabs(prodr)>MAXLONG)||(fabs(prodi)>MAXLONG))
	  { cerr << "Error: overflow in Quad * long multiplication!"<< endl; 
	    long d=0,e=0; d/=e;}
	return Quad((long)prodr,(long)prodi);
      }
    void operator*=(long m) { *this=(*this)*m; }

    Quad operator* (const Quad& b) const
      {
	double ar=(double)r, ai=(double)i, br=(double)b.r, bi=(double)b.i;
	double prodr = ar * br - ai * bi * Field::n;
	double prodi = ar * bi + ai * br;
	if (Field::t) prodi += ai * bi;
	if ((fabs(prodr)>MAXLONG)||(fabs(prodi)>MAXLONG))
	  { cerr << "Error: overflow in Quad * Quad multiplication!"<< endl;
	    long d=0,e=0; d/=e;}
	return Quad((long)prodr,(long)prodi);
      }
    void operator*=(const Quad& b) { *this=(*this)*b; }

//  obsolete versions of division with undiagnosed overflows
//  Quad operator/ (long b) const {return qdivi(*this,b);}
//  void operator/=(long b) {*this=qdivi(*this,b);}
//  Quad operator/ (const Quad& b) const 
//    {
//	if (b.i==0) return *this/b.r;
//	else return qdivi(multi(*this,b.conj()),b.norm());
//    }
//  void operator/=(const Quad& b) 
//    {
//	if (b.i==0) *this/=b.r;
//	else *this=qdivi(multi(*this,b.conj()),b.norm());
//    }
//  friend Quad operator/(long m, const Quad& a);

    Quad operator/ (long b) const 
      {
	if (Field::t)
	  { long ansi = roundover(i,b);
	    long ansr = roundover(2*r + i - b*ansi, 2*b);
	    return Quad(ansr,ansi);
	  }
	else
	  { return Quad(roundover(r,b), roundover(i,b)); }
      }

    void operator/=(long b) { *this = *this/b; }

    Quad operator/ (const Quad& b) const;
    void operator/=(const Quad& b) { *this = *this/b; }

    friend Quad operator/(long m, const Quad& a);

    int is_pos() const
      {
	if (Field::nunits > 2)
	  {                                                 // the old pos13
	    return ( ((i>=0)&&(r>0)) || ((r==0)&&(i==0)) );
	  }
	else
	  {                                                 // the old pos2
	    return ( (i>0) || ((i==0)&&(r>=0)) );
	  }
      }

    Quad pos_assoc() const                 // return the "positive" associate
      {                                    // replaces  Quad makepos(Quad)
	Quad ans = *this;
	while (!ans.is_pos()) ans *= Field::fundunit;
	return ans;
      }

// UGLY: make positive in situ
//
// void makepos() {while(!pos()) (*this)*=Field::fundunit; }     //ugly
//
		      
/*    operator complex() const;    */
    
    friend ostream& operator<<(ostream& s, const Quad& x);
    friend istream& operator>>(istream& s, Quad& x);
};

// for backward compatibility only - not recommended usage
/*
inline Quad quadconj(const Quad&a) { return a.conj(); }
inline long quadnorm(const Quad&a) { return a.norm(); }
*/

inline double realnorm(const Quad& z) {return sqrt(z.dnorm());} 

char* to_string(const Quad& a);  // outputs to a (new) string
inline Quad inverse(Quad& a)
{
    return Quad(1)/a;
}
inline Quad operator% (const Quad& a, const Quad& b) 
{ 
    return a-(b*(a/b));
}

inline long real(const Quad& a) {return a.r;}
inline long imag(const Quad& a) {return a.i;}


inline Quad operator+(long m, const Quad& a) 
{
    return Quad(m+a.r,a.i);
}
inline Quad operator-(long m, const Quad& a) 
{
    return Quad(m-a.r,-a.i);
}

inline Quad operator*(long m, const Quad& a) 
{ return a*m; }

inline Quad operator/(long m, const Quad& a) 
{ return Quad(m)/a; }

inline istream& operator>>(istream& s, Quad& x) 
{
  return s >> x.r >> x.i;
}

typedef Tlist<Quad> Quadlist;
typedef Tvar<Quad> Quadvar;

Quadlist residues(const Quad& a);

class Quadunits {
  static Quadlist ulist;
public:
  static void init()
    { ulist=Quadlist(Field::nunits);
      ulist[0]=1;
      for(long i=1; i<Field::nunits; i++) ulist[i]=Field::fundunit*ulist[i-1];
    }
  static Quad u(long i) { return ulist(i); }
};

#endif

// END OF FILE quadarith.h
