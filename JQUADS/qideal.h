// FILE qideal.h
// should only be included if ideal-methods is selected, i.e. MAX_CLASSNUM>1
#ifndef __QIDEAL_H__
#define __QIDEAL_H__

#include "quadarith.h"   // which includes "arith.h"

class Qideal;

class Qideal  {
  long a,b,c;   // the Z-basis in standard form  I=c[a,b+w]
  long iclass;  // 0 if principal, 1 if not (since h<=2), -1 if undetermined
  Quad g0,g1;   // if iclass!=-1 then I=<g0,g1>, g0 non-zero of minl norm,
                // thus  iclass=0 iff I=<g0>, when g1 is undefined
                // currently g1=0 if I princ a priori, o/w as found by fill()
  public:
//constructors 
  Qideal();                          // base initialization for class `Qideal'
  Qideal(const Qideal& );            // the copy constructor
  Qideal(const long&, const long&, const long&); // construct given Z-basis abc
  Qideal(const long&);               // princ ideal gen by long
  Qideal(const Quad&);               // princ ideal gen by Quad
/*  Qideal(const Quad&, const Quad&);  // ideal spanned by two Quads */
  ~Qideal()  {;}                     //destructor - currently nothing to do!
// member access
  long ay() const {return a;}
  long bee() const {return b;}
  long cee() const {return c;}
  Quad gee0() const {return g0;}
//
//operators
//
  int operator== (const Qideal& f) const {return (a==f.a)&&(b==f.b)&&(c==f.c);}
  int operator!= (const Qideal& f) const {return (a!=f.a)||(b!=f.b)||(c!=f.c);}
//
  Qideal operator+(const long&) const;
  Qideal operator+(const Quad&) const;
  Qideal operator+(const Qideal&) const;
  friend Qideal operator+(const long&, const Qideal&);
  friend Qideal operator+(const Quad&, const Qideal&);
  void operator+=(const long&);
  void operator+=(const Quad&);
  void operator+=(const Qideal&);
//
  Qideal operator*(const long&) const;
  Qideal operator*(const Quad&) const;
  Qideal operator*(const Qideal&) const;
  friend Qideal operator*(const long&, const Qideal&);
  friend Qideal operator*(const Quad&, const Qideal&);
  void operator*=(const long&);
  void operator*=(const Quad&);
  void operator*=(const Qideal&);

  Qideal ideal_prod_coeffs(const Qideal&, Quad&, Quad&, Quad&, Quad&) const;
  Qideal princprod(const Qideal&q, Quad&alpha, Quad&beta) const;

//
  Qideal operator/(const long&) const;
  Qideal operator/(const Quad&) const;
  Qideal operator/(const Qideal&) const;
  friend Qideal operator/(const long&, const Qideal&);
  friend Qideal operator/(const Quad&, const Qideal&);
  void operator/=(const long&);
  void operator/=(const Quad&);
  void operator/=(const Qideal&);
//
//functions defined in qideal.cc unless inline
  int isprincipal();         // fills iclass if necessary
  int divides(const long&) const ;
  int divides(const Quad&) const ;
  int divides(const Qideal&) const ;
  int div(const Qideal&f) const { return divides(f); }

  double dnorm() const { double ans=a; ans*=c; ans*=c; return ans; }
  long norm() const
    { double dn=dnorm();
      if (dn>MAXLONG) cerr << "Error: Qideal::norm overflow!" << endl;
      return (long)dn;
    }

  Qideal conj() const;            // returns the conjugate ideal

  Qideal pos_assoc() const { Qideal ans=*this; ans.fill(); return ans; }
  friend Qideal fill(const Qideal&f) { return f.pos_assoc(); }

  friend int comax(const Qideal&, const Qideal&, Quad&, Quad&);
  // returns 1 iff Qideals are comax, when returns one Quad from each, sum = 1

// i/o
  friend ostream& operator<<(ostream& s, const Qideal& x);
  friend istream& operator>>(istream& s, Qideal& x);
private:
  int ok() const;                 // checks that [a,b+w] *is* an ideal
  void fill();                    // determines iclass, g0, g1
  void abc_from_HNF(long*&);
  Quad princprod_coeff_alpha(long*&z) const;   // for use only in...
  Quad princprod_coeff_beta(long*&z) const;    // ... ideal_prod_coeffs
  Quad elt_spanned_by(const long&, const long&) const;
};

char* to_string(const Qideal& a);  // outputs to a (new) string

int div(const Qideal&, const Qideal&);  // returns 1 iff 1st divides 2nd

long val(const Qideal& factor, const Qideal& dividend);


#endif

// END OF FILE qideal.h
