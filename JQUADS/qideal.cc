// FILE qideal.cc

#include <iostream>
#include "qideal.h"

////////////////////////////////////
// implementation of class Qideal //
////////////////////////////////////

// private -- converts output from findzbasis to standard Qideal basis
void Qideal::abc_from_HNF(long*& basis)
{ c=basis[1];
  b=basis[0]/c;  a=basis[2]/c;  if (a<0) {a=-a;}
  b=b%a;
  if (b<0) {b+=a;}     // not sure what % does with negative numbers
  if (c<0) {c=-c;}
  if (!ok())
    { cerr <<"***Warning: "<< *this <<" not ok in abc_from_HNF***"<<endl;} 
}

// private -- fills other data fields given a,b,c
void Qideal::fill()
{ if (iclass==-1)
    { Quad alpha=Quad(a), beta=Quad(b,1);
      findminquads(alpha, beta, g0, g1);
      if (g0.div(alpha) && g0.div(beta)) 
	{ iclass=0;        // could set g1=0
	  while (!g0.is_pos()) { g0*=Field::fundunit; g1*=Field::fundunit;}
	}
      else
	{ iclass=1; }      // assumes h=2
      g0*=c; g1*=c;
//cout<<*this<<"   "<<iclass<<" "<<g0<<" "<<g1<<endl;
    }
}

Qideal::Qideal()
{ a=1; b=0; c=1; iclass=0; g0=1; g1=0;}   //default ideal is whole ring

Qideal::Qideal(const Qideal& i)   // the copy constructor
{ a=i.a; b=i.b; c=i.c; iclass=i.iclass;
  if(iclass!=-1){g0=i.g0; g1=i.g1;}
}

Qideal::Qideal(const long& aa, const long& bb, const long& cc)
{ if (cc==0) { a=1; b=0; c=0; iclass=0; g0=0; g1=0; }
else 
  { c=abs(cc); a=abs(aa); b=(bb%a);
    if (b<0) {b+=a;}     // not sure what % does with negative numbers
    if (!ok())           // legal to call member function from within ctor
      { cerr <<"***Warning: Invalid ideal "<< *this <<" ***"<<endl;} 
    iclass=-1;
  }
}
// If this constructor gets called a lot in situations in which cc[aa,bb + w]
// is guaranteed to be an ideal and in standard form, then it would be worth
// dispensing with standardisation in favour of { a=aa; b=bb; c=cc; iclass=-1;}
// There could be an optional 4th parameter with default 'standardise'
// and other option 'accept as is' - it would then be up to the programmer to
// ensure that the latter is only called when it should be.
// For the moment, it seems better to enforce standard form so that other
// functions can rely on it.

Qideal::Qideal(const long& n)       // princ ideal gen by long
{ a=1; b=0; c=abs(n); iclass=0; g0=c; g1=0; }
  
Qideal::Qideal(const Quad& alpha)
{ if (alpha.imag()==0) {*this=Qideal(alpha.real());}
else
  { 
    iclass=0; g0=alpha.pos_assoc(); g1=0;

    // try and make the numbers a bit smaller
    long g=gcd(alpha.real(),alpha.imag());
    Quad beta=alpha/g;
    long normb=beta.norm();
    
    long* rv = new long[3];
    long* iv = new long[3];
    long* basis = new long[3];

    rv[0] = normb;              iv[0] = 0;            // to help findzbasis !
    rv[1] = beta.real();        iv[1] = beta.imag();   
    rv[2] = -Field::n*iv[1];    iv[2] = rv[1] + Field::t*iv[1];
    //cout<<"About to call findzbasis with rv="<<rv<<", iv="<<iv<<", basis="<<basis<<endl;
    findzbasis(3,rv,iv,basis);
    abc_from_HNF(basis);
    c*=g;                        // g assumed positive !

#ifdef testbezout
    //CHECK:
    Qideal aldash(a,b,c);  // prevents copying of iclass! 
    if (aldash.isprincipal() && (aldash.g0==g0)) {;}  //OK 
    else 
      { cerr<<"Problem in constructor Qideal::Qideal(const Quad&) !!!"<<endl;
	cerr<<"isprincipal unable to confirm principal generator!"<<endl;
	cerr<<"alpha = "<<alpha<<" but g0 = "<<aldash.g0<<" for "<<*this<<endl;
      }
#endif
    delete[] basis; delete[] iv; delete[] rv;
  }
}

/*
// Naively coded: if we ever need this procedure, should throw in an integer,
// such as the gcd of normalpha, normbeta, to help findzbasis !
//
Qideal::Qideal(const Quad& alpha, const Quad& beta)
{ if (beta==0) {*this=Qideal(alpha);}
else
  {
    long* rv = new long[4];
    long* iv = new long[4];
    long* basis = new long[3];
    rv[0] = real(alpha); iv[0] = imag(alpha);   
    rv[1] = -Field::n*iv[0];    iv[1] = rv[0] + Field::t*iv[0];
    rv[2] = real(beta); iv[2] = imag(beta);   
    rv[3] = -Field::n*iv[2];    iv[3] = rv[2] + Field::t*iv[2];
    //cout<<"About to call findzbasis with rv="<<rv<<", iv="<<iv<<", basis="<<basis<<endl;
    findzbasis(4,rv,iv,basis);
    abc_from_HNF(basis);
    delete[] basis; delete[] iv; delete[] rv;
    iclass=-1;
  }
}
*/

////////////////////////////////////
// Operators for the class Qideal //
////////////////////////////////////

Qideal Qideal::operator+(const long& d) const
{ 
  Qideal ans=*this;
  ans += d;
  return ans;
}

/*
//More efficient short-cut version
Qideal Qideal::operator+(const long& d)
{ 
  if (divides(d)) { return *this; }
  else
    if (c==0) { return Qideal(d); }
    else
      {
	long e,f;
	long g = bezout(c,d,e,f);
	c/=g; d/=g;
	long h = gcd(a*c,d);
	ans = Qideal(h,e*b*c,g);
	ans.iclass=-1;
	return ans;
      }
}
*/

Qideal Qideal::operator+(const Quad& alpha) const
{ Qideal ans=*this;
  ans+=alpha;
  return ans;
}

Qideal Qideal::operator+(const Qideal& f) const
{ Qideal ans=*this;
  ans+=f;
  return ans;
}

Qideal operator+(const long&a, const Qideal& f)
{ return f+a; }

Qideal operator+(const Quad&a, const Qideal& f)
{ return f+a; }

void Qideal::operator+=(const long& aa)
{ 
  if (divides(aa))   {;}  //ideal remains unchanged
  else
    if (c==0) { *this = Qideal(aa); }
    else
      {
	long* rv = new long[4];
	long* iv = new long[4];
	long* basis = new long[3];
	rv[0] = a*c;   iv[0] = 0;
	rv[1] = b*c;   iv[1] = c;
	rv[2] = aa;    iv[2] = 0;
	rv[3] = 0;     iv[3] = aa;
	findzbasis(4,rv,iv,basis);
	abc_from_HNF(basis);
	delete[] basis; delete[] iv; delete[] rv;
	iclass=-1;
      }
}

void Qideal::operator+=(const Quad& alpha)
{ 
  if (divides(alpha)) {;}  //ideal remains unchanged
  else
    {
      if (c==0) {*this = Qideal(alpha);}
      else
	{
	  long* rv = new long[4];
	  long* iv = new long[4];
	  long* basis = new long[3];
	  rv[0] = a*c;              iv[0] = 0;
	  rv[1] = b*c;              iv[1] = c;
	  rv[2] = real(alpha);      iv[2] = imag(alpha);
	  rv[3] = -Field::n*iv[2];   iv[3] = rv[2] + Field::t*iv[2];
	  findzbasis(4,rv,iv,basis);
	  abc_from_HNF(basis);
	  delete[] basis; delete[] iv; delete[] rv;
	  iclass=-1;
	}
    }
}

void Qideal::operator+=(const Qideal& f)
{ 
  if (divides(f)) {;}  //ideal remains unchanged
  else
    {
      if (f.divides(*this)) { *this=f;}
      else
	{
	  long* rv = new long[4];
	  long* iv = new long[4];
	  long* basis = new long[3];
	  rv[0] = a*c;            iv[0] = 0;
	  rv[1] = b*c;            iv[1] = c;
	  rv[2] = (f.a)*(f.c);    iv[2] = 0;
	  rv[3] = (f.b)*(f.c);    iv[3] = f.c;
	  findzbasis(4,rv,iv,basis);
	  abc_from_HNF(basis);
	  delete[] basis; delete[] iv; delete[] rv;
	  iclass=-1;  
	}
    }
}

Qideal Qideal::operator*(const long& d) const
{
  if (d==0) { return Qideal(0); }
  else
    { Qideal ans = Qideal(a,b,c*d);
      if (iclass!=-1) { ans.iclass=iclass; ans.g0=g0*d; ans.g1=g1*d; }
      return ans;
    }
}

Qideal Qideal::operator*(const Quad& alpha) const
{ Qideal ans=*this;
  ans*=alpha;
  return ans;
}

Qideal Qideal::operator*(const Qideal& f) const
{ Qideal ans=*this;
  ans*=f;
  return ans;
}

Qideal operator*(const long& d, const Qideal& f)
{ return f*d; }

Qideal operator*(const Quad& a, const Qideal& f)
{ return f*a; }

void Qideal::operator*=(const long& d)
{
  if (d==0) { *this=Qideal(0); }
  else
    { 
      c*=abs(d);
      if (iclass!=-1) { g0*=d; g1*=d; }
    }
}

void Qideal::operator*=(const Quad& alpha)
{ if (c==0 || alpha.norm()==1) {;}    // unchanged
else if (alpha==0) { *this=Qideal(0); }
else if (iclass==0) { *this=Qideal(g0*alpha); } 
else
  { // ans = c[a,b+w][alpha,alpha w] = c[a*alpha, *=w, (b+w)alpha, *=w]
    cerr << "Hup!"<<endl;
    long* rv = new long[5];
    long* iv = new long[5];
    long* basis = new long[3];
    
    rv[0] = a*(alpha.norm());    iv[0] = 0;  // to help findzbasis !

    rv[1] = a*(alpha.real());    iv[1] = a*(alpha.imag());
    rv[2] = -Field::n*iv[1];     iv[2] = rv[1] + Field::t*iv[1];
    Quad temp=Quad(b,1)*alpha;
    rv[3] = temp.real();         iv[3] = temp.imag();
    rv[4] = -Field::n*iv[3];     iv[4] = rv[3] + Field::t*iv[3];

    findzbasis(5,rv,iv,basis);
    long savec=c;
    abc_from_HNF(basis);
    c*=savec;
    delete[] basis; delete[] iv; delete[] rv;
    if (iclass!=-1) {g0*=alpha; g1*=alpha;}
  }
}

void Qideal::operator*=(const Qideal& f)
{ if (c==0 || ((f.c==1)&&(f.a==1)))   {;}    // unchanged
else if ((f.c==0) || ((c==1)&&(a==1))) { *this=f; }  // f unchanged
else if ((iclass==0)&&(f.iclass==0)) { *this=Qideal(g0*(f.g0)); } 
else
  { 
    // ans = c[a,b+w]f.c[f.a, f.b + w]
    //     = c(f.c)[a(f.a), a(f.b + w), (b + w)(f.a), (b+w)(f.b + w)]
    long* rv = new long[4];
    long* iv = new long[4];
    long* basis = new long[3];
    rv[0] = a*(f.a);          iv[0] = 0;
    rv[1] = a*(f.b);          iv[1] = a;
    rv[2] = b*(f.a);          iv[2] = f.a;
    rv[3] = b*(f.b)-Field::n;  iv[3] = b + f.b + Field::t;
    findzbasis(4,rv,iv,basis);
    long savec=c;
    abc_from_HNF(basis);
    c*=savec*(f.c);
    delete[] basis; delete[] iv; delete[] rv;
    iclass=-1;
  }
}

// next two functions needed by ideal_prod_coeffs
Quad Qideal::princprod_coeff_alpha(long*&z) const
{
  Quad ans = Quad( z[3]*a + z[0]*b, z[0] + z[4]*a ) * c;
  if (ans != elt_spanned_by( z[3] - z[4]*b, z[0] + z[4]*a))
    cerr << "Warning: possible overflow in princprod_coeff_alpha!"<<endl;
  return ans;
}

Quad Qideal::princprod_coeff_beta(long*&z) const
{
  return elt_spanned_by( z[1], z[2]);
}

Qideal Qideal::ideal_prod_coeffs(const Qideal&q, Quad&alphax, Quad&betax,
				 Quad&alphay, Quad&betay) const
{
// INPUT: ideals *this = p = c*p0 = c[a,b+w],  q = c'*q0 = c'[a',b'+w]
// OUTPUT: the product ideal pq=cc'C[A,B+w]
//         together with alphax,betax,alphay,betay \in q such that
//         cc'AC = (alphay)ac + (betay)(b+w)c
//     cc'(B+w)C = (alphax)ac + (betax)(b+w)c

// STEP ONE: find the Z-basis of the product of the primitive ideals
// Spanning set: [a,b+w][a',b'+w]=[a(b'+w),a'(b+w),(b+w)(b'+w),aa',aa'w]
// (last term thrown in to help specialfindzbasiscoeffsmod)

  long* rv = new long[5];
  long* iv = new long[5];
  long*  x = new long[5];
  long*  y = new long[5];

  long* basis = new long[3];
    
  rv[0] =  a * q.b ;                iv[0] = a;
  rv[1] =  b * q.a ;                iv[1] = q.a;
  rv[2] =  b * q.b - Field::n;      iv[2] = b + q.b + Field::t;
  rv[3] =  a * q.a;                 iv[3] = 0;
  rv[4] =  0;                       iv[4] = a * q.a;

  specialfindzbasiscoeffsmod(5,rv,iv,basis,x,y);

  // STEP TWO: standardise basis, taking care to adjust x,y accordingly
  // (this would not be needed if the HNF procedure already standardised)

  if (basis[1]<0) { basis[0]*=-1; basis[1]*=-1;
		    for (long i=0; i<5; i++) x[i]*=-1;}
  if (basis[2]<0) { basis[2]*=-1; for(long i=0; i<5; i++) y[i]*=-1; }
  long newbasis0 = basis[0]%basis[2];           // not sure what % does ...
  if (newbasis0<0) { newbasis0 += basis[2]; }   // ... with negative arguments
  long tmpquot = (newbasis0 - basis[0])/basis[2];
  basis[0] = newbasis0;
  for (long i=0; i<5; i++) x[i] += tmpquot*y[i];

  alphax= q.princprod_coeff_alpha(x);
  betax = q.princprod_coeff_beta(x);
  alphay= q.princprod_coeff_alpha(y);
  betay = q.princprod_coeff_beta(y);

  Qideal p0q0;                         // ugly bugfix:
  p0q0.iclass=-1;                      // ... there should be a suitable ...
  p0q0.abc_from_HNF(basis);            // ... constructor call

  if ((p0q0.c != basis[1])||
      (p0q0.a * p0q0.c != basis[2])||
      (p0q0.b * p0q0.c != basis[0]))
    { cerr << "Error: inconsistent standardisation in princprod"<<endl;
      long e,f=0; e/=f;
    }

  // should now have:
  // pq=c[a,b+w]c'[a',b'+w]=cc'C[A,B+w]
  // (BC+Cw)cc' = (alphax)ac + betax(b+w)c
  //    (AC)cc' = (alphay)ac + betay(b+w)c

  delete[] basis; delete[] y; delete[] x; delete[] iv; delete[] rv;
  return p0q0*(c*q.c);
}

Qideal Qideal::princprod(const Qideal&q, Quad&alpha, Quad&beta) const
//INPUT: ideals *this = p, q whose product is principal, equal to <gen>, say
//OUTPUT: the product ideal, together with alpha, beta \in q s.t.
//        gen = ac * alpha + (b+w)c * beta
{
  Quad alphax,betax,alphay,betay,g0;
  Qideal ans = ideal_prod_coeffs(q, alphax,betax,alphay,betay);
  long* coeff = findminquadcoeffs( ans.a, Quad(ans.b,1), g0);
  // now g0 = coeff[1] * ans.a + coeff[0] * (ans.b+w)
  // gen= (ans.c)*g0;

  if ((!g0.div(ans.a))||(!g0.div(Quad(ans.b,1))))
    {
      cerr << "Error: princprod did not find principal generator!"<<endl;
      exit(1);
    }

  alpha= coeff[0]*alphax + coeff[1]*alphay;
  beta = coeff[0]*betax  + coeff[1]*betay;
  delete[] coeff;
  return ans;
}

Qideal Qideal::operator/(const long&n) const
{ Qideal ans=*this;
  ans/=n;
  return ans;
}

Qideal Qideal::operator/(const Quad&alpha) const
{ Qideal ans=*this;
  ans/=alpha;
  return ans;
}

Qideal Qideal::operator/(const Qideal&f) const
{ Qideal ans=*this;
  ans/=f;
  return ans;
}

// more efficient than first promoting dividend, though doubtful if ever used
Qideal operator/(const long&n, const Qideal&f)
{ Qideal ans=f.conj()*n;
  return ans/(f.norm());
}

// more efficient than first promoting dividend, though doubtful if ever used
Qideal operator/(const Quad&alpha, const Qideal&f)
{ Qideal ans=f.conj()*alpha;
  return ans/(f.norm());
}

void Qideal::operator/=(const long&n)
{ long na=abs(n);
  if (c%na)
    { cerr<<"***inexact division of "<<*this<<" by "<<n<<" ***"<<endl;}
  c/=na;
  if (iclass!=-1) { g0/=na; g1/=na; }
}

void Qideal::operator/=(const Quad&alpha)
{ if (imag(alpha)==0) { *this /= real(alpha); }
else
  {
    int errflag=0;

    // trying to divide c[a,b+w] by alpha
    // first divide out largest possible integer
    long g = gcd(alpha.real(),alpha.imag());
    if (c%g) {errflag=1;}
    c/=g;  // pretend alpha/=g also!!

    // After scaling we have  alpha=x+yw  with  gcd(x,y)=1
    // Moreover, the conjugate beta is  x-yw  or  x+y -yw,  with gcd=1
    Quad beta = (alpha/g).conj();
    long n = beta.norm();
    // Ans is Z-span of ca/alpha = (ca beta)/n,  c(b+w)/alpha = (c(b+w)beta)/n
    // The remark on gcds implies n|ca.
    g = gcd(n,c);
    // The norm of the answer is (ac^2)/n = a/(n/g) * (c/g)^2 * g
    c/=g;  n/=g;
    // Ans is Z-span of c((a/n)beta), c(((b+w)beta)/n) with norm c^2*(a/n)*g
    // We can throw in normans to help findzbasis
    if (a%n) {errflag=1;}
    long normans = (a/n)*g;
    Quad gen1 = beta * (a/n);
    Quad gen2 = beta * Quad(b,1);
    if (ndiv(n,gen2)) {errflag=1;}
    gen2/=n;
    // answer = c * ( Z-span of gen1,gen2, and optionally normans )
    if (errflag)                // shouldn't happen
      { cerr<<"Error: inexact division of "<<*this<<" by Quad "<<alpha<<endl;}
    long savec = c;

    long* rv = new long[3];
    long* iv = new long[3];
    long* basis = new long[3];
    rv[0] = normans;    iv[0] = 0;
    rv[1] = real(gen1); iv[1] = imag(gen1); 
    rv[2] = real(gen2); iv[2] = imag(gen2);   
    //cout<<"About to call findzbasis with rv="<<rv<<", iv="<<iv<<", basis="<<basis<<endl;
    findzbasis(3,rv,iv,basis);
    abc_from_HNF(basis);
    delete[] basis; delete[] iv; delete[] rv;
    
    c*= savec;
    if (iclass!=-1) { g0/=alpha; g1/=alpha; }
  }
}

// obsolete version
/*
void Qideal::operator/=(const Quad&alpha)
{ (*this)*=alpha.conj();
  long na=alpha.norm(), tmp=c/na;
  if (tmp*na!=c)                // shouldn't happen
    { cerr<<"***inexact division of "<<*this<<" by "<<alpha<<" ***"<<endl;}
  c=tmp;
  if (iclass!=-1) { g0/=na; g1/=na; }
}
*/

void Qideal::operator/=(const Qideal&f)
{ (*this) *= f.conj();
  long nf=f.norm(), tmp=c/nf;
  if (tmp*nf!=c)                // shouldn't happen
    { cerr << "***inexact ideal division of "<<*this<<" by "<<f<<" ***"<<endl;}
  c=tmp;
  if (iclass!=-1) { g0/=nf; g1/=nf; }
}

////////////////////////////////////////////////////////
// Qideal:: non-operator member functions and friends //
////////////////////////////////////////////////////////

int Qideal::ok() const
{ return ( (a!=0) && ((( xmodmul(b,b,a) + b*Field::t + Field::n)%a)==0) ) ;}
////////////////////////////////////////////////////////////////////////
// Naive version which first overflows for [47743,47434+a], d=(?)5  ////
//{ return ( (a!=0) && ((( b*b + b*Field::t + Field::n)%a)==0) ) ;} ////
////////////////////////////////////////////////////////////////////////

int Qideal::isprincipal()
{ if (iclass==-1) fill();  return (iclass==0); }


int Qideal::divides(const long& nn) const
{ if (c==0) {return nn==0;} else return (nn%(a*c)==0);}

int Qideal::divides(const Quad& alpha) const
{
  if (c==0) {return (alpha==0);}
  else
    {
      int ans=((real(alpha)%c==0) &&
	       (imag(alpha)%c==0) &&
	       ((real(alpha)-b*imag(alpha))%(c*a)==0) );
      return ans;
    }
}

int Qideal::divides(const Qideal& f) const
{return (divides(f.a*f.c)) && (divides(Quad(f.b*f.c,f.c))) ;}

Qideal Qideal::conj() const
{ Qideal ans=*this;
  ans.b= -ans.b - Field::t;
  if (ans.b<0) {ans.b += ans.a;}
  if (ans.iclass!=-1)
    { ans.g0 = ans.g0.conj();
      ans.g1 = ans.g1.conj(); }
  return ans;
}

Quad Qideal::elt_spanned_by(const long& x, const long& y) const
{
  return c*( x*a + y * Quad(b,1));
}

istream& operator>>(istream& s, Qideal& y)
{
  Qideal x;
  do { s >> x.a >> x.b >> x.c;  
       if (x.ok()) { y=Qideal(x.a,x.b,x.c); return s; }  // assign via c'tor
       else {cerr << "retry: ";}
     } while (1);
}

// less elegant version
/*
istream& operator>>(istream& s, Qideal& y)
{ Qideal x;
  int dud;
  do {s >> x.a >> x.b >> x.c;
      dud=!x.ok(); 
      if (dud) {cerr << "retry: "; }
    } 
   while (dud);
  y=Qideal(x.a, x.b, x.c);    // assign via c'tor
  return s;
}
*/

ostream& operator<<(ostream& s, const Qideal& x)
{
  long a=x.a, b=x.b, c=x.c;
  if (c!=1) s << c;
  s << "[" << a << ",";
  if (b!=0) s << b << "+";
  s << Field::name << "]";
  return s;
}

char* to_string(const Qideal& a)
{
  ostringstream ans;
  ans << a;
  return (char*)ans.str().c_str();
}

int div(const Qideal&fraka, const Qideal&frakb) // rtns 1 iff 1st divides 2nd
{
  return fraka.divides(frakb);
}
  
long val(const Qideal& factor, const Qideal& dividend)
{
  if ((dividend==0) || (factor.norm()<=1))
    {
      cerr<<"Warning: 9999 returned in Qideal version of val"<<endl;
      return 9999;
    }
  long e = 0; 
  Qideal n=dividend;
  while (factor.divides(n)) {n/=factor; e++; }
  return e;
}

int comax(const Qideal&ip, const Qideal&iq, Quad&p, Quad&q)
// returns 1 iff ip,iq are comax, when returns p,q elts of ip,iq s.t. p+q=1
{
  long g11,g12,g13,gg1,a3,g21,g22,g23;
  long pac = ip.a*ip.c;
  long qac = iq.a*iq.c;
  g11 = bezout ( pac, qac, g12, g13);
  gg1 = gcd ( ip.c, iq.c );
  a3 = ( (ip.c / gg1 ) * iq.c ) * ( iq.b - ip.b);
  g21 = bezout ( g11, a3, g22, g23);
  if (g21==1)
    {
      long gg2 = gcd(pac,qac);
      long n1 = xmodmul(g22,g12, qac/gg2);
      long n2 = -(iq.c/gg1)*g23;
      p = ip.elt_spanned_by(n1,n2);
      q = 1 - p;

// obsolete way of finding q
//      n3 = g22*g13;
//      n4 = (ip.c/gg1)*g23;
//      q = iq.elt_spanned_by(n3,n4);

      if (!iq.divides(q))
	cerr << "Error: comax returning q="<<q<<", not in "<<iq<<endl;
      return 1;
    }
  else {return 0;}
}

// END OF FILE qideal.cc

