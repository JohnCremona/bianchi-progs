// FILE qideal.cc

#include <iostream>
#include "intprocs.h"
#include "qideal.h"

////////////////////////////////////
// implementation of class Qideal //
////////////////////////////////////

// private -- converts output from findzbasis to standard Qideal basis
void Qideal::abc_from_HNF(vector<long>& basis)
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
{
  if (iclass!=-1) return;
  // cout<<"Filling in data for ideal"<<(*this)<<endl;
  Quad alpha=Quad(a), beta=Quad(b,1);
  findminquad(alpha, beta, g0, g1);
  // cout<<"findminquad for the primitive part returns g0="<<g0<<", g1="<<g1<<endl;
  if (div(g0,alpha) && div(g0,beta)) // principal!
    {
      iclass=0;
      g1 = 0;
    }
  else
    {
      iclass=1;
      g1 = makepos(c*g1);
    }
  g0 = makepos(c*g0);
}

Qideal::Qideal()
{ a=1; b=0; c=1; iclass=0; g0=1; g1=0; index=1;}   //default ideal is whole ring

Qideal::Qideal(const Qideal& i)   // the copy constructor
{ a=i.a; b=i.b; c=i.c; iclass=i.iclass; index=i.index;
  if(iclass!=-1)
    {
      g0=i.g0; g1=i.g1;
    }
}

Qideal::Qideal(const long& aa, const long& bb, const long& cc)
{
  if (cc==0)
    {
      a=1; b=0; c=0; iclass=0; g0=0; g1=0; index=1;
    }
  else
  {
    c = abs(cc);
    a = abs(aa);
    b = posmod(bb,a);
    if (!ok())
      {
        cerr <<"***Warning: Invalid ideal parameters ("<< aa << "," << bb << "," << cc << ") ***"<<endl;
      }
    iclass=-1;
    index=-1;
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
  : a(1), b(0), c(abs(n)), iclass(0), g0(abs(n)), g1(0), index(-1)
{ ; }

// Utility used to construct Qideals given two lists of integers rv,
// iv, or arbitrary but equal length, defining the Z-module spaned by
// all [rv[i], iv[i]]

Qideal::Qideal(const vector<long>& rv, const vector<long>& iv)  // ideal from Z-gens */
{
  vector<long> basis;
  findzbasis(rv,iv,basis);
  abc_from_HNF(basis);
  iclass=-1;
  index=-1;
}

Qideal::Qideal(const vector<Quad>& gens)       // ideal spanned by list of Quads */
{
  vector<long> rv, iv;
  for(vector<Quad>::const_iterator g = gens.begin(); g!=gens.end(); g++)
    {
      Quad a = *g;
      rv.push_back(a.re());
      iv.push_back(a.im());
      a *=Quad::w;
      rv.push_back(a.re());
      iv.push_back(a.im());
    }
  *this = Qideal(rv, iv);
}

Qideal::Qideal(const Quad& alpha) // principal ideal
{
  vector<Quad> gens = {alpha};
  *this = Qideal(gens);
  iclass=0;
  g0=makepos(alpha); g1=0;
  index=-1;
}


////////////////////////////////////
// Operators for the class Qideal //
////////////////////////////////////

Qideal Qideal::operator+(const long& d) const
{
  Qideal ans=*this;
  ans += d;
  return ans;
}

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
  if (divides(aa))
    return;  //ideal remains unchanged
  if (c==0)
    {
      *this = Qideal(aa);
      return;
    }
  vector<long> rv = get_rv(), iv = get_iv();
  rv.push_back(aa);
  rv.push_back(0);
  iv.push_back(0);
  iv.push_back(aa);
  *this = Qideal(rv, iv);
}

void Qideal::operator+=(const Quad& alpha)
{
  if (divides(alpha))
    {
      return;  //ideal remains unchanged
    }
  Qideal A(alpha);
  if (A.divides(*this))
    {
      *this = A;
      return;
    }

  vector<long> rv = get_rv(), iv = get_iv(), rva = A.get_rv(), iva = A.get_iv();
  rv.insert(rv.end(), rva.begin(), rva.end());
  iv.insert(iv.end(), iva.begin(), iva.end());
  *this = Qideal(rv, iv);
}

void Qideal::operator+=(const Qideal& f)
{
  if (divides(f))
    return;  //ideal remains unchanged
  if (f.divides(*this))
    {
      *this=f;
      return;
    }
  vector<long> rv1 = get_rv(), iv1 = get_iv();
  vector<long> rv2 = f.get_rv(), iv2 = f.get_iv();
  rv1.insert( rv1.end(), rv2.begin(), rv2.end() );
  iv1.insert( iv1.end(), iv2.begin(), iv2.end() );
  *this = Qideal(rv1, iv1);
}

Qideal Qideal::operator*(const long& d) const
{
  if (d==0)
    {
      return Qideal(0);
    }
  Qideal ans = Qideal(a,b,c*d);
  if (iclass!=-1)
    {
      ans.iclass=iclass; ans.g0=g0*d; ans.g1=g1*d;
    }
  return ans;
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
  if (d==0)
    {
      *this=Qideal(0);
      return;
    }
  c*=abs(d);
  if (iclass!=-1)
    { g0*=d; g1*=d; }
}

void Qideal::operator*=(const Quad& alpha)
{
  if (c==0 || alpha.norm()==1) return;              // this is unchanged
  if (alpha==0) { *this=Qideal(0); return;}         // this becomes 0
  if (alpha.im()==0) {*this *= alpha.re(); return;} // alpha in Z
  if (iclass==0) { *this=Qideal(g0*alpha); return;} // this is principal

  vector<Quad> gens = {alpha*Quad(a), alpha*Quad(b,1)}; // without the factor c
  long savec=c;
  *this = Qideal(gens);
  c*=savec;
  if (iclass!=-1) {g0*=alpha; g1*=alpha;}
}

void Qideal::operator*=(const Qideal& f)
{
  if (c==0 || ((f.c==1)&&(f.a==1))) return;    // unchanged
  if ((f.c==0) || ((c==1)&&(a==1))) { *this=f; return; }  // f unchanged
  if ((iclass==0)&&(f.iclass==0)) { *this=Qideal(g0*(f.g0)); return; }

  Quad g1(a), g2(b,1), h1(f.a), h2(f.b,1);
  vector<Quad> gens = {g1*h1, g1*h2, g2*h1, g2*h2};
  long savec=c;
  *this = Qideal(gens);
  c*=savec*(f.c);
  iclass=-1;
}

int Qideal::contains(const Quad& alpha) const
{
  std::ldiv_t qr1 = ldiv(alpha.re(), c);
  if (qr1.rem !=0) { return 0;}
  std::ldiv_t qr2 = ldiv(alpha.im(), c);
  return (qr2.rem==0) && ((qr1.quot - qr2.quot * b) % a == 0);
}

// next two functions needed by ideal_prod_coeffs
Quad Qideal::princprod_coeff_alpha(vector<long>&z) const
{
  Quad ans = Quad( z[3]*a + z[0]*b, z[0] + z[4]*a ) * c;
  if (ans != elt_spanned_by( z[3] - z[4]*b, z[0] + z[4]*a))
    cerr << "Warning: possible overflow in princprod_coeff_alpha!"<<endl;
  return ans;
}

Quad Qideal::princprod_coeff_beta(vector<long>&z) const
{
  return elt_spanned_by( z[1], z[2]);
}

Qideal Qideal::ideal_prod_coeffs(const Qideal&q,
                                 Quad&alphax, Quad&betax,
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

  long n = Quad::n, t=Quad::t;
  vector<long> rv = {a*q.b, b*q.a, b*q.b-n, a*q.a, 0};
  vector<long> iv = {    a,   q.a, b+q.b+t,     0, a*q.a};
  vector<long> x, y, basis;
  specialfindzbasiscoeffsmod(rv,iv,basis,x,y);

  // STEP TWO: standardise basis, taking care to adjust x,y accordingly
  // (this would not be needed if the HNF procedure already standardised)

  if (basis[1]<0)
    {
      basis[0]*=-1; basis[1]*=-1;
      for (long i=0; i<5; i++) x[i]*=-1;
    }
  if (basis[2]<0)
    {
      basis[2]*=-1;
      for(long i=0; i<5; i++) y[i]*=-1;
    }
  long newbasis0 = posmod(basis[0], basis[2]);
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

  return p0q0*(c*q.c);
}

Qideal Qideal::princprod(const Qideal&q, Quad&alpha, Quad&beta) const
//INPUT: ideals *this = p, q whose product is principal, equal to <gen>, say
//OUTPUT: the product ideal, together with alpha, beta \in q s.t.
//        gen = ac * alpha + (b+w)c * beta
{
  Quad alphax,betax,alphay,betay,g0;
  Qideal ans = ideal_prod_coeffs(q, alphax,betax,alphay,betay);
  vector<long> coeff = findminquadcoeffs( ans.a, Quad(ans.b,1), g0);
  // now g0 = coeff[1] * ans.a + coeff[0] * (ans.b+w)
  // gen= (ans.c)*g0;

  if ((ndiv(g0, Quad(ans.a,0)))||(ndiv(g0, Quad(ans.b,1))))
    {
      cerr << "Error: princprod did not find principal generator!"<<endl;
      exit(1);
    }

  alpha= coeff[0]*alphax + coeff[1]*alphay;
  beta = coeff[0]*betax  + coeff[1]*betay;
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
{
  (*this) *= alpha.conj();
  long na = alpha.norm();
  std::ldiv_t qr = ldiv(a, na);
  if (qr.rem!=0)                // shouldn't happen
    { cerr << "***inexact ideal division of "<<*this<<" by "<<alpha<<" ***"<<endl;}
  c = qr.quot;
  if (iclass!=-1) { g0/=na; g1/=na; }
}

void Qideal::operator/=(const Qideal&f)
{ Qideal keep = *this;
  (*this) *= f.conj();
  long nf = f.norm();
  std::ldiv_t qr = ldiv(c, nf);
  if (qr.rem!=0)                // shouldn't happen
    { cerr << "***inexact ideal division of "<<keep<<" by "<<f<<" ***"<<endl;}
  c = qr.quot;
  if (iclass!=-1) { g0/=nf; g1/=nf; }
}

////////////////////////////////////////////////////////
// Qideal:: non-operator member functions and friends //
////////////////////////////////////////////////////////

int Qideal::ok() const
{
 return quadnorm(Quad(b,1))%a==0;
 //return ( (a!=0) && ((( xmodmul(b,b,a) + b*Quad::t + Quad::n)%a)==0) ) ;
}

////////////////////////////////////////////////////////////////////////
// Naive version which first overflows for [47743,47434+a], d=(?)5  ////
//{ return ( (a!=0) && ((( b*b + b*Quad::t + Quad::n)%a)==0) ) ;} ////
////////////////////////////////////////////////////////////////////////

int Qideal::is_principal()
{
  //  cout<<"Testing "<<(*this)<< " for being principal"<<endl;
  if (iclass==-1)
    {
      //      cout<<" ... not yet known, calling fill()"<<endl;
      fill();
    }
  //  cout<<" ...iclass = "<<iclass<<", returning "<<(iclass==0)<<endl;
  return (iclass==0);
}

Qideal Qideal::conj() const
{ Qideal ans=*this;
  ans.b= -ans.b - Quad::t;
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

ostream& operator<<(ostream& s, const Qideal& x)
{
  long a=x.a, b=x.b, c=x.c;
  if (c!=1) s << c;
  s << "[" << a << ",";
  if (b!=0) s << b << "+";
  s << Quad::name << "]";
  return s;
}

char* to_string(const Qideal& a)
{
  ostringstream ans;
  ans << a;
  return (char*)ans.str().c_str();
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

      if (!iq.divides(q))
	cerr << "Error: comax returning q="<<q<<", not in "<<iq<<endl;
      return 1;
    }
  else {return 0;}
}

vector<Qideal> primitive_ideals_with_norm(long N, int both_conj)
// N is the norm of a primitive ideal iff it has no inert prime
// factors and ramified primes divid it at most once
{
  vector<Qideal> ans;
  vector<long> pdivs_norm = pdivs(N);
  for(vector<long>::const_iterator pi=pdivs_norm.begin(); pi!=pdivs_norm.end(); pi++)
    {
      long p=*pi;
      int s = Quad::chi(p);
      if ((s==-1) || ((s==0) && val(p,N)>1)) return ans;
    }
  // now the local tests pass so there are primitive ideals of norm N.
  // The have HNF [N,b+w] where b (mod N) satisfies
  //     b^2+t*b+n = 0 (mod N)
  //
  // Stupid implementation here: try all b mod N
  long t = Quad::t, n=Quad::n;
  long b, maxb = (both_conj? N-1: (N-t)/2);  // rounded down
  for (b=0; b<=maxb; b++)
    if (((b*(b + t) + n) % N) == 0)
      ans.push_back(Qideal(N,b,1));
  return ans;
}

vector<Qideal> ideals_with_norm(long N, int both_conj)
// I = c*I0 uniquely, where c^2|N and I0 is primitive of norm N/c^2
{
  vector<Qideal> ans;
  vector<long> clist = sqdivs(N); // list of c such that c^2|N
  for (vector<long>::const_iterator ci=clist.begin(); ci!=clist.end(); ci++)
    {
      long c = *ci;
      vector<Qideal> primitives = primitive_ideals_with_norm((N/c)/c, both_conj);
      for (vector<Qideal>::const_iterator Ji = primitives.begin(); Ji!=primitives.end(); Ji++)
        ans.push_back(c*(*Ji));
    }
  return ans;
}

vector<Qideal> ideals_with_bounded_norm(long maxnorm, int both_conj)
{
  vector<Qideal> ans;
  for (long N=1; N<=maxnorm; N++)
    {
      vector<Qideal> IN = ideals_with_norm(N, both_conj);
      ans.insert(ans.end(), IN.begin(), IN.end());
    }
  return ans;
}

void Qideal::set_index(int ind)
{
  if (index>0) // index already set, nothing to do
    return;
  if (ind!=0)  // we are given the index so use it
    {
      index=ind;
      return;
    }
  // compute the index by finding the sorted list of all ideals of this norm
  long n = norm();
  vector<Qideal> II = Qideal_lists::ideals_with_norm(n);
  vector<Qideal>::iterator i = find(II.begin(), II.end(), *this);
  if (i==II.end())
    {
      cerr<<"Error in finding index of "<<(*this)<<" (norm "<<n<<") in list of all ideals of norm "<<n<<": "
          <<II<<endl;
      exit(1);
    }
  index = std::distance(II.begin(), i);
  // cout<<"In set_index(), have just set index of ideal "<<(*this)<<" to "<<index<<endl;
}

string ideal_label(Qideal& I)  // returns label of ideal I
{
  stringstream s;
  s << I.norm() << "." << I.get_index();
  return s.str();
}

Qideal Qideal_from_norm_index(long N, int i) // i'th ideal of norm N
{
  if (i<1)
    {
      cerr<<"Qideal_from_norm_index(): i must be a positive integer"<<endl;
      exit(1);
    }
  vector<Qideal> II = Qideal_lists::ideals_with_norm(N);
  if (i>(int)II.size())
    {
      cerr<<"Qideal_from_norm_index(): i="<<i<<" is greater than "<<II.size()
          <<", the number of ideals of norm N="<<N<<endl;
      exit(1);
    }
  return II[i-1];
}

Qideal::Qideal(const string& s)           // ideal from label N.i
// need to parse s to obtain N and i as postive integers
{
  string Nstr, istr;
  stringstream ss(s);
  std::getline(ss, Nstr, '.');
  std::getline(ss, istr);
  long N;
  int i;
  stringstream(Nstr)>>N;
  stringstream(istr)>>i;
  *this = Qideal_from_norm_index(N,i);
}

// END OF FILE qideal.cc

