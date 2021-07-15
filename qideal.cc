// FILE qideal.cc

#include <iostream>
#include "intprocs.h"
#include "qideal.h"

////////////////////////////////////
// implementation of class Qideal //
////////////////////////////////////

// private -- converts output from findzbasis to standard Qideal basis
void Qideal::abc_from_HNF(vector<long>& basis)
{ c = abs(basis[1]);
  a = abs(basis[2]/c);
  b = posmod(basis[0]/c,a);
  ac = a*c;
  nm = ac*c;
  if (!ok())
    { cerr <<"***Warning: "<< *this <<" not ok in abc_from_HNF***"<<endl;}
}

// private -- fills other data fields given a,b,c
void Qideal::fill()
{
  if (iclass!=-1) return;
  //  cout<<"Filling in data for ideal"<<(*this)<<endl;
  g0=Quad(a); g1=Quad(b,1);
  unimod U;
  sl2z_reduce(g0,g1, U);
  if ((quadconj(g0)*g1).im()<0)
    cout<<"Badly oriented Z-basis in fill() 1"<<endl;
  //  cout<<"sl2z_reduce for the primitive part returns g0="<<g0<<", g1="<<g1<<endl;
  iclass = (div(g0,g1)? 0: 1); // 0 means principal
  // scale by content:
  g1 *= c;
  g0 *= c;
  //  cout<<" -- now gens are "<<gens()<<endl;
  if ((quadconj(g0)*g1).im()<0)
    cout<<"Badly oriented Z-basis in fill() 2"<<endl;
}

Qideal::Qideal()
{ a=c=ac=nm=1; b=0; iclass=0; g0=1; g1=Quad::w; index=1;}   //default ideal is whole ring

Qideal::Qideal(const Qideal& i)   // the copy constructor
{ a=i.a; b=i.b; c=i.c; ac=i.ac; nm=i.nm; iclass=i.iclass; index=i.index;
  if(iclass!=-1)
    {
      g0=i.g0; g1=i.g1;
    }
}

Qideal::Qideal(const long& aa, const long& bb, const long& cc)
{
  if (cc==0)
    {
      a=1; b=0; c=ac=nm=0; iclass=0; g0=0; g1=0; index=1;
    }
  else
  {
    c = abs(cc);
    a = abs(aa);
    b = posmod(bb,a);
    ac = a*c;
    nm = ac*c;
    if (!ok())
      {
        cerr <<"***Warning: Invalid ideal parameters ("<< aa << "," << bb << "," << cc << ") ***"<<endl;
      }
    iclass=-1;
    index=-1;
  }
}

Qideal::Qideal(const long& n)       // princ ideal gen by long
  : a(1), b(0), c(abs(n)), iclass(0), index(-1)
{
  ac=a*c; nm=ac*c;
  g0 = Quad(abs(n));
  g1 = Quad(0, abs(n));
}

// Utility used to construct Qideals given two lists of integers rv,
// iv, or arbitrary but equal length, defining the Z-module spaned by
// all [rv[i], iv[i]]

Qideal::Qideal(const vector<long>& rv, const vector<long>& iv)  // ideal from Z-gens
{
  vector<long> basis;
  findzbasis(rv,iv,basis);
  abc_from_HNF(basis);
  iclass=-1;
  index=-1;
}

Qideal::Qideal(const vector<Quad>& gens)       // ideal spanned by list of Quads
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
  g0=makepos(alpha); g1=Quad::w*g0;
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
  c *= abs(d);
  ac *= abs(d);
  nm *= (d*d);
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
  long fac = c;
  *this = Qideal(gens);
  c *= fac;
  ac *= fac;
  nm *= (fac*fac);
  if (iclass!=-1) {g0*=alpha; g1*=alpha;}
}

void Qideal::operator*=(const Qideal& f)
{
  if (c==0 || ((f.c==1)&&(f.a==1))) return;    // unchanged
  if ((f.c==0) || ((c==1)&&(a==1))) { *this=f; return; }  // f unchanged
  if ((iclass==0)&&(f.iclass==0)) { *this=Qideal(g0*(f.g0)); return; }

  Quad g1(a), g2(b,1), h1(f.a), h2(f.b,1);
  vector<Quad> gens = {g1*h1, g1*h2, g2*h1, g2*h2};
  long fac = c*(f.c);
  *this = Qideal(gens);
  c *= fac;
  ac *= fac;
  nm *= (fac*fac);
  iclass=index=-1;
}

int Qideal::contains(const Quad& alpha) const
{
  return (alpha.i%c==0) && ((alpha.r-b*alpha.i)%ac==0);
}

// return 1 iff this is coprime to I; if so, set r in this and s in I with r+s=1
int Qideal::is_coprime_to(const Qideal&J, Quad&r, Quad&s) const
{
  vector<long> v = {ac, J.ac, c*J.c*(b-J.b)}, w;
  int t = vecbezout(v, w);
  if (t==1)
    {
      r =   zcombo(w[0],  J.c * w[2]);
      s = J.zcombo(w[1],   -c * w[2]);
      assert (r+s==1);
    }
  return (t==1);
}

// return 1 iff this is coprime to alpha; if so, set inverse so an inverse of alpha modulo this
int Qideal::is_coprime_to(const Quad& alpha, Quad& inverse)
{
  if (alpha==0) return 0;
  Quad r, s;
  int t = is_coprime_to(Qideal(alpha), r, s);
  if(t==1) // then r is in this and s in (alpha) with r+s=1
    {
      inverse = reduce(s/alpha);
      assert (divides(alpha*inverse-1));
    }
  return t;
}

// reduction of alpha modulo this ideal
Quad Qideal::reduce(const Quad& alpha)
{
  if (iclass==-1)
    fill();
  //  cout<<"Reducing "<<alpha<<" mod "<<(*this)<<" with gens "<<gens()<<endl;
  if ((quadconj(g0)*g1).im()<0)
    cout<<"Badly oriented Z-basis "<<endl;
  return reduce_mod_zbasis(alpha, g0, g1);
}

// The i'th residue is Quad(x,y) for x mod a*c, y mod c; i = a*c*y+x

// Map from i to res (only depends on i mod norm)
Quad Qideal::resnum(int i) // the i'the residue mod this, in standard order (0'th is 0)
{
  std::ldiv_t qr = ldiv(posmod(i, nm), ac);
  return reduce(Quad(qr.rem, qr.quot));
}

int Qideal::numres(const Quad& alpha) // the index of a residue mod this, in standard order (0'th is 0)
{
  long y = posmod(alpha.im(), c);
  long x = posmod(alpha.re()-b*(alpha.im()-y), ac);
  return x + ac*y;
}

vector<Quad> Qideal::residues() // list of residues, sorted
{
  vector<Quad> res; res.reserve(nm);
  for (long i=0; i<nm; i++) res.push_back(resnum(i));
  return res;
}

pair<vector<Quad>, vector<Quad>> Qideal::invertible_residues() // lists of invertible residues, and their inverses
{
  vector<Quad> res, inv;
  res.reserve(nm);
  inv.reserve(nm);
  for (long i=0; i<nm; i++)
    {
      Quad r = resnum(i), s;
      if (is_coprime_to(r, s))
        {
          assert (divides(r*s-1));
          res.push_back(r);
          inv.push_back(s);
        }
    }
  return {res, inv};
}

// next two functions needed by ideal_prod_coeffs
Quad Qideal::princprod_coeff_alpha(vector<long>&z) const
{
  Quad ans = Quad( z[3]*a + z[0]*b, z[0] + z[4]*a ) * c;
  if (ans != zcombo( z[3] - z[4]*b, z[0] + z[4]*a))
    cerr << "Warning: possible overflow in princprod_coeff_alpha!"<<endl;
  return ans;
}

Quad Qideal::princprod_coeff_beta(vector<long>&z) const
{
  return zcombo( z[1], z[2]);
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
  std::ldiv_t qr = ldiv(c, na);
  if (qr.rem!=0)
    {
      cerr<<"***inexact division of "<<*this<<" by integer "<<n<<" ***"<<endl;
    }
  c = qr.quot;
  ac = a*c;
  nm = ac*c;
  if (iclass!=-1) { g0/=na; g1/=na; }
}

void Qideal::operator/=(const Quad&alpha)
{
  (*this) *= alpha.conj();
  long na = alpha.norm();
  std::ldiv_t qr = ldiv(c, na);
  if (qr.rem!=0)
    {
      cerr << "***inexact ideal division of "<<*this<<" by Quad "<<alpha<<" ***"<<endl;
    }
  c = qr.quot;
  ac = a*c;
  nm = ac*c;
  if (iclass!=-1) { g0/=na; g1/=na; }
}

void Qideal::operator/=(const Qideal&f)
{ Qideal keep = *this;
  (*this) *= f.conj();
  long nf = f.norm();
  std::ldiv_t qr = ldiv(c, nf);
  if (qr.rem!=0)                // shouldn't happen
    {
      cerr << "***inexact division of "<<keep<<" by ideal "<<f<<" of norm "<<nf<<" ***"<<endl;
    }
  c = qr.quot;
  ac = a*c;
  nm = ac*c;
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

int Qideal::is_principal()
{
  if (iclass==-1) fill();
  return (iclass==0);
}

Qideal Qideal::conj() const
{ Qideal ans=*this;
  ans.b= posmod(-b - Quad::t, a);
  if (iclass!=-1)
    {
      ans.g0 =  g0.conj();
      ans.g1 = -g1.conj(); // preserve orientation
    }
  return ans;
}

istream& operator>>(istream& s, Qideal& y)
{
  long a, b, c;
  do {
    s >> a >> b >> c;
    if ((a!=0) && (c!=0) && ((b*(b+Quad::t)+1)%a==0)) // i.e. N(b+w)=0 (mod a)
      {
        y = Qideal(a,b,c);
        return s;
      }  // assign via c'tor
    else
      {
        cerr << "retry: ";
      }
  }
  while (1);
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

vector<Qideal> primitive_ideals_with_norm(long N, int both_conj)
// N is the norm of a primitive ideal iff it has no inert prime
// factors and ramified primes divide it at most once
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
  // They have HNF [N,b+w] where b (mod N) satisfies
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
  stringstream ss(s);
  string Nstr, istr;
  std::getline(ss, Nstr, '.');
  std::getline(ss, istr);
  long N, i;
  stringstream(Nstr)>>N;
  stringstream(istr)>>i;
  *this = Qideal_from_norm_index(N,i);
}

// END OF FILE qideal.cc

