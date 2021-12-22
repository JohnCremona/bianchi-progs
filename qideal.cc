// FILE qideal.cc

#include "intprocs.h"
#include "mat22.h"
#include "qideal.h"

////////////////////////////////////
// implementation of class Qideal //
////////////////////////////////////

// private -- converts output from findzbasis to standard Qideal basis
// Here basis = {a,b,c} where the Z-basis is [a,b], [c,0]

void Qideal::abc_from_HNF(const vector<QUINT>& basis)
{ c = abs(basis[1]);
  a = abs(basis[2]/c);
  b = posmod(basis[0]/c,a);
  ac = a*c;
  nm = ac*c;
  if (!ok())
    { cerr <<"***Warning: "<< *this <<" not ok in abc_from_HNF***"<<endl;}
}

Quad Qideal::gen()  // smallest element, so a generator iff principal
{
  fill();
  if (!pos(g0))
    cout<<"After fill(), generator g0 = "<<g0<<" not normalised"<<endl;
  return g0;
}

vector<Quad> Qideal::gens() // reduced Z-module gens
{
  fill();
  return {g0, g1};
}

// private -- fills other data fields given a,b,c
void Qideal::fill()
{
  if (iclass!=-1) return;
  //  cout<<"Filling in data for ideal"<<(*this)<<endl;
  g0=Quad(a); g1=Quad(b,BIGINT(1));
  unimod U;
  sl2z_reduce(g0,g1, U);
  if ((quadconj(g0)*g1).im()<0)
    cout<<"Badly oriented Z-basis in fill() 1"<<endl;
  //  cout<<"sl2z_reduce for the primitive part returns g0="<<g0<<", g1="<<g1<<endl;
  iclass = (div(g0,g1)? 0: 1); // 0 means principal
  // scale by content:
  g1 *= c;
  g0 *= c;
  // normalize g0 by simultaneous scaling by units:
  while(!pos(g0))
    {
      g0*=fundunit;
      g1*=fundunit;
    }
  //  cout<<" -- now gens are "<<gens()<<endl;
  if ((quadconj(g0)*g1).im()<0)
    cout<<"Badly oriented Z-basis in fill() 2"<<endl;
  if (!pos(g0))
    cout<<"After fill(), generator g0 = "<<g0<<" not normalised"<<endl;
}

Qideal::Qideal()   //default ideal is whole ring
  : a(BIGINT(1)), b(BIGINT(0)), c(BIGINT(1)), ac(BIGINT(1)), nm(BIGINT(1)), iclass(0), g0(Quad::one), g1(Quad::w), index(1), F(0)
{ ; }

Qideal::Qideal(const Qideal& i)   // the copy constructor
  : a(i.a), b(i.b), c(i.c), ac(i.ac), nm(i.nm), iclass(i.iclass), index(i.index), F(0)
{ // don't copy the pointer F (though we could copy the underlying Factorizarion)
  if(iclass!=-1)
    {
      g0=i.g0; g1=i.g1;
    }
}

Qideal::~Qideal()
{
  // if (F!=0)
  //   {
  //     delete F;
  //     F = 0;
  //   }
}


Qideal::Qideal(const QUINT& aa, const QUINT& bb, const QUINT& cc)
{
  if (::is_zero(cc))
    {
      a=BIGINT(1); b=c=ac=nm=BIGINT(0); iclass=0; g0=g1=Quad::zero; index=1;
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
  F=0;
}

Qideal::Qideal(const QUINT& n)       // princ ideal gen by QUINT
  : a(BIGINT(1)), b(BIGINT(0)), c(abs(n)), iclass(0), g0(abs(n)), g1(BIGINT(0), abs(n)), index(-1), F(0)
{
  ac=a*c; nm=ac*c;
}

// Utility used to construct Qideals given two lists of integers rv,
// iv, or arbitrary but equal length, defining the Z-module spaned by
// all [rv[i], iv[i]]

Qideal::Qideal(const vector<QUINT>& rv, const vector<QUINT>& iv)  // ideal from Z-gens
{
  vector<QUINT> basis;
  findzbasis(rv,iv,basis);
  // cout<<"basis = "<<basis<<endl;
  abc_from_HNF(basis);
  iclass=-1;
  index=-1;
  F=0;
}

Qideal::Qideal(const vector<Quad>& gens)       // ideal spanned by list of Quads
{
  vector<QUINT> rv, iv;
  for(vector<Quad>::const_iterator g = gens.begin(); g!=gens.end(); ++g)
    {
      Quad a = *g;
      rv.push_back(a.re());
      iv.push_back(a.im());
      a *=Quad::w;
      rv.push_back(a.re());
      iv.push_back(a.im());
    }
  // cout<<"rv = "<<rv<<", iv = "<<iv<<endl;
  *this = Qideal(rv, iv);
}

Qideal::Qideal(const Quad& alpha) // principal ideal
{
  if (alpha.nm==0)
    *this = Qideal(0);
  else
    {
      vector<Quad> gens(1,alpha);
      *this = Qideal(gens);
    }
  iclass=0;
  g0=makepos(alpha); g1=Quad::w*g0;
  index=-1;
  F=0;
}


////////////////////////////////////
// Operators for the class Qideal //
////////////////////////////////////

Qideal Qideal::operator+(const QUINT& d) const
{
  if (is_zero())
    return Qideal(d);
  if (d==0)
    return *this;
  Qideal ans=*this;
  ans += d;
  return ans;
}

Qideal Qideal::operator+(const Quad& alpha) const
{
  if (is_zero())
    return Qideal(alpha);
  if (alpha.is_zero())
    return *this;
  Qideal ans=*this;
  ans+=alpha;
  return ans;
}

Qideal Qideal::operator+(const Qideal& f) const
{
  if (is_zero())
    return f;
  if (f.is_zero())
    return *this;
  Qideal ans=*this;
  ans+=f;
  return ans;
}

Qideal operator+(const QUINT&a, const Qideal& f)
{ return f+a; }

Qideal operator+(const Quad&a, const Qideal& f)
{ return f+a; }

void Qideal::operator+=(const QUINT& aa)
{
  if (is_zero())
    {
      *this = Qideal(aa);
      return;
    }
  if (divides(aa))
    return;  //ideal remains unchanged
  vector<QUINT> rv = get_rv(), iv = get_iv();
  bigint z(0);
  rv.push_back(aa);
  rv.push_back(z);
  iv.push_back(z);
  iv.push_back(aa);
  *this = Qideal(rv, iv);
}

void Qideal::operator+=(const Quad& alpha)
{
  if (is_zero())
    {
      *this = Qideal(alpha);
      return;
    }
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

  vector<QUINT> rv = get_rv(), iv = get_iv(), rva = A.get_rv(), iva = A.get_iv();
  rv.insert(rv.end(), rva.begin(), rva.end());
  iv.insert(iv.end(), iva.begin(), iva.end());
  *this = Qideal(rv, iv);
}

void Qideal::operator+=(const Qideal& f)
{
  if (is_zero())
    {
      *this = f;
      return;
    }
  if (divides(f))
    return;  //ideal remains unchanged
  if (f.divides(*this))
    {
      *this=f;
      return;
    }
  vector<QUINT> rv1 = get_rv(), iv1 = get_iv();
  vector<QUINT> rv2 = f.get_rv(), iv2 = f.get_iv();
  rv1.insert( rv1.end(), rv2.begin(), rv2.end() );
  iv1.insert( iv1.end(), iv2.begin(), iv2.end() );
  *this = Qideal(rv1, iv1);
}

Qideal Qideal::operator*(const QUINT& d) const
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

Qideal Qideal::operator*(Qideal& f) const
{
  Qideal ans(*this);
  ans*=f;
  return ans;
}

Qideal operator*(const QUINT& d, const Qideal& f)
{ return f*d; }

Qideal operator*(const Quad& a, const Qideal& f)
{ return f*a; }

void Qideal::operator*=(const QUINT& d)
{
  if (d==0)
    {
      *this=Qideal(0);
      return;
    }
  QUINT dd = abs(d);
  c *= dd;
  ac *= dd;
  nm *= (d*d);
  if (iclass!=-1)
    {
      g0*=dd;
      g1*=dd;
    }
}

void Qideal::operator*=(const Quad& alpha)
{
  if (c==0 || alpha.norm()==1) return;              // this is unchanged
  if (alpha.is_zero()) { *this=Qideal(alpha); return;}         // this becomes 0
  if (alpha.im()==0) {*this *= alpha.re(); return;} // alpha in Z
  if (iclass==0) {
    // cout<<"operator *= with this="<<(*this)<<" (which is principal with generator "<<g0<<") and "<<alpha<<endl;
    *this=Qideal(g0*alpha);
    // cout<<" product is "<<(*this)<<" with generator "<<g0<<endl;
    iclass=0;
    return;
  } // this is principal

  vector<Quad> gens = {alpha*Quad(a), alpha*Quad(b,BIGINT(1))}; // without the factor c
  QUINT fac = c;
  *this = Qideal(gens);
  c *= fac;
  ac *= fac;
  nm *= (fac*fac);
  if (iclass!=-1)
    {
      g0*=alpha;
      g1*=alpha;
      while(!pos(g0))
        {
          g0*=fundunit;
          g1*=fundunit;
        }
    }
}

void Qideal::operator*=(Qideal& f)
{
  if (c==0 || f.nm==1) return;    // unchanged
  if (f.c==0 || nm==1) { *this=f; return; }  // f unchanged

  if (is_principal())
    {
      if (f.is_principal())
        {
          *this = Qideal(g0*f.g0);
          return;
        }
      else
        {
          *this = Qideal({g0*f.g0, g0*f.g1});
          return;
        }
    }
  else
    {
      if (f.is_principal())
        {
          *this = Qideal({g0*f.g0, g1*f.g0});
          return;
        }
    }
  bigint z1(1);
  Quad x1(a), x2(b,z1), y1(f.a), y2(f.b,z1);
  vector<Quad> gens = {x1*y1, x1*y2, x2*y1, x2*y2};
  // cout<<"operator *= with this = "<<(*this)<<" and "<<f<<endl;
  // cout<<"primitive gens of this: "<<x1<<", "<<x2<<endl;
  // cout<<"primitive gens of that: "<<y1<<", "<<y2<<endl;
  // cout<<"primitive gens of product: "<<gens<<endl;
  QUINT fac = c*(f.c);
  *this = Qideal(gens);
  // cout<<" - primitive product is "<<(*this)<<endl;
  c *= fac;
  ac *= fac;
  nm *= (fac*fac);
  // cout<<" - full product is "<<(*this)<<endl;
  iclass=index=-1;
}

int Qideal::contains(const Quad& alpha) const
{
  if (is_zero())
    return alpha.is_zero();
  else
    return div(c,alpha.i) && div(ac,alpha.r-b*alpha.i);
}

// return 1 iff this is coprime to J; if so, set r in this and s in J with r+s=1
int Qideal::is_coprime_to(Qideal&J, Quad&r, Quad&s)
{
  vector<QUINT> v = {ac, J.ac, c*J.c*(b-J.b)}, w;
  if (vecbezout(v, w)!=1)
    return 0;
  r =   zcombo(w[0],  J.c * w[2]);
  s =   J.zcombo(w[1], -c * w[2]);
  assert (r+s==Quad::one);
  assert (contains(r));
  assert (J.contains(s));
  if (r.nm<0 || s.nm<0) // can only happen if there was overflow
    {
      cout<<"I="<<(*this)<<", J="<<J<<": initial r="<<r<<" (norm "<<r.nm<<"), s=1-r="<<s<<" (norm "<<s.nm<<")"<<"; N(IJ) = "<<r.nm*J.nm<<endl;
      cout<<"v = "<<v<<", w="<<w<<endl;
      exit(1);
    }
  Qideal IJ = (*this)*J;
  Quad r1 =   IJ.reduce(r);
  // cout << "I="<<(*this)<<", J="<<J<<", IJ="<<IJ<<" (norm "<<IJ.nm<<"): r="<<r<<" (norm "<<r.nm<<")";
  // cout << " reduces mod IJ to "<<r1<<endl;
  assert (contains(r1));
  assert (J.contains(Quad::one-r1));
  r = r1;
  s = Quad::one-r;
  // cout<<"I="<<(*this)<<", J="<<J<<": r="<<r<<", s=1-r="<<s<<endl;
  return 1;
}

// return 1 iff this is coprime to alpha; if so, set inverse so an inverse of alpha modulo this
int Qideal::is_coprime_to(const Quad& alpha, Quad& inverse)
{
  if (alpha.is_zero()) return 0;
  Quad r, s;
  Qideal A(alpha);
  int t = is_coprime_to(A, r, s);
  if(t==1) // then r is in this and s in (alpha) with r+s=1
    {
      inverse = reduce(s/alpha);
      assert (divides(alpha*inverse-Quad::one));
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
Quad Qideal::resnum(long i) // the i'the residue mod this, in standard order (0'th is 0)
{
  bigint quot, rem;
  ::divides(posmod(BIGINT(i), nm), ac, quot, rem);
  return reduce(Quad(rem, quot));
}

long Qideal::numres(const Quad& alpha) const // the index of a residue mod this, in standard order (0'th is 0)
{
  QUINT y = posmod(alpha.im(), c);
  QUINT x = posmod(alpha.re()-b*(alpha.im()-y), ac);
  return I2long(x + ac*y);
}

vector<Quad> Qideal::residues() // list of residues, sorted
{
  long i, n = I2long(nm);
  vector<Quad> res; res.reserve(n);
  for (i=0; i<n; i++) res.push_back(resnum(i));
  return res;
}

// lists of invertible residues, and their inverses
pair<vector<Quad>, vector<Quad>> Qideal::invertible_residues()
{
  long i, n = I2long(nm);
  vector<Quad> res, inv;
  res.reserve(n);
  inv.reserve(n);
  for (i=0; i<n; i++)
    {
      Quad r = resnum(i), s;
      if (is_coprime_to(r, s))
        {
          assert (divides(r*s-Quad::one));
          res.push_back(r);
          inv.push_back(s);
        }
    }
  return {res, inv};
}

// An AB-matrix with given first column
mat22 AB_matrix(const Quad& a, const Quad& c)
{
  Quad b, d;
  if (!quadbezout(a,c,d,b).is_zero()) // then (a,c)=(g) and a*d+b*c=g
    {
      return mat22(a,-b,c,d);
    }
  // otherwise (a,c) is not principal (and a,c are nonzero)
  Qideal I({a,c});
  Qideal I0 = Qideal(a)/I, I1 = Qideal(c)/I;
  Quad r0, r1;
  int t = I0.is_coprime_to(I1, r0, r1);
  assert (t && "I0, I1 coprime");
  QUINT g = I.norm();
  b = -(r1*g)/c;
  assert (b*c == -r1*g);
  d =  (r0*g)/a;
  assert (a*d == r0*g);
  if (Qideal({b,d}) != I.conj())
    {
      cerr<<"a = "<<a<<", Norm(a) = "<<quadnorm(a)<<endl;
      cerr<<"c = "<<c<<", Norm(c) = "<<quadnorm(c)<<endl;
      cerr<<"I = (a,c) = "<<I<<", g = Norm(I) = "<<g<<endl;
      cerr<<"a/I = "<<I0<<endl;
      cerr<<"c/I = "<<I1<<endl;
      cerr<<"r0 = "<<r0<<" (norm "<<quadnorm(r0)<<"), r1 = 1-r0 = "<<r1<<" (norm "<<quadnorm(r1)<<")"<<endl;
      cerr<<"b = -r1*g/c = "<<b<<endl;
      cerr<<"d =  r0*g/a = "<<d<<endl;
      cerr<<"conj(I) = "<<I.conj()<<endl;
      cerr<<"(b,d) = "<<Qideal({b,d})<<endl;
    }
  assert (Qideal({b,d}) == I.conj());
  mat22 M(a, b, c, d);
  assert (M.det()==g);
  return M;
}

// Assuming this*J is principal, sets g to a generator and returns a
// 2x2 matrix of determinant g whose columns generate this and J,
// the first column being (g0,g1)
mat22 Qideal::AB_matrix(Qideal&J, Quad&g)
{
  if (is_principal())
    {
      J.fill();
      g = g0.nm;
      return mat22(g0,Quad::zero,Quad::zero,J.g0);
    }
  Qideal IJ = (*this)*J, C = (*this).conj();
  if (!IJ.is_principal(g))
    cerr<<"ideals "<<(*this)<<" and "<<J<<" do not have principal product in AB_matrix()"<<endl;
  Qideal I0 = g0*C/nm, I1 = g1*C/nm; // = (g0)/this, (g1)/this
  Quad r0, r1;
  I0.is_coprime_to(I1, r0, r1);
  Quad h0 = -(r1*g)/g1;
  Quad h1 =  (r0*g)/g0;
  mat22 M(g0, h0, g1, h1);
  assert (M.det()==g);
  assert (Qideal({h0,h1}) == J);
  // cout << "M = "<<M<<" with det "<<M.det()<<endl;
  return M;
}

mat22 Qideal::AB_matrix()
{
  if (is_principal())
    {
      return mat22(g0,Quad::zero,Quad::zero,g0.conj());
    }
  Qideal C = (*this).conj();
  Qideal I0 = (g0*C)/nm, I1 = (g1*C)/nm; // = (g0)/this, (g1)/this
  Quad r0, r1;
  I0.is_coprime_to(I1, r0, r1);
  Quad h0 = -(r1*nm)/g1;
  Quad h1 =  (r0*nm)/g0;
  mat22 M(g0, h0, g1, h1);
  assert (M.det()==nm);
  assert (Qideal({h0,h1}) == C);
  return M;
}


Qideal Qideal::operator/(const QUINT&n) const
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
Qideal operator/(const QUINT&n, const Qideal&f)
{ Qideal ans=f.conj()*n;
  return ans/(f.norm());
}

// more efficient than first promoting dividend, though doubtful if ever used
Qideal operator/(const Quad&alpha, const Qideal&f)
{ Qideal ans=f.conj()*alpha;
  return ans/(f.norm());
}

void Qideal::operator/=(const QUINT&n)
{ QUINT na=abs(n);
  if (na==1) return;
  // cout<<"applying operator/= to ideal "<<(*this)<<" and "<<n<<endl;
  bigint quot, rem;
  if (!::divides(c, na, quot, rem))
    {
      cerr<<"***inexact division of "<<*this<<" by integer "<<n<<" ***"<<endl;
      exit(1);
    }
  c = quot;
  ac = a*c;
  nm = ac*c;
  if (iclass!=-1) { g0/=na; g1/=na; }
  // cout<<" -- after /= this = "<<(*this)<<endl;
}

void Qideal::operator/=(const Quad&alpha)
{
  if (alpha.nm==BIGINT(1)) return;
  (*this) *= alpha.conj();
  QUINT na = alpha.norm();
  bigint quot, rem;
  if (!::divides(c, na, quot, rem))
    {
      cerr << "***inexact ideal division of "<<*this<<" by Quad "<<alpha<<" ***"<<endl;
      exit(1);
    }
  c = quot;
  ac = a*c;
  nm = ac*c;
  if (iclass!=-1) { g0/=na; g1/=na; }
}

void Qideal::operator/=(const Qideal&f)
{
  QUINT nf = f.norm();
  if (nf==BIGINT(1)) return;
  Qideal keep = *this, fc = f.conj();
  // cout<<"dividing "<<(*this)<<" by "<<f<<endl;
  (*this) *= fc;
  //cout<<" - after multiplying by the conjugate: "<<(*this)<<endl;
  bigint quot, rem;
  if (!::divides(c, nf, quot, rem))
    {
      cerr << "***inexact division of "<<keep<<" by ideal "<<f<<" of norm "<<nf<<" ***"<<endl;
      exit(1);
    }
  c = quot;
  ac = a*c;
  nm = ac*c;
  if (iclass!=-1) { g0/=nf; g1/=nf; }
}

Qideal Qideal::intersection(const Qideal& I)
{
  if (contains(I)) return I;
  if (I.contains(*this)) return *this;
  Qideal G = I/((*this)+I);
  return (*this)*G;
}

// with nonzero a in this, return b such that this=(a,b)
Quad Qideal::second_generator(const Quad& a)
{
  if (a.is_zero() or not contains(a))
    {
      cerr << "method second_generator(a) called with a="<<a<<" not a nonzero element of "<<(*this)<<endl;
      return Quad::zero;
    }
  Quad b, d;
  // if this is principal ignore a and return a generator
  if (is_principal(b)) return b;
  // now a itself cannot be a generator, so (a)=I*J for some J
  Qideal J = Qideal(a)/(*this);
  // find an ideal in inverse class to this, coprime to J
  equivalent_coprime_to(J, b, d, 1); // this*J=(b) (and d=1)
  // now (a,b) = (a)+(b) = I*(J+P) = I with I=this, P=ideal in previous line (discarded)
  return b;
}

// return J such that J^2 is equivalent to this (or J^2*this is
// principal if anti==1), or if no such J exists (i.e., if the ideal
// class is not a square), return the zero ideal.  (implemented in
// primes.cc)
Qideal Qideal::sqrt_class(int anti)
{
  for (vector<Qideal>::const_iterator A = Quad::class_group.begin(); A!=Quad::class_group.end(); ++A)
    {
      Qideal Asq = *A;
      Asq *= Asq;
      if (anti)
        Asq = Asq.conj();
      if (is_equivalent(Asq))
        return *A;
    }
  return Qideal(Quad::zero);
}

mat22 Qideal::AB_matrix_of_level(const Qideal&N, Quad&g)
{
  if (is_principal())
    {
      g = g0.nm;
      return mat22(g0,Quad::zero,Quad::zero,g0.conj());
    }
  // find a small element of this intersection N
  Qideal L = intersection(N);
  L.fill();
  Quad a = L.g0, b, d, r, s;
  Qideal J = Qideal(a)/(*this); // I*J=(a)
  Qideal P = equivalent_coprime_to(J, b, d, 1); // I*P=(b) with J+P=(1)
  J.is_coprime_to(P, r, s); // r+s=1, r in J, s in P
                            // so r*b in I*J*P = (a)*P
  g = b;
  return mat22(b, -r*(b/a), a, s);
}

mat22 Qideal::AB_matrix_of_level(const Qideal&J, const Qideal&N, Quad&g)
{
  // find a small element of this intersection N:
  Qideal L = intersection(N);
  L.fill();
  Quad cc = L.g0; // smallest nonzero element of I also in N
  Quad aa = second_generator(cc);
  assert (Qideal({aa,cc}) == *this);

  // find the determinant of the matrix to be constructed:
  Qideal IJ = J*(*this);
  IJ.fill();
  g = IJ.g0;

  // Construct the second column:
  Qideal Iconj = conj();
  Qideal I0 = aa*Iconj/nm;
  Qideal I1 = cc*Iconj/nm;
  Quad r0, r1;
  I0.is_coprime_to(I1, r0, r1);
  Quad bb = -(r1*g)/cc;
  Quad dd =  (r0*g)/aa;
  mat22 M(aa, bb, cc, dd);
  assert (M.det()==g);
  assert (Qideal({bb,dd}) == J);
  assert (N.divides(M.entry(1,0)));
  return M;
}

////////////////////////////////////////////////////////
// Qideal:: non-operator member functions and friends //
////////////////////////////////////////////////////////

int Qideal::ok() const
{
  return div(a, b*(b + Quad::t) + Quad::n);
}

int Qideal::is_principal()
{
  if (iclass==-1) fill();
  return (iclass==0);
}

int Qideal::is_principal(Quad& g)
{
  if (iclass==-1) fill();
  if (iclass==0) {g=g0; return 1;}
  return 0;
}

Qideal Qideal::conj() const
{ Qideal ans=*this;
  ans.b= posmod(-b - Quad::t, a);
  ans.index = (b==ans.b? index: -1);
  if (iclass!=-1)
    {
      ans.g0 =  g0.conj();
      ans.g1 = -g1.conj(); // preserve orientation
      while(!pos(ans.g0))
        {
          ans.g0*=fundunit;
          ans.g1*=fundunit;
        }
    }
  return ans;
}

// Ideal input: either a label string N.i, or as an alternative for a
// principal ideal, two integers a b for (a+b*w)

istream& operator>>(istream& s, Qideal& y)
{
  string st;
  s >> st;
  if (st.find(".")==string::npos) // string contains no "."
    {
      QUINT r = BIGINT(st.c_str()), i;
      s >> i;
      Quad alpha(r,i);
      y = Qideal(alpha);
    }
  else
    y = Qideal(st);
  return s;
}

ostream& operator<<(ostream& s, const Qideal& x)
{
  QUINT a=x.a, b=x.b, c=x.c;
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
  if ((dividend.norm()==0) || (factor.norm()<=1))
    {
      cerr<<"Warning: 9999 returned in Qideal version of val"<<endl;
      return 9999;
    }
  long e = 0;
  Qideal n=dividend;
  while (factor.divides(n)) {n/=factor; e++; }
  return e;
}

vector<Qideal> primitive_ideals_with_norm(QUINT N, int both_conj)
// N is the norm of a primitive ideal iff it has no inert prime
// factors and ramified primes divide it at most once
{
  vector<Qideal> ans;
  vector<QUINT> pdivs_norm = pdivs(N);
  for(vector<QUINT>::const_iterator pi=pdivs_norm.begin(); pi!=pdivs_norm.end(); ++pi)
    {
      QUINT p=*pi;
      int s = Quad::chi(p);
      if ((s==-1) || ((s==0) && val(p,N)>1)) return ans;
    }
  // now the local tests pass so there are primitive ideals of norm N.
  // They have HNF [N,b+w] where b (mod N) satisfies
  //     b^2+t*b+n = 0 (mod N)
  //
  // Stupid implementation here: try all b mod N
  QUINT t(Quad::t), n(Quad::n), z1(1);
  QUINT b, maxb = (both_conj? N-1: (N-t)/2);  // rounded down
  for (b=0; b<=maxb; b++)
    if (div(N,(b*(b + t) + n)))
      ans.push_back(Qideal(N,b,z1));
  //cout<<" primitive ideals with norm "<<N<<" both_conj="<<both_conj<<"): "<<ans<<endl;
  return ans;
}

vector<Qideal> ideals_with_norm(QUINT N, int both_conj)
// I = c*I0 uniquely, where c^2|N and I0 is primitive of norm N/c^2
{
  vector<Qideal> ans;
  vector<QUINT> clist = sqdivs(N); // list of c such that c^2|N
  for (vector<QUINT>::const_iterator ci=clist.begin(); ci!=clist.end(); ++ci)
    {
      QUINT c = *ci;
      vector<Qideal> primitives = primitive_ideals_with_norm((N/c)/c, both_conj);
      for (vector<Qideal>::const_iterator Ji = primitives.begin(); Ji!=primitives.end(); ++Ji)
        ans.push_back(c*(*Ji));
    }
  //cout<<" ideals with norm "<<N<<" both_conj="<<both_conj<<"): "<<ans<<endl;
  return ans;
}

vector<Qideal> ideals_with_bounded_norm(QUINT maxnorm, int both_conj)
{
  vector<Qideal> ans;
  for (QUINT N=BIGINT(1); N<=maxnorm; N++)
    {
      vector<Qideal> IN = ideals_with_norm(N, both_conj);
      ans.insert(ans.end(), IN.begin(), IN.end());
    }
  return ans;
}

// NB The label of an ideal is a string of the form N.i where N is the
// norm and i its index *based at 1*.

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
  QUINT n = norm();
  vector<Qideal> II = Qideal_lists::ideals_with_norm(n);
  vector<Qideal>::iterator i = find(II.begin(), II.end(), *this);
  if (i==II.end())
    {
      cerr<<"Error in finding index of "<<(*this)<<" (norm "<<n<<") in list of all ideals of norm "<<n<<": "
          <<II<<endl;
      exit(1);
    }
  index = std::distance(II.begin(), i) + 1;
  // cout<<"In set_index(), have just set index of ideal "<<(*this)<<" to "<<index<<endl;
}

string ideal_label(Qideal& I)  // returns label of ideal I
{
  stringstream s;
  s << I.norm() << "." << I.get_index();
  return s.str();
}

string gens_string(Qideal& I)  // returns string of gens, of the form (x) if principal or (x,y) ideal I
{
  stringstream s;
  s << "(";
  if (I.is_principal())
    s << I.g0;
  else
    s << I.g0 <<","<< I.g1;
  s << ")";
  return s.str();
}

Qideal Qideal_from_norm_index(QUINT N, int i) // i'th ideal of norm N
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
  return II[i-1];  // offset since labels are indexed from 1
}

Qideal::Qideal(const string& s)           // ideal from label N.i
// need to parse s to obtain N and i as postive integers
{
  stringstream ss(s);
  string Nstr, istr;
  std::getline(ss, Nstr, '.');
  std::getline(ss, istr);
  QUINT N;
  stringstream(Nstr)>>N;
  if (N==0)
    *this = Qideal(0);
  else
    {
      long i;
      stringstream(istr)>>i;
      *this = Qideal_from_norm_index(N,i);
    }
}

// If a is in (c,d) return 1 and x,y such that a=c*x+d*y, else return 0
int express2gens(const Quad& a, const Quad& c, const Quad& d, Quad& x, Quad& y)
{
  Qideal I({c,d});
  if (!I.contains(a)) return 0;
  Qideal J = I.conj();
  Qideal I0 = (c*J)/I.norm(),  I1 = (d*J)/I.norm();
  I0.is_coprime_to(I1, x, y); // x+y=1
  // Now x*I subset (c) and y*I subset (d), so if a in I then a*x/c, a*y/d are integral
  x = (a*x)/c;
  y = (a*y)/d;
  assert (c*x+d*y==a);
  return 1;
}

// largest factor of this coprime to I
Qideal Qideal::divide_out(const Qideal& I)
{
  Qideal J = *this;
  Qideal G = J+I;
  while (G.nm!=1)
    {
      J /= G;
      G = J+I;
    }
  return J;
}

// Return 1 iff this is coprime to (c,d); if so, set x,y so c*x+d*y =1
// modulo this ideal.  If fix=1, ensure that y is coprime to this ideal.
int Qideal::is_coprime_to(const Quad& c, const Quad& d, Quad& x, Quad& y, int fix)
{
  Qideal I({c,d});
  Quad r,s,t,u;
  if (!is_coprime_to(I, r, s))
    return 0;

  // else r+s=1 with r in this, s in (c,d)
  express2gens(s, c, d, x, y);  // now c*x+d*y = s = 1-r
  x = reduce(x);
  y = reduce(y); // does not change c*x+d*y modulo this
  //  cout<<"Before adjustment, x="<<x<<" and y="<<y<<endl;
  if (fix && !is_coprime_to(y))
    {
      // Fixing is needed: replace (x,y) by (x-t*d,y+t*c) with y+t*c
      // coprime to this. Here, t must be divisible by all primes
      // dividing this not dividing y, but none of those dividing y.

      Qideal A(y);
      divide_out(A).equivalent_coprime_to(A, t, u, 1);
      y = reduce(y+t*c);
      x = reduce(x-t*d);
      assert (is_coprime_to(y));
      //      cout<<"After adjustment with t="<<t<<", x="<<x<<" and y="<<y<<endl;
    }
  assert (contains(c*x+d*y-Quad::one));
  return 1;
}

// brute force test whether a is a square of some element of reslist, mod M

int squaremod(const Quad& a, const Qideal& M, const vector<Quad>& reslist)
{
  if (M.contains(a)) return 0;
  vector<Quad>::const_iterator r=reslist.begin();
  while(r!=reslist.end())
    {
      Quad res=*r++;
      if(M.contains(res*res-a)) return +1;
    }
  return -1;
}

vector<int> makechitable(const Qideal& L, const vector<Quad>& reslist)
{
  vector<int> chi;
  if(reslist.size()==1)
    chi.push_back(1);
  else
    {
      vector<Quad>::const_iterator r=reslist.begin();
      while(r!=reslist.end())
	chi.push_back(squaremod(*r++,L,reslist));
    }
  return chi;
}


// END OF FILE qideal.cc
