// File FIELD.CC: class for working with number fields for Hecke eigenvalues
//////////////////////////////////////////////////////////////////////////

#include "field.h"

//#define DEBUG_ARITH

Field* FieldQQ = new Field();

Field::Field(const ZZX& p, string a, int verb)
{
  if (IsMonic(p) && IsIrreducible(p))
    {
      // Set A to be the companion matrix of p.  The function
      // CompanionMatrix(p) returns a mat_ZZ so we do this manually.
      d = deg(p);
      mat_m m(d,d);
      ZZ one(1);
      for(int i=1; i<d; i++)
        {
          m(i+1,i) = one;
          m(i,d) = -coeff(p, i-1);
        }
      m(d,d) = -coeff(p, d-1);
      // Finally call the other constructor
      *this = Field(m, one, a, verb);
    }
}

Field::Field() // defaults to Q
{
  d=1;
  ZZ one(1);
  denom=one;
  A.init(1,1); // zero matrix
  B.init(1,1); // zero matrix
  C.init(1,1); // zero matrix
  B(1,1) = Bdet = Bdet1 = Bdet2 = Bdet3 = one;
  Binv = B;
  var = ""; // will not be used anyway
}

Field::Field(const mat_m& m, const ZZ& den, string a, int verb)
  : var(a), d(m.nrows()), denom(den), A(m)
{
  if (verb)
    {
      cout << "----------------------------"<<endl;
      cout << "In Field constructor" << endl;
    }
  minpoly = scaled_charpoly(mat_to_mat_ZZ(A), denom);
  if (verb)
    cout << " - min poly = " << ::str(minpoly) << ", generator " << var << endl;


  // Compute change of basis matrix B, with column j equal to
  // denom^(n-j)*A^(j-1)v for j from 1 to d
  vec_m v(d);
  v[1] = pow(denom,d-1); // so v=[1,0,...,0]*denom^(d-1)
  // cout<<"v = "<<v<<endl;
  B.init(d,d);
  B.setcol(1,v);
  for(int i=2; i<=d; i++)
    {
      v = A*v / denom;
      B.setcol(i,v);
    }
  Bcontent = B.content();
  // cout << "Content of original B is " << Bcontent << endl;
  B /= Bcontent;
  Binv.init(d,d);
  Bdet = inverse(B,Binv); // so B*Binv = Bdet*identity
  Bdet1 = B(1,1);
  Bdet2 = Binv(1,1);
  Bdet3 = Bdet2*denom;
  mat_m I = mat_m::identity_matrix(d);
  assert (B*Binv == Bdet*I);
  // Now we should have Binv*A*B = denom*Bdet * C, where C = companion
  // matrix of minpoly. i.e. the B-conjugate of A can be divided by
  // denom to give the integral matrix C.
  C = (Binv*A*B) / (denom*Bdet);
  assert (CompanionMatrix(minpoly) == mat_to_mat_ZZ(C));

  Cpowers.resize(d);
  Cpowers[0] = I;
  for (int i=1; i<d; i++)
    Cpowers[i] = C*Cpowers[i-1];
  if(verb)
    {
      cout<<"basis  matrix = ";
      output_flat_matrix(Binv);
      cout<<endl;
      cout<<"inverse basis matrix = ";
      output_flat_matrix(B);
      cout<<endl;
      cout << "basis factor  = " << Bdet1 << "*" << Bdet2 << "=" << Bdet <<endl;
      if (verb>1)
        {
          cout<<"companion matrix  = ";
          output_flat_matrix(C);
          cout<<endl;
        }
      cout << "Leaving Field constructor" << endl;
      cout << "----------------------------"<<endl;
    }
}

void Field::display_bases(ostream&s) const
{
  if (isQ())
    return;
  mat_m I(mat_m::identity_matrix(d));

  s << "Powers of A (i.e. powers of " << var << " in A-embedding):\n";
  vector<mat_m> Apowers(d);
  Apowers[0] = I;
  for (int i=1; i<d; i++)
    Apowers[i] = A*Apowers[i-1];
  for (auto Apow: Apowers)
    {
      Apow.output(s);
      s<<endl;
      s<<endl;
    }
  ZZ fac = pow(denom, d-1) * Bdet;
  s << "Basis in A-embedding, scaled by "<< denom << "^" << d-1 << " * " << Bdet << " = " << fac <<":\n";
  fac = Bdet*Bcontent;
  s << "(first columns should be standard basis vectors * "<< fac <<")\n";
  for(int i=1; i<=d; i++)
    {
      vec_m coli = Binv.col(i);
      for(int j=1; j<=d; j++)
        coli[j] *= pow(denom, d-j);
      mat_m M = lin_comb_mats(coli, Apowers);
      M.output(s);
      s<<endl;
      s<<endl;
      if (!(M.col(1)==fac*I.col(i)))
        {
          cout<<"Scaling problem:  column(1) = "<<M.col(1)
              <<" but fac*e_"<<i<<" = "<<fac*I.col(i)<<endl;
        }
    }
  s << "Powers of C (i.e. powers of " << var << " in C-embedding):\n";
  for (auto Cpow: Cpowers)
    {
      Cpow.output(s);
      s<<endl;
      s<<endl;
    }
  s << "Basis in C-embedding, scaled by "<<Bdet<<":\n";
  s << "(first columns should be columns of basis())\n";
  for(int i=1; i<=d; i++)
    {
      mat_m M = lin_comb_mats(Binv.col(i), Cpowers);
      M.output(s);
      s<<endl;
      s<<endl;
      assert (M.col(1) == Binv.col(i));
    }
}

void Field::display(ostream&s, int raw)
{
  static const ZZ one(1);
  string fpol = ::str(minpoly);
  if (isQ() || (d==1))
    {
      s << "Q" << endl;
      return;
    }
  s << "Q(" << var << ") with defining polynomial "<< fpol <<" of degree "<<d;
  if (d==2)
    s << ", discriminant "<<discriminant(minpoly);
  s << endl;
  if(raw && !isQ())
    {
      s << "   Raw basis with respect to alpha-power basis:\n";
      for(int i=1; i<=d; i++)
        {
          FieldElement bi = element(vec_m::unit_vector(d,i), one, 1);
          s << "   #"<<i<<": "<< bi << endl;
          s << "   (with char poly = " << ::str(bi.charpoly())
            << " and minpoly = " << ::str(bi.minpoly()) << "\n";
        }
    }
}

////////////////////////////////////////////////////////////////////////

// raw means the given coords are w.r.t. the B-basis
FieldElement::FieldElement( Field* HF, const vec_m& c, const ZZ& d, int raw)
    :F(HF), coords(c), denom(d)
{
  // cout<<"Constructing a FieldElement in Q("<<F->var<<")\n";
  if (raw)
    {
      coords = F->Binv *coords;
      denom *= F->Bdet3;
    }
  cancel();
  // cout<<" - finished constructing " << (*this) << "\n";
}

FieldElement Field::element(const vec_m& c, const ZZ& d, int raw)
{
  return FieldElement(this, c, d, raw);
}

FieldElement Field::rational(const bigrational& x)
{
  if (isQ())
    return FieldElement(x);
  else
    return FieldElement(this, x.num(), x.den());
}

FieldElement Field::rational(const ZZ& x)
{
  return rational(bigrational(x));
}

FieldElement Field::rational(long x)
{
  return rational(bigrational(x));
}

FieldElement Field::rational(int x)
{
  return rational(bigrational(x));
}

FieldElement Field::zero()
{
  return rational(0);
}

FieldElement Field::one()
{
  return rational(1);
}

FieldElement Field::minus_one()
{
  return rational(-1);
}

FieldElement Field::two()
{
  return rational(2);
}

FieldElement Field::minus_two()
{
  return rational(-2);
}

FieldElement Field::gen()
{
  if (isQ())
    return rational(1);
  else
    return FieldElement(this, vec_m::unit_vector(d, 2));
}

void FieldElement::cancel() // divides through by gcd(content(coords, denom))
{
  if (IsOne(denom))
    return;
  if (denom<0)
    {
      denom = -denom;
      coords = -coords;
    }
  ZZ g = gcd(content(coords), denom);
  if (IsOne(g))
    return;
  denom /=g;
  coords /= g;
}

int FieldElement::is_zero() const
{
  if (F->isQ())
    return val.is_zero();
  return trivial(coords);
}

int FieldElement::is_one() const
{
  static const bigrational one(1);
  if (F->isQ())
    return val==one;
  return IsOne(denom) && coords == vec_m::unit_vector(F->d,1);
}

int FieldElement::is_minus_one() const
{
  static const bigrational minus_one(-1);
  if (F->isQ())
    return val==minus_one;
  return IsOne(denom) && coords == -vec_m::unit_vector(F->d,1);
}

string FieldElement::str() const
{
  // cout << "In FieldElement::str() with var = " << F->var << endl;
  // cout << "coords = " << coords << ", denom = " << denom << endl;
  ostringstream s;
  if (F->isQ())
    {
      s << val;
      return s.str();
    }
  string n = ::str(coords, F->var);
  if (n[0]=='+')
    n.erase(0,1);
  if (denom==1)
    return n;
  s << "(" << n << ")/" << denom;
  return s.str();
}

int FieldElement::operator==(const FieldElement& b) const
{
  if (F->isQ())
    {
      return val==b.val;
    }
  return F==b.F && denom==b.denom && coords==b.coords;
}

int FieldElement::operator!=(const FieldElement& b) const
{
  if (F->isQ())
    {
      return val!=b.val;
    }
  return (F!=b.F) || (denom!=b.denom) || (coords!=b.coords);
}

mat_m FieldElement::matrix() const // ignores denom, not used for Q
{
  if (F->isQ())
    return mat_m::scalar_matrix(1, num(val));
  return lin_comb_mats(coords, F->Cpowers);
}

// the charpoly is a power of the irreducible minpoly NB This monic
// integer polynomial is the min poly of the numerator, it ignores the
// denominator!
ZZX FieldElement::minpoly() const
{
  if (F->isQ())
    {
      ZZX mp;
      SetX(mp);
      SetCoeff(mp, 0, -num(val));
      return mp;
    }
  ZZX cp = charpoly();
  ZZX mp = factor(cp)[0].a;
  return mp;
}

// add b to this
void FieldElement::operator+=(const FieldElement& b)
{
  if (F!=b.F)
    {
      cerr << "Attempt to add elements of different fields!" << endl;
      exit(1);
    }
  if (b.is_zero())
    return;
  if (F->isQ())
    {
      val += b.val;
      return;
    }
  coords = b.denom*coords + denom*b.coords;
  denom *= b.denom;
  cancel();
}

FieldElement FieldElement::operator+(const FieldElement& b) const
{
  if (F!=b.F)
    {
      cerr << "Attempt to add elements of different fields!" << endl;
      exit(1);
    }
  FieldElement a = *this;
  if (!b.is_zero())
    a += b;
  return a;
}

FieldElement FieldElement::operator-() const
{
  if (F->isQ())
    return FieldElement(-val);
  return FieldElement(F, -coords, denom);
}

// subtract b
void FieldElement::operator-=(const FieldElement& b)
{
  if (F!=b.F)
    {
      cerr << "Attempt to subtract elements of different fields!" << endl;
      exit(1);
    }
  if (b.is_zero())
    return;
  if (F->isQ())
    {
      val -= b.val;
      return;
    }
  coords = b.denom*coords - denom*b.coords;
  denom *= b.denom;
  cancel();
}

FieldElement FieldElement::operator-(const FieldElement& b) const
{
  if (F!=b.F)
    {
      cerr << "Attempt to subtract elements of different fields!" << endl;
      exit(1);
    }
  FieldElement a = *this;
  if (!b.is_zero())
    a -= b;
  return a;
}

void FieldElement::operator*=(const FieldElement& b) // multiply by b
{
  if (F!=b.F)
    {
      cerr << "Attempt to multiply elements of different fields!" << endl;
      cerr << "LHS field: "; F->display(); cerr << endl;
      cerr << "RHS field: "; b.F->display(); cerr << endl;
      exit(1);
    }
  if (is_zero()) return;
  if (b.is_zero()) {*this = b; return;}
  if (F->isQ())
    {
      val *= b.val;
      return;
    }
  coords = (matrix()*b.matrix()).col(1);
  denom *= b.denom;
  cancel();
}

FieldElement FieldElement::operator*(const FieldElement& b) const
{
  if (F!=b.F)
    {
      cerr << "Attempt to multiply elements of different fields!" << endl;
      cerr << "LHS field: "; F->display(); cerr << endl;
      cerr << "RHS field: "; b.F->display(); cerr << endl;
      exit(1);
    }
  if (b.is_zero()) return b;
  FieldElement a = *this;
  if (!is_zero())
    a *= b;
  return a;
}

FieldElement FieldElement::inverse() const // raise error if zero
{
  if (is_zero())
    {
      cerr << "Attempt to invert zero!" << endl;
      exit(1);
    }
  if (is_one() || is_minus_one())
    return FieldElement(*this);
  if (F->isQ())
    return FieldElement(recip(val));

  mat_m M = matrix(), Minv;
  ZZ Mdet = ::inverse(M,Minv); // so M*Minv = Mdet*identity
  FieldElement ans = FieldElement(F, denom*Minv.col(1), Mdet);
  assert (operator*(ans).is_one());
  return ans;
}

void FieldElement::operator/=(const FieldElement& b)      // divide by b
{
  if (b.is_zero())
    {
      cerr << "Attempt to divide by zero!" << endl;
      exit(1);
    }
  if (F!=b.F)
    {
      cerr << "Attempt to divide elements of different fields!" << endl;
      exit(1);
    }
  if (is_zero())
    return;
  if (F->isQ())
    {
      val /= b.val;
      return;
    }
  operator*=(b.inverse());
}

FieldElement FieldElement::operator/(const FieldElement& b) const // raise error if b is zero
{
  if (b.is_zero())
    {
      cerr << "Attempt to divide by zero!" << endl;
      exit(1);
    }
  if (F!=b.F)
    {
      cerr << "Attempt to divide elements of different fields!" << endl;
      exit(1);
    }
  FieldElement a = *this;
  if (!is_zero())
    a /= b;
  return a;
}

FieldElement evaluate(const ZZX& f, const FieldElement a)
{
  FieldElement fa(a.F, to_ZZ(0));
  for(int i=deg(f); i>=0; i--)
    {
      fa += coeff(f,i);
      if(i)
        fa *= a;
    }
  return fa;

}

// NB for a in F, either [Q(sqrt(a))=Q(a)] or [Q(sqrt(a)):Q(a)]=2.
// The first function only applies when a has maximal degree:
// return 1 and r s.t. r^2=this, with deg(r)=degree(), else 0
int FieldElement::is_absolute_square(FieldElement& r) const
{
  if (F->isQ())
    return (val.is_square(r.val));
  // field not Q, reduce to integral case if necessary
  if (::is_one(denom))
    return is_absolute_integral_square(r);
  FieldElement x = FieldElement(F, denom*coords);
  int res = x.is_absolute_integral_square(r);
  if (res)
    r /= denom;
  return res;
}

// Same as above if the min poly is known and denom=1
int FieldElement::is_absolute_integral_square(FieldElement& r)  const
{
  if (F->isQ())
    return (val.is_square(r.val));
  assert (::is_one(denom));
  ZZX g, g0, g1, f = minpoly();
  if (::is_square(f, g)) // Tests if f(x^2) = (+/-) g(x)*g(-x)
    {
      parity_split(g, g0, g1); // g(x) = g0(x^2) + x*g1(x^2)
      // Now 0 = g(r) = g0(a)+r*g1(a) and g1(a)!=0
      r = - evaluate(g0,*this)/evaluate(g1,*this);
      if (r*r== *this)
        return 1;
      else
        {
          cout << (*this) << " has scaled min poly " << ::str(f)
               << " whose double has factor " << ::str(g)
               << " with even part g0 = " << ::str(g0)
               << " and odd part g1 = " << ::str(g1) << endl;
          cout << "These evaluate to " << evaluate(g0,*this) << " and " << evaluate(g1,*this)
               << " with negative quotient/denom r = " << r
               << " but r*r = " << r*r << endl;
          return 0;
        }
    }
  return 0;
}

// The second function applies in general:
// return 1 and r s.t. r^2=this, with deg(r)=degree(), else 0
//#define DEBUG_IS_SQUARE
int FieldElement::is_square(FieldElement& r, int ntries) const
{
#ifdef DEBUG_IS_SQUARE
  cout << "Testing whether " << (*this) << " is a square" << endl;
#endif
  if (is_zero() || is_one())
    {
      r = *this;
      return 1;
    }
  if (F->isQ())
    return (val.is_square(r.val));

  // If this has full degree, or if the relative degree is odd, just
  // call is_absolute_square()
  int d = F->d;
  ZZX m = minpoly();
  if ((d/deg(m))%2)
    {
#ifdef DEBUG_IS_SQUARE
      cout << "Is_square() works directly" << endl;
#endif
      int res = is_absolute_square(r);
#ifdef DEBUG_IS_SQUARE
      if (res)
        cout << (*this) << " is square with sqrt(" << (*this) << ") = " << r << endl;
      else
        cout << (*this) << " is not a square" << endl;
#endif
      return res;
    }
#ifdef DEBUG_IS_SQUARE
  cout << "Is_square() multipying by successive squares" << endl;
#endif

  // Otherwise we multiply this by a "random" square to have full degree first
  FieldElement b = F->gen();
  for (int i=0; i<ntries; i++, b+=to_ZZ(1))
    {
      FieldElement abb = (*this)*b*b, rb(F);
      ZZX mb = abb.minpoly();
      if ((d/deg(mb))%2)
        {
#ifdef DEBUG_IS_SQUARE
          cout << "Is_square() succeeds after " << i+1 << " tries, using b = " << b << endl;
#endif
          int res = abb.is_absolute_square(rb);
          if (res)
            r = rb/b;
          return res;
        }
      // else keep trying
    } // end of loop over shifts i
  cout << "is_square() fails on " << (*this) << " after " << ntries << " tries" << endl;
  return 0;
}

//#define DEBUG_SQUARES

unsigned int FieldModSq::get_index(const FieldElement& a, FieldElement& s)
{
  unsigned int i=0;
  for (auto x: elements)
    {
      if ((a*x).is_square(s))
        {
          assert (a*elements[i] == s*s);
          // Now a*x = s^2 so a = x*(s/x)^2
          // but we want a = x*s^2
          s /= x;
          if (F->isQ()&& s.val.num()<0)
            s=-s;
          assert (a == elements[i]*s*s);
          return i;
        }
      i++;
    }

  // We get here if a (mod squares) is not in the current group.  In
  // particular, it is not a square.

  // Increment the rank:
  unsigned int j = 1<<r; // This will be the new index (with a possible offset)
  r++;

  // If the field is Q we adjoin the squarefree part of a as the new generator:
  if (F->isQ())
    {
      bigrational g1, s1;
      sqfdecomp(a.val, g1, s1); // a = g1*s1^2 with g1 squarefree
      s = FieldElement(s1);
#ifdef DEBUG_SQUARES
      cout << "New generator " << g1 << " for Q^*/(Q^*)^2 from a = " << a << endl;
#endif
      gens.push_back(FieldElement(g1));
#ifdef DEBUG_SQUARES
      cout << "rank of k*/(k*)^2 grows to " << r << " after adding generator " << g1 << endl;
      cout << "elements were " << elements << endl;
#endif
      vector<FieldElement> new_elements(elements.size(), FieldElement(F));
      std::transform(elements.begin(), elements.end(), new_elements.begin(),
                     [g1](const FieldElement& x){return FieldElement(squarefree_product(x.val,g1));});
      elements.insert(elements.end(), new_elements.begin(), new_elements.end());
#ifdef DEBUG_SQUARES
  cout << "elements are now " << elements << endl;
#endif
      return j;
    }

  // Test whether a small integer times an existing rep will do

#ifdef DEBUG_SQUARES
  cout << "Trying to simplify a = " << a << " mod squares" << endl;
#endif
  for (auto u : {-1,2,-2,3,-3,5,-5})
    {
#ifdef DEBUG_SQUARES
      cout << "Trying u = " << u << endl;
#endif
      FieldElement U = F->rational(u);
      FieldElement au = a*U;
      i = 0;
      for (auto x: elements)
        {
          if ((au*x).is_square(s))
            {
#ifdef DEBUG_SQUARES
              cout << "Success with u = " << u << " and x = " << x << endl;
#endif
              // now a*u*x = s^2
              s /= (U*x);
              // now a = s^2 * u*x
              gens.push_back(U);
#ifdef DEBUG_SQUARES
              cout << "rank of k*/(k*)^2 grows to " << r << " after adding generator " << U << endl;
              cout << "elements were " << elements << endl;
#endif
              vector<FieldElement> new_elements(elements.size(), FieldElement(F));
              std::transform(elements.begin(), elements.end(), new_elements.begin(),
                             [U](const FieldElement& x){return U*x;});
              elements.insert(elements.end(), new_elements.begin(), new_elements.end());
#ifdef DEBUG_SQUARES
              cout << "elements are now " << elements << endl;
#endif
              return i + j;
            }
          i++;
        } // end of loop over elements x
    }  // end of loop over small u

  // If we reach here, we failed to find a small new generator, so we use a itself
#ifdef DEBUG_SQUARES
  cout << "No small u worked, so we take " << a << " as new generator" << endl;
#endif
  s = F->one();
  gens.push_back(a);
#ifdef DEBUG_SQUARES
  cout << "rank of k*/(k*)^2 grows to " << r << " after adding generator " << a << endl;
  cout << "elements were " << elements << endl;
#endif
  vector<FieldElement> new_elements(elements.size(), FieldElement(F));
      std::transform(elements.begin(), elements.end(), new_elements.begin(),
                     [a](const FieldElement& x){return a*x;});
  elements.insert(elements.end(), new_elements.begin(), new_elements.end());
#ifdef DEBUG_SQUARES
  cout << "elements are now " << elements << endl;
#endif
  return j;
}

string FieldModSq::elt_str(unsigned int i) const
{
  static const string v = "r";
  ostringstream s;
  if (i==0)
    {
      s << "1";
    }
  else
    {
      int first = 1;
      for (unsigned int j=0; j<r; j++)
        {
          if (bit(i,j))
            {
              if (!first) s << "*";
              s << v << (j+1);
              first = 0;
            }
        }
    }
  return s.str();
}

void FieldModSq::display()
{
  cout << "FieldModSq data:\n";
  cout << "Field  is " <<flush;
  F->display();
  cout << "rank = " << r;
  cout << ", gens = " << gens << endl;
  cout << "order = " << elements.size();
  cout << ", elements " << endl;
  // cout << "[ ";
  for (unsigned int i=0; i<elements.size(); i++)
    cout << i << ": " << elt_str(i) << " = " << elements[i] << endl;
  // cout << "]" << endl;
}

int Eigenvalue::operator==(const Eigenvalue& b) const
{
  return (SqCl==b.SqCl) && (root_index==b.root_index) && (xf==b.xf) && (a==b.a);
}

int Eigenvalue::operator!=(const Eigenvalue& b) const
{
  return (SqCl!=b.SqCl) || (root_index!=b.root_index) || (xf!=b.xf) || (a!=b.a);
}

// When i=sqrt(-1) is the first element of SqCl normalise using
// sqrt(-r)*(1+i)=-sqrt(r)*(1-i) and similar.
// sqrt(-r)*(1-i)=+sqrt(r)*(1+i) and similar.
// NB The xf field will only be non-zero when i is there.
void Eigenvalue::normalise()
{
  if (xf && root_index&1) // true when root_index is odd and xf nonzero
    {
      if (xf==1) a = -a;
      xf = -xf;
      root_index -= 1;
    }
}

Eigenvalue Eigenvalue::inverse() const // raise error if zero      // inverse
{
  if (is_zero())
    cerr << "Attempt to invert " << *(this) << endl;
  FieldElement b = a*root_part();
  if (xf) b *= ZZ(2);
  Eigenvalue ans(b.inverse(), SqCl, root_index, -xf);
#ifdef DEBUG_ARITH
  cout << "Inverse of " << (*this) << " is " << ans << endl;
  assert (((*this)*ans).is_one());
#endif
  return ans;
}

Eigenvalue Eigenvalue::times_i() const
{
  Eigenvalue ans = *this;
  if (ans.xf==0)
    {
      if (ans.root_index&1)
        {
          // we already have a factor i, so remove it and negate
          ans.root_index -= 1;
          ans.a *= ZZ(-1);
        }
      else
        {
          // we do not, so add it
          ans.root_index += 1;
        }
#ifdef DEBUG_ARITH
      cout << "times_i(" << (*this) << ") = " << ans << endl;
#endif
      return ans;
    }
  // else (ans.xf==-1) // factor (1-i)
  else if (ans.xf==1)
    {
      if (ans.root_index&1) // factor i
        {
          // we already have a factor i, so remove it and negate
          ans.root_index -= 1;
        }
      else
        {
          // we do not, so add it: change 1+i into 1-i and negate
          ans.xf = -ans.xf;
        }
      ans.a *= ZZ(-1);
#ifdef DEBUG_ARITH
      cout << "times_i(" << (*this) << ") = " << ans << endl;
#endif
      return ans;
    }
  // else (ans.xf==-1) // factor (1-i)
  if (ans.root_index&1) // factor i
    {
      // we already have a factor i, so remove it and negate
      ans.root_index -= 1;
      ans.a *= ZZ(-1);
    }
  else
    {
      // we do not, so add it: change 1-i into 1+i
      ans.xf = -ans.xf;
    }
#ifdef DEBUG_ARITH
      cout << "times_i(" << (*this) << ") = " << ans << endl;
#endif
  return ans;
}

Eigenvalue Eigenvalue::times_minus_i() const
{
  return - times_i();
}

Eigenvalue Eigenvalue::operator*(Eigenvalue b) const
{
  if (is_zero()) return Eigenvalue(*this);
  if (b.is_zero()) return b;
#ifdef DEBUG_ARITH
  cout << "Multiplying " << (*this) << " by " << b << endl;
  cout << "[" << a << "*sqrt(" << root_part() << ")*" << extra_factor() << "]";
  cout << " * ";
  cout << "[" << b.a << "*sqrt(" << b.root_part() << ")*" << b.extra_factor() << "]";
  cout << endl;
#endif

  // Multiply the coefficients:
  FieldElement c = a * b.a;
#ifdef DEBUG_ARITH
  cout << "Product of coefficients = " << c << endl;
#endif

  // Multiply the root parts:
  FieldElement r, s(a.F->one());
  unsigned int j;
  if (root_index==0)
    j = b.root_index;
  else if (b.root_index==0)
    j = root_index;
  else if (b.root_index==root_index)
    {
      j = 0;
      s = root_part();
      c *= s;
    }
  else
    // NB In this case s is only determined up to sign, hence so is c,
    // hence so is the final product
    {
      r = root_part() * b.root_part();
      j = SqCl->get_index(r, s);
      // Now r = s^2 * elt(j)
      c *= s;
      assert (r == s*s*SqCl->elt(j));
    }
#ifdef DEBUG_ARITH
  cout << "Product of root parts = " << s << " * " << "sqrt(" << SqCl->elt(j) << ")" << endl;
#endif

  // Form the product without the last factors:
  Eigenvalue ans(c, SqCl, j); // sets ans.xf to 0

#ifdef DEBUG_ARITH
  cout << "Before setting last factor, ans = " << ans << endl;
#endif

  if (xf==0)
    ans.xf = b.xf;
  else
    {
      if (b.xf==0)
        ans.xf = xf;
      else
        {
          // now both are nonzero
          if (xf!=b.xf) // (1+i)*(1-i)=2
            ans.a *= ZZ(2);
          else
            {
              if (xf==1) // b.xf=1 too, (1+i)^2=2i
                ans = ans * Eigenvalue(a.F->two(), SqCl, 1);
              else // now xf=b.xf=-1, (1-i)^2=-2i
                ans = ans * Eigenvalue(a.F->minus_two(), SqCl, 1);
            }
        }
    }
  if (ans.xf) // else normalise does nothing
    {
#ifdef DEBUG_ARITH
      cout << "Before normalising, product " << ans << endl;
#endif
      ans.normalise();
    }
#ifdef DEBUG_ARITH
  cout << "Returning product " << ans << endl;
#endif
  return ans;
}

Eigenvalue Eigenvalue::operator/(Eigenvalue b) const
{
  if (is_zero()) return Eigenvalue(*this);
#ifdef DEBUG_ARITH
  cout << "Dividing " << (*this) << " by " << b << endl;
  cout << "[" << a << "*sqrt(" << root_part() << ")*" << extra_factor() << "]";
  cout << " / ";
  cout << "[" << b.a << "*sqrt(" << b.root_part() << ")*" << b.extra_factor() << "]";
  cout << endl;
#endif

  // Divide coefficients:
  FieldElement c = a/b.a;
#ifdef DEBUG_ARITH
  cout << "Quotient of coefficients = " << c << endl;
#endif

  // Divide root parts
  FieldElement r, s(a.F->one());
  unsigned int j;
  if (b.root_index==0)
    {
      j = root_index;
#ifdef DEBUG_ARITH
      cout << "Quotient of root parts = " << "sqrt(" << SqCl->elt(j) << ")" << endl;
#endif
    }
  else if (root_index==0)
    {
      // 1/sqrt(r) = (1/r)*sqrt(r)
      c /= b.root_part();
      j = b.root_index;
#ifdef DEBUG_ARITH
      cout << "Quotient of root parts = (1/" << c << ") * " << "sqrt(" << SqCl->elt(j) << ")" << endl;
#endif
    }
  else
    {
      // sqrt(r1)/sqrt(r2) = s*sqrt(r3) where r1/r2 = s^2*r3
      r = root_part() / b.root_part();
      j = SqCl->get_index(r, s);
      c *= s;
#ifdef DEBUG_ARITH
      cout << "Quotient of root parts = " << s << " * " << "sqrt(" << SqCl->elt(j) << ")" << endl;
#endif
      assert (r == s*s*SqCl->elt(j));
    }

  Eigenvalue ans(c, SqCl, j);
#ifdef DEBUG_ARITH
  cout << "Before adjusting last factor, ans = " << ans << endl;
#endif
  if (b.xf==0)
    ans = Eigenvalue(c, SqCl, j, xf);
  else if (xf==b.xf)
    ans = Eigenvalue(c, SqCl, j);
  else if (xf==0)
    ans = Eigenvalue(c/ZZ(2), SqCl, j, -b.xf);
  else if (xf==1)
    ans = Eigenvalue(c, SqCl, j) * Eigenvalue(a.F->one(), SqCl, 1);
  else
    ans = Eigenvalue(-c, SqCl, j) * Eigenvalue(a.F->one(), SqCl, 1);

  if (ans.xf) // else normalise does nothing
    {
#ifdef DEBUG_ARITH
      cout << "Before normalising, quotient " << ans << endl;
#endif
      ans.normalise();
    }
#ifdef DEBUG_ARITH
  Eigenvalue check = ans*b;
  if (!(check == (*this)))
    {
      cerr << "**************\n"
           << "Quotient ("<<(*this)<<")/("<<b<<") returns " << ans
           << " but "<<b<<"*"<<ans<<" = "<< check
           << " -- wrong!"
           <<endl;
      exit(1);
    }
  check = b.inverse()*(*this);
  if (check != ans)
    {
      cerr << "**************\n"
           << "Quotient ("<<(*this)<<")/("<<b<<") returns " << ans
           << " but "<<(b.inverse())<<"*("<<(*this)<<") = "<< check
           << " -- wrong!"
           <<endl;
      exit(1);
    }
  cout << "Returning quotient " << ans << endl;
#endif
  return ans;
}

string Eigenvalue::str() const
{
  if (is_zero())
    return "0";
  if (is_one())
    return "1";
  if (is_minus_one())
    return "-1";
  if (root_index==0 && xf==0)
    return a.str();

  ostringstream s;
  int QQ = a.F->isQ();

  // output the first factor
  if (a.is_one()) {;}
  else if (a.is_minus_one()) {s<<"-";}
  else if(QQ) {s<<a<<"*";}
  else {s<<"("<<a<<")*";}

  // output the second (sqrt) factor if nontrivial
  if (root_index) // then a involves sqrts
    {
      FieldElement r = SqCl->elt(root_index);
      if (r.is_minus_one())
        s << "i";
      else
        {
          s << "sqrt(";
          if (QQ)
            s << r;
          else
            s << SqCl->elt_str(root_index);
          s  << ")";
        }
      if (xf!=0)
        s << "*";
    }
  // output extra factor (1+i or 1-i) if present
  if (xf!=0)
    s << extra_factor();
  return s.str();
}
