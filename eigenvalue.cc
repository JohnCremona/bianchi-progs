// File EIGENVALUE.CC: classes for working with Hecke eigenvalues
/////////////////////////////////////////////////////////////////

#include "eigenvalue.h"

string FieldModSq::str(int raw) const
{
  ostringstream s;
  if (raw)
    {
      s << r;
      for (auto g:gens)
        s << " " << g.str(1);
    }
  else
    {
      s << "Base field " << F->str()
        << ", rank = " << r
        << ", gens = " << gens
        << ", order = " << elements.size()
        << ", elements: " << elements;
      // for (unsigned int i=0; i<elements.size(); i++)
      //   s << i << ": " << elt_str(i) << " = " << elements[i] << "\n";
    }
  return s.str();
}

// from r gens make the list of 2^r elements
void FieldModSq::make_elements()
{
  elements = {F->one()};
  for (auto g: gens)
    {
      vector<FieldElement> new_elements(elements.size(), FieldElement(F));
      std::transform(elements.begin(), elements.end(), new_elements.begin(),
                     [g](const FieldElement& x){return g*x;});
      elements.insert(elements.end(), new_elements.begin(), new_elements.end());
    }
}

// Input from a raw string format; x's field must be already set
istream& operator>>(istream& s, FieldModSq& x)
{
  // cout << "Reading FieldModSq data, base field is " << *(x.F) << endl;
  s >> x.r;
  // cout << "rank = " << x.r << endl;
  x.gens.resize(x.r, x.F->zero());
  // s >> x.gens; // does not work
  for (auto& g: x.gens)
    s >> g;
  // cout << "gens = " << x.gens << endl;
  x.real_flag = !(x.r>0 && x.gens[0]==x.F->minus_one());
  // cout << "real flag = " << x.real_flag << endl;
  x.make_elements();
  // cout << "elements are " << x.elements <<endl;
  return s;
}

//#define DEBUG_SQUARES

// Compute the index of a nonzero element.

// If a belongs to the current group return i and set s, where a =
// elements[i]*s^2.

// If a does not belong to the subgroup (mod squares):
//   if update (default):
//      append a to gens, increment r, set s=1 return the new r;
//   else:
//      do not change the group, return -1.
unsigned int FieldModSq::get_index(const FieldElement& a, FieldElement& s, int update)
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
          if (F->isQ()&& s.get_val().num()<0)
            s=-s;
          assert (a == elements[i]*s*s);
          return i;
        }
      i++;
    }

  // We get here if a (mod squares) is not in the current group.  In
  // particular, it is not a square.

  if (!update)
    return -1;

  // Now we update the group.  First increment the rank:
  unsigned int j = 1<<r; // This will be the new index (with a possible offset)
  r++;

  // unset the real flag if a=-1
  if (a.is_minus_one())
    real_flag = 0;

  // If the field is Q we adjoin the squarefree part of a as the new generator:
  if (F->isQ())
    {
      bigrational g1, s1, ar(a.get_val());
      int flip = (!real_flag && (!a.is_minus_one()) && (ar.num()<0));
      if (flip)
        {
#ifdef DEBUG_SQUARES
          cout << "Replacing " << ar << " by " << -ar << endl;
#endif
          ar = -ar;
        }
      sqfdecomp(ar, g1, s1); // a = g1*s1^2 with g1 squarefree
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
                     [g1](const FieldElement& x){return FieldElement(squarefree_product(x.get_val(),g1));});
      elements.insert(elements.end(), new_elements.begin(), new_elements.end());
#ifdef DEBUG_SQUARES
  cout << "elements are now " << elements << endl;
#endif
      return j+flip;
    }

  // Test whether a small integer times an existing rep will do

#ifdef DEBUG_SQUARES
  cout << "Trying to simplify a = " << a << " mod squares" << endl;
#endif
  for (auto u : {-1,2,-2,3,-3,5,-5,6,-6,7,-7})
    {
      if (is_complex() && u<0)
        continue;
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
  // (or -a in the complex case if a is rational and negative)
#ifdef DEBUG_SQUARES
  cout << "No small u worked, so we take " << a << " as new generator" << endl;
#endif
  s = F->one();
  FieldElement b(a);
  bigrational ra;
  int flip = (is_complex() && (!b.is_minus_one()) && b.is_rational(ra) && (ra.num()<0));
  if (flip) // adjoin -a not a if a is negative rational
    {
      b = -a;
#ifdef DEBUG_SQUARES
      cout << "Replacing " << a << " by " << b << endl;
#endif
    }
  gens.push_back(b);
#ifdef DEBUG_SQUARES
  cout << "rank of k*/(k*)^2 grows to " << r << " after adding generator " << b << endl;
  cout << "elements were " << elements << endl;
#endif
  vector<FieldElement> new_elements(elements.size(), FieldElement(F));
      std::transform(elements.begin(), elements.end(), new_elements.begin(),
                     [b](const FieldElement& x){return b*x;});
  elements.insert(elements.end(), new_elements.begin(), new_elements.end());
#ifdef DEBUG_SQUARES
  cout << "elements are now " << elements << endl;
#endif
  return j+flip;
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


bigrational Eigenvalue::norm() const
{
  int d1 = SqCl->order(); // a power of 2
  bigrational anorm(a.norm());
  if (d1==1)
    return anorm;
  anorm = bigrational(pow(anorm.num(), d1), pow(anorm.den(), d1));
  int half_d1 = d1/2; // only used when d1>1 so is even
  if (root_index>1)   // then d1>1 so is even
    {
      bigrational rnorm((-root_part()).norm());
      anorm *= bigrational(pow(rnorm.num(), half_d1), pow(rnorm.den(), half_d1));
    }
  if (xf)
    {
      anorm *= bigrational(pow(2,half_d1 * SqCl->field()->degree()));
    }
  return anorm;
}

bigrational Eigenvalue::trace() const
{
  bigrational atrace(a.trace());
  if (SqCl->rank()==0)
    return atrace;
  atrace *= ZZ(SqCl->order());
  bigrational zero;
  if (xf==0)
    {
      return (root_index==0? atrace: zero);
    }
  // now xf =+-1 so there's a factor of 1+-i
  if (root_index>1)
    return zero;
  if ((root_index==1) && (xf==1)) // i*(1+i)=i-1, i*(1-i)=i+1
    atrace = -atrace;
  return atrace;
}

// integer multiple of i, assuming not real
Eigenvalue eye(FieldModSq* S, const ZZ& n)
{
  assert (S->is_complex());
  return Eigenvalue(FieldElement(S->field(), n), S, 1, 0);
}

// Return an embedding into an absolute field (optionally
// polredabs'ed) together with a list of images of the gens.  If the
// rank is 0 return the identity.
//#define DEBUG_ABS_FIELD
FieldIso FieldModSq::absolute_field_embedding(vector<FieldElement>& im_gens, string newvar, int reduce) const
{
#ifdef DEBUG_ABS_FIELD
  cout << " - in absolute_field_embedding() for ";
  display();
  cout << endl;
  cout << " : base field is " << *F << endl;
#endif
  FieldIso emb(F);                    // starting with the identity,

  // case of trivial extension: do nothing (ignore newvar and reduce parameters)
  if (r==0)
    {
#ifdef DEBUG_ABS_FIELD
      cout << " : returning trivial embedding (identity)" << endl;
#endif
      return emb;
    }

  // Field* Fext = (Field*)emb.codom();  // emb maps F to Fext
  Field* Fext = (Field*)(emb.codom());  // emb maps F to Fext
  im_gens.clear();
  int i = 0;
  FieldElement x, sqrt_x;
  for (auto g: gens)
    {
      i++;
#ifdef DEBUG_ABS_FIELD
      cout << i << ": adjoining sqrt(" << g << ")" << flush;
#endif
      x = emb(g); // = r in current field Fext
#ifdef DEBUG_ABS_FIELD
      cout << " = sqrt(" << x << ")" << endl;
#endif
      // create the next iso in the chain
      newvar = F->get_var() + std::to_string(i);
      FieldIso emb1(Fext->sqrt_embedding(x, newvar, sqrt_x, 0)); // no reduction now
#ifdef DEBUG_ABS_FIELD
      cout << " : next simple embedding is \n" << emb1 << endl;
#endif
      // update the field extension
      Fext = (Field*)emb1.codom();
#ifdef DEBUG_ABS_FIELD
      cout << " : next field extension is " << *Fext << endl;
#endif
      // update the embedding of F
      emb.postcompose(emb1);
#ifdef DEBUG_ABS_FIELD
      cout << " : next embedding is " << emb << endl;
#endif
      // map existing im_gens into new Fext
      im_gens = emb1(im_gens);
      //std::for_each(im_gens.begin(), im_gens.end(), [emb1](FieldElement& a){a = emb1(a);});
      // append the new sqrt
      im_gens.push_back(sqrt_x);
#ifdef DEBUG_ABS_FIELD
      cout << " : im_gens is now " << im_gens << endl;
      cout << " in fields\n";
      for (auto z: im_gens) cout << *z.field() << endl;
#endif
    }
  // Final reduction (if requested) and seeting of variable name provided
  if (reduce)
    {
      FieldIso emb1(Fext->reduction_isomorphism(newvar));
#ifdef DEBUG_ABS_FIELD
      cout << " : reduction iso is " << emb1 << endl;
#endif
      // map existing im_gens into new Fext
      im_gens = emb1(im_gens);
      //std::for_each(im_gens.begin(), im_gens.end(), [emb1](FieldElement& x){x = emb1(x);});
      // update the embedding of F
      emb.postcompose(emb1);
#ifdef DEBUG_ABS_FIELD
      cout << " : final embedding is " << emb << endl;
#endif
    }
  else
    {
      Fext->set_var(newvar);
#ifdef DEBUG_ABS_FIELD
      cout << " : final embedding is " << emb << endl;
#endif
    }
  return emb;
}

//#define DEBUG_CONJ
Eigenvalue Eigenvalue::conj() const
{
  Eigenvalue ans = *this;
#ifdef DEBUG_CONJ
  cout << "** Conjugating " << ans << endl;
#endif
  if (SqCl->is_real())
    {
      // cout << "** field is real, returning " << ans << endl;
      return ans;
    }
  // We assume that the first gen mod squares is -1, and that when xf
  // is nonzero, root_index is even

#ifdef DEBUG_CONJ
  if (xf==0)
    {
      if (root_index&1)
        cout << "** xf=0 and index is odd, returning " << -ans << endl;
      else
        cout << "** xf=0 and index is even, returning " << ans << endl;
    }
  else
    {
      assert(root_index%2==0);
      cout << "** xf!=0 and index is even, returning " << Eigenvalue(a, SqCl, root_index, -xf) << endl;
    }
#endif

  return (xf==0?
          (root_index&1? -ans : ans) // negate a iff root_index is odd
          :
          Eigenvalue(a, SqCl, root_index, -xf)// flip the sign of xf
          );
}

Eigenvalue Eigenvalue::operator*(const Eigenvalue& b) const
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
  FieldElement r, s(a.field()->one());
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
                ans = ans * Eigenvalue(a.field()->two(), SqCl, 1);
              else // now xf=b.xf=-1, (1-i)^2=-2i
                ans = ans * Eigenvalue(a.field()->minus_two(), SqCl, 1);
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

Eigenvalue Eigenvalue::operator/(const Eigenvalue& b) const
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
  FieldElement r, s(a.field()->one());
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
    ans = Eigenvalue(c, SqCl, j) * Eigenvalue(a.field()->one(), SqCl, 1);
  else
    ans = Eigenvalue(-c, SqCl, j) * Eigenvalue(a.field()->one(), SqCl, 1);

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

// Input from a raw string format; x's SqCl (and hence its field) must
// be already set
istream& operator>>(istream& s, Eigenvalue& x)
{
  // cout << "Reading an eigenvalue, F = " << *(x.a.field()) << ", Fmodsq = " << *(x.SqCl) << endl;
  s >> x.a;
  // cout << "Base value = " << x.a << endl;
  if (x.SqCl->rank())
    {
      s >> x.root_index;
      // cout << "root_index = " << x.root_index << endl;
      if (x.SqCl->is_complex())
        {
          s >> x.xf;
          // cout << "extra factor = " << x.xf << endl;
        }
    }
  // cout << "Eigenvalue read: " << x << endl;
  return s;
}

string Eigenvalue::str(int raw) const
{
  ostringstream s;
  if (raw)
    {
      s << a.str(1);
      if (SqCl->rank())
        {
          s << " " << root_index;
          if (SqCl->is_complex())
            s << " " << xf;
        }
      return s.str();
    }

  if (is_zero())
    return "0";
  if (is_one())
    return "1";
  if (is_minus_one())
    return "-1";
  if (root_index==0 && xf==0)
    return a.str();

  int QQ = a.field()->isQ();

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

// embed an Eigenvalue into the absolute field Fabs, given an
// embedding of F into Fabs and images of the FieldModSq gens in Fabs
//#define DEBUG_EMBED_EIGS
FieldElement embed_eigenvalue(const Eigenvalue& ap, const FieldIso& emb, const vector<FieldElement>& im_gens)
{
#ifdef DEBUG_EMBED_EIGS
  cout << "Embedding Eigenvalue " << ap << endl;
  cout << " via embedding " << emb << endl;
#endif
  FieldElement a(emb(ap.base()));
#ifdef DEBUG_EMBED_EIGS
  cout << "Base = " << ap.base() << " --> " << a << endl;
#endif
  if (ap.is_zero())
    return a;

  unsigned int s = ap.parent()->order();
  unsigned int apri = ap.index();
  int xf = ap.xfac();
#ifdef DEBUG_EMBED_EIGS
  cout << "index = " << apri << ", root part = " << ap.root_part() << endl;
  cout << "extra factor code = " << xf << endl;
#endif
  if ((apri==0)&&(xf==0)) // quick return when field extension is trivial or ap has no extra factors
    {
#ifdef DEBUG_EMBED_EIGS
      cout << "No extra factors, returning " << a << endl;
#endif
      return a;
    }

  // Multiply by those im_gens for which the corresponding bit of
  // ap.root_index is 1
  for (unsigned int i=0; i<s; i++)
    if (bit(apri, i))
      {
#ifdef DEBUG_EMBED_EIGS
        cout << "About to multiply " << a << " in field " << *(a.field()) << "(" << a.field() << ")"
             << " by " << im_gens[i] << " in field " << *(im_gens[i].field()) << "(" << im_gens[i].field()
             << ")" << endl;
#endif
        a *= im_gens[i];
#ifdef DEBUG_EMBED_EIGS
        cout << "Multiplying by " << im_gens[i] << " gives " << a << endl;
#endif
      }
#ifdef DEBUG_EMBED_EIGS
  cout << "Multiplying by all sqrt factors gives " << a << endl;
#endif

  // multiply by 1+i or 1-1 if required
  if (ap.xfac())
    {
      FieldElement one = a.field()->rational(1);
      if (ap.xfac()==1)
        {
#ifdef DEBUG_EMBED_EIGS
          cout << "Multiplying by (1 + " << im_gens[0] << ") gives " << a << endl;
#endif
          a *= (one+im_gens[0]);
        }
      if (ap.xfac()==-1)
        {
          a *= (one-im_gens[0]);
#ifdef DEBUG_EMBED_EIGS
          cout << "Multiplying by (1 - " << im_gens[0] << ") gives " << a << endl;
#endif
        }
    }
#ifdef DEBUG_EMBED_EIGS
  cout << "Final embedded value is " << a << endl;
#endif
  return a;
}
