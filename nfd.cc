// FILE nfd.cc: implementation of class Newspace for newforms of any dimension
//////////////////////////////////////////////////////////////////////////

#include "nfd.h"
#include "field.h"
#include "oldforms.h"
#include "qidloop.h"

newform_comparison newform_cmp;

#define genus_class_triviality_bound 10 // if aP=0 for this number of
                                        // good primes in a genus
                                        // class then we assume that a
                                        // newform is self-twist and
                                        // all aP for P in this class
                                        // are 0.

// When this is set, use the old spaces info read in from file, with
// the correct multiplicities even when there are self-twist newforms
// at lower levels, to compute the new char poly of any linear
// combination of principal Hecke operators.

// Otherwise, use recursion to find the new polys. This does not
// always get the multiplicities right when there are self-twist
// newforms at lower levels (which can only happen when the class
// number is even).

#define USE_OLD_SPACES

// For class group C4 only, return v where v[i] is the index of ideal
// class c^i for one generator class c
vector<int> C4classes()
{
  if (is_C4())
    {
      // an ideal in class c which generates the class group:
      Qideal Pc = Quad::class_group_2_cotorsion[1];
      // The ideal classes are labelled 0,1,2,3 with 0 for principal but the order varies:
      return {0, Pc.ideal_class(), (Pc*Pc).ideal_class(), Pc.conj().ideal_class()};
    }
  else
    {
      cerr << "Cannot call C4classes() unless the class group is C4" << endl;
      return {};
    }
}

Newform::Newform(Newspace* x, int ind, const ZZX& f, int verbose)
  :nf(x), index(ind), lab(codeletter(ind-1))
{
  d = deg(f);
  string fstring = str(f);
  if (verbose)
    cout << "Constructing Newform from factor f = "<< fstring <<" of degree "<<d<<endl;

  // Compute f(T); since T is scaled by dH and f(X) is not, we
  // evaluate dH^d*f(X/dH) at T; that is, we scale the coefficient of
  // X^i by dH^(d-i):
  if (verbose)
    cout << "Finding kernel of f(T)..."<<endl;
  S = kernel(to_mat(evaluate(scale_poly_up(f, to_ZZ(nf->dH)), nf->T_mat)));
  if(dim(S)!=d)
    {
      cout<<"Problem: eigenspace has wrong dimension "<<dim(S)<<", not "<<d<<endl;
      exit(1);
    }
  denom_abs=(nf->dH)*denom(S);
  if (verbose)
    {
      cout<<"Finished constructing subspace S of dimension "<<d
          <<", absolute denom = "<<denom_abs<<endl;
      cout<<"Computing A, the restriction of T..." <<flush;
    }
  mat A = transpose(restrict_mat(to_mat(nf->T_mat),S));
  if(verbose)
    cout<<"done. Checking its char poly..."<<endl;

  // Check that (scaled) charpoly(A) = f
  ZZX cpA = scaled_charpoly(mat_to_mat_ZZ(A), to_ZZ(denom_abs), to_ZZ(nf->H1->hmod));
  if (cpA != f)
    {
      cout<<endl;
      cout<<"Error: f(X) =            "<<fstring<<endl;
      cout<<"but scaled_charpoly(A) = "<<str(cpA)<<endl;
      exit(1);
    }

  if(verbose)
    {
      cout<<"done."<<endl;
      cout<<"A (the matrix of T restricted to S) = ";
      output_flat_matrix(A);
      if(denom_abs>1)
        cout<<" / " << denom_abs;
      cout<<endl;
      cout<<"f(X) is the min poly of A"<<endl;
    }

  // compute projcoord, precomputed projections the basis of S
  projcoord = nf->H1->FR.get_coord() * basis(S);

  // Compute Hecke field basis (expressing the basis on which we will
  // express eigenvalues w.r.t. the power basis on the roots of f)

  if (d==1)
    {
      F = F0 = FieldQQ;
      Fiso = FieldIso(F);
    }
  else
    {
      string var = codeletter(index-1);
      F0 = new Field(to_mat_m(A), to_ZZ(denom_abs), var, verbose>1);
      Fiso = F0->reduction_isomorphism(var);
      F = (Field*)Fiso.codom();
      if (Fiso.is_nontrivial() && verbose)
        {
          cout << "[replacing original Hecke field with polynomial " << ::str(F0->poly())
               << " with polredabs reduced field with polynomial " << ::str(F->poly()) << "]" << endl;
        }
    }
  if (verbose)
    {
      cout <<"Principal Hecke field data:" << endl;
      cout << *F << endl;
    }

  // Compute character values
  int n2r = Quad::class_group_2_rank;
  if (n2r)
    {
      epsvec.resize(nf->eps_ops.size());
      std::transform(nf->eps_ops.begin(), nf->eps_ops.end(), epsvec.begin(),
                     [this](const matop& op){return (eps(op)>0?+1:-1);});
      if (verbose)
        cout<<"genus char values: "<<epsvec<<endl;
      triv_char = std::all_of(epsvec.cbegin(), epsvec.cend(), [](int i) { return i == +1; });
    }
  else
    triv_char=1;

  bc = (nf->N.is_Galois_stable()? -1 : 0);
  bct = -1;  // means unknown
  cm = 1;    // means unknown
  sfe = 0;   // means unknown

  // Initialise book-keeping data for eigenvalue computation
  Fmodsq = new FieldModSq(F);
  possible_self_twists = nf->possible_self_twists;
  int self_twist_possible = !possible_self_twists.empty();
  if (verbose && n2r)
    {
      if (self_twist_possible)
        cout << "Possible self-twist discriminants: " << possible_self_twists << endl;
      else
        cout << "Self-twist not possible at this level" << endl;
    }
  self_twist_flag = (self_twist_possible ? -1 : 0); // -1 means not yet decided
  genus_class_trivial_counter.resize(1<<n2r, 0);
  genus_classes_nonzero.resize(1,0);
  genus_classes_no_ext.resize(1,0);
  genus_class_ideals.resize(1,Qideal(ONE));
  genus_class_no_ext_ideals.resize(1,Qideal(ONE));
  genus_class_aP.resize(1,Eigenvalue(F->one(), Fmodsq, 0));
}

// Constructor which will read from file
Newform::Newform(Newspace* x, int i, int verbose)
  :nf(x), index(i), lab(codeletter(i-1))
{
  if (verbose)
    cout << "Constructing Newform from file data "<< endl;
  if (!input_from_file(verbose))
    cerr << "Unable to read Newform " << lab << endl;
}

string Newspace::short_label()
{
  return level_label;
}

string Newform::short_label() const
{
  return nf->short_label() + string("-") + lab;
}

string Newspace::long_label()
{
  return field_label() + string("-") + level_label;
}

string Newform::long_label() const
{
  return nf->long_label() + string("-") + lab;
}

string Newspace::conj_label() const
{
  Qideal Nbar = N.conj();
  return label(Nbar);
}

string Newform::conj_label() const
{
  return nf->conj_label() + string("-") + lab;
}

string Newspace::long_conj_label() const
{
  return field_label() + string("-") + conj_label();
}

string Newform::long_conj_label() const
{
  return nf->long_conj_label() + string("-") + lab;
}

// eigenvalue of a general principal operator:
FieldElement Newform::eig(const matop& T)
{
  nf->H1->projcoord = projcoord;
  //      cout << "Matrix of "<<T.name()<<" is\n" << nf->H1->calcop_restricted(T, S, 0, 0) << endl;
  vec_m apv = to_vec_m(nf->H1->applyop(T, nf->H1->freemods[pivots(S)[1] -1], 1)); // 1: proj to S
  //      cout << "ap vector = " << apv <<endl;
  static const ZZ one(1);
  if (F0->isQ())
    {
      FieldElement ap(bigrational(apv[1], to_ZZ(denom_abs)));
      //            cout << "ap = " << ap << endl;
      return ap;
    }
  else
    {
      FieldElement ap(F0, apv, one, 1); // raw=1
      if (!Fiso.is_identity())
        {
          // cout << "ap0 = " << ap << " --> ";
          ap = Fiso(ap);
          // cout << "ap = " << ap << endl;
        }
      return ap;
    }
}

// eigenvalue of AutoHeckeOp(P):
FieldElement Newform::ap(const Quadprime& P)
{
  return eig(AutoHeckeOp(P,nf->N));
}

// eigenvalue +-1 of a scalar involution operator
int Newform::eps(const matop& T) // T should be a scalar +- identity
{
  //   nf->H1->projcoord = projcoord;
  //   cout << "Matrix of "<<T.name()<<" is\n" << nf->H1->calcop_restricted(T, S, 0, 0) << endl;
  FieldElement e = eig(T);
  //   cout << "Computed e = " << e << endl;
  //   cout << "(should be +1 or -1)" << endl;
  if (e.is_one())
    return 1;
  if (e.is_minus_one())
    return -1;
  cerr << "eps(" << T.name() << ") returns " << e << ", not +1 or -1" << endl;
  exit(1);
}

// eigenvalue of a (good) prime from aPmap if P is in there;
// otherwise either raise an error (if stored_only=1) or (not yet
// implemented) compute it.
Eigenvalue Newform::eig(const Quadprime& P, int stored_only)
{
  auto it = aPmap.find(P);
  if (it!=aPmap.end())
    return it->second;
  if (stored_only && triv_char)
    {
      cerr << "Eigenvalue for P="<<P<<" not yet computed" << endl;
    }
  return Eigenvalue();
  // add code here to compute one eigenvalue (for trivial char case
  // only?)
}

// coefficient in Fabs of integral ideal M from aMmap or computed (and
// stored in aMmap) using multiplicative relations. Trivial character
// only.
FieldElement Newform::aM(Qideal& M) // not const Qideal& as we factor it
{
  if (!triv_char)
    {
      cout << "General Fourier coefficients only implemented for forms with trivial character!" << endl;
      return Fabs->zero();
    }

  // return cached value if we have it
  auto it = aMmap.find(M);
  if (it!=aMmap.end())
    return it->second;

  // otherwise compute it using multiplicative relations from
  // precomputed prime coefficients
  FieldElement a;
  if (M.norm().is_one())
    a = Fabs->one();
  else
    {
      Factorization F = M.factorization();
      if (F.is_prime_power())
        {
          Quadprime P = F.prime(0);
          int good = !P.divides(nf->N);
          int e = F.exponent(0);
          // aPmap has the aP for bad P as well as good, from compute_AL_eigs()
          FieldElement aP = embed_eigenvalue(aPmap[P], abs_emb, im_gens);
          FieldElement nP = Fabs->rational(I2long(P.norm()));
          switch (e) {
          case 1: // prime
            a = aP;
            break;
          case 2: // a(P^2) = a(P)^2+N(P)
            a = aP*aP;
            if (good)
              a += nP;
            break;
          default: // general case a(P^e) = a(P)*a(P^(e-1)) + N(P)*a(P^(e-2))
            Qideal MoverP = M/P;
            a = aP*aM(MoverP);
            if (good)
              {
                MoverP /= P;
                a += nP*aM(MoverP);
              }
          }
        } // end of prime power case
      else
        {
          Qideal Q = F.prime_power(0);
          Qideal MoverQ = M/Q;
          a = aM(Q) * aM(MoverQ);
        }
    }
  aMmap[M] = a;
  return a;
}

// Assuming aPmap filled, fill aMmap (Fourier coefficients)
void Newform::compute_coefficients(int ntp, int verbose)
{
  if (!triv_char && verbose)
    {
      cout << "General Fourier coefficients only implemented for forms with trivial character!" << endl;
      return;
    }

  if (aPmap.empty())
    compute_eigs(ntp, verbose);

  // Compute a(M) for M of norm up to maxP (the max N(P) for which we
  // have aP).  NB trace_list[0] is the trace of a(1) which should be
  // the absolute Hecke field degree
  Qidealooper Iloop(1, I2long(maxP), 1, 1); // 1: both conjugates; 1: sorted
  while (Iloop.not_finished())
    {
      Qideal I = Iloop.next();
      FieldElement a = aM(I); // computes and caches aI
      trace_list.push_back(a.trace().num()); // traces of coefficients are integral
    }
}

// Principal eigenvalue of AutoHeckeOp(P) for a good prime P, from
// stored aP in aPmap.  Only implemented for trivial character
// (where this is the eigenvalue of P or P^2) or C4 class group.
// 'biglevel' is a multiple of the current level, auxiliary ideals A
// must be coprime to this, not just to the current level
//#define DEBUG_EIGPAUTO
FieldElement Newform::eigPauto(Quadprime& P, const Qideal& biglevel, int verb)
{
#ifdef DEBUG_EIGPAUTO
  verb = 1;
#endif
  int n2r = Quad::class_group_2_rank;
  if (verb)
    cout << "In eigPauto(P) with P = " << P << endl;
  Eigenvalue aP = eig(P);
  if (verb)
    cout << " - stored eig for P is " << aP << endl;
  FieldElement a = aP.base();
  if (n2r==0) // odd class number case is easy
    {
      if (verb)
        cout << " - odd class number, returning aP = " << a << endl;
      return a;
    }
  if (triv_char)
    {
      Quadprime PP=P;
      if (PP.genus_class()==0) // then there is no sqrt factor
        {
          if (verb)
            cout << " - trivial genus class, returning aP = " << a << endl;
        }
      else
        {
          ZZ p(I2long(P.norm())); // stupid conversion from INT to ZZ
          a = (aP*aP).base(); // (a_P)^2
          if (verb)
            cout << " - non-trivial genus class, aP^2 = " << a << endl;
          a = a - p;
          if (verb)
            cout << " - hence a_{P^2} = " << a << endl;
        }
      return a;
    }
  if (is_C4())
    // C4 class group, character chi_1, values on classes 0,1,2,3 are 1,i,-1,-i
    // The ideal classes are labelled 0,1,2,3 with 0 for principal but the order varies:
    {
      vector<int> C4C = C4classes();
#ifdef DEBUG_EIGPAUTO
      cout << "classes 1,c,c^2,c^3 are " << C4C << endl;
#endif
      int ic1 = C4C[1], ic2 = C4C[2], ic3 = C4C[3];
      int cP = P.ideal_class();
      if (cP==0) // principal
        {
          if (verb)
            cout << " - principal, returning aP = " << a << endl;
          return a;
        }
      if (cP==ic2) // Return the eig of T(P)T(A,A)
        {
          Qideal A = P.sqrt_coprime_to(biglevel);
          int cA = A.ideal_class(); // c or c^3
          assert ((cA==ic1) || (cA==ic3));
          // Here [A]^2=[P] implies [A] = c or c^3, so chi(A)=i or -i;
          // so the eig of T(P)T(A,A) is aP*i or -aP*i;
          // and aP=a*i, so aP*i=-a, so we return -a or +a respectively.
          if (cA==ic1)
            a = -a;
          return a;
        }
      // now cP = ic1 or ic3 and we return the eig of T(P^2)T(A,A)
#ifdef DEBUG_EIGPAUTO
      assert ((cP==ic1) || (cP==ic3));
      cout << "[P] = "  << cP << endl;
#endif
      Qideal A = P.equivalent_mod_2_coprime_to(biglevel,1);
      int cA = A.ideal_class(); // c or c^3
#ifdef DEBUG_EIGPAUTO
      assert ((cA==ic1) || (cA==ic3));
      cout << "A = "  << A << endl;
      cout << "[A] = "  << cA << endl;
#endif
      // a_P     = a*sqrt(r)*(1+/-i)
      // (a_P)^2 = +/- 2*a^2*r*i
#ifdef DEBUG_EIGPAUTO
      cout << "(aP)^2 = "  << aP*aP << endl;
#endif
      a = (aP*aP).base(); // =  (a_P)^2/i
#ifdef DEBUG_EIGPAUTO
      cout << "(aP)^2 / i = "  << a << endl;
#endif
      ZZ p(I2long(P.norm())); // stupid conversion from INT to ZZ
      if (cP==ic1)
        a -= p;
      else
        a += p;
#ifdef DEBUG_EIGPAUTO
      cout << "(a_{P^2}) / i = "  << a << endl;
#endif
      // Now a =  a_{P^2} / i
      if (cA==ic1)
        a = -a;
#ifdef DEBUG_EIGPAUTO
      cout << "(a_{P^2}) *chi(A) = "  << a << endl;
#endif
      // Now a = a_{P^2}*chi(A)
      return a;
    }
  cerr << "Newform::eigPauto() only implemented for forms with trivial char or class group C4" << endl;
  return FieldElement(F);
}

// Principal eigenvalue of a linear combination of the above:
FieldElement Newform::eig_lin_comb(const vector<Quadprime>& Plist, const vector<scalar>& coeffs,
                                   const Qideal& biglevel, int verb)
{
  FieldElement eig(F->zero());
  auto Pi = Plist.begin();
  auto ci = coeffs.begin();
  while (Pi!=Plist.end())
    {
      scalar c = *ci++;
      Quadprime P = *Pi++;
      if (c!=0)
        eig += eigPauto(P, biglevel, verb) * to_ZZ(c);
    }
  return eig;
}

// Characteristic polynomial of such a linear combination:
ZZX Newform::char_pol_lin_comb(const vector<Quadprime>& Plist, const vector<scalar>& coeffs,
                               const Qideal& biglevel, int verb)
{
  return eig_lin_comb(Plist, coeffs, biglevel, verb).charpoly();
}

//#define DEBUG_BC

void Newform::check_base_change(void)
{
  if (bct!=-1) // already set
    return;
  Qideal N(nf->N);
  if(!(N.is_Galois_stable()))
    {
      bc = 0;
      // may still be bct
    }
  if (aPmap.empty())
    {
      cerr << "Cannot check for base change: no aP computed!" << endl;
      return;
    }
  bc = 1;       // will set to 0 if not bc
  bct = 1;      // will set to 0 if not bct (NB bc implies bct)
  for (auto PaP : aPmap)
    {
      if (!bct)
        break; // then alo bc==0
      Quadprime P = PaP.first;
      Eigenvalue aP = PaP.second;
#ifdef DEBUG_BC
      cout<<"P = "<<P<<" has aP="<<aP<<endl;
#endif
      if(P.is_Galois_stable()) // this prime is inert or ramified
        continue;
      if(P.divides(N))
        continue;
      Quadprime Pbar = P.conj();
      if(Pbar.divides(N))
        continue;
      auto it = aPmap.find(Pbar);
      if (it==aPmap.end()) // the conjugate aP is not known
        continue;
      Eigenvalue aPbar = it->second;
#ifdef DEBUG_BC
          cout<<"Pbar = "<<Pbar<<" has aPbar = "<<aPbar<<endl;
#endif
          if(aP!=aPbar) // aP, aP differ
            {
#ifdef DEBUG_BC
              cout<<"Different -- not base-change"<<endl;
#endif
              bc = 0;
              if (aP!=-aPbar) // aP, aP differ even up to sign
                {
#ifdef DEBUG_BC
                  cout<<"Different up to sign -- not base-change"<<endl;
#endif
                  bct = 0;
                }
            }
    }
#ifdef DEBUG_BC
  if (bc)
    cout<<"Base-change"<<endl;
  else if (bct)
    cout<<"Not base-change, but base-change up to twist"<<endl;
  else
    cout<<"Not base-change even up to twist"<<endl;
#endif
  return;
}

Newspace::Newspace(homspace* h1, int maxnp, int maxc, int triv_char_only, int verb)
  : verbose(verb), H1(h1)
{
  N = H1->N;
  level_label = label(N);
  if(verbose)
    cout << "In Newspace constructor at level " << level_label << " = " << N << endl;
  dimH = H1->h1dim();
  cdimH = H1->h1cuspdim();
  cdimH_tc = H1->trivial_character_subspace_dimension(1); // cuspidal=1
  dH = H1->h1cdenom();
  hmod = H1->hmod;
  Ndivs = alldivs(N, 1); // 1 for proper, i.e. excude N itself
  badprimes = N.factorization().sorted_primes();
  possible_self_twists = N.possible_unramified_twists();

  if(verbose && dH>1)
    cout<<"H has dimension "<<dimH<<", cuspidal dimension "<<cdimH<<", denominator "<<dH<<endl;

  if (cdimH==0)
    {
      if (verbose)
        cout << "Cuspidal dimension is 0, so newspace is trivial" << endl;
      return;
    }

  if (triv_char_only && cdimH_tc==0)
    {
      if (verbose)
        cout << "Trivial character cuspidal dimension is 0, so trivial character newspace is trivial" << endl;
      return;
    }

  // Make the unramified quadratic character operators
  t2ideals = make_nulist(N);
  eps_ops.resize(t2ideals.size());
  std::transform(t2ideals.begin(), t2ideals.end(), eps_ops.begin(),
                 [this](Qideal& A){return CharOp(A, N);});
  if (verbose&& t2ideals.size())
    {
      cout << "Unramified quadratic character operators: ";
      for (auto T: eps_ops) cout << T.name() << " ";
      cout << endl;
    }
  // Find the splitting operator
  find_T(maxnp, maxc, triv_char_only);
  // Construct the newforms if that succeeded
  if (split_ok)
    {
      if (verbose)
        {
          cout << "Using splitting operator T = " << T_name << endl;
          cout << "Number of irreducible components is " << factors.size() << endl;
        }
      int i=1;
      newforms.reserve(factors.size());
      for (const auto& f: factors)
        {
          // The index i is stored in the newform on construction, but
          // these will change when we sort them
          newforms.push_back(Newform(this, i, f, verbose));
          ++i;
        }
    }
  else
    // abort if not
    cout << "Unable to find a suitable splitting operator!" << endl;

  // Sort the newforms (by character, dimension, polynomial; we have no traces yet)
  sort_newforms();
}

// sort the list of newforms using newform_cmp
void Newspace::sort_newforms()
{
  // Sort the newforms (by character, dimension, polynomial)
  std::sort(newforms.begin(), newforms.end(), newform_cmp);
  int i=1;
  for (Newform& f: newforms)
    f.set_index(i++);
}

// constructor from file
Newspace::Newspace(const Qideal& level, int verb)
{
  input_from_file(level, verb);
}

// Compute the char poly of T on the new cuspidal subspace using the
// oldspaces to obtain the old factors with correct multiplicities.

// If triv_char=0: requires oldspace data for forms with all genus
// characters, so will only work over fields where this is
// implemented.

// If triv_char=1: only requires oldspace data for forms with
// trivial genus characters, so works over all fields.
ZZX Newspace::new_cuspidal_poly(const vector<Quadprime>& Plist, const vector<scalar>& coeffs,
                                const gmatop &T, int triv_char)
{
  // This will use the caches
  if (verbose>1)
    cout << "In new_cuspidal_poly() with T = " <<T.name()
         << ", triv_char = " << triv_char << endl;

  // Get/compute the char poly of T on the cuspidal subspace /
  // cuspidal, trivial character, subspace
  ZZX f_all = get_poly(N, T, 1, triv_char, H1->modulus); // cuspidal=1
  if (verbose>1)
    {
      cout << "f_all = " << str(f_all) << endl;
      display_factors(f_all);
    }
  ZZX f_old;
  set(f_old); // set = 1
  for (auto D: Ndivs) // Ndivs contains all *proper* D|N
    {
      Newspace* NSD = get_Newspace(D, verbose);
      if (NSD->nforms()==0)
        continue;
      Qideal M = N/D;
      if (verbose>1)
        cout << "Divisor D = " << label(D) << " has " << NSD->nforms() << " newforms" <<endl;
      vector<Qideal> divs = alldivs(M);
      int m_default = divs.size();
      if (verbose>1)
        cout << "# divisors of N/D is " << m_default <<endl;
      for (auto form: NSD->newforms)
        {
          if (form.triv_char==0 && triv_char)
            continue;
          ZZX f_D = form.char_pol_lin_comb(Plist, coeffs, N, verbose>1);
          if (verbose>1)
            {
              cout << "T's poly for " << form.label_suffix() << " is f_D = " << str(f_D) << endl;
              display_factors(f_D);
            }
          INT CMD = form.self_twist_discriminant();
          int m = (CMD.is_zero()? m_default : old_multiplicity(CMD,  divs));
          if (verbose>1)
            cout << "multiplicity m = " << m << endl;
          for (int i=0; i<m; i++)
            f_old *= f_D;
        }
    }
  if (verbose>1)
    {
      cout << "The product of all old polys is " << str(f_old) << endl;
      display_factors(f_old);
    }
  //essentially f_new = f_all / f_old; but checking divisibility
  ZZX f_new, rem;
  DivRem(f_new, rem, f_all, f_old);
  if (!IsZero(rem))
    {
      cout << "Problem in Newspace::new_cuspidal_poly() at level "<<level_label<<endl;
      cout << "Operator " << T.name() << " has full poly" << endl;
      display_factors(f_all);
      cout << "and old poly" << endl;
      display_factors(f_old);
      cout << " which does not divide the full poly. " << endl;
      exit(1);
    }

  // cache this new poly
  string NT = NTkey(N,T);
  if (verbose)
    cout<<"Caching new cuspidal poly for " << NT << ": " << str(f_new) << endl;
  if (triv_char)
    tc_new_cuspidal_poly_dict[NT] = f_new;
  else
    new_cuspidal_poly_dict[NT] = f_new;
  return f_new;
}

// Return true iff this combo of ops T has char poly on the new
// cuspidal subspace which is squarefree and coprime to both the old
// cuspidal poly and the full Eisenstein poly. f_new is set to the new
// cuspidal poly. If triv_char=1 then same for the char poly of T on
// the new cuspidal trivial-character subspace, in which case f_new
// must also be coprime to the char poly of T on the new cuspidal
// nontrivial char subspace.
int Newspace::valid_splitting_combo(const vector<Quadprime>& Plist, const vector<scalar>& coeffs,
                                    const gmatop &T, int triv_char, ZZX& f_new)
{
  // f_new is the char poly of T on either the new cuspidal subspace,
  // or the trivial character subspace of that:
  f_new = new_cuspidal_poly(Plist, coeffs, T, triv_char);
  if (!IsSquareFree(f_new))
    {
      if (verbose>1)
        cout << "\n NO: new Hecke polynomial "<<str(f_new)
             << " for " << T.name() << " is not squarefree" << endl;
      return 0;
    }
  // f_full is the char poly of T on the full space (not just the
  // cuspidal subspace, or the new subspace, or (if triv_char==1) just
  // the trivial character subspace:
  ZZX f_full = get_poly(N, T, 0, 0, H1->modulus); // cuspidal=0, triv_char=0

  // Hence the quotient f_other is the char poly of T on the
  // complement of the new cuspidal (or new cuspidal trivial char0
  // subspace.  We want f_new to be coprime to this:
  ZZX f_other = f_full / f_new;
  if (!AreCoprime(f_new, f_other))
    {
      if (verbose>1)
        cout << "\n NO: new Hecke polynomial "<<str(f_new)
             << " for " << T.name()
             <<" is not coprime to old Hecke polynomial * Eisenstein polynomial"<<str(f_other)<<endl
             <<" (full polynomial is "<<str(f_full)<<")"<<endl;
      return 0;
    }
  if (verbose>1)
    cout << "\n YES: new Hecke polynomial for " << T.name()
         << " is squarefree and coprime to old * Eisenstein Hecke polynomial" << endl;
  return 1;
}

// Find a linear combination T of up to maxnp operators (T_{A,A} or
// T_P) with coefficients up to maxc, whose char poly on the new
// cuspidal subspace is squarefree and coprime to its char polys on
// the oldspace and non-cuspidal subspace.  If triv_char_only==1, T
// must have squarefree char poly on the trivial character new
// cuspidal subspace and also be coprime to its char poly on the
// nontrivial char part of the newspace.
//
// This function manages the linear combinations, withe the validity
// testing done by valid_splitting_combo().
//
// Set split_ok=1 if successful else 0.

void Newspace::find_T(int maxnp, int maxc, int triv_char_only)
{
  split_ok = 0;
  vector<Quadprime> Plist = make_goodprimes1(N, maxnp, 1); // only_one_conj=1
  vector<matop> TPlist(maxnp);
  std::transform(Plist.begin(), Plist.end(), TPlist.begin(),
                 [this](Quadprime P){return AutoHeckeOp(P,N);});
  vector<matop> ops = TPlist;
  // ops.insert(ops.end(), TPlist.begin(), TPlist.end());
  vector<vector<int>> lincombs = all_linear_combinations(ops.size(), maxc);
  if (verbose)
    {
      cout << "Trying linear combinations with coefficients up to "<<maxc
           <<" of " << ops.size() << " operators";
      if (verbose>1)
        {
          cout << " (";
          for (auto T: ops)
            cout << " " << T.name();
          cout << ")";
        }
      cout << endl;
    }
  gmatop T_op;
  ZZX f;
  for (auto lc: lincombs)
    {
      vector<scalar> ilc(lc.size());
      std::transform(lc.begin(), lc.end(), ilc.begin(), [](int c){return scalar(c);});
      T_op = gmatop(ops, ilc);
      T_name = T_op.name();
      if (verbose)
        {
          cout << "Trying "<<lc<<": "<<T_name<<"..."<<flush;
        }
#ifdef USE_OLD_SPACES
      split_ok = valid_splitting_combo(Plist, ilc, T_op, triv_char_only, f);
#else
      split_ok = test_splitting_operator(N, T_op, H1->modulus, verbose>1);
#endif
      if (split_ok)
        {
          if (verbose)
            cout<<"OK!"<<endl;
          break;
        }
      if (verbose)
        cout << " no good, continuing..." << endl;
    }

  if (split_ok)
    {
      if (verbose)
        cout << " OK: using operator " << T_name << " to split off newforms" << endl;
      T_mat = heckeop(T_op, 0, 1); // not cuspidal,  dual
      if (verbose)
        cout << " Getting new cuspidal poly for " << T_name << " from cache" << endl;
      f = get_new_poly(N, T_op, 1, triv_char_only, H1->modulus); // cuspidal=1, triv_char=0 (cached)
      if (verbose)
        cout << "  New cuspidal poly for " << T_name << " from cache is " << str(f) << endl;
    }
  else
    return;

  factors.clear();
  if (verbose)
    {
      cout << " New cuspidal Hecke polynomial for operator " << T_name
           <<" is "<<str(f)<<endl;
    }
  NTL::vec_ZZX NTL_factors= SFFactor(f);
  ::sort(NTL_factors.begin(), NTL_factors.end(), poly_cmp);
  int nfactors = NTL_factors.length();
  if (verbose)
    cout<<"Irreducible factors:"<<endl;
  for(int i=0; i<nfactors; i++)
    {
      ZZX fi = NTL_factors[i];
      if (verbose)
        cout<<(i+1)<<":\t"<<str(fi)<<"\t(degree "<<deg(fi)<<")"<<endl;
      factors.push_back(fi);
    }
  return;
}

Eigenvalue Newform::compute_one_principal_eig(int i, const matop& T, int store, int verbose)
{
  string Tname = T.name();
  if (verbose)
    cout << "Computing eigenvalue of " << Tname << "..." << flush;
  Eigenvalue a(eig(T), Fmodsq);
  if (store)
    eigmap[{i,Tname}] = a;
  if (verbose)
    cout << ":\t" << a << endl;
  return a;
}

// Fill dict eigmap of eigenvalues of first ntp principal operators
void Newform::compute_principal_eigs(int nap, int verbose)
{
  Qideal N(nf->N);
  int nop = 0; // used to index the dict
  int ip = 0;  // counts the primes used
  for ( auto& P : Quadprimes::list)
    {
      ip++;
      nop++;
      if (ip>nap)
        break;
      if (P.divides(N))
        continue;
      compute_one_principal_eig(nop, AutoHeckeOp(P, N), 1, verbose);

      int jp = 0;
      for ( auto& Q : Quadprimes::list)
        {
          if (Q.divides(N) || (P==Q))
            continue;
          jp++;
          if (jp>=ip)
            break;
          // Compute an eigenvalue for P*Q if that has square class
          Qideal PQ  = P*Q;
          if (!(PQ.has_square_class()))
            continue;
          nop++;
          if (PQ.is_principal())
            compute_one_principal_eig(nop, HeckePQOp(P, Q, N), 1, verbose);
          else
            {
              Qideal A = PQ.sqrt_coprime_to(N);
              compute_one_principal_eig(nop, HeckePQChiOp(P, Q, A, N), 1, verbose);
            }
        }
    } // end of prime loop
}

// Fill aPmap, dict of eigenvalues of first ntp good primes
void Newform::compute_eigs(int ntp, int verbose)
{
  if (triv_char)
    {
      compute_eigs_triv_char(ntp, verbose); // a(P) for good P
      compute_AL_eigs(ntp, verbose);        // e(Q) and a(Q) for bad Q
    }
  else
    {
      if (is_C4())
        compute_eigs_C4(ntp, verbose);      // a(P) for good P
      else
        compute_principal_eigs(ntp, verbose);
    }
  if (Quad::class_group_2_rank)
    {
      string var = F->get_var() + string("1");
      abs_emb = Fmodsq->absolute_field_embedding(im_gens, var, 1); // reduce=1
    }
  else // F = Fabs and abs_emb = identity
    {
      abs_emb = FieldIso(F);
    }
  Fabs = (Field*)abs_emb.codom();
  if (verbose && Quad::class_group_2_rank)
    {
      cout << "absolute embedding:\n" << abs_emb << endl
           << "absolute full Hecke field:\n" << *Fabs << endl <<endl;
    }
  check_base_change();
  if (triv_char) // must be done after Fabs is computed
    compute_coefficients(ntp, verbose);   // a(M) and traces for N(M)<=max N(P)
}

// Fill dict aPmap of eigenvalues of first ntp good primes, class group C4 only
void Newform::compute_eigs_C4(int ntp, int verbose)
{
  // If the character (restricted to Cl[2]) is trivial we take chi=1
  // and use the general code:
  if (triv_char)
    {
      compute_eigs_triv_char(ntp, verbose);
      return;
    }

  Qideal N = nf->N;
  FieldElement b(F), b1(F);

  // add -1 to the generators of F^*/(F*)^2 (b1=1 but is ignored)
  int j = Fmodsq->get_index(F->minus_one(),b1); // -1 = elements[1]*b1^2
  assert (j==1);

  // Otherwise the character (restricted to Cl[2]) is non-trivial so
  // is quartic. We pick a generator c of the class group and take
  // chi(c)=+i.  The ideal classes are labelled 0,1,2,3 with 0 for
  // principal but the order varies:
  vector<int> C4C = C4classes();
  int ic1 = C4C[1], ic2 = C4C[2], ic3 = C4C[3];
  Eigenvalue eig_1(FieldElement(F,ZZ(1)),Fmodsq);
  Eigenvalue eig_i(FieldElement(F,ZZ(1)),Fmodsq,1);
  vector<Eigenvalue> chi(4, eig_1), chi_inv(4, eig_1);
  chi[ic1] = chi_inv[ic3] = eig_i;
  chi[ic2] = chi_inv[ic2] = -eig_1;
  chi[ic3] = chi_inv[ic1] = -eig_i;
  // so chi(A) = chi[cA]  and 1/chi(A) = chi_inv[cA] where cA = ideal_class(A)

  // We will store a prime P0 in class c or c^3 with nonzero aP0 and
  // set this flag to 1:
  Quadprime P0;
  Eigenvalue aP0, aP0_inverse;
  int cP0 = -1; // to avoid a compiler warning about not being set
  int P0_set = 0;

  int ip = 0;  // counts the primes used (including bad primes which are skipped here)
  maxP = INT(0);
  for ( auto& P : Quadprimes::list)
    {
      ip++;
      if (ip>ntp)
        break;
      if (P.divides(N))
        continue; // we compute AL eigs separately, later
      maxP = max(maxP, P.norm());
      int cP = P.ideal_class();
      string Pname = prime_label(P);
      if (cP==0) // P principal, compute T(P) directly
        {
          matop T = HeckePOp(P, N);
          aPmap[P] = compute_one_principal_eig(ip, T, 1, verbose);
          continue;
        }
      if (cP==ic2) // P in class c^2, compute T(P)*T(A,A) where [A]=c or c^3
        {
          Qideal A = P.sqrt_coprime_to(N);
          int cA = A.ideal_class(); // c or c^3
          assert ((cA==ic1) || (cA==ic3));
          matop T = HeckePChiOp(P, A, N);
          FieldElement aPchiP = compute_one_principal_eig(ip, T, 1, verbose).base();
          Eigenvalue aP = Eigenvalue(aPchiP, Fmodsq) * chi_inv[cA];
          if (verbose)
            cout << " a("<<P<<") = " << aP << endl;
          aPmap[P] = aP;
          continue;
        }
      // Now P is in class c or c^3
      if (verbose>1)
        cout << "P = " << P << ", chi(P) = " << chi[cP] << endl;

      if (P0_set)
        {
          // If P and P0 are in opposite classes it is simpler: compute a(P*P0), then set a(P)=a(P*P0)/aP0
          if (cP!=cP0) // then P*P0 is principal
            {
              matop T = HeckePQOp(P, P0, N);
              Eigenvalue aPaP0 = compute_one_principal_eig(ip, T, 1, verbose);
              Eigenvalue aP = aPaP0 * aP0_inverse;
              assert (aP*aP0==aPaP0);
              if (verbose)
                {
                  if (verbose>1)
                    cout << "a("<<P<<") = a("<<P<<"*"<<P0<<")/a("<<P0<<") = "<<aPaP0<<"/"<<aP0<<" = "<<aP<<endl;
                  else
                    cout << " a("<<P<<") = " << aP << endl;
                }
              aPmap[P] = aP;
              continue;
            }
          // Now P and P0 have the same class (c or c^3) so we compute
          // T(P*P0)*T(P0,P0) and divide the eigenvalue by aP0 and
          // chi(P0):
          Qideal B = P0*P; // so B is square-free and of square class
          Qideal C = B.sqrt_coprime_to(N); // so B*C^2 = P0*P*C^2 is principal
          matop T = HeckeBChiOp(B, C, N);
          Eigenvalue aPaP0chiC = compute_one_principal_eig(ip, T, 1, verbose);
          Eigenvalue aP = aPaP0chiC*aP0_inverse * chi_inv[C.ideal_class()];
          aPmap[P] = aP;
          continue;
        }
      // Now P is the first prime in one of these classes. We compute
      // T(P)*T(P,P) to find a(P)^2, pick a square root (if nonzero),
      // store the P in P0, aP in aP0, and set the flag.

      FieldElement p(F->rational(I2long(P.norm()))); // stupid conversion from INT to ZZ
      matop T =  HeckeP2ChiOp(P, P, nf->N);
      FieldElement a = compute_one_principal_eig(ip, T, 1, verbose).base();
      if (verbose>1)
        cout << "a(P^2)*chi(P) = " << a << endl;
      // The following uses chi(P)^2=-1
      Eigenvalue aP_2 = Eigenvalue(p-a, Fmodsq) * chi[cP]; // a(P)^2
      if (verbose>1)
        {
          Eigenvalue aP2chiP(a,Fmodsq);
          Eigenvalue a_P2 = aP2chiP * chi_inv[cP]; // a(P^2)
          cout << "a(P^2) = " << a_P2 << endl;
          cout << "a(P)^2 = " << aP_2 << endl;
        }
      b = aP_2.base() * ZZ(2);
      if (b.is_zero())
        {
          if (verbose)
            cout << "a(" <<P<<") = 0"<<endl;
          aPmap[P] = aP_2;
          continue; // without setting the flag
        }

      // Now a(P)   = (1/2) * sqrt(b) * (1+i).
      // We adjoin sqrt(b) to F(Fmodsq) if necessary.

      // Case (1) If b=b1^2 with b1 in F then a(P) = (b1/2)*(1+i).

      // Case (2) If b=-b1^2 with b1 in F then a(P) = (b1/2)*i*(1+i) =
      // (-b1/2)*(1-i).

      // Case (3) In general, b=b1^2*r and a(P) =
      // (b1/2)*sqrt(r)*(1+i), where r may be a new generator.

      j = Fmodsq->get_index(b,b1);      // +b = elements[j]*b1^2
      // All cases are covered by the following: if j is odd then
      // elements[j] involves i=sqrt(-1), so we: replace j by j-1,
      // negate b1, and use factor 1-i instead of 1+i
      b1 = b1/ZZ(2);
      int xf = 1;
      if (j&1) // absorb factor i into (1+i) factor: i*(1+i)=-(1-i)
        {
          b1 = -b1;
          xf = -1;
          j -= 1;
        }
      Eigenvalue aP(b1, Fmodsq, j, xf);
      aPmap[P] = aP;
      if (verbose)
        cout << "a(" <<P<<") = "<<aP<<endl;

      // Store P and aP and set the flag:
      P0 = P;
      aP0 = aP;
      aP0_inverse = aP.inverse();
      cP0 = cP;
      P0_set = 1;
      if (verbose)
        cout << "Stored eigenvalue " << aP0 << " with inverse " << aP0_inverse << " for prime " << P0 << endl;
    } // end of loop over P

  if (P0_set)
    self_twist_flag = 0;
  else
    {
      if (verbose)
        cout << "All primes in the nontrivial genus class have eigenvalue 0.\n"
             << "Flagging this newform as having self-twist.\n";
      self_twist_flag = +1;
      CMD = possible_self_twists[0];
    }

  if (verbose)
    {
      cout << "The full Hecke field is k_f(i";
      if (Fmodsq->rank()>1)
        cout << ", sqrt("<<Fmodsq->gen(1)<<")";
      cout << ")" << endl;
    }
}

// Fill aPmap, dict of eigenvalues of good primes in the first ntp primes
void Newform::compute_eigs_triv_char(int ntp, int verbose)
{
  if (!triv_char)
    return;
  int nap = 0;
  Qideal N(nf->N);
  maxP = INT(0);
  auto pr = Quadprimes::list.begin();
  while((pr!=Quadprimes::list.end()) && (nap<ntp))
    {
      Quadprime P = *pr++;
      nap++;
      if (P.divides(N))
        continue; // compute AL eigs separately, later
      maxP = max(maxP, P.norm());
      Eigenvalue aP;
      long c = P.genus_class(1); // 1 means reduce mod Quad::class_group_2_rank
      if (c==0) // P has trivial genus class
        {
          aP = compute_one_principal_eig(nap, AutoHeckeOp(P, N), 1, verbose);
          aPmap[P] = aP;
          if (verbose)
            cout<<" -- P="<<P<<" has trivial genus class, aP = " << aP << endl;
          if (aP.is_zero())
            genus_class_trivial_counter[c] +=1;
          continue;
        }

      // Now P does not have square class, its genus class is c>0
      // See if we already have an eigenvalue for this genus class

      if (verbose)
        cout<<" -- P="<<P<<" has genus class "<<c<<", genus_classes covered so far: "<<genus_classes_nonzero<<endl;
      auto ci = std::find(genus_classes_nonzero.begin(), genus_classes_nonzero.end(), c);
      if (ci != genus_classes_nonzero.end()) // P is in a known genus class
        {
          int i = ci - genus_classes_nonzero.begin();
          Qideal A = genus_class_ideals[i];
          Qideal B = A*P; // so B is square-free and of square class
          if (verbose)
            cout << " -- computing T("<<P<<") using T("<<label(B)<<") = T(P)*T(A) with A = "
                 <<label(A)<<endl;
          if (B.is_principal())         // compute T(B)
            {
              Eigenvalue APeig = compute_one_principal_eig(nap, HeckeBOp(B, N), 1, verbose);
              if (APeig.is_zero())
                {
                  aP = APeig;
                  if (verbose)
                    cout << "a(P*A) = 0 so a(P) = 0" << endl;
                }
              else
                {
                  aP = APeig / genus_class_aP[i];
                  if (verbose)
                    cout << "P*A has eigenvalue " << APeig
                         << "; dividing by " << genus_class_aP[i]
                         << " gives " << aP << endl;
                }
            }
          else                          // compute T(B)*T(C,C)
            {
              Qideal C = B.sqrt_coprime_to(N); // so A*P*C^2 is principal
              Eigenvalue APCeig = compute_one_principal_eig(nap, HeckeBChiOp(B,C, N), 1, verbose);
              if (APCeig.is_zero())
                {
                  aP = APCeig;
                  if (verbose)
                    cout << "a(P*A)chi(C) = 0 so a(P) = 0" << endl;
                }
              else
                {
                  aP = APCeig / genus_class_aP[i];
                  if (verbose)
                    cout << "P*A has eigenvalue " << APCeig
                         << "; dividing by " << genus_class_aP[i]
                         << " gives " << aP << endl;
                }
            }
          aPmap[P] = aP;
          // See whether P is a better genus class rep than the one we have:
          if ((P.norm()<A.norm()) && (!aP.is_zero()))
            {
              genus_class_ideals[i] = P;
              genus_class_aP[i] = aP;
            }
          continue;
        }

      // Now we have a new genus class, compute a_{P}^2 unless we
      // already have at least genus_class_triviality_bound zeros in
      // this class and the level admits nontrivial self twists.  NB 5
      // is not enough for field 299, level 100.2
      if (verbose)
        cout << " -- P=" <<P<<" has genus class "<<c<<", genus_class_trivial_counter = "<<genus_class_trivial_counter<<endl;
      if ((possible_self_twists.size()>0) && (genus_class_trivial_counter[c] >= genus_class_triviality_bound))
        {
          if (verbose)
            cout << " -- P = " <<P<<": genus class "<<c<<" has "<<genus_class_trivial_counter[c]
                 <<" zero eigenvalues, so assuming self-twist, and taking aP=0"<<endl;
          aPmap[P] = Eigenvalue(F->zero(), Fmodsq, 0);
          continue;
        }

      if (verbose)
        cout << " -- P="<<P<<": computing T(P^2) to get a(P^2) and hence a(P)^2" << endl;
      Qideal P2 = P*P;
      FieldElement aP2; // not an Eigenvalue yet as we'll be adding N(P)
      if (P2.is_principal())  // compute T(P^2)
        {
          aP2 = compute_one_principal_eig(nap, HeckeP2Op(P, N), 1, verbose).base();
        }
      else // T(P^2)*T(A,A) with (A*P)^2 principal
        {
          Qideal A = P.equivalent_mod_2_coprime_to(N, 1);
          aP2 = compute_one_principal_eig(nap, HeckeP2ChiOp(P,A, N), 1, verbose).base();
        }
      // Now aP2 is the eigenvalue of T(P^2)
      if (verbose)
        cout << " -- a(P^2) = " << aP2 << endl;
      ZZ normP(I2long(P.norm())); // stupid conversion from INT to ZZ
      aP2 += normP;
      // Now aP2 is the eigenvalue of T(P)^2
      if (verbose)
        cout << " -- a(P)^2 = " << aP2 << endl;

      // Check if this eigenvalue is 0
      if (aP2.is_zero())
        {
          aPmap[P] = Eigenvalue(aP2, Fmodsq, 0);
          genus_class_trivial_counter[c] +=1;
          if (verbose)
            cout << " -- genus_class_trivial_counter for class "<<c
                 <<" is now "<<genus_class_trivial_counter[c]<<endl;
          continue;
        }

      // Look up the square class of aP2 in Fmodsq to see if we need to extend the field
      unsigned int old_order = Fmodsq->order();
      FieldElement s = F->one();
      unsigned int j = Fmodsq->get_index(aP2, s);
      int Hecke_field_extended = j>=old_order;
      if (verbose)
        {
          cout<<" -- P="<<P<<", a(P)^2 = "<<aP2<<" is in square class #"<<j;
          if (Hecke_field_extended)
            cout<<" (new class, rank mod squares is now " << Fmodsq->rank() << ")";
          cout << endl;
        }
      aP = Eigenvalue(s, Fmodsq, j);
      aPmap[P] = aP;
      if (verbose)
        cout << " -- taking a(P) = " << aP << endl;

      // If this new genus class did not require an extension of the
      // Hecke field, record that in genus_classes_no_ext, doubling
      // the size of that subgroup by appending the binary sum of each
      // element and c

      if (!Hecke_field_extended)
        {
          long oldsize = genus_classes_no_ext.size();
          if (verbose)
            cout << "Doubling the size of genus_classes_no_ext from " << oldsize
                 << " to " << 2*oldsize << endl;
          genus_classes_no_ext.resize(2*oldsize);
          genus_class_no_ext_ideals.resize(2*oldsize);
          for (int i = 0; i<oldsize; i++)
            {
              genus_classes_no_ext[oldsize+i] = genus_classes_no_ext[i]^c;
              genus_class_no_ext_ideals[oldsize+i] = genus_class_no_ext_ideals[i]*P;
            }
        }

      // update genus_classes_nonzero (append binary sum of each and c)
      // update genus_class_ideals (append product of each and P)
      // update genus_class_aP (append product of each and aP)
      long oldsize = genus_classes_nonzero.size();
      if (verbose)
        cout<<" -- doubling number of genus classes covered by P with nonzero aP from "<<oldsize
            <<" to "<<2*oldsize<<" using P="<<P<<endl;
      genus_classes_nonzero.resize(2*oldsize);
      genus_class_ideals.resize(2*oldsize);
      genus_class_aP.resize(2*oldsize);
      genus_class_trivial_counter[c] = 0;
      for (int i = 0; i<oldsize; i++)
        {
          genus_classes_nonzero[oldsize+i] = genus_classes_nonzero[i]^c;
          genus_class_ideals[oldsize+i] = genus_class_ideals[i]*P;
          genus_class_aP[oldsize+i] = genus_class_aP[i]*aP;
        }
      if (genus_classes_nonzero.size() == genus_class_trivial_counter.size())
        {
          if (self_twist_flag!=0) // otherwise we already know it is not self-twist
            {
              self_twist_flag = 0;
              if (verbose)
                cout << "All genus classes now have a nonzero eigenvalue, so this form is *not* self-twist" << endl;
            }
        }
      else
        {
          // we can possibly eliminate some of possible_self_twists now,
          // namely those whose characters chi have chi(P)=-1, since such
          // a newform would have to have aP=0:
          int n_before = possible_self_twists.size();
          // We can now eliminate any self-twist discriminants D
          // such that P.genus_character(D)=-1
          possible_self_twists.erase(std::remove_if(possible_self_twists.begin(),
                                                    possible_self_twists.end(),
                                                    [&P](INT D) { return P.genus_character(D)==-1;}),
                                     possible_self_twists.end());
          int n_after = possible_self_twists.size();
          if ((n_before>n_after) && (verbose))
            cout<<" - after erasing "<<(n_before-n_after)
                <<" possible self-twist discriminants, these remain: "<<possible_self_twists << endl;
        }
    } // end of primes loop

  // Report on possible self-twist, if not all the genus classes have
  // a nonzero eigenvalue:

  int n2r = Quad::class_group_2_rank;
  unsigned int ngenera = 1<<n2r;

  if (self_twist_flag == -1)
    {
      if (verbose)
        cout << "Only " << genus_classes_nonzero.size() << " out of " << ngenera
             << " genus classes have non-zero eigenvalues!\n"
             << "Flagging this newform as having self-twist.\n";
      int nst = possible_self_twists.size();
      if (nst==0)
        cout << "*** problem! all possible unramified self-twist discriminants have been eliminiated!" << endl;
      else if (nst==1)
        {
          CMD = possible_self_twists[0];
          if (verbose)
            cout << "Self-twist discriminant is " << CMD << endl;
        }
      else
        if (verbose)
          cout << "*** Self-twist discriminant is one of " << possible_self_twists << endl;
      self_twist_flag = +1;
      CMD = possible_self_twists[0]; // NB needs code to detect which
                                     // one if more than one, but this
                                     // is OK for C4 class group
    }

  // Report on how many genus classes have eigenvalues in the
  // principal Hecke field

  unsigned int ngenera2 = ngenera;
  if (!CMD.is_zero())
    ngenera2 = ngenera/2;

  if (verbose && n2r)
    {
      cout << "Out of " << ngenera << " genus classes, ";
      if (!CMD.is_zero())
        cout << "of which " << ngenera2 << " have nonzero aP, ";
      cout << genus_classes_no_ext.size()
           << " have eigenvalues in the principal Hecke field." << endl;
      cout << "These classes are " << genus_classes_no_ext
           << ", with ideals       " << genus_class_no_ext_ideals << endl;
      cout << "This means that the number of unramified quadratic twists of this newform is "
           << genus_classes_no_ext.size() << endl;
    }
  if (ngenera2 != genus_classes_no_ext.size() * Fmodsq->order())
    {
      cerr << "ngenera = " << ngenera << endl;
      cerr << "ngenera2 = " << ngenera << endl;
      cerr << "Fmodsq->order() = " << Fmodsq->order() << endl;
      cerr << "genus_classes_no_ext = " << genus_classes_no_ext << endl;
    }
  assert (ngenera2 == genus_classes_no_ext.size() * Fmodsq->order());

  set_unramified_twist_discriminants(); // fill the unramified_twist_discriminants listl
} // end of compute_eigs_triv_char

// Fill dict eQmap *after* aPmap, and set sfe, ONLY in triv_char case
void Newform::compute_AL_eigs(int ntp, int verbose)
{
  // Compute AL-eigs to fill map<Quadprime, Eigenvalue> eQmap and also
  // the entries in aPmap indexed by bad primes.

  if (!triv_char)  // we don't deal with pseudo-eigenvalues
    return;

  if (aPmap.empty())
    compute_eigs(ntp, verbose);

  // Except in the trivial char case or if all Q^e|N are principal,
  // this requires aPmap to be filled with a nonzero eigenvalue for
  // primes P in each class inverse to each [Q^e].

  sfe = -1; // since the sign is *minus* the product of the AL eigs

  if (nf->badprimes.empty())  // nothing to do
    return;

  if (verbose)
    cout << "Computing W(Q) eigenvalues for Q in " << nf->badprimes << endl;

  Eigenvalue zero(F->zero(), Fmodsq);
  Qideal N(nf->N);
  for (auto Q: nf->badprimes)
    {
      if (verbose)
        cout << "Computing W(Q) eigenvalue for Q = " << Q << endl;
      int e = val(Q,N);
      assert (e>0);
      Qideal Qe = Q;
      for (int i=1; i<e; i++)
        Qe*=Q;

      int eQ;

      if (Qe.is_principal()) // direct computation
        {
          eQ = eps(AtkinLehnerQOp(Q,N));
          if (verbose)
            cout << "Direct computation gives eigenvalue " << eQ << endl;
        }
      else
        {
          if (Qe.has_square_class())
            {
              Qideal A = Qe.sqrt_coprime_to(N);
              eQ = eps(AtkinLehnerQChiOp(Q,A,N));
              if (verbose)
                cout << "Indirect computation gives eigenvalue " << eQ << endl;
            }
          else
            {
              Quadprime P;
              Eigenvalue aP = zero;
              eQ = 0; // 0 as place-holder

              // Use W(Q)*T(P) where P*Q is principal and a(P) is nonzero
              for (auto Pi: Quadprimes::list)
                {
                  if (Pi.divides(N))
                    continue;
                  if (!(Pi*Qe).is_principal())
                    continue;
                  auto it = aPmap.find(Pi);
                  if (it==aPmap.end())
                    continue;
                  P = Pi;
                  aP = it->second;   // = aPmap[P];
                  if (!aP.is_zero()) // found a good P with aP nonzero
                    {
                      if (verbose)
                        cout << "Indirect computation using T(P)W(Q) with P="<<P<<", aP = " << aP << endl;
                      Eigenvalue aPeQ(eig(HeckePALQOp(P, Q, N)), Fmodsq);
                      eQ = (aPeQ==aP? 1 : -1);
                      break;
                    }
                } // end of loop over P
              if (eQ==0) // then we never found a good P with aP known and nonzero
                cout << "Unable to compute W("<<Q<<"), need to compute more aP" << endl;
            }
        }
      // Store the AL eigenvalue in eQmap and update sfe:
      eQmap[Q] = eQ;
      sfe *= eQ;

      // aQ is minus the AL eigenvalue if Q||N else 0
      if (e>1)
        aPmap[Q] = zero;
      else
        aPmap[Q] = Eigenvalue(FieldElement(F, ZZ(-eQ)), Fmodsq);

    } // end of loop over bad primes Q
} // end of Newform::compute_AL_eigs()

// output basis for the Principal Hecke field and character of one newform
// If full, also output multiplicative basis for the full Hecke field
// Optionally aP and AL (if trivial char) data too
// Optionally traces too
// Optionally principal eigs too
void Newform::display(int aP, int AL, int principal_eigs, int traces) const
{
  cout << "Newform #" << index << " (" << long_label() << ")" << endl;
  // Information about the character:

  int n2r = Quad::class_group_2_rank;
  if (n2r>0)
    {
      if (n2r==1)
        cout << " - Genus character value: " <<  (epsvec[0]>0? "+1 (trivial)": "-1");
      else
        {
          cout << " - Genus character values: " <<  epsvec;
          if (triv_char)
            cout << " (trivial)";
        }
      cout <<endl;
    }

  // Information about the dimension (degree of principal Hecke field):

  cout << " - Dimension: "<<d<<endl;

  // Information about base-change:

  if (bct==-1)
    {
      cout << " - Newform's base-change status not determined" << endl;
    }
  else
    {
      if (bc==1)
        cout << " - Newform is base-change" << endl;
      else if (bct==1)
        cout << " - Newform is twisted base-change but not base-change" << endl;
      else if (bct==0)
        cout << " - Newform is not base-change or twisted base-change" << endl;
      else
        {
          if (nf->N.is_Galois_stable())
            cout << " - Newform's base-change status not determined" << endl;
          else
            cout << " - Newform is not base-change" << endl;
        }
    }

  // Information about unramified self-twist:

  if ((n2r>0) && (self_twist_flag!=-1))
    {
      if (self_twist_flag==0)
        cout << " - Newform is not unramified self-twist" << endl;
      else
        {
          if (CMD.is_zero())
            cout << " - Newform appears to be unramified self-twist, but discriminant not determined" << endl;
          else
            cout << " - Newform is unramified self-twist by discriminant "<< CMD << endl;
        }
    }

  // Information about sign of functional equation:

  if (sfe!=0)
    {
      cout << " - Sign of functional equation = " << (sfe>0? "+1" : "-1") << endl;
    }


  // Information about Hecke field(s):

  if (n2r>0)
    cout << " - Principal Hecke field k_f = ";
  else
    cout << " - Hecke field k_f = ";
  cout << *F;
  if (n2r>0) // display full Hecke field
    {
      int FisQ = F->isQ();
      int FisCM = Fmodsq->is_complex();
      int r = Fmodsq->rank();
      if (triv_char)
        {
          if (r==0)
            cout << " - Full Hecke field is the same as the Principal Hecke field" << endl;
          else
            {
              // If the base field is Q we output something like
              // "Q(sqrt(2), sqrt(5))", otherwise something like
              // "k_f(sqrt(r1), sqrt(r2)) where r1 = ..., r2 = ...".
              cout << " - Full Hecke field k_F = "
                   << (FisQ? "Q" : "k_f")
                   << "(";
              for (int i=0; i<r; i++)
                {
                  if (i) cout << ", ";
                  cout << "sqrt(";
                  if (FisQ)
                    cout << Fmodsq->gen(i);
                  else
                    cout << "r" << (i+1);
                  cout << ")";
                }
              cout << ")";
              if (!(FisQ))
                {
                  cout << " where ";
                  for (int i=0; i<r; i++)
                    {
                      if (i) cout << ", ";
                      cout << "r" << (i+1) << " = " << Fmodsq->gen(i);
                    }
                }
              cout << endl;
              cout << "                        = " << *Fabs <<endl;
              if (!FisQ)
                {
                  cout << "        [where " << F->gen() << " = " << abs_emb(F->gen());
                  for (int i=0; i<r; i++)
                    {
                      if (FisCM && i==0)
                        cout << ", i = " << im_gens[i];
                      else
                        cout << ", sqrt(r" << (i+1) << ") = " << im_gens[i];
                    }
                  cout << "]" << endl;
                }
            }
        } // end of trivial char case
      else // special case code for C4 class group
        {
          if (is_C4())
            {
              cout << " - Full Hecke field k_F = "
                   << (FisQ? "Q" : "k_f")
                   << "(i";
              if (Fmodsq->rank()>1)
                cout << ", sqrt(r2)), where r2 = " << Fmodsq->gen(1);
              else
                cout << ")";
              cout << endl;
              if (nf->verbose>1)
                cout << *Fmodsq << endl;
              cout << "                        = " << *Fabs <<endl;
              if (!FisQ)
                {
                  cout << "        [where " << F->gen() << " = " << abs_emb(F->gen());
                  for (int i=0; i<r; i++)
                    {
                      if (FisCM && i==0)
                        cout << ", i = " << im_gens[i];
                      else
                        cout << ", sqrt(r" << (i+1) << ") = " << im_gens[i];
                    }
                  cout << "]" << endl;
                }
            }
          else
            {
              cout << " - Full Hecke field not determined except for trivial character or class group C4" << endl;
            }
        } // end of non-trivial char case
    }
  if (triv_char && AL)
    {
      cout << endl;
      display_AL();
    }
  if (aP)
    {
      cout << endl;
      display_aP();
    }
  if (principal_eigs)
    {
      cout << endl;
      display_principal_eigs();
    }
  if (triv_char && traces)
    {
      cout << endl;
      cout << "First " << trace_list.size() << " traces: " << trace_list << endl;
    }
}

// Display aP data (trivial char or C4 fields)
//#define CHECK_TRACES
void Newform::display_aP() const
{
  if (aPmap.empty())
    {
      cout << "No aP known" << endl;
      return;
    }
  cout << "aP for first " << aPmap.size() << " primes:" << endl;
  for (auto x : aPmap)
    {
      Eigenvalue aP = x.second;
      cout << x.first << ":\t" << aP;
#ifdef CHECK_TRACES
      bigrational taP = aP.trace();
      cout << " (trace = " << taP << ")";
#endif
      if (Quad::class_group_2_rank && Fmodsq->rank()>0)
        {
          FieldElement bP = embed_eigenvalue(aP, abs_emb, im_gens);
          cout << " \t= " << bP;
#ifdef CHECK_TRACES
          bigrational tbP = bP.trace();
          cout << " (trace = " << tbP << ")";
          assert (taP==tbP);
#endif
        }
      auto it = eQmap.find(x.first);
      if (it != eQmap.end())
        {
          cout << "\t(AL eigenvalue = " << it->second << ")";
        }
      cout << endl;
    }
}

// Display A-L eigenvalues (trivial char or C4 fields)
void Newform::display_AL() const
{
  if (nf->N.norm()==1)
    {
      cout << "No Atkin-Lehner eigenvalues as level is (1)" << endl;
      return;
    }
  cout << "Atkin-Lehner eigenvalues:" << endl;
  for (auto x: eQmap)
    cout << x.first << ":\t" << x.second << endl;
}

// Display principal eigenvalues
void Newform::display_principal_eigs() const
{
  cout << "Eigenvalues of principal operators:" << endl;
  for (auto x: eigmap)
    cout << x.first.second << ":\t" << x.second << endl;
}

// return +1 for base-change, -1 for twisted bc, 0 for neither, 2 for don't know
int Newform::base_change_code(void) const
{
  if (bc==1) return 1;
  else if (bct==1) return -1;
  else if (bct==0) return 0;
  else return 2;
}

map<Quadprime, Eigenvalue> Newform::TP_eigs(int ntp, int verbose)
{
  if (aPmap.empty())
    compute_eigs(ntp, verbose);
  return aPmap;
}

map<Quadprime, int> Newform::AL_eigs(int ntp, int verbose)
{
  if (eQmap.empty())
    compute_AL_eigs(ntp, verbose);
  return eQmap;
}

map<pair<int,string>, Eigenvalue> Newform::principal_eigs(int nap, int verbose)
{
  if (eigmap.empty())
    compute_principal_eigs(nap, verbose);
  return eigmap;
}

// filename
string Newform::filename(int conj) const
{
  stringstream s;
  s << newspaces_directory() << "/" << (conj? conj_label():short_label());
  return s.str();
}

string Newspace::filename(int conj)
{
  stringstream s;
  s << newspaces_directory() << "/" << (conj? conj_label(): short_label());
  return s.str();
}

// newform file output only implemented for forms with trivial
// character or over C4 fields
void Newform::output_to_file(int conj) const
{
  if (!(is_C4() || triv_char))
    return;
  ofstream out;
  out.open(filename(conj).c_str());

  // Field, level, letter-code:
  out << field_label() << " "
      << (conj? nf->conj_label(): nf->short_label())
      << " " << lab << endl;
  int n2r = Quad::class_group_2_rank;

  // Principal dimension
  out << dimension(0);
  // only if class number even..
  if (n2r)
    {
      // full dimension:
      out << " " << dimension(1);
      // character
      out << " ";
      vec_out(out, epsvec, "", ""); // no pre-/post []
    }
  out << endl;

  // Principal Hecke field:
  out << F->str(1) << endl;  // raw=1

  // Full Hecke field (even class number only)

  if (n2r)
    {
      // Multiquadratic extension data:
      out << Fmodsq->str(1) << endl;  // raw=1

      // Absolute field:
      out << Fabs->str(1) << endl;  // raw=1

      // Embedding of F into Fabs: matrix entries and denom
      out << abs_emb.str(1) <<endl;

      // Images of sqrt gens in Fabs:
      for (auto g: im_gens)
        out << g.str(1) << " ";

      // Self-twist discriminant or 0 (if class number even):
      out << "\n" << CMD << endl;
    }

  // Base-change code  +1 for base-change, -1 for twisted bc, 0 for neither, 2 for don't know
  out << base_change_code() << endl;

  out << endl;

  if (triv_char)
    {
      // Output A-L eigenvalues if computed, else output nothing
      vector<Quadprime> bads = nf->badprimes;
      if (conj)
        {
          Qideal Nbar = (nf->N).conj();
          bads = Nbar.factorization().sorted_primes();
        }
      for (auto& Q: bads)
        {
          int eQ=0;
          if (!eQmap.empty())
            eQ = eQmap.at(conj? Q.conj() : Q);
          out << prime_label(Q) << " " << eQ << endl;
        }
      out << endl;
    }

  // Output aP, preserving the order of the primes in thee list
  for (auto x: aPmap)
    {
      if (conj)
        {
          // NB1 For conjugate primes which are both good or both bad,
          // we will have made sure that if one is in the map then so
          // is the other, but it may happen that one is bad and its
          // conjugate is good and missing. e.g. field 17, level
          // 106.1, nap=15.

          // NB2 We swap the eigenvalues of P and conj(P); if the
          // character is nontrivial (then we must have C4 class
          // group, or there would not be any a(P) anyway), we also
          // swap the factors (1+i) and (1-i), as otherwise we would
          // be storing a newform with character chi_3 not chi_1.
          auto it = aPmap.find(x.first.conj());
          if (it!=aPmap.end())
            {
              Eigenvalue aP = aPmap.at(x.first.conj());
              // cout << "Conjugate of " << aP << " is ";
              aP = aP.conj();
              // cout << aP << endl;
              out << x.first << " " << aP.str(1) << endl;
            }
        }
      else
        {
          Eigenvalue aP = x.second;
          out << x.first << " " << aP.str(1) << endl;
        }
    }
}

//#define DEBUG_TWIST

// Construct another newform which is the unramified quadratic twist
// of this one by D, where D is a discriminant divisor
Newform Newform::unram_quadratic_twist(const INT& D) const
{
#ifdef DEBUG_TWIST
  cout << "About to twist one newform by " << D << endl;
#endif
  Newform f = *this; // copy

  // Twist aP (change a copy then reassign):
  auto aPmapD(aPmap);
  for (const auto& x: aPmap)
    {
      auto P = x.first;
      if (P.genus_character(D) == -1)
        aPmapD[P].negate();
    }
  f.aPmap = aPmapD;

  // Twist eQ (change a copy then reassign):
  auto eQmapD(eQmap);
  for (const auto& x: eQmap)
    {
      auto Q = x.first;
      if (Q.genus_character(D) == -1)
        eQmapD[Q] *= -1;
    }
  f.eQmap = eQmapD;

  // Twist aM and traces:
  Qidealooper Mloop(1, I2long(maxP), 1, 1); // 1: both conjugates; 1: sorted
  int iM=0;
  while (Mloop.not_finished())
    {
      Qideal M = Mloop.next();
#ifdef DEBUG_TWIST
      cout << "Ideal M = " << M << " has character " << M.genus_character(D) <<endl;
#endif
      if (M.genus_character(D) == -1)
        {
#ifdef DEBUG_TWIST
          cout << "Before twisting, aM = " << f.aMmap[M] << " with trace " << f.trace_list[iM] << endl;
#endif
          f.aMmap[M].negate();
          f.trace_list[iM] *= -1;
#ifdef DEBUG_TWIST
          cout << "After twisting, aM = " << f.aMmap[M] << " with trace " << f.trace_list[iM] << endl;
#endif
        }
      iM++;
    }
  return f;
}

//#define DEBUG_UNRAM_TWISTS

void Newform::set_unramified_twist_discriminants()
{
  if (Quad::class_group_2_rank==0)
    return;
  if (!triv_char)
    return;
  unramified_twist_discriminants.clear();

#ifdef DEBUG_UNRAM_TWISTS
  cout << "Determining all unramified quadratic twists up to Galois conjugacy" <<endl;
  cout << "  genus_classes_no_ext =    " << genus_classes_no_ext << endl;
  cout << "  genus_class_no_ext_ideals " << genus_class_no_ext_ideals << endl;
#endif
  // We use the subgroup H of the genus class group whose element list
  // is genus_classes_no_ext with representative ideals in
  // genus_class_no_ext_ideals.

  // Special case 1: genus_classes_no_ext = [0]. Then all unramified
  // quadratic twists are conjugates, so we return an empty list.
  // This will always be the case if the class number is odd, so the
  // class group's 2-rank is 1.

  int nclasses = genus_classes_no_ext.size();
  if (nclasses == 1)
    {
#ifdef DEBUG_UNRAM_TWISTS
      cout << "No nontrivial unramified quadratic twists up to Galois conjugacy" << endl;
#endif
    return;
    }

  // Special case 2: genus_classes_no_ext = all genus classes. Then
  // all unramified quadratic twists are pairwise non-conjugate, so we
  // return all of them (except the trivial one).  This is the only
  // other case if the class group's 2-rank is 1.

  int n2 = Quad::class_group_2_torsion.size();

  if (nclasses == n2)
    {
#ifdef DEBUG_UNRAM_TWISTS
      cout << "All nontrivial unramified quadratic twists are distinct up to Galois conjugacy" << endl;
#endif
      unramified_twist_discriminants.reserve(Quad::all_disc_factors.size()-1);
      std::copy(Quad::all_disc_factors.begin()+1, Quad::all_disc_factors.end(),
                unramified_twist_discriminants.end());
      return;
    }

  // Now the class group's 2-rank must be at least 2.

#ifdef DEBUG_UNRAM_TWISTS
  cout << "Finding D such that the unramified quadratic twists by D are nontrivial and distinct up to Galois conjugacy" << endl;
#endif
  // genus_class_no_ext_ideals is a list of ideals in the classes in
  // genus_classes_no_ext.size().  To find a list of the twists D we
  // need, we evaluate the genus character for D at all of these and
  // only keep D if the list of values is new.
  vector<vector<int>> D_values_list(1, vector<int>(genus_class_no_ext_ideals.size(), 1));
#ifdef DEBUG_UNRAM_TWISTS
  cout << "Initial D_values_list" << endl;
  for (auto L: D_values_list) cout << L << endl;
#endif
  vector<INT> unramified_twist_discriminants;
  for (auto& D : Quad::all_disc_factors)
    {
      if (!D.is_one())
        {
#ifdef DEBUG_UNRAM_TWISTS
          cout << "\nTesting D = " << D << endl;
#endif
          vector<int> D_values;
          for (auto I: genus_class_no_ext_ideals)
            D_values.push_back(I.genus_character(D));
#ifdef DEBUG_UNRAM_TWISTS
          cout << "D_values = [I.genus_character(D) for I in " << genus_class_no_ext_ideals << "] = "
               << D_values << endl;
#endif
          auto it = std::find(D_values_list.begin(), D_values_list.end(), D_values);
          if (it == D_values_list.end())
            {
#ifdef DEBUG_UNRAM_TWISTS
              cout << "New D_values vector, adding to list and keeping D = " << D << endl;
#endif
              D_values_list.push_back(D_values);
              unramified_twist_discriminants.push_back(D);
#ifdef DEBUG_UNRAM_TWISTS
              cout << "D_values_list is now" << endl;
              for (auto L: D_values_list) cout << L << endl;
              cout << "unramified_twist_discriminants is now " << unramified_twist_discriminants << endl;
#endif
            }
#ifdef DEBUG_UNRAM_TWISTS
          else
            cout << "Repeat D_values vector, not keeping D = " << D << endl;
#endif
        }
    }
#ifdef DEBUG_UNRAM_TWISTS
  cout << "The nontrivial unramified quadratic twists up to Galois conjugacy are "
       << " by D in " << unramified_twist_discriminants << endl;
#endif
}

// Construct all newforms which are nontrivial unramified quadratic
// twists of this one, up to Galois conjugacy
vector<Newform> Newform::all_unram_quadratic_twists() const
{
  vector<Newform> twisted_newforms; // will hold the result
  if ((Quad::class_group_2_rank==0) ||  (!triv_char))
    return twisted_newforms;

  twisted_newforms.reserve(unramified_twist_discriminants.size());
  for (auto& D : unramified_twist_discriminants)
    {
#ifdef DEBUG_UNRAM_TWISTS
      cout << "Applying unramified quadratic twist  by " <<  D << endl;
#endif
      Newform fD = unram_quadratic_twist(D);
      twisted_newforms.push_back(fD);
    }
  return twisted_newforms;
}

// Input newform data (needs lab to be set to construct the filename).
// Returns 0 if data file missing, else 1
int Newform::input_from_file(int verb)
{
  string fname = filename();
  if (verb)
    cout << "Reading newform " << lab << " from " << fname << " (verb="<<verb<<")"<<endl;
  ifstream fdata(fname.c_str());
  if (!fdata.is_open())
    {
      cerr << "Newform file " << fname << " missing" << endl;
      return 0;
    }
  // Check field, level, letter-code:
  string dat;
  fdata >> dat;
  assert (dat==field_label());
  fdata >> dat;
  assert (dat==nf->short_label());
  fdata >> dat;
  assert (dat==lab);

  int n2r = Quad::class_group_2_rank;

  // Principal dimension and (if class number even) full dimension:
  fdata >> d;
  int fd = d;
  if (n2r)
    fdata >> fd;
  if (verb>1)
    cout << "--> Principal dim = " << d << " full dim = " << fd << endl;
  // character vector epsvec
  triv_char = 1;
  epsvec.resize(n2r);
  if (n2r)
    {
      for (auto&e: epsvec)
        fdata >> e;
      // set triv_char flag
      triv_char = std::all_of(epsvec.cbegin(), epsvec.cend(), [](int i) { return i == +1; });
      if (verb>1)
        {
          cout << "--> Character = " << epsvec;
          if (triv_char) cout << " (trivial)";
          cout << endl;
        }
    }

  // Principal Hecke field:
  F = new Field();
  fdata >> &F;
  if (verb>1)
    cout << "--> Principal field is " << *F << endl;

  // Full Hecke field data (if class number even)
  Fmodsq = new FieldModSq(F);
  Fabs = new Field();

  if (n2r)
    {
      // Multiquadratic extension data:
      fdata >> *Fmodsq;
      if (verb>1)
        cout << "--> Field-mod-squares data: " << *Fmodsq << endl;

      // Absolute field:
      fdata >> &Fabs;
      if (verb>1)
        cout << "--> Absolute field: " << *Fabs << endl;

      // Embedding of F into Fabs: matrix entries and denom
      abs_emb = FieldIso(F, Fabs);
      fdata >> abs_emb;
      if (verb>1)
        cout << "--> Embedding into absolute field: " << abs_emb << endl;

      // Images of sqrt gens in Fabs:
      im_gens.resize(Fmodsq->rank(), FieldElement(Fabs));
      for (auto& g: im_gens)
        fdata >> g;
      if (verb>1)
        cout << "--> sqrt gens in absolute field: " << im_gens << endl;

      // Self-twist discriminant or 0 (if class number even):
      fdata >> CMD;
      self_twist_flag = (CMD<0);
      if (verb>1)
        {
          cout << "--> CMD = " << CMD;
          if (self_twist_flag)
            cout << " (self-twist)";
          cout  << endl;
        }
    }

  // Base-change code  +1 for base-change, -1 for twisted bc, 0 for neither, 2 for don't know
  int bcc;
  fdata >> bcc;
  bc=bct=-1;
  if (bcc==1)       {bc = bct = 1;}
  else if (bcc==-1) {bc = 0; bct = 1;}
  else if (bcc==0)  {bc = bct = 0;}

  if (verb>1)
    {
      cout << "Before reading eigenvalues:" << endl;
      display();
    }

  sfe = 0; // dummy value in case of nontrivial character
  if (triv_char)
    {
      sfe = -1;
      //cout << "Now reading eQ for Q in " << nf->badprimes <<endl;
      // A-L eigenvalues (will read nothing if level is (1))
      for (auto Q: nf->badprimes)
        {
          string Qlab = prime_label(Q);
          if (verb>1)
            cout << "reading AL eig for Q = " << Q << " = " << Qlab << endl;
          string lab;
          fdata >> lab;
          if (verb>1)
            cout << "--> Q has label " << Qlab << ", file label " << lab << endl;
          if (lab!=Qlab)
            cerr << "!!! Q has label " << Qlab << " but read label " << lab << endl;
          assert (lab==Qlab);

          int eQ;
          fdata >> eQ;
          eQmap[Q] = eQ;
          sfe *= eQ;
          if (verb>1)
            cout << "AL eigenvalue for " << Qlab << " is " << eQ << endl;
        }
      if (verb>1)
        {
          cout << "After reading eQ, before reading aP:" << endl;
          display();
        }
    }

  // aP
  string Plab;
  Quadprime P;
  Eigenvalue aP(F, Fmodsq);
  // read whitespace, so if there are no aP on file it does not try to read any
  fdata >> ws;
  // keep reading lines until end of file
  while (!fdata.eof())
    {
      fdata >> Plab // prime label
            >> aP   // eigenvalue data
            >> ws;  // eat whitespace, including newline
      if (verb>1)
        cout << "read prime label " << Plab << endl;

      P = Quadprime(Plab);
      aPmap[P] = aP;
      if (verb>1)
        cout << "--> P = " << Plab
             << ": a_P = "<<aPmap[P]<<endl;
    }
  if (verb)
    {
      cout << "After reading everything from " << fname <<":" << endl;
      display();
    }
  return 1;
}

void Newspace::output_to_file(int conj)
{
  //  int echo=0; // Set to 1 to echo what is written to the file for debugging
  ofstream out;
  out.open(filename(conj).c_str());

  // Field, level, number of newforms:
  out << field_label() << " "
      << (conj? conj_label() : short_label())
      << " " << newforms.size() << endl;
  if (newforms.empty())
    {
      out.close();
      return;
    }

  int n2r = Quad::class_group_2_rank;

  // Homological dimensions:
  vector<int> dims = dimensions(0);
  for (auto d: dims) out << " " << d;
  out << endl;

  // Full dimensions and trivial character flag (for even class number):
  if (n2r)
    {
      dims = dimensions(1);
      for (auto d: dims)
        out << " " << d;
      out << endl;
      for (const auto& f: newforms)
        out << " " << f.triv_char;
      out << endl;
    }

  out.close();
  for (const auto& f: newforms)
    {
      if (is_C4() || f.triv_char)
        f.output_to_file(conj);
    }
}

int Newspace::input_from_file(const Qideal& level, int verb)
{
  N = level;
  level_label = label(N);
  if (verb)
    cout << "In Newspace::input_from_file() with N="<<level_label << " = " << N << endl;
  Ndivs = alldivs(N, 1); // 1 for proper, i.e. excude N itself
  badprimes = N.factorization().sorted_primes();

  string fname = filename();
  ifstream fdata(fname.c_str());
  if (!fdata.is_open())
    {
      cerr << "File " << fname << " not available for Newspace input" << endl;
      return 0;
    }
  string dat;
  fdata >> dat;
  assert (dat==field_label());
  fdata >> dat;
  assert (dat==level_label);
  int nnf;
  fdata >> nnf;
  if (verb)
    cout << "-> Level "<<level_label<<" has " << nnf << " newforms in "<<fname<<endl;
  if (nnf==0)
    {
      fdata.close();
      return 1;
    }

  vector<int> pdims(nnf);
  for (auto& d: pdims)
    fdata >> d;
  if (verb>1)
    cout << "-> dims: " << pdims <<endl;
  vector<int> fdims = pdims;
  vector<int> triv_char_flags(nnf, 1); // set all to 1 as default if !n2r
  int n2r = Quad::class_group_2_rank;
  if (n2r)
    {
      for (auto& d: fdims)
        fdata >> d;
      if (verb>1)
        cout << "-> full dims: " << fdims <<endl;
      for (auto& tc: triv_char_flags)
        fdata >> tc;
      if (verb>1)
        cout << "-> trivial characater flags: " << triv_char_flags <<endl;
    }

  // This Newform constructor will read its data from file
  int C4 = is_C4();
  for (int i=1; i<=nnf; i++)
    {
      if (C4 || triv_char_flags[i-1])
        {
          if (verb)
            cout << "About to read newform #" << i << " from file" << endl;
          newforms.push_back(Newform(this, i, verb));
        }
    }

  if (verb>1)
    {
      cout << "Finished reading Newspace data with " << nnf << " newforms from " << fname << endl;
      for (int i=0; i<nnf; i++)
        {
          int tc = triv_char_flags[i];
          if (C4 || tc)
            {
              cout << "#" << (i+1) << ": " << "dim = " << pdims[i];
              if (n2r)
                cout << ", full dim = " << fdims[i];
              if (C4)
                cout << (tc? ", trivial" : ", non-trivial") << " character";
              cout << endl;
            }
        }
    }
  fdata.close();
  return 1;
}


// For each newform f, create and append all its unramified
// quadratic twists and resort.  This does nothing if the class
// number is odd.  We only add twists which are not Galois
// conjugates: that is, we twist by D in nontrivial cosets of a
// subgroup of the group of all discriminant divisors.
void Newspace::add_unram_quadratic_twists()
{
  if (Quad::class_group_2_rank == 0)
    return;

#ifdef DEBUG_TWIST
  cout << "Applying unramified quadratic twists " << endl;
#endif
  vector<Newform> all_extra_newforms;
  for ( const auto& F : newforms)
    if (F.triv_char)
      {
        vector<Newform> extra_newforms = F.all_unram_quadratic_twists();
#ifdef DEBUG_TWIST
        cout << "Adding " << extra_newforms.size() << " extra newforms, ";
#endif
        all_extra_newforms.insert(all_extra_newforms.end(), extra_newforms.begin(), extra_newforms.end());
#ifdef DEBUG_TWIST
        cout << "number is now " << newforms.size() << endl;
#endif
      }
  newforms.insert(newforms.end(), all_extra_newforms.begin(), all_extra_newforms.end());
  sort_newforms();
}

// output basis for the Hecke field and character of all newforms
void Newspace::display_newforms(int aP, int AL, int principal_eigs, int traces, int triv_char_only) const
{
  for ( auto& F : newforms)
    {
      if ((!triv_char_only) || F.triv_char)
        {
          F.display(aP, AL, principal_eigs, traces);
          cout<<endl;
        }
    }
}

// Return a list of the degrees of the principal or full Hecke fields
vector<int> Newspace::dimensions(int full) const
{
  vector<int> dims(newforms.size());
  std::transform(newforms.begin(), newforms.end(), dims.begin(),
                 [full](const Newform& F){return F.dimension(full);});
  return dims;
}

mat_m Newspace::heckeop(const gmatop& T, int cuspidal, int dual) const
{
  return to_mat_m(H1->calcop(T, cuspidal, dual, 0)); // 0 display
}

mat_m Newspace::heckeop(const matop& T, int cuspidal, int dual) const
{
  return to_mat_m(H1->calcop(T, cuspidal, dual, 0)); // 0 display
}

mat_m Newspace::heckeop(Quadprime& P, int cuspidal, int dual)
{
  return heckeop(AutoHeckeOp(P, N), cuspidal, dual);
}

// dict of Newspaces read from file
map<string,Newspace*> Newspace_dict;  // Key: label(N)

Newspace* get_Newspace(const Qideal& N, int verb)
{
  string Nlabel = label(N);
  if (Newspace_dict.find(Nlabel) != Newspace_dict.end())
    {
      if (verb)
        cout << "Newspace at level " << Nlabel << " retrieved from cache" << endl;
      return Newspace_dict[Nlabel];
    }
  if (verb)
    cout << "Newspace at level " << Nlabel << " not in cache, reading from file..." << endl;
  Newspace* NSP = new Newspace(N, verb);
  Newspace_dict[Nlabel] = NSP;
  if (verb)
    cout << "Newspace at level " << Nlabel << " read from file and cached" << endl;
  return NSP;
}
