// FILE HOMSPACE.CC: Implemention of class homspace

#include "cusp.h"
#include "homspace.h"
#include <assert.h>

homspace::homspace(const Qideal& I, scalar mod, int hp, int verb, scalar ch)
  :verbose(verb), plusflag(hp), N(I), P1(I), characteristic(ch), modulus(mod), hmod(0),
   triv_char_subdim(-1), triv_char_dual_subdim(-1), triv_char_cuspidal_subdim(-1)
{
  Quad::setup_geometry("geodata", 0); // will do nothing after the first time called
  nsymb = P1.size();

  if (verbose)
    {
      cout << nsymb << " symbols";
      if (!Quad::is_Euclidean)
        {
          cout << ", " << Quad::SD.n_alph();
          if (Quad::class_number>1)
            cout << " + " << Quad::SD.n_sig(1);
          cout << " types, total " << (Quad::SD.n_alph()+Quad::SD.n_sig(1))*nsymb << " edges" << endl;
        }
      else
        cout << endl;
    }

  ER = edge_relations(&P1, hp, verb, characteristic);
  ngens = ER.get_ngens();

  FR = face_relations(&ER, modulus, hp, verb, characteristic); // fills relmat with the relations and solves
  denom1 = FR.get_denom();
  dimension = FR.get_rank();
  hmod = FR.get_hmod();

  make_freemods();

  kernel_delta();

  if(verbose)
    {
      cout << "number of cusps = " << ncusps << endl;
      cout << "dimension = " << dimension << endl;
      cout << "cuspidal dimension = " << cuspidal_dimension << endl;
      cout<<"denom1 = "<<denom1<<endl;
      cout<<"denom2 = "<<denom2<<endl;
      cout<<"denom3 = "<<denom3<<endl;
      cout << "Finished constructing homspace.\n";
    }
}

// for i >=0, the i'th edge as a modsym
modsym homspace::edge_generator(int i)
{
  pair<int, int> st = ER.symbol_number_and_type(i);
  mat22 U(P1.lift_to_SL2(st.first));
  int type = st.second;
  RatQuad a(type>=0? Quad::SD.alist[type]: Quad::SD.slist[-type]);
  return modsym(U(a), U.image_oo());
}

void homspace::make_freemods()
{
  if (dimension==0) return;

  int i;
  modsym m;

  freegens.resize(dimension);
  for (i=0; i<dimension; i++)
    freegens[i] = ER.gen(FR.gen(i+1));
  if (verbose)
    {
      cout << "freegens: ";
      for (i=0; i<dimension; i++) cout << freegens[i] << " ";
      cout << endl;
      cout << "Freemods:\n";
    }

  for (i=0; i<dimension; i++)
    {
      m = edge_generator(freegens[i]);
      freemods.push_back(m);
      if (verbose)
        {
          int j = freegens[i];
          pair<int, int> st = ER.symbol_number_and_type(j);
          int s = st.first;  // (c:d) symbol number
          int t = st.second; // symbol type (negative for singular edges)
          mat22 U = P1.lift_to_SL2(s);
          cout<<"--lifting symbol #"<<s<<" to SL2: "<<U
              <<", type "<<t<<" --> "<<m<<"\n"
              <<i<<": "<<m<<endl;
        }
    }

  if (verbose)
    {
      cout<<"Checking that freemods convert back to unit vectors:"<<endl;
    }
  vec ei(dimension); // initially 0
  for (i=0; i<dimension; i++)
    {
      m = freemods[i];
      if (verbose)
        cout<< m << " --> " << flush;
      vec v = chain(m);
      ei[i+1] = denom1;
      if (v!=ei && v!=-ei)
        {
          cerr<<endl;
          if (!verbose) cerr<< m << " --> ";
          cerr<< v<<endl;
          cerr<<" *** WRONG, should be "<<ei<<endl;
          exit(1);
        }
      else
        {
          if (verbose) cout << " OK"<<endl;
        }
      ei[i+1] = 0;
    }
}

void homspace::kernel_delta()
{
  if (verbose)
    cout<<"Computing boundary map"<<endl;
  cusplist cusps(N, plusflag);
  deltamat.init(2*dimension,dimension); // 2*dimension is an upper bound
  int nr = 0; // # distinct cusps encountered, = # rows of deltamat used
  for (int i=0; i<dimension; i++)
    {
      modsym m = freemods[i];
      RatQuad a = m.alpha(), b = m.beta();
      if (verbose>1)
        cout<<"Finding index of cusp "<<b<<"..."<<flush;
      int ib = cusps.index(b);
      if (verbose>1)
        cout<<" index is "<<ib<<endl;
      deltamat(ib+1, i+1) += 1;  // N.B. offset of 1
      nr = max(nr, ib+1);
      if (verbose>1)
        cout<<"Finding index of cusp "<<a<<"..."<<flush;
      int ia = cusps.index(a);
      if (verbose>1)
        cout<<" index is "<<ia<<endl;
      deltamat(ia+1, i+1) -= 1;
      nr = max(nr, ia+1);
      if (verbose)
        cout << "#"<<i<<" -->  C"<<(ib+1)<<" - C"<<(ia+1)<<endl;
    }
  ncusps=cusps.count();
  assert(nr==ncusps);
  deltamat = deltamat.slice(nr, dimension);
  if(verbose)
    {
      cout<<ncusps<<" inequivalent cusps encountered "
          <<" (deltamat has "<<nr<<"rows)"<<endl;
      if(verbose>1)
        cout<<"Matrix of boundary map = "<<deltamat<<endl;
    }
  if (characteristic!=0) modulus = scalar(characteristic);
  denom2 = scalar(1);  // relative denominator of ker(delta)
  smat sdeltamat(deltamat);
  kern = kernel(sdeltamat, modulus);
  // kern is a now subspace modulo modulus.  We now try to lift it to
  // char 0 by lifting the basis matrix and allowing a denominator d2
  const smat& basiskern = basis(kern);
  // Lift basis(kern) to char 0 to get the denominator (but kern
  // itself remains in characteristic modulus).
  if (characteristic==0)
    {
      smat sk;
      int ok = liftmat(basiskern,modulus,sk,denom2);
      if (ok) // replace basis(kern) by sk
        kern = ssubspace(sk, kern.pivs(), modulus);
      else
        cout << "**!!!** failed to lift modular kernel to char 0\n" << endl;
    }

  tkernbas = transpose(basiskern);         // dim(kern) x rank
  if(verbose>1)
    cout<<"tkernbas = "<<tkernbas.as_mat()<<endl;

  cuspidal_dimension = dim(kern); // "h1cuspdim()"
  denom3 = denom1 * denom2; // absolute denominator of ker(delta) "h1cdenom()"
  if (verbose)
    {
      cout << "Basis of ker(delta):\n";
      cout << basiskern.as_mat();
      cout << "pivots: " << pivots(kern) << endl;
      for (int i=0; i<dimension; i++)
        cout << "generator "<< i << ": " << freemods[i] << endl;
      cout << "homspace denom = " << denom1 <<endl;
      cout << "cuspidal relative denom = " << denom2 <<endl;
      cout << "cuspidal absolute denom = " << denom3 <<endl;
    }
}

// This method constructs the conjugate homspace (which might be the
// same again, if N is Galois-stable), and the maps between them
// induced by the obvious conjugation map on modular symbols.
// Theoretically, both spaces should have the same dimension, and the
// maps between them should be isomorphisms.

int homspace::check_conjugate(int verb)
{
  Qideal Nconj = N.conj();
  if (verb) cout<<"Constructing conjugate homspace for level "<<label(Nconj)<<"..."<<endl;
  homspace H1conj(Nconj, modulus, plusflag, 0);
  if (H1conj.dimension!=dimension)
    {
      if (verb) cout<<"Error: this space and the conjugate space have differeent dimensions"<<endl;
      return 0;
    }
  if (verb)
    cout<<" - dimensions equal ("<<dimension<<"), constructing matrix of conjugation map..."<<endl;
  mat conjmat1(dimension,dimension);
  for (int i=0; i<dimension; i++)
    {
      conjmat1.setcol(i+1, H1conj.chain(freemods[i].conj()));
    }
  int conjmatrank1 = smat(conjmat1).rank(modulus);
  if (verb) cout<<" - conjugation map has rank "<<conjmatrank1<<endl;

  // Now the reverse map
  mat conjmat2(dimension,dimension);
  for (int i=0; i<dimension; i++)
    {
      conjmat2.setcol(i+1,chain(H1conj.freemods[i].conj()));
    }
  int conjmatrank2 = smat(conjmat2).rank(modulus);
  if (verb) cout<<" - reverse conjugation map has rank "<<conjmatrank2<<endl;
  if (conjmatrank1==dimension && conjmatrank2==dimension)
    {
      if (verb)
        cout<<" - spaces are isomorphic, as expected."<<endl;
      return 1;
    }
  else
    {
      if (verb)
        {
          cout << " - the forwards  conjugation map has rank "<<conjmatrank1<<endl;
          cout << " - the backwards conjugation map has rank "<<conjmatrank2<<endl;
          cout << " products are\n"<<conjmat1*conjmat2<<"\nand"<<conjmat2*conjmat1<<endl;
        }
      return 0;
    }
}

//#define DEBUG_CHAIN

vec homspace::chaincd(const Quad& c, const Quad& d, int type, int proj)
{
  int ind = P1.index(c,d);
  int i= ER.coords(ind, type);
#ifdef DEBUG_CHAIN
  cout<<"Symbol ("<<c<<":"<<d<<") has index "<<ind<<" plus offset "<< ER.offset(type) <<" = "<<ind+ER.offset(type)
       <<", giving coordindex "<<i;
#endif
  if (i)
    {
      vec ans = reduce_mod_p(scalar(sign(i)) * (proj? projcoord.row(abs(i)) : coords(abs(i))), hmod);
#ifdef DEBUG_CHAIN
      cout << ": coordinate vector "<<ans<<endl;
#endif
      return ans;
    }
  else
    return vec((proj? projcoord.ncols(): dimension)); // zero vector
}
#undef DEBUG_CHAIN

vec homspace::chain(const RatQuad& alpha, const RatQuad& beta, int proj)
// Instead of just {return chain(beta, proj) - chain(alpha, proj);},
// we could apply a version of "Karim's trick" when either alpha or
// beta is principal.  But experiment showed that this is actually a
// bit slower.
{
  return reduce_mod_p(chain(beta, proj) - chain(alpha, proj),hmod);
}
#if(0)
{
  if (alpha.is_principal())
    {
      Quad a(alpha.num()), b(alpha.den()), x, y;
      Quad g(quadbezout(a,b, x,y)); // g=ax+by=1
#ifdef DEBUG_CHAIN
      cout<<"alpha = "<<alpha<<" = a/b with a="<<a<<", b="<<b<<", gcd="<<g<<endl;
#endif
      mat22 M(b,-a, x,y);    // det(M)=g=1 and M(alpha) = 0
      assert (M.is_unimodular());
      Quad c = N.reduce(x), d = N.reduce(-b);
#ifdef DEBUG_CHAIN
      cout<<"Computing alpha->beta chain {"<<alpha<<","<<beta<<"}\n";
      cout<<"   translated to {0, "<<M(beta)<<"} with c="<<c<<", d="<<d<<"\n";
#endif
      return chain(M(beta), proj, c, d);
    }
  if (beta.is_principal())
    {
      return -chain(beta, alpha, proj);
    }
#ifdef DEBUG_CHAIN
  cerr<<"chain(alpha,beta) with alpha="<<alpha<<" and beta="<<beta<<" non-principal"<<endl;
#endif
}
#endif

vec homspace::chain(const Quad& aa, const Quad& bb, int proj, const Quad& cc, const Quad& dd)
{
  Quad e, a(aa), b(bb), c(cc), d(dd), q, f;
  vec ans = chaincd(c,d,0,proj); // this is the path {0,oo} when (c:d)=(0:1) (the default)
  int t=0, u;
#ifdef DEBUG_CHAIN
  //   if (!Quad::is_Euclidean)
  cout<<" INIT (c:d)_0=("<<c<<":"<<d<<")_0 = "<< modsym(lift_to_SL2(N,c,d),0)<<") AT "<< RatQuad(a,b,1) << endl;
#endif
  while (!b.is_zero())
     {
       Quad::SD.pseudo_euclidean_step(a,b, t, c,d);
       //c = N.reduce(c); d = N.reduce(d); // reduce modulo the level

       // either t>=0 and we have a standard edge:
       if (t>=0)
         {
           u = Quad::SD.a_inv[t];
#ifdef DEBUG_CHAIN
           cout<<" STEP (t="<<t<<", t'="<<u<<", (c:d)_t'=("<<c<<":"<<d<<")_"<<u<<" = "<< modsym(lift_to_SL2(N,c,d),u)<<") TO "<<RatQuad(a,b,1) << endl;
#endif
           // Look up this symbol, convert to a vector w.r.t. homology basis
           ans = reduce_mod_p(ans + chaincd(c, d, u, proj), hmod);
#ifdef DEBUG_CHAIN
           cout<<" partial coordinate vector = "<<ans<<endl;
#endif
         }
       else // t<0 means that the last step took us to sigma[|t|] via
            // a translation only.  We will not reach b=0; instead we
            // finish off by subtracting M{sigma[|t|],oo} where M has
            // second row (c,d). [See Lingham's thesis, p.77]
         {
           ans = reduce_mod_p(ans - chaincd(c, d, t, proj), hmod);
#ifdef DEBUG_CHAIN
           cout<<" full coordinate vector = "<<ans<<endl;
#endif
           return ans;
         }
     }

   // We get here when b=0, so no singular edge was used
#ifdef DEBUG_CHAIN
   cout<<" full coordinate vector = "<<ans<<endl;
#endif
   return ans;
}

vec homspace::applyop(const matop& T, const RatQuad& alpha, int proj)
{ vec ans(dimension);
  if (proj) ans.init(projcoord.ncols());
  std::for_each(T.mats.begin(), T.mats.end(),
                [this, alpha, proj, &ans] (const mat22& M)
                {ans = reduce_mod_p(ans + chain(M(alpha), proj), hmod);});
  return ans;
}

vec homspace::applyop(const matop& T, const modsym& m, int proj)
{ vec ans(dimension);
  if (proj) ans.init(projcoord.ncols());
  std::for_each(T.mats.begin(), T.mats.end(),
                [this, m, proj, &ans] (const mat22& M)
                {ans = reduce_mod_p(ans + chain(M(m.alpha()), M(m.beta()), proj), hmod);});
  return ans;
}

mat homspace::calcop(const matop& T, int cuspidal, int dual, int display)
{
  if(display)
    cout<<"Computing " << T.name() <<"...";
  mat m(dimension,dimension);
  for (int j=0; j<dimension; j++)
     { vec colj = applyop(T,freemods[j]);
       m.setcol(j+1,colj);
     }
  if(cuspidal)
    m = restrict_mat(smat(m),kern).as_mat();
  if(dual)
    m = transpose(m);
  if (display)
    {
      cout<<"done."<<endl;
      cout << "Matrix of " << T.name() << " = " << m;
      if (dimension>1) cout << endl;
    }
  return m;
}

mat homspace::calcop(const gmatop& T, int cuspidal, int dual, int display)
{
  if(display)
    cout<<"Computing " << T.name() <<"...";
  mat m(dimension,dimension);
  auto ci = T.coeffs.begin();
  auto Ti = T.ops.begin();
  while (ci!=T.coeffs.end())
    {
      scalar c = *ci++;
      if (c !=0 )
        {
          mat mi = calcop(*Ti, cuspidal, dual, display);
          if (c!=1)
            mi *= c;
          m += mi;
        }
      ++Ti;
    }
  if (display)
    {
      cout<<"done."<<endl;
      cout << "Matrix of " << T.name() << " = " << m;
      if (dimension>1) cout << endl;
    }
  return m;
}

ZZX homspace::charpoly(const matop& T, int cuspidal)
{
  ZZ den = to_ZZ(cuspidal? denom3: denom1);
  return scaled_charpoly(mat_to_mat_ZZ(calcop(T,cuspidal,0,0)), den, to_ZZ(hmod));
}

ZZX homspace::charpoly(const gmatop& T, int cuspidal)
{
  ZZ den = to_ZZ(cuspidal? denom3: denom1);
  return scaled_charpoly(mat_to_mat_ZZ(calcop(T,cuspidal,0)), den, to_ZZ(hmod));
}

mat homspace::calcop_cols(const matop& T, const vec_i& jlist, int verb)
{
  if(verb)
    cout<<"Computing " << T.name() <<"...";
  int i, d = dim(jlist);
  mat m(d,dimension);
  for (i=1; i<=d; i++)
    {
      int j = jlist[i];
      vec colj = applyop(T,freemods[j-1]);
      m.setcol(i,colj);
     }
  if (verb)
    cout<<"done."<<endl;
  return m;
}

smat homspace::s_calcop_cols(const matop& T, const vec_i& jlist, int verb)
{
  int i, d = dim(jlist);
  if(verb)
    cout<<"Computing " << T.name() <<"..."<<flush;//" in s_calcop_cols() with d="<<d<<", jlist="<<jlist<<"...";
  smat m(d,dimension);
  for (i=1; i<=d; i++)
    {
      int j = jlist[i];
      svec colj(applyop(T,freemods[j-1]));
      m.setrow(i,colj);
     }
  if (verb)
    cout<<"done."<<endl;
  // if (verb)
  //   {
  //     cout << "Matrix of " << T.name() << " = ";
  //     if (dimension>1) cout << "\n";
  //     cout<<m.as_mat();
  //   }
  return m;
}

smat homspace::s_calcop(const matop& T, int cuspidal, int dual, int display)
{
  if(display)
    cout<<"Computing " << T.name() <<"..."<<flush;//" in s_calcop()..."<<flush;

  // First compute the matrix of T on the full space (dimension x
  // dimension).  Since smats have setrow() but not setcol() we do
  // this by rows, giving the transpose of the actual matrix:
  smat m(dimension,dimension);
  for (int j=0; j<dimension; j++)
    {
      svec colj(applyop(T,freemods[j]));
      m.setrow(j+1, colj);
     }

  // If we want the matrix of the restriction to the cuspidal subspace
  // we must apply restrict_mat to the untransposed matrix:
  if(cuspidal)
    {
      m = restrict_mat(transpose(m),kern);
      if(dual) // then re-transpose
        {
          m = transpose(m);
        }
    }
  else
    {
      if(!dual) // then un-transpose
        {
          m=transpose(m);
        }
    }
  if(display)
    cout<<"done."<<endl;
  return m;
}

mat homspace::calcop_restricted(const matop& T, const subspace& s, int dual, int display)
{
  int d=dim(s);
  if(display)
    cout<<"Computing " << T.name()
        <<" restricted to subspace of dimension "<<d<<" ..."<<flush;
  mat m(d,dimension);
  for (int j=0; j<d; j++)
     {
       int jj = pivots(s)[j+1]-1;
       vec colj = applyop(T,freemods[jj]);
       m.setrow(j+1,colj);
     }
  if(hmod!=0)
    m = matmulmodp(m,basis(s),hmod);
  else
    m = m*basis(s);
  if(!dual) m=transpose(m); // as above code computes the transpose
  // if (display) cout << "Matrix of " << T.name() << " = " << m;
  // if (display && (dimension>1)) cout << endl;
  if (display)
    cout<<"done."<<endl;
  return m;
}

smat homspace::s_calcop_restricted(const matop& T, const ssubspace& s, int dual, int display)
{
  int d=dim(s);
  if(display)
    cout<<"Computing " << T.name()// <<" in s_calcop_restricted()"
        <<" restricted to subspace of dimension "<<d<<" ..."<<flush;
  smat m(d,dimension);
  for (int j=1; j<=d; j++)
     {
       int jj = pivots(s)[j];
       svec colj(applyop(T,freemods[jj-1]));
       m.setrow(j,colj);
     }
  m = mult_mod_p(m,basis(s), modulus);
  if(!dual) m=transpose(m); // as above code computes the transpose
  if(display)
    cout<<"done."<<endl;
  return m;
}

vec homspace::maninvector(Quadprime& P, int proj)
{
  int t = P.is_principal();
  assert (t && "P should be principal in maninvector()");
  vector<Quad> resmodp=P.residues();
  Quad p = P.gen();
  vec ans = chain(Quad::zero,p, proj), part;             // =0, but sets the right length.
  int first = 1;
  for (const auto& r : resmodp)
    {
      if (first) {first=0; continue;} // skip resmodp[0]
      ans = reduce_mod_p(ans + chain(r,p, proj), hmod);
    }
  return ans;
}

vec homspace::manintwist(const Quad& lambda, const vector<Quad>& res, vector<int> chitable, int proj)
{
  vec ans = chain(Quad::zero,lambda, proj), part;          // =0, but sets the right length.
  auto chi=chitable.begin();
  auto r=res.begin();
  while(r!=res.end())
    ans = reduce_mod_p(ans + scalar(*chi++)*chain(*r++,lambda, proj), hmod);
 return ans;
}

// Return the (dual) subspace of the full space cut out by the given
// eigenvalues for the basis of Quad::class_group_2_rank unramified
// characters.  If dual=0 (only) and cuspidal=1, then return the
// subspace of the cuspidal subspace similarly cut out.
ssubspace homspace::unramified_character_subspace(const vector<int>& eigs, int cuspidal, int dual)
{
  if (Quad::class_group_2_rank==0) // no characters, so return cuspidal or full space
    return (cuspidal? kern: ssubspace(dimension));

  if (cuspidal && dual)
    cout << "options cuspidal=1, dual=1 invalid in unramified_character_subspace()" <<endl;

  vector<Qideal> nulist = make_nulist(N);
  auto nui = nulist.begin();
  auto ei = eigs.begin();
  scalar den = denom1;

  if (!dual) // start with the full or cuspidal subspace and compute
             // successive eigenspaces:
    {
      ssubspace s = (cuspidal? kern : ssubspace(dimension));
      int subdim = dim(s);
      while (nui!=nulist.end() && subdim>0)
        {
          smat m = restrict_mat(s_calcop(CharOp(*nui++, N), 0, 0, 0), s); // cuspidal=0, dual=0, display=0
          s = combine(s, eigenspace(m, (*ei++)*den, modulus));
          subdim = dim(s);
        }
      return s;
    }

  // Now dual=1 and cuspidal=0.  Compute eigenspace for first
  // character, then compute successive dual subspaces, using
  // calcop_restricted, which is more efficient.

  smat m = s_calcop(CharOp(*nui++, N), cuspidal, dual, /*display*/ 0);
  scalar eig = (*ei++)*den;
  ssubspace s = eigenspace(m, eig, modulus);
  int subdim = dim(s);

  while (nui!=nulist.end() && subdim>0)
    {
      m = s_calcop_restricted(CharOp(*nui++, N), s, 1, 0); // dual=1, display=0
      s = combine(s, eigenspace(m, (*ei++)*den, modulus));
      subdim = dim(s);
    }
  return s;
}

pair<int,int> homspace::unramified_character_subspace_dimensions(const vector<int>& eigs)
{
  ssubspace s = unramified_character_subspace(eigs, 0, 1); // cuspidal=0, dual=1
  return {dim(s), (mult_mod_p(tkernbas, s.bas(), modulus)).rank(modulus)};
}

// return triv_char_subspace, after computing if necessary
// NB cannot have cuspidal=dual=1
ssubspace homspace::trivial_character_subspace(int cuspidal, int dual)
{
  int& subd = (dual? triv_char_dual_subdim:
               (cuspidal? triv_char_cuspidal_subdim: triv_char_subdim));
  ssubspace& subs = (dual? triv_char_dual_subspace:
                     (cuspidal? triv_char_cuspidal_subspace: triv_char_subspace));
  if (subd==-1)
    {
      auto all_ones = vector<int>(Quad::class_group_2_rank, +1);
      subs = unramified_character_subspace(all_ones, cuspidal, dual);
      subd = dim(subs);
    }
  return subs;
}

// list of (total,cuspidal) dimensions of subspaces on which all T(A,A)
// act trivially with self-twist by unramified quadratic char D for
// each D (including D=1, meaning no self-twist)
vector<pair<int,int>> homspace::trivial_character_subspace_dimensions_by_twist(int use_lower_bounds, int use_cuspidal_lower_bounds, vector<int> lower_bounds, vector<int> cuspidal_lower_bounds)
{ //verbose=2;
  scalar den = h1denom();
  pair<int,int> subdims = {dimension, cuspidal_dimension};
  vector<pair<int,int>> dimlist;
  if (Quad::class_group_2_rank==0)
    {
      dimlist.push_back(subdims);
      return dimlist;
    }

  if(verbose>1)
    {
      cout<<"\nFull dimension = "<<dimension<<", cuspidal dimension = "<<cuspidal_dimension<<endl;
      cout<<"Finding trivial character subspace..."<<flush;
    }

  ssubspace s = trivial_character_subspace(0, 1); // cuspidal=0, dual=1

  pair<int,int> subdims0 = {dim(s), (mult_mod_p(tkernbas, s.bas(), modulus)).rank(modulus)};
  // we'll subtract dimensions of nontrivial self-twist spaces from dimlist[0]
  dimlist.push_back(subdims0);

  if(verbose>1)
    {
      cout<<"...done.  Dimension = "<<subdims0.first<<", cuspidal dimension = "<<subdims0.second<<endl;
      cout<<"Pushing these onto dimlist, which is now ";
      cout<<"[";
      for(auto di=dimlist.begin(); di!=dimlist.end(); ++di)
        {
          if (di!=dimlist.begin()) cout<<", ";
          cout << "(" << (di->first) << "," << (di->second) << ")";
        }
      cout<<"]"<<endl;
    }

  // In the loop, s does not change, but the subspace sD depends on D

  auto Di = Quad::all_disc_factors.begin()+1;
  auto lbds = lower_bounds.begin()+1;
  auto clbds = cuspidal_lower_bounds.begin()+1;
  int lbd=0, clbd=0;
  while(Di!=Quad::all_disc_factors.end())
    {
      INT D = *Di++;
      if(verbose>1)
        cout<<"D = "<<D<<":"<<endl;
      if (dimlist[0].first==0) // then previous D have exhausted the space
        {
          subdims = {0,0};
          dimlist.push_back(subdims);
          if(verbose>1)
            {
              cout << " whole space accounted for by previous D" << endl;
            }
          continue;
        }
      ssubspace sD = s;
      if (use_lower_bounds) lbd = *lbds++;
      if (use_cuspidal_lower_bounds) clbd = *clbds++;
      subdims = subdims0;
      int subdim = subdims.first;
      int MAXNREPEATS = 5;   // 4 not enough for d=165, level 1.1, D=-4
      int nrepeats = 0;      // stop when dimension has not changed MAXNREPEATS times
      QuadprimeLooper Pi(N); // loop over primes not dividing N
      int ip = 0, np = 10;   // only use first few non-square-class primes

      while (ip<np && subdims.first>0 && Pi.ok() && nrepeats<MAXNREPEATS)
        {
          if (use_lower_bounds && subdims.first<= lbd) break;
          if (use_cuspidal_lower_bounds && subdims.second<= clbd) break;
          Quadprime P = Pi;
          if (P.genus_character(D) == -1)
            {
              if(verbose>1)
                cout<<"Forcing aP=0 for P = "<<P<<": current dimension is "<<dim(sD)<<endl;
              ip++;
              long Pnorm = I2long(P.norm());
              scalar eig = -den*Pnorm;
              Qideal P2 = P*P, A;
              matop op;
              if (P2.is_principal())
                op = HeckeP2Op(P, N);
              else   // compute T(P^2)*T(A,A) where (A*P)^2 is principal
                {
                  A = P.equivalent_mod_2_coprime_to(N, 1);
                  op = HeckeP2ChiOp(P,A,N);
                }
              if(verbose>2)
                {
                  if (P2.is_principal())
                    cout << " - computed "<<op.length()<<" op matrices for T(P^2)" << endl;
                  else
                    cout << " - computed "<<op.length()<<" op matrices for T(P^2)*T(A,A) for A = " << label(A) << endl;
                }
              smat m = s_calcop_restricted(op, sD, 1, (verbose>1)); // dual, no display
              if(verbose>1)
                {
                  cout << " - computed matrix of this op restricted to current subspace" << endl;
                  cout << " - computing subeigenspace for eigenvalue " << eig << endl;
                }
              ssubspace newsD = combine(sD, eigenspace(m, eig, modulus));
              int newsubdim = dim(newsD);
              if(verbose>1)
                {
                  cout << " - subeigenspace has dimension " << newsubdim << ": ";
                  if (newsubdim==subdim)
                    cout << "repeat #"<<(nrepeats+1);
                  else
                    cout << "reduced by "<<(subdim-newsubdim);
                  cout << endl;
                }
              if (newsubdim==subdim)
                {
                  nrepeats++;
                }
              else
                {
                  nrepeats=0;
                  sD = newsD;
                  subdim = newsubdim;
                  subdims = {subdim,(mult_mod_p(tkernbas, sD.bas(), modulus)).rank(modulus)};
                }
              ++Pi; // extra increment, so we don't use both conjugates
            }
          ++Pi; // increment prime
        }
      // Now, sD is the D-self-twist subspace.  We record either its
      // dimension, or the cuspidal subdimension:
      dimlist.push_back(subdims);
      dimlist[0].first -= subdims.first;
      dimlist[0].second -= subdims.second;
      if(verbose>1)
        {
          cout<<"Pushing onto dimlist, which is now ";
          cout<<"[";
          for(auto di=dimlist.begin(); di!=dimlist.end(); ++di)
            {
              if (di!=dimlist.begin()) cout<<", ";
              cout << "(" << (di->first) << "," << (di->second) << ")";
            }
          cout<<"]"<<endl;
        }
    }
  return dimlist;
}

// list of total or cuspidal dimensions of subspaces on which all T(A,A)
// act trivially with self-twist by unramified quadratic char D for
// each D (including D=1, meaning no self-twist)
vector<int> homspace::trivial_character_subspace_dimensions_by_twist(int cuspidal, int use_lower_bounds,
                                                                     vector<int> cuspidal_lower_bounds)
{
  if (cuspidal)
    {
      auto dims = trivial_character_subspace_dimensions_by_twist(0, use_lower_bounds, {}, cuspidal_lower_bounds);
      vector<int> dims1(dims.size());
      std::transform(dims.begin(), dims.end(), dims1.begin(),
                     [] (const pair<int,int>& di) {return di.second;});
      return dims1;
    }
  else
    {
      auto dims = trivial_character_subspace_dimensions_by_twist(use_lower_bounds, 0, cuspidal_lower_bounds, {});
      vector<int> dims1(dims.size());
      std::transform(dims.begin(), dims.end(), dims1.begin(),
                     [] (const pair<int,int>& di) {return di.first;});
      return dims1;
    }
}

// test for cuspidality (of a non-dual subspace only)
int homspace::is_cuspidal(const subspace& s) const
{
  return (deltamat*s.bas()).is_zero();
}

// Functions for caching homspaces, full Hecke polynomials and new Hecke polynomials
// Keys are strings of the form Nlabel (for homspace) or Nlabel-Plabel (for Hecke polynomials)

string Nkey(const Qideal& N)
{
  return label(N);
}

string Nmodpkey(const Qideal& N, const scalar p)
{
  stringstream s;
  s << label(N) << "-mod-"<<p;
  return s.str();
}

string NPkey(const Qideal& N, const Qideal& P)
{
  stringstream s;
  s << label(N) << "-" << label(P);
  return s.str();
}

string NTkey(const Qideal& N, const matop& T)
{
  stringstream s;
  s << label(N) << "-" << T.name();
  return s.str();
}

// identical code to previous
string NTkey(const Qideal& N, const gmatop& T)
{
  stringstream s;
  s << label(N) << "-" << T.name();
  return s.str();
}

string NPmodpkey(const Qideal& N, const Quadprime& P, scalar p)
{
  stringstream s;
  s << label(N) << "-" << label(P) << "-mod-"<<p;
  return s.str();
}

// cache of homspaces keyed by level
map<string,homspace*> H1_dict;

// cache of operator matrices keyed by level and opname, not restricted to cuspidal
map<string, mat> full_mat_dict;

// cache of operator charpolys keyed by level and opname, not restricted to cuspidal
map<string, ZZX> poly_dict;
// cache of operator charpolys keyed by level and opname, restricted to cuspidal
map<string, ZZX> cuspidal_poly_dict;
// cache of operator new charpolys keyed by level and opname, not restricted to cuspidal
map<string, ZZX> new_poly_dict;
// cache of operator new charpolys keyed by level and opname, restricted to cuspidal
map<string, ZZX> new_cuspidal_poly_dict;
// char polys on trivial char subspace
map<string, ZZX> tc_poly_dict;
// char polys on trivial char cuspidal subspace
map<string, ZZX> tc_cuspidal_poly_dict;
// char polys on trivial char new subspace
map<string, ZZX> tc_new_poly_dict;
// char polys on trivial char new cuspidal subspace
map<string, ZZX> tc_new_cuspidal_poly_dict;

// mod-p versions
map<string,homspace*> H1_modp_dict;
map<string, ZZ_pX> full_poly_modp_dict;
map<string, ZZ_pX> new_poly_modp_dict;

void output_poly_dict(ostream& os, map<string, ZZX> D)
{
  for (auto key_pol: D)
    os<<key_pol.first<<" "<<key_pol.second<<endl;
}

map<string, ZZX> input_poly_dict(istream& is)
{
  map<string, ZZX> D;
  string key;
  ZZX poly;
  while (!is.eof())
    {
      is >> key >> poly;
      D[key] = poly;
    }
  return D;
}

void output_poly_dict(ostream& os, map<string, ZZ_pX> D)
{
  for (auto key_pol: D)
    os<<key_pol.first<<" "<<key_pol.second<<endl;
}

map<string, ZZ_pX> input_poly_dict(istream& is,  const ZZ& p)
{
  cout << "Reading poly dict..."<<flush;
  ZZ_p::init(ZZ(p));
  map<string, ZZ_pX> D;
  string key;
  ZZ_pX poly;
  while (!is.eof())
    {
      is >> key >> poly;
      D[key] = poly;
    }
  cout << "done.  Read " << D.size() << " items" << endl;
  return D;
}

homspace* get_homspace(const Qideal& N, scalar mod)
{
  string Nlabel = label(N);
  if (H1_dict.find(Nlabel) != H1_dict.end())
    return H1_dict[Nlabel];
  homspace* H = new homspace(N, mod, 1); // cuspidal=1
  H1_dict[Nlabel] = H;
  return H;
}

//#define DEBUG_GET_FULL_MAT

// Key is label(N)-T.name()
// Value is matrix of T on the full space (not restricted to cuspidal subspace)
mat get_full_mat(const Qideal& N,  const matop& T, const scalar& mod)
{
  string NT = NTkey(N,T);
#ifdef DEBUG_GET_FULL_MAT
  cout << "In get_full_mat() with matop key " << NT << endl;
#endif
  if (full_mat_dict.find(NT) != full_mat_dict.end())
    {
#ifdef DEBUG_GET_FULL_MAT
      cout << "key " << NT << " is in full_mat_dict, returning cached matrix" << endl;
#endif
      return full_mat_dict[NT];
    }
#ifdef DEBUG_GET_FULL_MAT
  cout << "key " << NT << " not in full_mat_dict, computing matrix" << endl;
#endif
  homspace* H = get_homspace(N, mod);
  mat M = H->calcop(T,0,0,0); // cuspidal=0, dual=0, display=0
  full_mat_dict[NT] = M;
#ifdef DEBUG_GET_FULL_MAT
  cout << "caching and returning matrix of " << NT << endl;
#endif
  return M;
}

// Key is label(N)-T.name()
// Value is matrix of T on the full space (not restricted to cuspidal subspace)
mat get_full_mat(const Qideal& N,  const gmatop& T, const scalar& mod)
{
  string NT = NTkey(N,T);
#ifdef DEBUG_GET_FULL_MAT
  cout << "In get_full_mat() with gmatop key " << NT << endl;
#endif
  if (full_mat_dict.find(NT) != full_mat_dict.end())
    {
#ifdef DEBUG_GET_FULL_MAT
      cout << "key " << NT << " is in full_mat_dict, returning cached matrix" << endl;
#endif
      return full_mat_dict[NT];
    }
#ifdef DEBUG_GET_FULL_MAT
  cout << "key " << NT << " not in full_mat_dict, computing matrix" << endl;
#endif
  int d = get_homspace(N, mod)->h1dim();
  mat M(d,d);
  if (d)
    {
      auto ci = T.coeffs.begin();
      auto Ti = T.ops.begin();
      while (ci!=T.coeffs.end())
        {
          scalar c = *ci++;
          if (c !=0 )
            {
              mat Mi = get_full_mat(N, *Ti, mod);
              if (c!=1)
                Mi *= c;
              M += Mi;
            }
          ++Ti;
        }
    }
  // We do not need to add to the dict if this gmatop consist of a
  // single matop, since that will have been done in the loop.
  if (full_mat_dict.find(NT) == full_mat_dict.end())
    full_mat_dict[NT] = M;
#ifdef DEBUG_GET_FULL_MAT
  cout << "caching and returning matrix of " << NT << endl;
#endif
  return M;
}

//#define DEBUG_GET_POLY

// from one of poly_dict, tc_poly_dict, cuspidal_poly_dict, tc_cuspidal_poly_dict
// depending on flags cuspidal & triv_char
ZZX get_poly(const Qideal& N,  const gmatop& T, int cuspidal, int triv_char, const scalar& mod)
{
  string NT = NTkey(N,T);
#ifdef DEBUG_GET_POLY
  cout << "In get_poly(), N = " << label(N) << ", T = " << T.name()
       << ", cuspidal = " << cuspidal
       << ", triv_char = " << triv_char
       << ", mod = " << mod
       <<endl;
  cout << "key = " << NT << endl;
#endif
  auto poly_cache = (cuspidal?
                     (triv_char? tc_cuspidal_poly_dict: cuspidal_poly_dict):
                     (triv_char? tc_poly_dict: poly_dict));
  auto new_poly_cache = (cuspidal?
                         (triv_char? tc_new_cuspidal_poly_dict: new_cuspidal_poly_dict):
                         (triv_char? tc_new_poly_dict: new_poly_dict));
  if (poly_cache.find(NT) != poly_cache.end())
    {
#ifdef DEBUG_GET_POLY
      cout << "key is in cache, returning " << str(poly_cache[NT]) << endl;
#endif
      return poly_cache[NT];
    }

  homspace* H = get_homspace(N, mod);
  scalar den = (cuspidal? H->h1cdenom() :H->h1denom());
#ifdef DEBUG_GET_POLY
  cout << "Homspace for level " << N << " obtained from get_homspace()" << endl;
  cout << "den =  " << den << endl;
  cout << "About to call get_full_mat()"<< endl;
#endif

  mat M = get_full_mat(N, T, mod);  // dimension x dimension
#ifdef DEBUG_GET_POLY
  cout << "Full matrix M of size " << M.nrows() << " obtained from get_full_mat()" << endl;
  cout << "den =  " << den << endl;
  //  output_flat_matrix(M);  cout << endl;
#endif
  if (cuspidal || triv_char)
    {
      // H->kern is the subspace ker(delta) of the full space
      // H->trivial_character_subspace([cuspidal]) is the [cuspidal] subspace on which all
      // unrmaified characters act trivially
      smat sM = smat(M);
      if (cuspidal)
        {
#ifdef DEBUG_GET_POLY
          cout << "Cuspidal case" << endl;
#endif
          if (Quad::class_group_2_rank && triv_char)
            {
              ssubspace tcsub = H->trivial_character_subspace(1, 0); // cuspidal=1, dual=0
              smat s;
              scalar tcden(1);
              int ok = liftmat(basis(tcsub),mod,s,tcden);
              if (ok) // replace basis(tcsub) by s
                tcsub = ssubspace(s, tcsub.pivs(), mod);
              else
                cout << "**!!!** failed to lift modular kernel to char 0\n" << endl;
              den *= tcden;
              sM = restrict_mat(sM, tcsub);
#ifdef DEBUG_GET_POLY
              cout << "Restricted M to cuspidal trivial char subspace" << endl;
              cout << "(which has relative denominator " << tcden
                   << ", absolute denominator " << den << ")" << endl;
#endif
            }
          else
            {
              sM = restrict_mat(sM,H->kern);  // cuspidal_dimension x cuspidal_dimension
#ifdef DEBUG_GET_POLY
              cout << "Restricted M to cuspidal subspace" << endl;
              cout << "(which has denominator " << den << ")" << endl;
#endif
            }
        }
      else // not cuspidal
        {
#ifdef DEBUG_GET_POLY
          cout << "Non-cuspidal case" << endl;
#endif
          if (Quad::class_group_2_rank && triv_char)
            {
              ssubspace tcsub = H->trivial_character_subspace(0, 0); // cuspidal=0, dual=0
              smat s;
              scalar tcden(1);
              int ok = liftmat(basis(tcsub),mod,s,tcden);
              if (ok) // replace basis(tcsub) by s
                tcsub = ssubspace(s, tcsub.pivs(), mod);
              else
                cout << "**!!!** failed to lift modular kernel to char 0\n" << endl;
              den *= tcden;
              sM = restrict_mat(sM, tcsub);
#ifdef DEBUG_GET_POLY
              cout << "Restricted M to non-cuspidal trivial char subspace" << endl;
              cout << "(which has relative denominator " << tcden
                   << ", absolute denominator " << den << ")" << endl;
#endif
            }
        }
      M = sM.as_mat();
    }
  ZZX full_poly =  scaled_charpoly(mat_to_mat_ZZ(M), to_ZZ(den), to_ZZ(H->hmod));
  poly_cache[NT] = full_poly;
  if (deg(full_poly)==0)
    new_poly_cache[NT] = full_poly;
  return full_poly;
}

// from one of new_poly_dict, tc_new_poly_dict, new_cuspidal_poly_dict, tc_new_cuspidal_poly_dict
// depending on flags cuspidal & triv_char

// NB In even class number this may raise an error when newspaces at
// lower levels have self-twist since the effect of this on lowering
// oldspace dimensions is ignored.
ZZX get_new_poly(const Qideal& N, const gmatop& T, int cuspidal, int triv_char, const scalar& mod)
{
  string NT = NTkey(N,T);
  auto new_poly_cache = (cuspidal?
                         (triv_char? tc_new_cuspidal_poly_dict: new_cuspidal_poly_dict):
                         (triv_char? tc_new_poly_dict: new_poly_dict));
  if (new_poly_cache.find(NT) != new_poly_cache.end())
    return new_poly_cache[NT];

  ZZX new_poly = get_poly(N, T, cuspidal, triv_char, mod);
  if (deg(new_poly)==0)
    {
      new_poly_cache[NT] = new_poly;
      return new_poly;
    }
  vector<Qideal> DD = alldivs(N);
  for( auto D : DD)
    {
      if (D==N)
        continue;
      ZZX new_poly_D = get_new_poly(D, T, cuspidal, triv_char, mod);
      if (deg(new_poly_D)==0)
        continue;
      Qideal M = N/D;
      int mult = alldivs(M).size();
      // The actual multiplicity may be less than this.  If an
      // irreducible factor of new_poly_D corresponds to a
      // newforms with self-twist by the unramified quadratic
      // character chi, then the old multiplicity of this factor
      // is the number of divisors D1 of N/D for which chi(D1)=+1.
      // BUT here we do not know which factors are self-twist.

      for (int i=0; i<mult; i++)
        {
          //essentially new_poly /= new_poly_D // but checking divisibility
          ZZX quo, rem;
          DivRem(quo, rem, new_poly, new_poly_D);
          if (IsZero(rem))
            new_poly = quo;
          else
            {
              cout << "Problem in get_new_poly("<<NT<<"), D="<<label(D)<<endl;
              cout << "Dividing " << str(new_poly) << " by " << str(new_poly_D)
                   << " gives quotient " << str(quo) <<", remainder "<< str(rem) << endl;
              cout << "Old multiplicities are smaller than expected."<<endl;
              cout << "This should be because a newform at a dividing level has inner twist"<<endl;
              if (Quad::class_group_2_rank==0)
                cout << "*** so it should not happen, as the class number "
                     << Quad::class_number<<" is odd!"<<endl;
            }
        }
      if (deg(new_poly)==0) // nothing left, new dimension must be 0
        break;
    } // end of loop over divisors
  new_poly_cache[NT] = new_poly;
  return new_poly;
}

// Return true iff T's new poly is squarefree and coprime to its old poly
int test_splitting_operator(const Qideal& N, const gmatop& T, const scalar& mod, int verbose)
{
  if (verbose)
    cout << "Testing " << T.name() << "..." << flush;
  ZZX f_new = get_new_poly(N, T, 1, 0, mod); // cuspidal=1, triv_char=0
  if (!IsSquareFree(f_new))
    {
      if (verbose>1)
        cout << "\n NO: new Hecke polynomial "<<str(f_new)<<" is not squarefree" << endl;
      return 0;
    }
  ZZX f_full = get_poly(N, T, 0, 0, mod); // cuspidal=0, triv_char=0
  ZZX f_old = f_full / f_new;
  if (!AreCoprime(f_new, f_old))
    {
      if (verbose>1)
        cout << "\n NO: new Hecke polynomial "<<str(f_new)
             <<" is not coprime to old Hecke polynomial "<<str(f_old)<<endl
             <<" (full polynomial is "<<str(f_full)<<")"<<endl;
      return 0;
    }
  if (verbose>1)
    cout << "\n YES: new Hecke polynomial is squarefree and coprime to old Hecke polynomial" << endl;
  return 1;
}

homspace* get_homspace_modp(const Qideal& N, scalar p)
{
  string Nlabel = Nmodpkey(N, p);
  auto res = H1_modp_dict.find(Nlabel);
  if (res==H1_modp_dict.end())
    {
      homspace* H = new homspace(N, p, 1, 0, p); // plus=1, verb=0
      H1_modp_dict[Nlabel] = H;
      return H;
    }
  else
    return H1_modp_dict[Nlabel];
}

ZZ_pX get_full_poly_modp(const Qideal& N,  const Quadprime& P, scalar p)
{
  string NP = NPmodpkey(N,P,p);
  auto res = full_poly_modp_dict.find(NP);
  if (res==full_poly_modp_dict.end())
    {
      homspace* H = get_homspace_modp(N, p);
      ZZ_pX full_poly = to_ZZ_pX(H->charpoly(HeckePOp(P, N), 1)); // 1 for cuspidal
      full_poly_modp_dict[NP] = full_poly;

      if (deg(full_poly)==0)
        {
          new_poly_modp_dict[NP] = full_poly;
        }
      return full_poly;
    }
  else
    return full_poly_modp_dict[NP];
}

ZZ_pX get_new_poly_modp(const Qideal& N, const Quadprime& P, scalar p)
{
  string NP = NPmodpkey(N,P,p);
  auto res = new_poly_modp_dict.find(NP);
  if (res==new_poly_modp_dict.end())
    {
      ZZ_pX new_poly = get_full_poly_modp(N, P, p);
      if (deg(new_poly)==0)
        {
          new_poly_modp_dict[NP] = new_poly;
          return new_poly;
        }
      vector<Qideal> DD = alldivs(N);
      for( auto D : DD)
        {
          if (D==N)
            continue;
          ZZ_pX new_poly_D = get_new_poly_modp(D, P, p);
          if (deg(new_poly_D)==0)
            continue;
          Qideal M = N/D;
          int mult = alldivs(M).size();
          // The actual multiplicity may be less than this.  If an
          // irreducible factor of new_poly_D corresponds to a
          // newforms with self-twist by the unramified quadratic
          // character chi, then the old multiplicity of this factor
          // is the number of divisors D1 of N/D for which chi(D1)=+1.
          // BUT here we do not know which factors are self-twist.

          for (int i=0; i<mult; i++)
            {
              //new_poly /= new_poly_D;
              ZZ_pX quo, rem;
              DivRem(quo, rem, new_poly, new_poly_D);
              if (IsZero(rem))
                new_poly = quo;
              else
                {
                  cout << "Problem in get_new_poly("<<NP<<"), D="<<label(D)<<endl;
                  cout << "Dividing " << new_poly << " by " << new_poly_D
                       << " gives quotient " << quo <<", remainder "<< rem << endl;
                  cout << "Old multiplicities are smaller than expected."<<endl;
                  cout << "This should be because a newform at a dividing level has inner twist"<<endl;
                  if (Quad::class_group_2_rank==0)
                    cout << "*** so it should not happen, as the class number "
                         << Quad::class_number<<" is odd!"<<endl;
                }
            }
          if (deg(new_poly)==0) // nothing left, new dimension must be 0
            {
              new_poly_modp_dict[NP] = new_poly;
              return new_poly;
            }
        }
      new_poly_modp_dict[NP] = new_poly;
      return new_poly;
    }
  else
    return new_poly_modp_dict[NP];
}
