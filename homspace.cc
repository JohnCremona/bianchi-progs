// FILE HOMSPACE.CC: Implemention of class homspace

#include "cusp.h"
#include "homspace.h"
#include "matprocs.h"
#include <assert.h>

homspace::homspace(const Qideal& I, scalar mod, int hp, int verb, scalar ch)
  :verbose(verb), plusflag(hp), N(I), P1(I), characteristic(ch), modulus(mod), hmod(0)
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
  vec_i pivs, npivs;
  scalar d2;
  smat sdeltamat(deltamat);
  kern = kernel(sdeltamat, modulus);
  if (characteristic==0)
    {
      smat sk;
      int ok = liftmat(smat_elim(sdeltamat,modulus).kernel(npivs,pivs),modulus,sk,d2);
      if (!ok)
        cout << "**!!!** failed to lift modular kernel to char 0\n" << endl;
    }
  else
    {
      d2 = 1;
    }

  const smat& basiskern = basis(kern);
  tkernbas = transpose(basiskern);         // dim(kern) x rank
  if(verbose>1)
    cout<<"tkernbas = "<<tkernbas.as_mat()<<endl;

  cuspidal_dimension = dim(kern);
  denom2 = d2;
  denom3 = denom1 * denom2;

  if (verbose)
    {
      cout << "Basis of ker(delta):\n";
      cout << basiskern.as_mat();
      cout << "pivots: " << pivots(kern) << endl;
      for (int i=0; i<dimension; i++)
        cout << "generator "<< i << ": " << freemods[i] << endl;
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
  if (verb) cout<<"Constructing conjugate homspace for level "<<ideal_label(Nconj)<<"..."<<endl;
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
  if(cuspidal) m = restrict_mat(smat(m),kern).as_mat();
  if(dual) m = transpose(m);
  if (display) cout << "Matrix of " << T.name() << " = " << m;
  if (display && (dimension>1)) cout << endl;
  if (display)
    cout<<"done."<<endl;
  return m;
}

ZZX homspace::charpoly(const matop& T, int cuspidal)
{
  ZZ den = to_ZZ(cuspidal? denom3: denom1);
  return scaled_charpoly(mat_to_mat_ZZ(calcop(T,cuspidal,0,0)), den, hmod);
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
  smat m(dimension,dimension);
  for (int j=0; j<dimension; j++)
    { svec colj(applyop(T,freemods[j]));
       m.setrow(j+1,colj);
     }
  if(cuspidal) // as above code computes the transpose
    {
      m = restrict_mat(transpose(m),kern);
      if(dual) m = transpose(m);
    }
  else
    if(!dual) {m=transpose(m);}
  // if (display)
  //   {
  //     cout << "Matrix of " << T.name() << " = ";
  //     if (dimension>1) cout << "\n";
  //     cout<<m.as_mat();
  //   }
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
  // if (display)
  //   {
  //     cout << "Matrix of " << T.name() << " = " << m.as_mat();
  //     if (dimension>1) cout << endl;
  //   }
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

// The subspace cut out by the given eigenvalues for the basis of
// Quad::class_group_2_rank unramified characters:
ssubspace homspace::unramified_character_subspace(const vector<int>& eigs)
{
  int subdim = h1dim();
  scalar den = h1denom();

  if (Quad::class_group_2_rank==0) // no characters, so return full (resp. cuspidal) space
    return ssubspace(subdim);

  vector<Qideal> nulist = make_nulist(N);
  auto nui = nulist.begin();
  auto ei = eigs.begin();

  // Compute eigenspace for first character, then compute successive
  // dual subspaces.

  int dual = 1;
  smat m = s_calcop(CharOp(*nui++, N), /*cuspidal*/ 0, dual, /*display*/ 0);
  scalar eig = (*ei++)*den;
  ssubspace s = eigenspace(m, eig, modulus);
  subdim = dim(s);

  for (; nui!=nulist.end() && subdim>0; ++ei)
    {
      m = s_calcop_restricted(CharOp(*nui++, N), s, dual, 0);
      eig = (*ei)*den;
      s = combine(s, eigenspace(m, eig, modulus));
      subdim = dim(s);
    }
  return s;
}

pair<int,int> homspace::unramified_character_subspace_dimensions(const vector<int>& eigs)
{
  ssubspace s = unramified_character_subspace(eigs);
  return {dim(s), (mult_mod_p(tkernbas, s.bas(), modulus)).rank(modulus)};
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

  ssubspace s = trivial_character_subspace();

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
                    cout << " - computed "<<op.length()<<" op matrices for T(P^2)*T(A,A) for A = " << ideal_label(A) << endl;
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

// test for cupidality (of a non-dual subspace only)
int homspace::is_cuspidal(const subspace& s) const
{
  return (deltamat*s.bas()).is_zero();
}

// functions for caching homspaces, full Hecke polynomials and new Hecke polynomials

map<Qideal,homspace*> H1_dict;
map<pair<Qideal,Quadprime>, ZZX> full_poly_dict;
map<pair<Qideal,Quadprime>, ZZX> new_poly_dict;

homspace* get_homspace(const Qideal& N, scalar mod)
{
  auto res = H1_dict.find(N);
  if (res==H1_dict.end())
    {
      homspace* H = new homspace(N, mod, 1); // cuspidal=1
      H1_dict[N] = H;
      return H;
    }
  else
    return H1_dict[N];
}

ZZX get_full_poly(const Qideal& N,  Quadprime& P, const scalar& mod)
{
  pair<Qideal,Quadprime> NP = {N,P};
  auto res = full_poly_dict.find(NP);
  if (res==full_poly_dict.end())
    {
      homspace* H = get_homspace(N, mod);
      ZZX full_poly = H->charpoly(HeckePOp(P, N), 1); // 1 for cuspidal
      full_poly_dict[NP] = full_poly;
      return full_poly;
    }
  else
    return full_poly_dict[NP];
}

ZZX get_new_poly(Qideal& N, Quadprime& P, const scalar& mod)
{
  pair<Qideal,Quadprime> NP = {N,P};
  auto res = new_poly_dict.find(NP);
  if (res==new_poly_dict.end())
    {
      ZZX new_poly = get_full_poly(N, P, mod);
      vector<Qideal> DD = alldivs(N);
      for( auto D : DD)
        {
          if (D==N)
            continue;
          ZZX new_poly_D = get_new_poly(D, P, mod);
          Qideal M = N/D;
          int mult = alldivs(M).size();
          for (int i=0; i<mult; i++)
            new_poly /= new_poly_D;
        }
      new_poly_dict[NP] = new_poly;
      return new_poly;
    }
  else
    return new_poly_dict[NP];
}
