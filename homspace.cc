// FILE HOMSPACE.CC: Implemention of class homspace

#include "euclid.h"
#include "cusp.h"
#include "homspace.h"
#include <assert.h>

homspace::homspace(const Qideal& I, int hp, int cuspid, int verb, long ch)
  :verbose(verb), cuspidal(cuspid), plusflag(hp), N(I), characteristic(ch), hmod(0)
{
  P1 = P1N(N);
  nsymb = P1.size();

  if (verbose)
    {
      cout << nsymb << " symbols";
      if (!Quad::is_Euclidean)
        {
          cout << ", " << n_alphas;
          if (Quad::class_number>1)
            cout << " + " << n_sigmas-1;
          cout << " types, total " << (n_alphas+n_sigmas-1)*nsymb << " edges" << endl;
        }
      else
        cout << endl;
    }

  ER = edge_relations(&P1, hp, verb, characteristic);
  ngens = ER.get_ngens();

  FR = face_relations(&ER, hp, verb, characteristic); // fills relmat with the relations and solves
  denom1 = FR.get_denom();
  rk = FR.get_rank();
  hmod = FR.get_hmod();

  make_freemods();

  kernel_delta();

  if(verbose)
    {
      cout << "number of cusps = " << ncusps << endl;
      if (cuspidal)
        cout << "dimension = " << dimension << endl;
      cout<<"denom1 = "<<denom1<<endl;
      cout<<"denom2 = "<<denom2<<endl;
      cout<<"denom3 = "<<denom3<<endl;
      cout << "Finished constructing homspace.\n";
    }
}

// for i >=0, the i'th edge as a modsym
modsym homspace::edge_generator(long i)
{
  pair<long, int> st = ER.symbol_number_and_type(i);
  return modsym(P1.lift_to_SL2(st.first), st.second);
}
 
void homspace::make_freemods()
{
  if (rk==0) return;

  long i;
  modsym m;

  freegens.resize(rk);
  for (i=0; i<rk; i++)
    freegens[i] = ER.gen(FR.gen(i+1));
  if (verbose)
    {
      cout << "freegens: ";
      for (i=0; i<rk; i++) cout << freegens[i] << " ";
      cout << endl;
      cout << "Freemods:\n";
    }

  for (i=0; i<rk; i++)
    {
      m = edge_generator(freegens[i]);
      freemods.push_back(m);
      if (verbose)
        {
          long j = freegens[i];
          pair<long, int> st = ER.symbol_number_and_type(j);
          long s = st.first;  // (c:d) symbol number
          long t = st.second; // symbol type (negative for singular edges)
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
  vec ei(rk); // initially 0
  for (i=0; i<rk; i++)
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
  mat deltamat(2*rk,rk);
  int i, ia, ib;
  for (i=0; i<rk; i++)
    {
      modsym m = freemods[i];
      RatQuad a = m.alpha(), b = m.beta();
      if (verbose>1)
        cout<<"Finding index of cusp "<<b<<"..."<<flush;
      ib = cusps.index(b);
      if (verbose>1)
        cout<<" index is "<<ib<<endl;
      deltamat(ib+1, i+1) += 1;  // N.B. offset of 1
      if (verbose>1)
        cout<<"Finding index of cusp "<<a<<"..."<<flush;
      ia = cusps.index(a);
      if (verbose>1)
        cout<<" index is "<<ia<<endl;
      deltamat(ia+1, i+1) -= 1;
      if (verbose)
        cout << "#"<<i<<" -->  C"<<(ib+1)<<" - C"<<(ia+1)<<endl;
    }
  ncusps=cusps.count();
  if(verbose)
    {
      cout<<ncusps<<" inequivalent cusps encountered "<<endl;
      if(verbose>1)
        cout<<"Matrix of boundary map = "<<deltamat<<endl;
    }
  scalar modulus = (characteristic==0? DEFAULT_MODULUS: characteristic);
  vec pivs, npivs;
  int d2;
  kern = kernel(smat(deltamat), modulus);
  int ok = 1;
  if (characteristic==0)
    {
      smat sk;
      ok = liftmat(smat_elim(deltamat, modulus).kernel(npivs,pivs),MODULUS,sk,d2);
      if (!ok)
        cout << "**!!!** failed to lift modular kernel to char 0\n" << endl;
    }
  else
    {
      d2 = 1;
    }

  tkernbas = transpose(kern.bas());         // dim(kern) x rank
  if(verbose>1)
    cout<<"tkernbas = "<<tkernbas.as_mat()<<endl;

  dimension = (cuspidal? dim(kern): rk);
  denom2 = d2;
  denom3 = denom1 * denom2;

  const smat& basiskern = basis(kern);
  if (verbose)
    {
      cout << "Basis of ker(delta):\n";
      cout << basiskern.as_mat();
      cout << "pivots: " << pivots(kern) << endl;
    }
  for (i=0; i<rk; i++)
    {
      int n = (cuspidal? ! trivial(basiskern.row(i+1).as_vec()) : 1);
      needed.push_back(n);
      if (verbose)
        {
          cout << "generator "<< i << ": " << freemods[i];
          if (!n) cout << " (not needed)";
          cout << endl;
        }
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
  homspace H1conj(Nconj, plusflag, cuspidal, 0);
  if (H1conj.rk!=rk)
    {
      if (verb) cout<<"Error: this space and the conjugate space have differeent dimensions"<<endl;
      return 0;
    }
  if (verb)
    cout<<" - dimensions equal ("<<rk<<"), constructing matrix of conjugation map..."<<endl;
  mat conjmat1(rk,rk);
  for (int i=0; i<rk; i++)
    {
      conjmat1.setcol(i+1, H1conj.chain(freemods[i].conj()));
    }
  long conjmatrank1 = smat(conjmat1).rank();
  if (verb) cout<<" - conjugation map has rank "<<conjmatrank1<<endl;

  // Now the reverse map
  mat conjmat2(rk,rk);
  for (int i=0; i<rk; i++)
    {
      conjmat2.setcol(i+1,chain(H1conj.freemods[i].conj()));
    }
  long conjmatrank2 = smat(conjmat2).rank();
  if (verb) cout<<" - reverse conjugation map has rank "<<conjmatrank2<<endl;
  if (conjmatrank1==rk && conjmatrank2==rk)
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
  long ind = P1.index(c,d);
  long i= ER.coords(ind, type);
#ifdef DEBUG_CHAIN
  cout<<"Symbol ("<<c<<":"<<d<<") has index "<<ind<<" plus offset "<< ER.offset(type) <<" = "<<ind+ER.offset(type)
       <<", giving coordindex "<<i;
#endif
  if (i)
    {
      vec ans = reduce_modp(sign(i) * (proj? projcoord.row(abs(i)) : coords(abs(i))), hmod);
#ifdef DEBUG_CHAIN
      cout << ": coordinate vector "<<ans<<endl;
#endif
      return ans;
    }
  else
    return vec((proj? projcoord.ncols(): rk)); // zero vector
}

vec homspace::chain(const RatQuad& alpha, const RatQuad& beta, int proj)
// Instead of just {return chain(beta, proj) - chain(alpha, proj);},
// we apply a version of "Karim's trick" -- though only when alpha is
// a principal cusp.
{
  Quad a(alpha.num()), b(alpha.den()), x, y;
  Quad g = quadbezout(a,b, x,y);
  //  cout<<"gcd("<<a<<","<<b<<") = " << g <<endl;
  if (g==Quad::one)
    {
      mat22 M(b,-a, x,y);    // det(M)=1 and M(alpha) = 0
      assert (M.is_unimodular());
      Quad c = N.reduce(x), d = N.reduce(-b);
#ifdef DEBUG_CHAIN
      cout<<"Computing alpha->beta chain {"<<alpha<<","<<beta<<"}\n";
      cout<<"   translated to {0, "<<M(beta)<<"} with c="<<c<<", d="<<d<<"\n";
#endif
      return chain(M(beta), proj, c, d);
    }
  else
    {
#ifdef DEBUG_CHAIN
      cerr<<"chain(alpha,beta) with alpha="<<alpha<<" non-principal"<<endl;
#endif
      return reduce_modp(chain(beta, proj) - chain(alpha, proj),hmod);
    }
}

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
       pseudo_euclidean_step(a,b, t, c,d);
       //c = N.reduce(c); d = N.reduce(d); // reduce modulo the level

       // either t>=0 and we have a standard edge:
       if (t>=0)
         {
           u = alpha_inv[t];
#ifdef DEBUG_CHAIN
           cout<<" STEP (t="<<t<<", t'="<<u<<", (c:d)_t'=("<<c<<":"<<d<<")_"<<u<<" = "<< modsym(lift_to_SL2(N,c,d),u)<<") TO "<<RatQuad(a,b,1) << endl;
#endif
           // Look up this symbol, convert to a vector w.r.t. homology basis
           ans = reduce_modp(ans + chaincd(c, d, u, proj), hmod);
#ifdef DEBUG_CHAIN
           cout<<" partial coordinate vector = "<<ans<<endl;
#endif
         }
       else // t<0 means that the last step took us to sigma[|t|] via
            // a translation only.  We will not reach b=0; instead we
            // finish off by subtracting M{sigma[|t|],oo} where M has
            // second row (c,d). [See Lingham's thesis, p.77]
         {
           ans = reduce_modp(ans - chaincd(c, d, t, proj), hmod);
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
{ vec ans(rk);
  if (proj) ans.init(projcoord.ncols());
  for (vector<mat22>::const_iterator mi = T.mats.begin(); mi!=T.mats.end(); ++mi)
    ans = reduce_modp(ans + chain((*mi)(alpha), proj), hmod);
  return ans;
}

vec homspace::applyop(const matop& T, const modsym& m, int proj)
{ vec ans(rk);
  if (proj) ans.init(projcoord.ncols());
  for (vector<mat22>::const_iterator mi = T.mats.begin(); mi!=T.mats.end(); ++mi)
    {
      mat22 M = *mi;
      ans = reduce_modp(ans + chain(M(m.alpha()), M(m.beta()), proj), hmod);
    }
  return ans;
}

mat homspace::calcop(const matop& T, int dual, int display)
{
  if(display)
    cout<<"Computing " << T.name() <<"...";
  mat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     { vec colj = applyop(T,freemods[j]);
       m.setcol(j+1,colj);
     }
  if(cuspidal) m = restrict_mat(smat(m),kern).as_mat();
  if(dual) m = transpose(m);
  // if (display) cout << "Matrix of " << T.name() << " = " << m;
  // if (display && (dimension>1)) cout << endl;
  if (display)
    cout<<"done."<<endl;
  return m;
}

mat homspace::calcop_cols(const matop& T, const vec& jlist, int verb)
{
  if(verb)
    cout<<"Computing " << T.name() <<"...";
  int i, d = dim(jlist);
  mat m(d,rk);
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

smat homspace::s_calcop_cols(const matop& T, const vec& jlist, int verb)
{
  int i, d = dim(jlist);
  if(verb)
    cout<<"Computing " << T.name() <<"..."<<flush;//" in s_calcop_cols() with d="<<d<<", jlist="<<jlist<<"...";
  smat m(d,rk);
  for (i=1; i<=d; i++)
    {
      int j = jlist[i];
      svec colj = applyop(T,freemods[j-1]);
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

smat homspace::s_calcop(const matop& T, int dual, int display)
{
  if(display)
    cout<<"Computing " << T.name() <<"..."<<flush;//" in s_calcop()..."<<flush;
  smat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     { svec colj = applyop(T,freemods[j]);
       m.setrow(j+1,colj);
     }
  if(cuspidal)
    {
      m = restrict_mat(transpose(m),kern);
      if(dual) m = transpose(m);
    }
  else if(!dual) {m=transpose(m);}
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
  if(display)
    cout<<"Computing " << T.name()
        <<" restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
  long d=dim(s);
  mat m(d,rk);
  for (long j=0; j<d; j++)
     {
       long jj = pivots(s)[j+1]-1;
       vec colj = applyop(T,freemods[jj]);
       m.setrow(j+1,colj);
     }
  if(hmod)
    m = matmulmodp(m,basis(s),hmod);
  else
    m = m*basis(s);
  if(!dual) m=transpose(m); // dual is default for restricted ops
  // if (display) cout << "Matrix of " << T.name() << " = " << m;
  // if (display && (dimension>1)) cout << endl;
  if (display)
    cout<<"done."<<endl;
  return m;
}

smat homspace::s_calcop_restricted(const matop& T, const ssubspace& s, int dual, int display)
{
  if(display)
    cout<<"Computing " << T.name()// <<" in s_calcop_restricted()"
        <<" restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
  long d=dim(s);
  smat m(d,rk);
  for (long j=1; j<=d; j++)
     {
       long jj = pivots(s)[j];
       svec colj = applyop(T,freemods[jj-1]);
       m.setrow(j,colj);
     }
  m = mult_mod_p(m,basis(s),MODULUS);
  m.reduce_mod_p();
  if(!dual) m=transpose(m); // dual is default for restricted ops
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
  vector<Quad>::const_iterator res=resmodp.begin()+1;
  while(res!=resmodp.end())
    ans = reduce_modp(ans + chain(*res++,p, proj), hmod);
  return ans;
}

vec homspace::manintwist(const Quad& lambda, const vector<Quad>& res, vector<int> chitable, int proj)
{
  vec ans = chain(Quad::zero,lambda, proj), part;          // =0, but sets the right length.
  vector<int>::const_iterator chi=chitable.begin();
  vector<Quad>::const_iterator r=res.begin();
  while(r!=res.end())
    ans = reduce_modp(ans + (*chi++)*chain(*r++,lambda, proj), hmod);
 return ans;
}

// The subspace cut out by the given eigenvalues for the basis of
// Quad::class_group_2_rank unramified characters:
ssubspace homspace::unramified_character_subspace(const vector<int>& eigs, int c, int dual)
{
  if (Quad::class_group_2_rank==0)
    {
      return ssubspace(c? h1cuspdim(): h1dim());
    }
  vector<Qideal> nulist = make_nulist(N);
  vector<Qideal>::iterator nui = nulist.begin();
  vector<int>::const_iterator ei = eigs.begin();
  smat m = s_calcop(CharOp(*nui++, N), dual, 0);
  // cout<<"Finding common eigenspace of "<<nulist.size()<<" character involutions"<<endl;
  if (c && !cuspidal)
    m = restrict_mat(m,kern);
  long den = (c? h1cdenom(): h1denom());
  ssubspace s = eigenspace(m, (*ei++)*den);
  int subdim = dim(s);

  for (; nui!=nulist.end() && subdim>0; ++ei)
    {
      m = s_calcop(CharOp(*nui++, N), dual, 0);
      if (c && !cuspidal)
        m = restrict_mat(m,kern);
      s = subeigenspace(m, (*ei)*den, s);
      subdim = dim(s);
    }
  return s;
}

// list of (cuspidal) dimensions of subspaces on which all T(A,A)
// act trivially with self-twist by unramified quadratic char D for
// each D (including D=1, meaning no self-twist)
vector<int> homspace::trivial_character_subspace_dimension_by_twist(int c)
{
  ssubspace s = trivial_character_subspace(c, 0); // cuspidal, not dual

  vector<int> dimlist;
  dimlist.push_back(dim(s));
  int n2r = Quad::class_group_2_rank;
  if (n2r==0)
    return dimlist;

  int stdim = 0; // total dim of nontrivial self-twist subspaces
  int den = (c? h1cdenom(): h1denom());

  for(auto Di = Quad::all_disc_factors.begin()+1; Di!=Quad::all_disc_factors.end(); ++Di)
    {
      QUINT D = *Di;
      ssubspace sD = s;
      int subdim = dim(s);

      QuadprimeLooper Pi(N); // loop over primes not dividing N
      int ip = 0, np = 10;   // only use 10 non-square-class primes
      while (ip<np && subdim>0 && Pi.ok())
        {
          Quadprime P = Pi;
          if (P.genus_character(D) == -1)
            {
              ip++;
              long Pnorm = I2long(P.norm());
              long eig = -den*Pnorm;
              Qideal P2 = P*P;
              matop op;
              if (P2.is_principal())
                op = HeckeP2Op(P, N);
              else   // compute T(P^2)*T(A,A)
                {
                  Qideal A = P.equivalent_coprime_to(N, 1);
                  op = HeckeP2ChiOp(P,A,N);
                }
              smat m = s_calcop(op, 0, 0); // not dual, no display
              if (c && !cuspidal)
                m = restrict_mat(m,kern);
              sD = subeigenspace(m, eig, sD);
              subdim = dim(sD);
            }
          ++Pi;
        }
      stdim += subdim;
      dimlist.push_back(subdim);
    }
  dimlist[0] = dim(s) - stdim;
  return dimlist;
}
