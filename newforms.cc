//  newforms.cc

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <functional>   // std::multiplies
#include <numeric>   // std::multiplies
#include "looper.h"
#include "newforms.h"

#define MAXDEPTH 20 // maximum depth for splitting off eigenspaces

#define genus_class_triviality_bound 10 // if aP=0 for this number of
                                        // good primes in a genus
                                        // class then we assume that a
                                        // newform is self-twist and
                                        // all aP for P in this class
                                        // are 0.

// Notes on scaling:
//
// For a newform F (however scaled / normalised), the *integral
// periods* of F are the integrals of (the differential associated to)
// F along paths between Gamma_0(N)-equivalent cusps. The differential
// is scaled by a normalising factor a_K, which depends only on the
// field: the factor is a_K=(4\pi)^2/(w*|D|), where w is the number of
// units and D the discriminant; this is chosen so that the integral
// of F over {0,oo}, denoted I_F({0,oo}) = L(F,1).

// [The factor of w accounts for F being a sum over nonzero elements of
// the ring of integers, with the coefficients of \alpha and u*\alpha
// being the same for units u, while L(F,a) is a sum over integral
// ideals.]

// These integral periods are all integral multiples of a minimal
// positive real period P.  The map from integral homology to the
// associated period multiple takes gamma in H_1(-,Z) to n_F(gamma) =
// (scaled integral of F over gamma)/P, and the image of this map is
// (by definition) the full integer ring Z.

// Similarly the *cuspidal periods* of F are its integrals along *all*
// paths gamma between cusps, whether or not they are
// Gamma_0(N)-equivalent, and these are the integral multiples of a
// least positive real cuspidal period P_cusp.  In particular, P is
// such a multiple, say P=c*P_cusp, so P*Z is the subgroup of P_cusp*Z
// of index c; c is called the "cuspidal factor", equal to c=P/P_cusp.

// Let n_cusp(gamma) = (integral of F over gamma)/P_cusp for gamma in
// H_1(-,Z; cusps).

// The map from gamma to n_cusp(gamma) is modular and an eigenfunction
// for Hecke, so is a primitive basis vector for the associated dual
// eigenspace.  Since we do linear algebra on the dual (by transposing
// Hecke matrices), we can compute this map, up to scaling, by simply
// taking the dot product of our dual basis vector v with the coords
// of gamma with respect to the homology basis.  Regardless of whether
// our "freegens" generates the whole integral homology (w.r.t. cusps)
// or only a sublattice of finite index, i.e. whether denom1=1 or >1,
// this does not matter since the map n_cusp is primitive, by
// definition, so to get the correct values (up to sign) we just have
// to divide out these dot products by the gcd of all of them.

// In more detail: each edge is a (possibly rational) linear
// combination of the "freegens" generating cuspidal homology; we have
// a common denominator of all these (called "denom1") which we can
// ignore.  We take the dot product of our dual basis vector v with
// each of these, and divide out by the content of the result, to get
// a new longer primitive vector w, whose coordindates w_i, one for
// each edge e_i (modulo edge relations), are such that the integral
// of F over e_i is w_i*P_cusp.  The vectors w_i are stored in each
// newform as 'basis', and form the columns of the matrix projcoord.

// Thus, by definition, the values of n_cusp(gamma) as gamma ranges
// over Gamma_0(N) are coprime integers.  The vector of these integers
// is the 'basis' component of each newform; it has length 'ngens'
// with coordinates indexed by edges modulo edge-relations.  We do not
// compute the longer vector of length 'ngens'; each of its entries is
// 0 or +/- one of the previous vector's entries, so is still a
// primitive integer vector.

// Hence the map n_cusp(gamma) is almost encoded in the matrix
// projcoord (created by calling make_projcoord()), whose rows are
// indexed by edges modulo edge relations. To get the value for the
// j'th newform on an edge (c:d)_t we find the edge number i =
// offset(t)+index(c,d)), use that to look up
// k = ER.coords(i) = ER.coords(index(c,d), t), and the value is
// sign(k)*projcoord[abs(k),j].

// For example, n_cusp({0,oo}) is obtained this way with
// (c,d,t)=(0,1,0).  Now L(F,1) is the (scaled) integral of F over
// {0,oo} and hence we obtain the integer L(F,1)/P_cusp =
// n_cusp({0,oo}).  For L(F,1)/P we divide by the cuspidal factor c
// for F, obtaining the rational n_cusp({0,oo})/c.  [P=c*P_cusp so L/P
// = L/(c*P_cusp) = (L/P_cusp)/c.]

// A second way to compute L(F,1)/P_cusp which only involves integral
// periods is to use a "Manin vector" mvp for some good prime p, which
// is the sum over x mod p of {0,x/p}, since (1+N(p)-a(p))*n_cusp({0,oo})
// = n_cusp(mvp), hence L/P_cusp = n_cusp(mvp)/(1+N(p)-a(p)).

// NB in the newform *constructor* we cannot use projcoord since that
// is only computed after all the newforms are found.

// Sorting functions

// Compare two integers, using the natural order ...-2,-1,0,1,2,...
int less_ap(long a, long b)
{
  return sign(a-b);
}

// Compare two integer vectors lexicographically, using less_ap():
int less_apvec(const vector<long>& v, const vector<long>& w)
{
  auto vi=v.begin(), wi=w.begin();
  while(vi!=v.end())
    {
      int s = less_ap(*vi++,*wi++);
      if(s) return s;
    }
  return 0;
}

struct newform_eigs_comparer {
  bool operator()(const newform& f, const newform& g)
  {
    return less_apvec(f.eigs,g.eigs)==-1;
  }
}
  less_newform_eigs;

struct newform_aplist_comparer {
  bool operator()(const newform& f, const newform& g)
  {
    return less_apvec(f.aplist,g.aplist)==-1;
  }
}
  less_newform_lmfdb;

newform::newform(newforms* nfs, const vec& v, const vector<long>& eigs)
  : eigs(eigs)
{
  //cout<<"Constructing newform with eigs "<<eigs<<" and basis "<<v<<endl;
  nf=nfs;

  // convert basis vector from coords w.r.t. homology-gens (edges mod
  // face-relations) to coords w.r.t. edge-gens:
  basis = nf->lengthen_basis(v);

  if (nf->verbose)
  {
    cout << "denom = " << nf->h1->denom1 << endl;
    cout << "short newform basis = "<<v<<endl;
    cout << "long  newform basis = "<<basis<<endl;
  }

  genus_classes.resize(1,0);
  genus_class_ideals.resize(1,Qideal(ONE));
  genus_class_aP.resize(1,1);
  fake = 0;
  // to be set later
  CMD = 0;  // will be set to an unramified negative discriminant if self-twist
  genus_classes_filled = 0; // will be set to 1 when all/half genus classes have a nonzero aP
  cm = 1;   // will be set to 0 or a negative square-free integer if CM
  bc = 4;   // will be set to 0 or a square-free integer if base-change or twist of b.c.
  sfe = pdot = dp0 = lambdadot = matdot = 0;
  genus_class_trivial_counter.resize(nf->nchi, 0);
  possible_self_twists = nf->possible_self_twists; // may be cut down on computing aP later
}

// Compute AL eigs and SFE: fills aqlist and sets sfe
// Requires h1
void newform::compute_AL()
{
  aqlist.resize(nf->badprimes.size());
  std::transform(nf->badprimes.begin(), nf->badprimes.end(), aqlist.begin(),
                 [this] ( Quadprime& Q ) {return eigenvalueAtkinLehner(Q);});
  sfe = std::accumulate(aqlist.begin(),aqlist.end(),-1,std::multiplies<long>());
}

// Compute cuspidalfactor (needs long basis and bigtkernbas)
void newform::compute_cuspidalfactor()
{
  if(nf->characteristic || Quad::class_number>1)
    {
      cuspidalfactor=1; // place-holder
    }
  else
    {
      cuspidalfactor = I2long(content((nf->h1->bigtkernbas)*basis));
      if(nf->verbose>1)
        cout<<"cuspidalfactor = "<<cuspidalfactor<<endl;
    }
}

// Compute L/P ratio (needs cuspidalfactor)
void newform::compute_loverp()
{
  int pdot0 = I2long(abs(nf->zero_infinity[index]));
  dp0  =  1 + (nf->nP0) - nf->aP0[index-1];  // aP0 is based at 0
  pdot = I2long(abs(nf->mvp[index]));

  // compute L/P as n_F({0,oo}) = n_cusp({0,oo})/c
  loverp =  rational(pdot0, cuspidalfactor);

  // compute L/P again using Manin vector
  rational loverp_mvp(pdot, dp0 * cuspidalfactor);

  // Check they agree:
  if (pdot != dp0*pdot0)
    {
      cout << "Inconsistent values for L/P computed two ways!"<<endl;
      cout << "pdot  = " <<pdot <<endl;
      cout << "pdot0 = " <<pdot0 <<endl;
      cout << "dp0   = " <<dp0 <<endl;
      cout << "cuspidalfactor = "<<cuspidalfactor<<endl;
    }
#ifdef DEBUG_LoverP
  else
    {
      cout << "Consistent values for L/P computed two ways!"<<endl;
    }
  {
      cout << "from {0,oo} directly: " << loverp <<endl;
      cout << "pdot0 = "<<pdot0<<endl;
      cout << "cuspidalfactor = "<<cuspidalfactor<<endl;
      cout << "from Manin vector:    " << loverp_mvp <<endl;
      cout << "pdot = "<<pdot<<endl;
      cout << "nP0 = "<<nf->nP0<<endl;
      cout << "iP0 = "<<nf->iP0<<endl;
      cout << "eigs (size "<<eigs.size()<<") = ";
      vec_out(cout, eigs, 10);
      cout<<endl;
      cout << "ap0 = "<<nf->aP0[index-1]<<endl;
      cout << "dp0 = "<<dp0<<endl;
      cout << "cuspidalfactor = "<<cuspidalfactor<<endl;
    }
#endif
 }

// When a newform has been read from file, we have the aqlist and
// aplist but not the sequence of eigs in order.  This is needed
// both for recovering the basis vector from the h1 (in case we want
// to compute more ap), and for computing oldform multiplcities.

void newform::eigs_from_data()
  // recreate eigs list (in case we need to recover basis vector):
  // start with unramified char eigs (all +1), then Tp eigs ap for
  // good p
{
  // cout<<"In eigs_from_data (level "<<ideal_label(nf->N)<<"), aplist = "<<aplist<<endl;
  int ch(nf->characteristic);
  if (ch == 0)      // the first n2r eigs are all +1
    eigs.resize(nf->n2r, +1);
  else
    eigs.resize(0, +1);

  // Get a(P) or a(P^2) from the a(P) in aplist, for good P
  auto pr=Quadprimes::list.begin();
  auto api=aplist.begin();
  while (((int)eigs.size() < nf->nap+nf->n2r) && (api!=aplist.end()))
    {
      Quadprime P = *pr;
      INT normP = P.norm();
      while ((P.divides(nf->N)) || (ch>0 && (normP%ch==0)))
        {
          // cout<<" - P = "<<ideal_label(P)<<": bad prime, skipping"<<endl;
          ++pr;
          ++api;
          P = *pr;
          normP = P.norm();
        }
      long ap = *api;
      if (!P.has_square_class()) // then we need the eigenvalue of T(P^2)
        {
          ap = ap*ap - I2long(normP);
        }
      eigs.push_back(ap);
      //cout<<" - P = "<<ideal_label(P)<<": eig = "<<ap<<endl;
      ++pr; ++api;
      if (pr == Quadprimes::list.end())
        break;
    }
  // cout<<" eigs_from_data produced eigs = "<<eigs<<endl;
  fill_in_genus_class_data();
}

// When a newform has been read from file, when the class number is
// even,before computing more ap, we need to fill in the genus class
// data for each newform.
void newform::fill_in_genus_class_data()
{
  genus_classes.resize(1,0);
  genus_class_ideals.resize(1,Qideal(ONE));
  genus_class_aP.resize(1,1);

  // Now we fill the genus classes with known ideals and eigenvalues,
  // one in each genus class:
  int m2r = 0; // 2-rank of genus classes so far filled
  int iP = -1;  // index of old prime P begin looked at
  int nap = aplist.size();
  auto pr = Quadprimes::list.begin();
  while((pr!=Quadprimes::list.end()) && (m2r<nf->n2r) && (iP<nap-1))
    {
      Quadprime P = *pr++;
      iP++;
      if (P.divides(nf->N))
        continue;
      long aP = aplist[iP];
      if (aP==0)
        continue;
      long c = P.genus_class(1); // 1 means reduce mod Quad::class_group_2_rank
      if (c==0)
        continue;
      // if (nf->verbose)
      //   cout<<"form #"<<i<<" has eigenvalue "<<aP<<" and genus class "<<c<<endl;
      // See if we already have an eigenvalue for this genus class
      auto ci = std::find(genus_classes.begin(), genus_classes.end(), c);
      if (ci == genus_classes.end()) // then we do not
        {
          long oldsize = genus_classes.size();
          genus_classes.resize(2*oldsize);
          genus_class_ideals.resize(2*oldsize);
          genus_class_aP.resize(2*oldsize);
          for (int j = 0; j<oldsize; j++)
            {
              genus_classes[oldsize+j] = genus_classes[j]^c;
              genus_class_ideals[oldsize+j] = genus_class_ideals[j]*P;
              genus_class_aP[oldsize+j] = genus_class_aP[j]*aP;
            }
          m2r++;
        }
    } // loop on primes
  if (nf->verbose>1)
    {
      cout<<"Finished filling in genus class data for form #"<<index<<endl;
      cout<<"genus classes: "<<genus_classes<<endl;
      cout<<"genus class ideals: "<<genus_class_ideals<<endl;
      cout<<"genus class eigenvalues: "<<genus_class_aP<<endl;
    }
}

// For M a *multiple* of this level N, make the list of eigs
// appropriate for the higher level, taking into account the primes
// P (if any) dividing M but not N. For such P we delete the a(P)
// from the sublist of T(P) eigenvalues.

vector<long> newform::oldform_eigs(Qideal& M)
{
  assert (nf->N.divides(M));

  eigs_from_data();
  vector<long> M_eigs;

  if (nf->verbose)
    {
      cout<<"Making oldform eigs at level "<<ideal_label(M)<<" from eigs at level "<<ideal_label(nf->N)<<endl;
      cout<<" - input eigs: "<<eigs<<endl;
    }
  // insert eigs for central characters:
  if (nf->characteristic == 0)  // the first n2r eigs are all +1
    M_eigs.resize(nf->n2r, +1);

  auto ei = eigs.begin() + (nf->n2r);
  for ( const auto& P : Quadprimes::list)
    {
      if (ei==eigs.end()) break;
      if (!P.divides(nf->N))
        {
          if (!P.divides(M)) // else this T(P) eigenvalue is ignored
            {
              if (nf->verbose>1)
                cout << " keeping eigenvalue "<<(*ei)<< " at "<<P<<endl;
              M_eigs.push_back(*ei);
            }
          ++ei;
        }
    }
  if (nf->verbose)
    {
      cout<<" - output eigs: "<<M_eigs<<endl;
    }
  return M_eigs;
}

newform::newform(newforms* nfs, int ind,
                 const vector<int>& intdata, const vector<Quad>& Quaddata,
                 const vector<long>& aq, const vector<long>& ap)
{
  nf=nfs;
  index = ind;
  Qideal N(nf->N);
  int ch(nf->characteristic);

  if (ch == 0)
    {
      sfe = intdata[0];
      //  cout<<"sfe from file = "<<sfe<<endl;
      pdot = intdata[1];
      dp0 = intdata[2];
      cuspidalfactor = intdata[3];
      loverp = rational(abs(pdot),dp0*cuspidalfactor);
      lambda = Quaddata[0];
      lambdadot = intdata[4];
      a = Quaddata[1];
      b = Quaddata[2];
      c = Quaddata[3];
      d = Quaddata[4];
      matdot = intdata[5];
      bc = intdata[6];
      cm = intdata[7];
      CMD = intdata[8];
      aqlist = aq;
      // Recompute sign of functional equation = minus product of all A-L eigenvalues
      int newsfe = std::accumulate(aqlist.begin(),aqlist.end(),-1,std::multiplies<long>());
      if (newsfe!=sfe)
        cout<<"Problem in data on file for level "<<ideal_label(N)<<": sfe = "<<sfe<<" and aqlist = "<<aqlist<<", but minus product of latter is "<<newsfe<<endl;
    }

  aplist = ap;
  genus_classes.resize(1,0);
  genus_class_ideals.resize(1,Qideal(ONE));
  genus_class_aP.resize(1,1);
  genus_class_trivial_counter.resize(nf->nchi, 0);
  fake = 0;

  // recreate eigs list (used to recover basis vector, and for oldclasses at higher levels):

  eigs_from_data();
}

// find (a,b,c,d) such that cusp b/d is equivalent to 0 and the
// integral over {0,M(0)} = {0,b/d} with M = [a,b;N*c,d] is a
// nonzero multiple "matdot" of the period P.

void newform::find_matrix()
{
  if(nf->verbose>1)
    cout<<"computing integration matrix for newform "<<index<<"..."<<flush;
  matdot=0;
  Qideal N(nf->N);
  for (Quadlooper dl(2, 1000, 1); dl.ok()&&!matdot; ++dl)
    { d=(Quad)dl;
      Qideal D(d);
      if (N.is_coprime_to(D))
        {
          vector<Quad> reslist = residues(d);
          for( const auto& bi : reslist)
            {
              b = bi;
              Qideal bN = b*N;
              if (D.is_coprime_to(bN, a, c))
                // found a candidate q=b/d: a+c=1 with d|a and b|c and c/b in N
                {
                  c /= -b;
                  a /= d; // now a*d-b*c=1 with c in N
                  assert (a*d-b*c==Quad::one);
                  RatQuad q(b,d);
                  matdot = I2long(abs((nf->h1->chain(q, 1))[index]));
                  //cout << "Period from {0,"<<q<<"}, unscaled multiple "<<matdot<<endl;
                  if (matdot)
                    {
                      if (divides(cuspidalfactor, matdot))
                        matdot /= cuspidalfactor;
                      else
                        cout << "Error: unscaled matdot = " << matdot
                             << " is not divisible by cuspidalfactor " << cuspidalfactor << endl;
                      if (matdot)
                        {
                          //cout << "Using period from {0,"<<q<<"}, multiple "<<matdot<<endl;
                          break;
                        }
                    }
                } // b coprime to d test
            } // loop over b
        } // d coprime to N test
    } // loop over d
  if(nf->verbose>1)
    cout<<"M = ["<<a<<","<<b<<";"<<c<<","<<d<<"] with factor "<<matdot<<endl;
}

//#define DEBUG_BC

int newform::base_change_code(void)
{
  if (bc==4) // not yet set
    {
#ifdef DEBUG_BC
      cout<<"bc not set, computing code..."<<endl;
#endif
      bc = 0;
      if (is_base_change())
        {
#ifdef DEBUG_BC
          cout<<" - form is bc..."<<flush;
#endif
          bc = base_change_discriminant();
#ifdef DEBUG_BC
          cout<<" with disc "<<bc<<endl;
#endif
        }
      else if (is_base_change_twist())
        {
#ifdef DEBUG_BC
          cout<<" - form is twist of bc..."<<flush;
#endif
          bc = -base_change_twist_discriminant();
#ifdef DEBUG_BC
          cout<<" with disc "<<-bc<<endl;
#endif
        }
    }
#ifdef DEBUG_BC
  else
    {
      cout << "bc already set to "<<bc<<endl;
    }
#endif
  return bc;
}

int newform::is_base_change(void)
{
  if(!(nf->N.is_Galois_stable()))
    return 0;
  auto ap = aplist.begin();
  auto pr=Quadprimes::list.begin();
  Qideal N(nf->N);
  while(ap!=aplist.end() && pr!=Quadprimes::list.end())
    {
      long api = *ap++;
      Quadprime p0 = *pr++;
#ifdef DEBUG_BC
      cout<<"p="<<p0<<" has ap="<<api<<endl;
#endif
      if(!p0.is_Galois_stable()) // this prime not inert or ramified
        {
          if (ap==aplist.end()) // the conjugate ap is not known
            return 1;
          long apj = *ap++;
          Quadprime P1 =  *pr++;
          // skip if either divides level:
          if(p0.divides(N))
            continue;
          if(P1.divides(N))
            continue;
#ifdef DEBUG_BC
          cout<<"Next prime "<<P1<<" has ap="<<apj<<endl;
#endif
          if(api!=apj) // ap mismatch
            {
#ifdef DEBUG_BC
              cout<<"Mismatch -- not base-change"<<endl;
#endif
              return 0;
            }
        }
    }
#ifdef DEBUG_BC
  cout<<"All OK -- base-change"<<endl;
#endif
  return 1;
}

int newform::is_base_change_twist(void)
{
  auto ap = aplist.begin();
  auto pr=Quadprimes::list.begin();
  Qideal N(nf->N);
  while(ap!=aplist.end() && pr!=Quadprimes::list.end())
    {
      long api = *ap++;
      Quadprime p0 = *pr++;
      //cout<<"p="<<p0<<" has ap="<<api<<endl;
      if(!p0.is_Galois_stable()) // this prime not inert or ramified
        {
          if (ap==aplist.end()) // the conjugate ap is not known
            {
              //cout<<"All OK -- base-change up to twist"<<endl;
              return 1;
            }
            // read next (conjugate) prime and eigenvalue:
          long apj = *ap++;
          Quadprime P1 =  *pr++;
          // skip if either divides level:
          if(p0.divides(N))
            continue;
          if(P1.divides(N))
            continue;
          // Check the ap agree up to sign:
          if(abs(api)!=abs(apj)) // ap mismatch
            {
              //cout<<"Mismatch -- not base-change-twist"<<endl;
              return 0;
            }
        }
    }
  //cout<<"All OK -- base-change up to twist"<<endl;
  return 1;
}

// if form is base-change, find the d s.t. the bc has eigenvalues in Q(sqrt(d))
int newform::base_change_discriminant(void)
{
  if (is_base_change()==0) return 0;
  int bcd = 1;
  Qideal N(nf->N);
  auto api = aplist.begin();
  auto pr=Quadprimes::list.begin();
  while(api!=aplist.end() && pr!=Quadprimes::list.end())
    {
      long ap = *api++;
      Quadprime P = *pr++;
      if(!P.is_inert())
        continue;
      if(P.divides(N)) // this prime is bad
        continue;
      long dp = ap+2*P.prime();
      //cout<<"p="<<p<<" has ap="<<ap<<", disc = "<<dp;
      dp = squarefree_part(dp);
      //cout<<" with squarefree part "<<dp<<endl;
      if (dp==0) continue;
      if (bcd==1) // first one
        {
          bcd = dp;
        }
      else
        {
          if (dp!=bcd) // mismatch: not possible?
            {
              cout<<"\nWarning from base_change_discriminant(): bcd="<<bcd<<", dp="<<dp<<"; returning 0"<<endl;
              //cout<<"mismatch: bcd=0"<<endl;
              return 0;
            }
        }
    }
  return bcd;
}

// if form is twist of base-change, find the d s.t. the bc has eigenvalues in Q(sqrt(d)), if any
//
// Method: for each split P, d will be the squarefree part of either
// 2*p+a_P or 2*p-a_P, for some choice of sign depending on P.  So we
// return d if there is a choice of signs giving the same d for all
// split P, otherwise 0.

//#define DEBUG_BCTD

int newform::base_change_twist_discriminant(void)
{
  if (is_base_change_twist()==0) return 0;
  int cmd = is_CM();
  if (cmd!=0)
    {
      int bcd = -squarefree_part((Quad::d)/cmd);
#ifdef DEBUG_BCTD
      cout << "base change twist and CM("<<cmd<<") so bcd = "<<bcd<<endl;
#endif
      return bcd;
    }
  int bcd1 = 0, bcd2 = 0; // flags that they have not yet been set
  int n_candidates = 0;
  Qideal N(nf->N);
  auto api = aplist.begin();
  auto pr=Quadprimes::list.begin();
  Quadprime P;
  long ap;
  while(pr!=Quadprimes::list.end() && n_candidates !=1)
    {
      P = *pr++;
      int use_P = P.is_inert() && !P.divides(N);
      if (!use_P)
        {
          if (api!=aplist.end())
            {
              ++api; // must increment even if though this P is not being used
            }
          continue;
        }
      if (api!=aplist.end())
        {
          ap = *api++; // must increment even if this P is not being used
#ifdef DEBUG_BCTD
          cout<<" - using stored aP = "<<ap<<" for P = "<<P<<"..."<<endl;
#endif
        }
      else
        {
          nf->makebases(); // does nothing if already made
#ifdef DEBUG_BCTD
          cout<<" - computing aP for P = "<<P<<"..."<<flush;
#endif
          ap = eigenvalueHecke(P);
#ifdef DEBUG_BCTD
          cout<<" done,  aP = "<<ap<<" for P = "<<P<<"..."<<endl;
#endif
        }
      long dp1 = ap  + 2*P.prime();
      long dp2 = dp1 - 2*ap;
#ifdef DEBUG_BCTD
      cout<<"P="<<P<<" has aP="<<ap<<", discs "<<dp1<<", "<<dp2;
#endif
      dp1 = squarefree_part(dp1);
      dp2 = squarefree_part(dp2);
#ifdef DEBUG_BCTD
      cout<<" with squarefree parts "<<dp1<<", "<<dp2<<endl;
#endif
      if (dp1*dp2==0) continue;
      if (n_candidates==0)  // first pair, store
        {
          bcd1 = dp1;
          bcd2 = dp2;
          n_candidates = (bcd1==bcd2? 1: 2);
#ifdef DEBUG_BCTD
          cout<<" possible d: "<<bcd1<<", "<<bcd2<<"; "<<n_candidates<<" candidate(s)"<<endl;
#endif
        }
      else // see if only one is a repeat, if so it's the value we want
        {
          if ((bcd1!=-1)&&(dp1!=bcd1)&&(dp2!=bcd1)) // then bcd1 is bogus
            {
#ifdef DEBUG_BCTD
              cout<<" eliminating d="<<bcd1<<endl;
#endif
              bcd1 = -1;
              n_candidates = int(bcd2>0);
            }
          if ((bcd2!=-1)&&(dp1!=bcd2)&&(dp2!=bcd2)) // then bcd2 is bogus
            {
#ifdef DEBUG_BCTD
              cout<<" eliminating d="<<bcd2<<endl;
#endif
              bcd2 = -1;
              n_candidates = int(bcd1>0);
            }
          if ((bcd1==-1)&&(bcd2==-1)) // then we have no remaining candidates
            {
#ifndef DEBUG_BCTD
              if (nf->verbose)
#endif
                cout<<"\nbase_change_twist_discriminant() finds no suitable d!"<<endl;
              return 0;
            }
          // otherwise we still have at least one candidate, and continue
        }
    }
  // At this point, we have a valid d if and only if exactly one of
  // bcd1, bcd2 is positive (or both if they are equal)
  if ((bcd1>0)&&(bcd2>0)&&(bcd1==bcd2))
    {
#ifdef DEBUG_BCTD
      cout<<" base_change_twist_discriminant() returns d = "<<bcd1<<endl;
#endif
      return bcd1;
    }
  if ((bcd1>0)&&(bcd2<=0))
    {
#ifdef DEBUG_BCTD
      cout<<" base_change_twist_discriminant() returns d = "<<bcd1<<endl;
#endif
      return bcd1;
    }
  if ((bcd2>0)&&(bcd1<=0))
    {
#ifdef DEBUG_BCTD
      cout<<" base_change_twist_discriminant() returns d = "<<bcd2<<endl;
#endif
      return bcd2;
    }
  cout<<"\nWarning from base_change_twist_discriminant(): bcd1="<<bcd1<<", bcd2="<<bcd2<<"; returning 0"<<endl;
  return 0;
}

// Test if form is CM, return 0 or the CM disc

int newform::is_CM(void)
{
  if (cm==1) // not already set
    {
      auto api = aplist.begin();
      auto pr=Quadprimes::list.begin();
      while(api!=aplist.end() && pr!=Quadprimes::list.end())
        {
          long ap = *api++;
          Quadprime P = *pr++;
          if (ap==0) continue;
          long dp = ap*ap-4*I2long(P.norm());
          //cout<<"P="<<P<<" has aP="<<ap<<", disc = "<<dp;
          if (dp==0) continue;
          dp = squarefree_part(dp);
          //cout<<" with squarefree part "<<dp<<endl;
          if (cm==1) // first one
            {
              cm = dp;
              //cout << "Setting cm to "<<cm<<endl;
              continue;
            }
          if (dp!=cm) // mismatch: not CM
            {
              //cout<<"mismatch: CM=0"<<endl;
              cm = 0;
              break;
            }
        }
    }
  //cout << "Returning cm = "<<cm<<endl;
  return cm;
}

// Return this twisted by the genus character associated to D
newform newform::twist(const INT& D)
{
  newform f = *this; //copy constructor
  if (D==ONE)
    return f;

  //cout<<"Twisting newform by discriminant "<<D<<endl;
  //cout<<"aP before: "<<f.aplist<<endl;

  // Twist the ap:
  auto Pi = Quadprimes::list.begin();
  auto aPi=f.aplist.begin();
  for (;
       aPi!=f.aplist.end();
       ++Pi, ++aPi)
    (*aPi) *= Pi->genus_character(D);

  //cout<<"aP after: "<<f.aplist<<endl;

  // Twist the AL eigenvalues:
  auto Qi = nf->badprimes.begin();
  auto aQi=f.aqlist.begin();
  for (;
       aQi!=f.aqlist.end();
       ++Qi, ++aQi)
    (*aQi) *= Qi->genus_character(D);
  // Twist the sfe:
  f.sfe *= nf->N.genus_character(D);
  return f;
}

void newforms::makeh1(void)
{
  if(!h1)
    {
      h1 = new homspace(N, /*plus*/ 1, /*verbose*/ 0, characteristic);
      nfhmod=hmod = h1->h1hmod();
    }
}

newforms::newforms(const Qideal& iN, int disp, long ch)
  : N(iN), verbose(disp), n2r(Quad::class_group_2_rank), characteristic(ch), have_bases(0)
{
  modulus = default_modulus<scalar>();
  //  cout<<"In newforms constructor with level = "<<ideal_label(N)<<endl;
  nchi = 1<<n2r;
  level_is_square = N.is_square();

  // nulist is a list of n2r ideals coprime to N whose classes generate the 2-torsion
  if ((characteristic==0) && (n2r > 0))
    {
      nulist = make_nulist(N);
      if (verbose>1)
        cout<<"nulist: "<<nulist<<endl;
      possible_self_twists = N.possible_unramified_twists();
      if (verbose>1)
        cout<<"possible unramified self twist discriminants at this level: "<<possible_self_twists<<endl;
    }

  // badprimes is a list of all primes Q|N
  badprimes = N.factorization().sorted_primes();
  nwq = 0; // prevents any W_Q being used for splitting

  // goodprimes is a list of at least nap good primes (excluding those
  // dividing characteristic if >0), including at least one principal
  // one which has index iP0;

  nap = MAXDEPTH;
  goodprimes = make_goodprimes(N, nap, iP0, characteristic);
  nap = goodprimes.size(); // it may be > original nap
  if (nap!=MAXDEPTH && verbose)
    cout<<" nap changed to "<<nap<<" since goodprimes = "<<goodprimes<<endl;
  P0 = goodprimes[iP0];
  nP0 = I2long(P0.norm());

// P0 is the smallest good principal prime: and iP0 its index (in
// plist, which starts with the bad primes and then the good
// primes in order).  P0 must be principal since we have only
// implemented maninvector() for principal primes.

  if (verbose>1)
    {
      cout << "good primes used: "<<goodprimes<<endl;
    }

  h1=0;
  of=0;
  nfhmod=0;
}

// instantiations of virtual functions required by the splitter_base class:
mat newforms::opmat(int i, int dual, int verb)
{
  return h1->calcop(h1matop(i),dual,verb);
}

vec newforms::opmat_col(int i, int j, int verb)
{
  return h1->calcop_col(h1matop(i),j, verb);
}

mat newforms::opmat_cols(int i, const vec_i& jlist, int verb)
{
  return h1->calcop_cols(h1matop(i),jlist, verb);
}

mat newforms::opmat_restricted(int i, const subspace& s, int dual, int verb)
{
  return h1->calcop_restricted(h1matop(i),s,dual,verb);
}

smat newforms::s_opmat(int i, int dual, int)
{
  return h1->s_calcop(h1matop(i),0, dual, verbose);
}

smat newforms::s_opmat_cols(int i, const vec_i& jlist, int)
{
  return h1->s_calcop_cols(h1matop(i),jlist, verbose);
}

smat newforms::s_opmat_restricted(int i, const ssubspace& s, int dual, int)
{
  return h1->s_calcop_restricted(h1matop(i),s,dual,0);
}

//#define DEBUG_LAMBDA
void newforms::find_lambdas()
{
  vector<int> gotlambda(n1ds);
  int i, nfound=0;
#ifdef DEBUG_LAMBDA
  if(verbose) cout<<"Looking for twisting primes.\n";
#endif
  for (i=0; i<n1ds; i++)
    if(nflist[i].pdot!=0)
      {
#ifdef DEBUG_LAMBDA
        if(verbose) cout<<"Newform "<<i<<": lambda=1 will do.\n";
#endif
        nflist[i].lambda=Quad::one;
        nflist[i].lambdadot=nflist[i].pdot;
        gotlambda[i]=1;
        nfound++;
      }
    else
      {
        nflist[i].lambda=Quad::zero; // indicates no lambda exists (yet)
        nflist[i].lambdadot=0;
	gotlambda[i]=0;
      }
#ifdef DEBUG_LAMBDA
  if(verbose)cout<<nfound<<" easy cases out of "<<n1ds<<endl;
#endif
  if (level_is_square || !N.is_principal() || n2r>0)
    {
      return;
    }

  for( auto& L : Quadprimes::list)
    {
      if (nfound==n1ds) break;
      if (L.divides(TWO)) continue;
      if (L.divides(N)) continue;
      if (!L.is_principal()) continue;
      Quad lam = L.gen();
#ifdef DEBUG_LAMBDA
      if(verbose)cout << "Testing lambda = " << L << " = ("<<lam<<"): principal, odd, good"<<endl;
#endif
      vector<Quad> lamres = L.residues();
      if(squaremod(fundunit,lam,lamres)==1)
        {
#ifdef DEBUG_LAMBDA
          if(verbose)cout<<"passed second test: fundamental unit is a square"<<endl;
#endif
          // TODO:  work out what to do here if N is not principal!
          int chimod  = squaremod(N.gen(),lam,lamres);
          vector<int> chitab = makechitable(lam,lamres);
          vec mvtw = h1->manintwist(lam,lamres,chitab, 1);
          for(int j=0; (j<n1ds)&&(nfound<n1ds); j++)
            {
              if(gotlambda[j]==0)
                {
#ifdef DEBUG_LAMBDA
                  if(verbose)cout<<"Newform # "<<j<<": ";
#endif
                  newform& nfj = nflist[j];
                  int ldot = I2long(abs(mvtw[j+1]));  // j based at 0 but vec mvtw based at 1
                  if(ldot&&((chimod*nfj.sfe)==+1))
                    {
#ifdef DEBUG_LAMBDA
                      if(verbose)cout<<"Success! ";
#endif
                      nfj.loverp = rational(ldot, nfj.cuspidalfactor);
                      nfj.lambda = lam;
                      nfj.lambdadot = ldot;
                      gotlambda[j] = 1;
                      nfound++;
                    }
                }
#ifdef DEBUG_LAMBDA
              if(verbose)cout<<endl;
#endif
            }
        }
    }
}

void newforms::find()
{
  if(verbose)
    cout<<"Constructing homspace at level "<<ideal_label(N)<<" ...\n";
  makeh1();

  // fill the h1matops and eigranges lists with empties; the functions
  // h1matop(i) and eigrange(i) will compute and store these the first
  // time they are needed.
  int maxdepth(MAXDEPTH);
  for (int i=0; i<maxdepth; i++)
    {
      h1matops.push_back(matop());
      eigranges.push_back(vector<long>());
    }

  nfhmod=hmod = h1->h1hmod();
  int dimcusp = h1->h1cuspdim();
  int dimall = h1->h1dim();

  if(verbose)
    {
      cout<<"Dimension = "<<dimall<<endl;
      cout<<" - cuspidal dimension = "<<dimcusp<<endl;
    }

  int mindepth, olddim1, olddim2;

  // find oldform dimensions (all, rational, non-rational) and their split by characters:
  if(characteristic==0)
    {
      if (verbose)
        cout<<"Retrieving oldform data for level "<<ideal_label(N)<<"...\n";
      of = new oldforms(this);
      if (verbose)
        of->display();

      dimtrivcuspold = of->olddimall;
      olddim1 = of->olddim1;
      olddim2 = of->olddim2;

      old1dims = of->old1dims;
      old2dims = of->old2dims;
      olddims = of->olddims;

      mindepth = n2r+3;  // if too small we may get fake rational newforms (eg d=22, level 121.1)
    }
  else
    {
      dimtrivcuspold=0;
      olddim1 = olddim2 = 0;
      mindepth = nap;
    }

  // find dimension of trivial character subspace and its split by characters:
  alldims = h1->trivial_character_subspace_dimensions_by_twist(1, 1, olddims);
  dimtrivcusp = std::accumulate(alldims.begin(), alldims.end(), 0, std::plus<int>());

  if(verbose && (n2r>0))
    {
      cout<<" - cuspidal trivial character subspace dimension = "<<dimtrivcusp;
      cout<<" = "<<alldims<< " split by self-twist character"<<endl;
    }

  // deduce newform dimensions and their split by characters (but we not yet know how much is rational):
  dimtrivcuspnew = dimtrivcusp - dimtrivcuspold; // total new dimension at this level
  assert (dimtrivcuspnew>=0 && "new dimension cannot be negative");
  newdims.resize(nchi);
  for (int i=0; i<nchi; i++)
    {
      newdims[i] = alldims[i] - olddims[i];
      // cout<<"i = "<<i<<", alldims[i]="<<alldims[i]<<", olddims[i]="<<olddims[i]<<endl;
      assert (newdims[i]>=0 && "components of new dimensions cannot be negative");
    }

  if(verbose && (characteristic==0))
    {
      // cout<<" - cuspidal trivial character subspace dimension = "<<dimtrivcusp;
      // cout<<" = "<<alldims<< " split by self-twist character"<<endl;
      cout<<"Trivial character cuspidal subspace:"<<endl;
      cout<<" - total dimension = "<<dimtrivcusp;
      if (n2r>0)
        cout<<" = "<< alldims<<" split by self-twist character";
      cout << endl;
      cout<<" - old dimension   = "<<dimtrivcuspold
          <<" (rational:"<<olddim1<<"; non-rational: "<<olddim2<<")";
      if (n2r>0)
        cout<<" = "<< olddims <<" ("<<old1dims<<"; "<<old2dims
            <<") split by self-twist character"<<endl;
      cout<<" - new dimension   = "<<dimtrivcuspnew;
      if (n2r>0)
        {
          cout<<" = "<<newdims<<" split by self-twist character";
          cout<<"\n unramified character discriminants: " << Quad::all_disc_factors;
        }
      cout<<endl;
      if(verbose>1)
        {
          cout<<"maxdepth = "<<maxdepth<<endl;
          cout<<"mindepth = "<<mindepth<<endl;
        }
    }

  n1ds = 0; // number of rational newforms found (will be incremented by the finder)
  if(dimtrivcuspnew>0)  // Else no newforms certainly so do no work!
    {
      if(verbose)
        cout<<"Finding rational newforms...\n";
      use_nf_number=-1; // flags to use() that the nfs found are new
      form_finder ff(this,modulus,1,maxdepth,mindepth,1,0,verbose);
      ff.find();
     }
  n2ds=dimtrivcuspnew-n1ds; // dimension of new, non-rational forms
  if(verbose>1)
    {
      cout<<"provisional n1ds = "<<n1ds<<endl;
      cout<<"provisional n2ds = "<<n2ds<<endl;
    }
  if (n2ds<0 && verbose)
    {
      cout << "found more possibly rational newforms than the total new dimension!\n";
      cout << "This means that at least "<<(-n2ds)<< " must be fake rationals"<<endl;
    }

  // We cannot yet split n1ds, and hence n2ds, by character since we
  // cannot detect self-twist until we have computed more ap -- unless
  // n1ds=0.  [It would be possible to test whether the basis vector
  // for each newform lies in the relevant subspaces intead.] The
  // splitting of n1ds is done in the getap() method.

  if (n1ds==0)
    {
      new1dims.resize(nchi, 0);
      new2dims = newdims;
    }

  if(verbose)
    {
      cout << "Total new cuspidal dimension " << dimtrivcuspnew << " made up as follows:\n";
      cout << "Rational: "<<n1ds<<"; non-rational: "<<n2ds<<endl;
      if (n1ds>0 && n2r>0)
        cout << " (this number of rational newforms might include fake rationals)" <<endl;
      // if (n2r>0)
      //   {
      //     cout << " rational non-self-twist: "<<new1dims[0]<<endl;
      //     for (int j=1; j<nchi; j++)
      //       {
      //         if (new1dims[j]>0)
      //           cout << " rational, self-twist by "<<Quad::all_disc_factors[j]<<": "<<new1dims[j]<<endl;
      //       }
      //     cout << " non-rational non-self-twist: "<<new2dims[0]<<endl;
      //     for (int j=1; j<nchi; j++)
      //       {
      //         if (new2dims[j]>0)
      //           cout << " non-rational, self-twist by "<<Quad::all_disc_factors[j]<<": "<<new2dims[j]<<endl;
      //       }
      //   }
    }

  // At this point the newforms in the list only have eigs and long basis vector
  fill_in_newform_data();
  nap=0;
  have_bases=1;
}

// fill in extra data in each newform:
void newforms::fill_in_newform_data(int AL, int CF, int LP, int M)
{
  if(n1ds==0) return; // no work to do
  make_projcoord();    // Compute homspace::projcoord and homspace::bigtkernbas
  //cout<<" - projcoord computed:\n"<<h1->projcoord<<endl;
  if (Quad::class_number == 1)
    {
      make_bigtkernbas();  //  before filling in newform data
      //cout<<" - bigtkernbas computed:\n"<<h1->bigtkernbas<<endl;
    }
  find_jlist();
  vec zero_infinity_coords = h1->chaincd(Quad::zero, Quad::one, 0, 0);
  zero_infinity = h1->chaincd(Quad::zero, Quad::one, 0, 1); // last 1 means use projcoord
  mvp = h1->maninvector(P0, 1);              // last 1 means use projcoord
#ifdef DEBUG_LoverP
  cout << "coords of {0,oo} is " << zero_infinity_coords << endl;
  cout << "proj   of {0,oo} is " << zero_infinity << endl;
  cout << "proj   of mvp    is " << mvp << endl;
#endif
  aP0 = apvec(P0);                         // vector of ap for first good principal prime
  if (verbose>1) cout << "found eigenvalues for P0="<<P0<<": "<<aP0<<endl;
  // Fill in data for each newform.
  for (int j=0; j<n1ds; j++)
    {
      newform& nfj = nflist[j];
      nfj.index = j+1;

      // compute A-L eigenvalues now in odd class number, else they are
      // computed in getap()
      if (AL && Quad::class_group_2_rank==0)
        nfj.compute_AL();

      if (characteristic>0)
        return;

      // compute cusidalfactor
      if (CF)
        {
          nfj.compute_cuspidalfactor();

          // compute L/P ratio (uses cuspidalfactor)
          if (LP)
            nfj.compute_loverp();

          // find M=[a,b;c,d] in Gamma_0(N) (i.e. N|c, det=1), so cusp b/d is
          // equivalent to 0 and the integral over {0,M(0)} = {0,b/d} is a
          // nonzero multiple "matdot" of the base period P0. (uses cuspidalfactor)
          if (M)
            nfj.find_matrix();
        }
    }

  // Find the twisting primes for each newform (more efficient to do
  // this here instead of within the newform constructors, as one lambda
  // might work for more than one newform). NB If the level is square
  // SFE=-1 then no such lambda will exist.
  if (characteristic==0)
    find_lambdas();
}

// Compute a long basis vector from a short one
vec newforms::lengthen_basis(const vec& sbasis)
{
  // convert short basis vector (coords w.r.t. face-gens) to a long one (coords w.r.t.edge-gens):
  if(hmod!=0)
    { // we don't have a mod p mat*vec
      mat vcol(dim(sbasis),1);
      vcol.setcol(1,sbasis);
      return reduce_mod_p(matmulmodp(h1->FR.coord, vcol, hmod).col(1), default_modulus<scalar>());
    }
  vec basis = (h1->FR.coord)*sbasis;
  if (characteristic==0)
    basis /= content(basis);  // basis is now independent of h1's denom1
  return basis;
}

// Create or update a newform from a short basis vector and
// eigenvalues.  The second vec argument is not used.  If
// use_nf_number is -1 we have found a new eigenvector, and create a
// newform from its short basis vector and eigenvalues; otherwise we
// just update the newform's basis vector.

void newforms::use(const vec& b1, const vec&, const vector<long> eigs)
{
  if (use_nf_number==-1)
    {
      //cout<<"Constructing newform with eigs "<<eigs <<" and short basis vector " << b1<<endl;
      // The constructor sets the basis field from the given short basis vector
      nflist.push_back(newform(this,b1,eigs));
      n1ds++;
      if (n1ds>dimtrivcuspnew)
        {
          cout << "*** Warning: in splitting eigenspaces (level "<<ideal_label(N)<<"): apparently found more ";
          cout << "1D rational newforms ("<< n1ds
               <<", possibly including fake rationals) than the total new cuspidal dimension ("
               <<dimtrivcusp<<") ***"<<endl;
          cout<<"Extra newform has eigs "<<eigs<<endl;
        }
    }
  else // store eigs and basis
    {
      // Update the newform from the short basis vector and eigenvalues
      nflist[use_nf_number].eigs = eigs;
      nflist[use_nf_number].basis = lengthen_basis(b1);
    }
}

void newforms::display(int detail)
{
 if (n1ds==0) {cout<<"No newforms."<<endl; return;}
 cout << "\n"<<n1ds<<" newform(s) at level " << ideal_label(N) << " = " << gens_string(N);
 if (n2r>0)
   cout << " (up to twist by unramified character)";
 cout  << ":"<< endl;
 if (detail)
   for(int i=0; i<n1ds; i++)
     {
       cout<<i+1<<":\t";
       nflist[i].display();
     }
}

void newform::display(void) const
{
  cout << "basis = " << basis;
  if (nf->characteristic==0)
    cout << ";\taqlist = " << aqlist;
  cout << ";\taplist = " << aplist << endl;
  if (nf->characteristic==0)
    {
      cout << "Sign of F.E. = " << sfe << endl;
      // stop outputting these as we do not use them
      // cout << "Twisting prime lambda = " << lambda << ", factor = " << lambdadot << endl;
#ifdef DEBUG_LoverP
      cout << "L/P ratio    = " << loverp << ", cuspidal factor = " << cuspidalfactor << endl;
      cout << "Integration matrix = [" << a << "," << b << ";" << c << "," << d << "], factor   = " << matdot << endl;
#endif
      if (CMD!=0)
        cout << "Unramified self-twist by discriminant "<<CMD<<endl;
    }
}

void newforms::list(long nap)
{
  //  string idlabel = (Quad::class_number==1? ideal_code(N.gen()): ideal_label(N));
  string idlabel = ideal_label(N);
  string idgens = gens_string(N), flabel = field_label();
  string s1 = flabel + " " + idlabel + " ";
  string s2 = " " + idgens + " 2 ";
  for(int i=0; i<n1ds; i++) // make sure these are set *before* making the unramified twists
    {
      nflist[i].base_change_code();
      nflist[i].is_CM();
    }
  if (n2r==0)
    for(int i=0; i<n1ds; i++)
      nflist[i].list(s1 + codeletter(i) + s2, nap);
  else
    {
      vector<newform> twisted_newforms;
      for(int i=0; i<n1ds; i++)
        {
          INT D = nflist[i].CMD;
          vector<INT> twists = disc_factors_mod_D((D==ZERO?ONE:D));
          int ntwists = twists.size();
          // if(D!=0)
          //   cout<<i<<": "<<ntwists<<" twists (D="<<D<<"), by "<<twists<<endl;
          for (int j = 0; j<ntwists; j++)
            twisted_newforms.push_back(nflist[i].twist(twists[j]));
        }
      ::sort(twisted_newforms.begin(), twisted_newforms.end(), less_newform_lmfdb);
      int nform = 0;
      for (auto nf=twisted_newforms.begin(); nf!=twisted_newforms.end(); ++nf)
        nf->list(s1 + codeletter(nform++) + s2, nap);
    }
}

void newform::list(string prefix, long nap)
{
  if(nap==-1) nap=aplist.size();

  cout << prefix;

  if (nf->characteristic>0)
    cout << nf->characteristic << " ";
  else
    {
      cout << base_change_code() << " ";
      cout << is_CM() << " ";
    }

  if (nf->characteristic==0)
    {
      cout << sfe << " " << loverp << " ";
      cout << "[";
      int first = 1;
      for( const auto& aq : aqlist)
        {
          if (!first) cout<<",";
          cout<<aq;
          first = 0;
        }
      cout << "] ";
    }
  // The x here is essentially a place-holder representing the Hecke
  // field defining polynomial so that the output here will be
  // consistent with newforms whose Hecke field has degree >1
  cout << "x ";
  cout << "[";
  int n = 0;
  for ( const auto& ap : aplist)
    {
      if (n==nap) break;
      if (n>0) cout<<",";
      if ((nf->characteristic>0) && (ap<0)) // we use -999 for omitted eigs
        cout << "*";
      else
        cout << ap;
      n++;
    }
  cout <<"]" <<endl;
}

// compute the eigenvalue for a single operator on this newform
long newform::eigenvalue(const matop& op, pair<long,long> apbounds, long factor)
{
  vec image = nf->h1->applyop(op, m0, 1);
  int top = I2long(image[index]); // where this is the i'th newform
  long ap;
  // The eigenvalue is now top/fac (which should divide exactly)
  int nfhmod = I2long(nf->nfhmod);
  if(nfhmod!=0)
    ap=mod(xmodmul(top,facinv,nfhmod), nfhmod);
  else
    {
      if (top%fac !=0)
        {
          cout<<"Problem in eigenvalue: image = "<<image<< " has "<< index <<" entry "
              <<top<<" which is not divisible by pivot "<<fac<<endl;
          cout<<flush;
        }
      ap = top/fac;
    }
  if (nf->characteristic>0)
    ap = posmod(ap, nf->characteristic);
  else // check it is in range (in characteristic 0 only)
    {
      long absfac = abs(factor);
      if (divides(factor,ap))
        {
          ap /= factor;
          if ((ap<absfac*apbounds.first) || (ap>absfac*apbounds.second))
            {
              cout<<"Error:  eigenvalue "<<ap<<" for operator "<<op.name()
                  <<" for form # "<< index <<" is outside valid range "
                  << factor << "*" << apbounds.first<<"..."<< factor << "*" <<apbounds.second<<endl;
              exit(1);
            }
        }
      else
        {
          cout<<"Error:  eigenvalue "<<ap<<" for operator "<<op.name()
              <<" for form # "<< index <<" is not a multiple of "
              << factor << endl;
          exit(1);
        }
    }
  return ap;
}

long newform::eigenvalueHecke(Quadprime& P, int verbose)
{
  if (fake)
    return 0; // No need to compute more for a fake rational
  Qideal N = nf->N;
  if (P.divides(N))
    return 0; // we'll deal with A-L eigs later
  long c = P.genus_class(1);  // 1 means reduce mod Quad::class_group_2_rank
  //cout << "P=" <<P<<" has genus character "<<P.genus_character()<<" and genus class "<<c<<endl;
  if (c==0) // then P has square class, compute aP directly via T(P) or T(P)*T(A,A)
    {
      if (P.is_principal())  // compute T(P)
        {
          if (verbose>1)
            cout << "form "<<index<<", computing T("<<P<<") directly as "<<P<<" is principal"<<endl;
          return eigenvalue(HeckePOp(P, N), eigenvalue_range(P));
        }
      else   // P has square class, compute T(P)*T(A,A)
        {
          Qideal A = P.sqrt_coprime_to(N);
          if (verbose>1)
            cout << "form "<<index<<", computing T("<<P<<") using T(P)*T(A,A) with A = "<<ideal_label(A)<<endl;
          return eigenvalue(HeckePChiOp(P,A,N), eigenvalue_range(P));
        }
    }
  else // P does not have square class, its genus class is c>0
    {
      // See if we already have an eigenvalue for this genus class
      if (verbose>1)
        cout<<"P="<<P<<" has genus class "<<c<<", genus_classes covered so far: "<<genus_classes<<endl;
      auto ci = std::find(genus_classes.begin(), genus_classes.end(), c);
      if (ci != genus_classes.end()) // then we do
        {
          int i = ci - genus_classes.begin();
          long factor = genus_class_aP[i];
          Qideal A = genus_class_ideals[i];
          Qideal B = A*P; // so B is square-free and of square class
          if (verbose>1)
            cout << "form "<<index<<", computing T("<<P<<") using T("<<ideal_label(B)<<") = T(P)*T(A) with A = "
                 <<ideal_label(A)
                 <<", with T(A) eigenvalue "<<factor<<endl;
          long aP;
          if (B.is_principal())               // compute T(B)
            aP = eigenvalue(HeckeBOp(B, N), eigenvalue_range(P), factor);
          else                                // compute T(B)*T(A,A)
            {
              Qideal C = B.sqrt_coprime_to(N);
              aP = eigenvalue(HeckeBChiOp(B,C,N), eigenvalue_range(P), factor);
            }
          // See whether P is a better genus class rep than the one we have:
          if ((P.norm()<A.norm()) && (aP!=0))
            {
              genus_class_ideals[i] = P;
              genus_class_aP[i] = aP;
            }
          return aP;
        }
      else // we have a new genus class, compute a_{P}^2 unless we already have at least 5 zeros in this class
           // and the level admits nontrivial self twists.
           // NB 5 would not be enough for field 299, level 100.2, without checking for possible self twists.
        {
          if (verbose>1)
            cout << "P=" <<P<<" has genus class "<<c<<", genus_class_trivial_counter = "<<genus_class_trivial_counter<<endl;
          if ((possible_self_twists.size()>0) && (genus_class_trivial_counter[c] >= genus_class_triviality_bound))
            {
              if (verbose>0)
                cout << "form "<<index<<", P = " <<P<<": genus class "<<c<<" has "<<genus_class_trivial_counter[c]
                     <<" zero eigenvalues, so assuming self-twist, and taking aP=0"<<endl;
              return 0;
            }
          if (verbose>1)
            cout << "form "<<index<<", P = "<<P<<": computing T(P^2) to get a(P^2) and hence a(P)^2" << endl;
          Qideal P2 = P*P;
          long aP, aP2;
          long normP = I2long(P.norm());
          if (P2.is_principal())  // compute T(P^2)
            {
              aP2 = eigenvalue(HeckeP2Op(P,N), eigenvalue_sq_range(P));
            }
          else // T(P^2)*T(A,A) with (A*P)^2 principal
            {
              Qideal A = P.equivalent_mod_2_coprime_to(N, 1);
              aP2 = eigenvalue(HeckeP2ChiOp(P,A,N), eigenvalue_sq_range(P));
            }
          // Now aP2 is the eigenvalue of T(P^2)
          aP2 += normP;
          // Now aP2 is the eigenvalue of T(P)^2
          if (verbose>1)
            cout << " - a(P)^2 = " << aP2 << endl;
          if (isqrt(aP2, aP))
            {
              if (aP!=0) // else we cannot use this as a new genus pivot
                {
                  if (verbose>1)
                    cout << " - taking a(P) = " << aP << endl;
                  // update genus_classes (append binary sum of each and c)
                  // update genus_class_ideals (append product of each and P)
                  // update genus_class_aP (append product of each and aP)
                  long oldsize = genus_classes.size();
                  if (verbose>1)
                    cout<<" -- form "<<index
                        <<": doubling number of genus classes covered by P with nonzero aP from "<<oldsize
                        <<" to "<<2*oldsize<<" using P="<<P<<endl;
                  genus_classes.resize(2*oldsize);
                  genus_class_ideals.resize(2*oldsize);
                  genus_class_aP.resize(2*oldsize);
                  genus_class_trivial_counter[c] = 0;
                  for (int i = 0; i<oldsize; i++)
                    {
                      genus_classes[oldsize+i] = genus_classes[i]^c;
                      genus_class_ideals[oldsize+i] = genus_class_ideals[i]*P;
                      genus_class_aP[oldsize+i] = genus_class_aP[i]*aP;
                    }
                  // we can possibly eliminate some of
                  // possible_self_twists now, namely those whose
                  // characters chi have chi(P)=-1, since such a
                  // newform would have to have aP=0:
                  int n_before = possible_self_twists.size();
                  if (0)//(n_before)
                    {
                      possible_self_twists.erase(std::remove_if(possible_self_twists.begin(),
                                                                possible_self_twists.end(),
                                                                [&P](INT D) { return P.genus_character(D) == -1;}),
                                                 possible_self_twists.end());
                      int n_after = possible_self_twists.size();
                      if ((n_before>n_after) && (verbose>1))
                        cout<<" - after erasing "<<(n_before-n_after)
                            <<" possible self-twist discriminants, these remain: "<<possible_self_twists << endl;
                    }
                  // We can now eliminate any self-twist discriminants
                  // not matching the square-free part of aP^2-4*N(P):
                  int d1 = squarefree_part(aP2-4*normP);
                  possible_self_twists.erase(std::remove_if(possible_self_twists.begin(),
                                                            possible_self_twists.end(),
                                                            [&d1](INT D) { return D!=d1;}),
                                             possible_self_twists.end());
                  int n_after = possible_self_twists.size();
                  if ((n_before>n_after) && (verbose>1))
                    cout<<" - after erasing "<<(n_before-n_after)
                        <<" possible self-twist discriminants, these remain: "<<possible_self_twists << endl;
                }
              else
                {
                  genus_class_trivial_counter[c] +=1;
                  if (verbose>1)
                    cout << "form "<<index<<", genus_class_trivial_counter for class "<<c
                        <<" is now "<<genus_class_trivial_counter[c]<<endl;
                }
            }
          else
            {
              if (verbose)
                {
                  cout<<" -- form "<<index<<": P="<<P<<", a(P)^2 = "<<aP2<<" is not a square!";
                  cout<<" -- tagging this as a fake rational, to be discarded"<<endl;
                }
              fake = 1;
            }
          return aP;
        }
    }
}

long newform::eigenvalueAtkinLehner(Quadprime& Q, int verbose)
{
  if (fake) return 0;
  Qideal N = nf->N;
  int e = val(Q,N);
  if (e==0)
    return 1;
  Qideal Qe = Q;
  while (--e) Qe*=Q;
  if (Qe.is_principal())
    {
      if (verbose)
        cout << "computing W("<<Q<<") directly"<<endl;
      return eigenvalue(AtkinLehnerQOp(Q,N), {-1,1});
    }
  if (Qe.has_square_class())
    {
      Qideal A = Qe.sqrt_coprime_to(N);
      if (verbose)
        cout << "computing W("<<Q<<") using W(Q)*T(A,A) with A = "<<ideal_label(A)<<endl;
      return eigenvalue(AtkinLehnerQChiOp(Q,A,N), {-1,1});
    }

  // Now we rely on already having aP for some P in the same genus
  // class as Q^e, and compute W(Q)T(P)

  // So far we have only implemented W(Q)T(P) when Q^e*P is
  // principal, not when it only has square class, so we must find a
  // good prime in the opposite ideal class to Q^e for which we know
  // aP (and it is nonzero).

  if (verbose)
    cout << "computing W("<<Q<<") using W(A)*T(P) for suitable P (opposite class to Q^e="
         <<ideal_label(Qe) <<", aP nonzero)"<<endl;
  for (auto Pi = Quadprimes::list.begin(); Pi!=Quadprimes::list.end(); ++Pi)
    {
      Quadprime P = *Pi;
      if (P.divides(N))
        continue;
      if (!Qe.is_anti_equivalent(P))
        continue;
      if (verbose)
        cout << " - trying P = "<<P<<" which is in the right ideal class"<<endl;
      vector<long int>::size_type i = Pi-Quadprimes::list.begin();
      long aP;
      if (i<aplist.size())
        {
          aP = aplist[i];
          if (verbose)
            cout << " - stored a(P) = "<<aP<<endl;
        }
      else
        {
          aP = eigenvalueHecke(P, verbose);
          if (verbose)
            cout << " - computed a(P) = "<<aP<<endl;
        }
      if (aP==0)
        {
          if (verbose)
            cout << " - no good as aP = "<<aP<<endl;
          continue;
        }
      if (verbose)
        cout << " - OK, aP = "<<aP<<endl;
      long eQ = eigenvalue(HeckePALQOp(P, Q, N), {-1,+1}, aP);
      if (verbose)
        cout << " - eigenvalue of T(P)W(Q) is "<<aP*eQ<<" so A-L eigenvalue = "<<eQ<<endl;
      return eQ;
    }
  cout << "*** Warning: unable to find a suitable auxiliary prime P for which "
       << "the eigenvalue of T(P)W(Q) is known and nonzero!  Returning A-L eigenvalue of 0"<<endl;
  return 0;
}

void newforms::sort_eigs(void)
{
  ::sort(nflist.begin(),nflist.end(),less_newform_eigs);
}

void newforms::sort_lmfdb(void)
{
  ::sort(nflist.begin(),nflist.end(),less_newform_lmfdb);
}

// Replaces matrix "coord" in homspace with "projcoord"
//
// Each has ngens rows (ngens = number of edges modulo
// edge-relations).  coord has rk columns, and the i'th row of coord
// gives the coordinates of the i'th generating edge with respect to
// the basis modulo face-relations, these begin implicitly scaled by
// denom1.  projcoord has n1ds columns, and the i'th row gives the
// coords with respect to the (partial) basis of eigenforms; its
// columns are primitive and uniquely determined (up to sign),
// independent of the homology basis, and in particular of denom1.

void newforms::make_projcoord()
{
  h1->projcoord.init(h1->ngens,n1ds);
  for (int j=1; j<=n1ds; j++)
    {
      h1->projcoord.setcol(j, nflist[j-1].basis);
    }
#ifdef DEBUG_LoverP
  cout << "projcoord (transpose):\n" << transpose(h1->projcoord) << endl;
  cout << "basis of ker(delta) (rows):\n" << h1->tkernbas.as_mat() <<endl;
#endif
}

// Set bigtkernbas member of homspace

// The rows of h1->tkernbas are a basis for ker(delta) with respect to
// the freegens basis for homology (rel cusps), but this is not a
// basis for the integral homology when h1->denom>1.  Instead we need
// a basis with respect to the gens.

// We only call this when the class number is 1
void newforms::make_bigtkernbas(void)
{
  // 'big' version bigtkernbas
  scalar modulus = default_modulus<scalar>();
  if (characteristic) modulus = scalar(characteristic);
  mat tcoord = transpose(h1->FR.get_coord());
  //cout<<" *** computed tcoord: "<<tcoord<<endl;
  smat bigdeltamat = mult_mod_p(h1->sdeltamat, smat(tcoord), modulus);
  //cout<<" *** computed bigdeltamat: "<<bigdeltamat<<endl;
  h1->bigtkernbas = transpose(basis(kernel(bigdeltamat, modulus)));
}

// try to read from file, and if no data file exists, finds from scratch and stores
void newforms::read_from_file_or_find()
{
  if (verbose>1)
    cout << " - reading newform data for level "<<ideal_label(N)<<endl;
  int ok = read_from_file();
  if (ok)
    {
      if (verbose>1)
        cout << " - successfully read newform data for level "<<ideal_label(N)<<endl;
      return;
    }
  if (verbose)
    cout << " - no newform data for level "<<ideal_label(N)<<" exists, finding newforms..."<<endl;
  find();
  if (verbose)
    cout << " - found "<<n1ds<<" newforms for level "<<ideal_label(N)<<endl;
  if (n1ds>0)
    {
      if (verbose)
        cout << " - computing eigenvalues numbers 1 to "<<max(nap,25)<<"... "<<endl;
      getap(1, max(nap,25), 0);
    }
  string eigfilename = (Quad::class_number==1? eigfile(N.gen(), characteristic): eigfile(N, characteristic));
  output_to_file(eigfilename);
  if(verbose)
    {
      cout << "  finished creating and storing newforms at level " << N << endl;
      if (verbose>1)
        display();
    }
}

int newforms::read_from_file()
{
  if(verbose)
    cout << "Retrieving newform data for N = " << ideal_label(N) << endl;

// Read newform data from file

  if(verbose>1) cout << "Getting newform data for " << N << endl;
  string eigfilename = (Quad::class_number==1? eigfile(N.gen(), characteristic): eigfile(N, characteristic));
  ifstream data(eigfilename.c_str());
  if (!data)
    {
      if(verbose)
        {
          cout << "No data file for level " << ideal_label(N);
          if (characteristic) cout << " mod " << characteristic;
          cout << endl;
        }
      return 0;
    }
  data >> n1ds >> n2ds;
  new1dims.resize(nchi);
  new2dims.resize(nchi);
  if (n2r>0)
    {
      for (int i=0; i<nchi; i++)
        data >> new1dims[i];
      for (int i=0; i<nchi; i++)
        data >> new2dims[i];
    }
  else
    {
      new1dims[0] = n1ds;
      new2dims[0] = n2ds;
    }
  data >> nap;
  if(verbose>1)
    {
      cout<<" read data for "<<n1ds<<" newforms at level "<<N
          <<", total new dimension = "<<(n1ds+n2ds)<<", nap = "<<nap<<endl;
    }
  if (n1ds==0)
    return 1;

  vector<vector<long> > aqs(n1ds), aps(n1ds), eigs(n1ds);
  vector<vector<int> > intdata(n1ds);   // sfe, pdot, dp0, cuspidalfactor, lambdadot, matdot
  vector<vector<Quad> > Quaddata(n1ds); // lambda, a, b, c, d

  int i;
  for(i=0; i<n1ds; i++)
    {
      eigs[i].resize(nap);
      aps[i].resize(nap);
      if (characteristic==0)
        {
          aqs[i].resize(badprimes.size());
          intdata[i].resize(9);
          Quaddata[i].resize(5);
        }
    }

  if (characteristic==0)
    {
      // Read the auxiliary data (unless in positive characteristic):
      for (i=0; i<n1ds; i++) data>>intdata[i][0];  // sfe
      for (i=0; i<n1ds; i++) data>>intdata[i][1];  // pdot
      for (i=0; i<n1ds; i++) data>>intdata[i][2];  // dp0
      for (i=0; i<n1ds; i++) data>>intdata[i][3];  // cuspidalfactor
      for (i=0; i<n1ds; i++) data>>Quaddata[i][0]; // lambda
      for (i=0; i<n1ds; i++) data>>intdata[i][4];  // lambdadot
      for (i=0; i<n1ds; i++) data>>Quaddata[i][1]; // a
      for (i=0; i<n1ds; i++) data>>Quaddata[i][2]; // b
      for (i=0; i<n1ds; i++) data>>Quaddata[i][3]; // c
      for (i=0; i<n1ds; i++) data>>Quaddata[i][4]; // d
      for (i=0; i<n1ds; i++) data>>intdata[i][5];  // matdot
      int iCMD=6;
      for (i=0; i<n1ds; i++) data>>intdata[i][6];  // bc
      //for (i=0; i<n1ds; i++) intdata[i][6]=4;  // temp fix to to recompute bc
      for (i=0; i<n1ds; i++) data>>intdata[i][7];  // cm
      iCMD=8;
      if (n2r>0)
        for (i=0; i<n1ds; i++) data>>intdata[i][iCMD];  // CMD
      else
        for (i=0; i<n1ds; i++) intdata[i][iCMD]=0;  // not used

      //  Read the W-eigenvalues at level M into aqs:
      for(i=0; i<(int)badprimes.size(); i++)
        for( auto& f : aqs)
          data >> f[i];
    }

  // Next read the coefficients at level M into aps:
  for(i=0; i<nap; i++)
    for( auto& f : aps)
      data >> f[i];

  data.close();

  if(verbose>1)
    {
      cout << "Finished reading newform data for level " << N << endl;
      if (characteristic==0)
        {
          cout << "aqs = " << endl;
          for(i=0; i<n1ds; i++) cout<<i<<": "<<aqs[i]<<endl;
        }
      cout << "aps = " << endl;
      for(i=0; i<n1ds; i++) cout<<i<<": "<<aps[i]<<endl;
      if (characteristic==0)
        {
          cout << "intdata = " << endl;
          for(i=0; i<n1ds; i++) cout<<i<<": "<<intdata[i]<<endl;
          cout << "Quaddata = " << endl;
          for(i=0; i<n1ds; i++) cout<<i<<": "<<Quaddata[i]<<endl;
        }
    }

// Extract number of newforms and their eigenvalues from this.

  if(verbose>1) cout << " read "<<n1ds << " newforms for N = " << ideal_label(N) << endl;

 // construct the newforms from this data
  for(i=0; i<n1ds; i++)
    {
      if(verbose>1)
        {
          cout << " constructing newform # " << i << endl;
          cout << " intdata  = "<< intdata[i] <<endl;
          cout << " Quaddata = "<< Quaddata[i] <<endl;
          cout << " aqs = "<< aqs[i] <<endl;
          cout << " aps = "<< aps[i] <<endl;
        }
      // the constructor here calls eigs_from_data() so these newforms have their eigs lists
      nflist.push_back(newform(this,i+1, intdata[i],Quaddata[i], aqs[i],aps[i]));
    }
  return 1;
}

void newforms::makebases(int extra_data)
{
  if(have_bases) return;
  makeh1();  // create the homology space if not yet
  sort_eigs();   // sort the newforms by their eigs list for efficient basis recovery
  // fill the h1matops and eigranges lists with empties; the functions
  // h1matop(i) and eigrange(i) will compute and store these the first
  // time they are needed.
  int maxdepth(MAXDEPTH);
  for (int i=0; i<maxdepth; i++)
    {
      h1matops.push_back(matop());
      eigranges.push_back(vector<long>());
    }
  form_finder splitspace(this, modulus, 1, maxdepth, 0, 1, 0, verbose);
  if(verbose) cout<<"About to recover "<<n1ds<<" newform bases (nap="<<nap<<")"<<endl;
  for (use_nf_number=0; use_nf_number<n1ds; use_nf_number++)
    {
      if (verbose)
        {
          cout<<"Recovering newform #"<<(use_nf_number+1) <<", eigs ";
          vec_out(cout, nflist[use_nf_number].eigs, 10);
          cout<< "...\n";
        }
      splitspace.splitoff(nflist[use_nf_number].eigs);
    }
  if(verbose) cout<<"Finished recovering newform bases, resorting back into lmfdb order..."<<endl;
  sort_lmfdb();
  if (extra_data)
    {
      if(verbose>1) cout<<"Filling in newform data..."<<endl;
      // no need to recompute ALs, but do recompute cuspidalfactor, loverp and integration  matrix
      fill_in_newform_data(0, 1, 1, 1);
    }
  have_bases=1;
  if(verbose)
    {
      cout<<"Finished makebases()";
      //if (extra_data) cout<<" with extra data";
      cout<<endl;
    }
}

// getap() calls apvec(P) to compute eigenvalues e for each newform,
// for the primes P in the given range, as specified above.

void newforms::getap(int first, int last, int verbose)
{
  if (n1ds==0)
    return;
  int nQP = Quadprimes::list.size();
  if(last>nQP)
    {
      last=nQP;
      cout<<"Cannot compute more than "<<nQP
          <<" ap since we only have that many primes precomputed"<<endl;
    }
  if(last<=nap)
    {
      cout<<"Already have "<<nap <<" ap " << "at level "<<N<<" so no need to compute more"<<endl;
    }
  // now nap < last <= nQP

  vector<int> nonsquarebadprimes;

  if (first==1)
    {
      if (n2r>0)
        {
          // In the case of odd class number we will have already computed
          // all the W(Q) eigenvalues.  Otherwise, we compute all those
          // eigenvalues we can now (i.e. those for which Q^e has square
          // class); for any which we cannot yet compute the stored values
          // will be 0 and will be computed later, once we have enough a_P
          // for P covering the genus classes. We store the list of
          // indices of bad Q for which this will be needed.

          for (auto Qi = badprimes.begin(); Qi!=badprimes.end(); ++Qi)
            {
              vector<long> aq = apvec(*Qi);
              if (aq[0]==0)
                nonsquarebadprimes.push_back(Qi-badprimes.begin());
              for (int j=0; j<n1ds; j++)
                nflist[j].aqlist.push_back(aq[j]);
            }
          // Sign of functional equation = minus product of all A-L eigenvalues
          if (characteristic==0 && nonsquarebadprimes.empty())
            for (int j=0; j<n1ds; j++)
              {
                nflist[j].sfe = std::accumulate(nflist[j].aqlist.begin(),nflist[j].aqlist.end(),-1,std::multiplies<long>());
              }
          else
            if (verbose)
              {
                cout<<" * computed A-L eigenvalues for "<<badprimes.size()-nonsquarebadprimes.size()<<" out of "<<badprimes.size()<<" bad primes"<<endl;
                cout<<" * missing:";
                for (int j=0; j<(int)nonsquarebadprimes.size(); j++)
                  cout<<badprimes[nonsquarebadprimes[j]]<<" ";
                cout<<endl;
              }
        }
    }
  else // first>1, i.e. we are computing more ap
    {
      if (verbose>1)
        cout << "We have "<<nap<<" eigenvalues out of "<<last<<" already, so we need to compute "<<(last-nap)<<" more."<<endl;
    }

  auto pr = Quadprimes::list.begin()+first-1;
  while((pr!=Quadprimes::list.end()) && (nap<last))
    // We should not stop when nap==last if for some newform its
    // genus_classes.size() is less than nchi, or nchi/2 if it has
    // extra twist.
    {
      Quadprime P = *pr++;
      long vp = val(P, N);
      vector<long> apv;

      // For bad primes P: if they have square class, or if first>1,
      // then we have already computed the eigenvalues and stored them
      // in nflist[*].aqlist. Here we only need to recover those
      // values when vp==1 to store their negatives in the
      // nflist[*].aplist.  If first==1 and they don't have square
      // class, we do nothing now except store 0's as placeholders in
      // the nflist[*].aplist.

      if (vp==0)
        {
          apv=apvec(P); // list of all T(P) eigenvalues
        }
      else // bad prime
        {
          int k = std::find(badprimes.begin(), badprimes.end(), P) - badprimes.begin();
          if (std::find(nonsquarebadprimes.begin(), nonsquarebadprimes.end(), k) == nonsquarebadprimes.end())
            {
              for (int j=0; j<n1ds; j++)
                apv.push_back(nflist[j].aqlist[k]);
            }
          else // we don't know the eigenvalues yet
            {
              apv.resize(n1ds,0);
            }
        }

      if(verbose)
        {
          string PorQ = (P.divides(N)? "Q": "P");
          cout<<PorQ<<" = "<<P<<" = "<<gens_string(P)<<"\tN("<<PorQ<<") = "<<P.norm()<<"\t";

          for (int i=0; i<n1ds; i++)
            {
              cout<<setw(5);
              long ap = apv[i];
              if ((characteristic>0) && (ap==-999))
                cout<<"?";
              else
                {
                  if (vp==0 || ap!=0)
                    cout << ap;
                  else
                    cout << "?";
                }
              cout <<" ";
            }
          cout << endl;
        }
      for (int i=0; i<n1ds; i++)
        {
          long ap = apv[i];
          long cp = ((vp==0)||(characteristic>0)? ap : (vp==1? -ap : 0));
          nflist[i].aplist.push_back(cp);
        }
      nap++;
    }

  // fill in A-L eigenvalues for any missing Q (only relevant for even class number)

  if (first==1 && !nonsquarebadprimes.empty())
    {
      if (verbose)
        cout<<"Computing missing A-L eigenvalues"<<endl;
      for (auto i=nonsquarebadprimes.begin(); i!=nonsquarebadprimes.end(); ++i)
        {
          Quadprime Q = badprimes[*i];
          if (verbose)
            cout<<"Computing W("<<Q<<")"<<endl;
          int e = val(Q,N);
          vector<long> aqv=apvec(Q); // list of all W(Q) eigenvalues

          for (int j=0; j<n1ds; j++)
            {
              //long aq = aqv[j]; //nflist[j].eigenvalueAtkinLehner(Q, verbose);
              long aq = nflist[j].eigenvalueAtkinLehner(Q, verbose);
              nflist[j].aqlist[*i] = aq;
              int k = std::find(Quadprimes::list.begin(), Quadprimes::list.end(), Q) - Quadprimes::list.begin();
              nflist[j].aplist[k] = (e>1? 0: -aq);
            }
        }
      // Sign of functional equation = minus product of all A-L eigenvalues
      if (characteristic==0)
        for (int j=0; j<n1ds; j++)
          {
            nflist[j].sfe = std::accumulate(nflist[j].aqlist.begin(),nflist[j].aqlist.end(),-1,std::multiplies<long>());
        }
    }

  if (first>1)
    return;

  int nfakes=0; // counts number of "fake rationals" to be deleted
  for (int i=0; i<n1ds; i++)
    {
      if (nflist[i].fake)
        nfakes++;
    }
  if (nfakes)
    {
      if (verbose)
        cout<<"The "<<n1ds<<" rational eigenspaces found originally include "
            <<nfakes<<" fake rational(s) which will now be deleted." << endl;
      int newn1ds=0;
      for (int i=0; i<n1ds; i++)
        {
          while (i<n1ds && nflist[i].fake)
            i++;
          if (i==n1ds)
            break; // for when we run off the end
          if (newn1ds<i)
            {
              if (verbose)
                cout<<"original eigenspace "<<i<<" being renumbered "<<newn1ds<<endl;
              nflist[newn1ds] = nflist[i];
            }
          newn1ds++;
        }
      assert(newn1ds+nfakes==n1ds);
      n1ds = newn1ds;
      nflist.resize(n1ds);
      n2ds += nfakes;
      assert(n1ds+n2ds==dimtrivcuspnew);
      if (verbose)
        {
          cout<<"Revised n1ds = "<<n1ds<<"; re-filling in newform data"<<endl;
        }
      fill_in_newform_data();
    }

  // Test each newform for being base-change (or twist of) or CM.

  // Also (only relevant for even class number), test each newform for
  // being unramified self-twist.  If so, then the number of
  // genus_classes should be 2^{n2r-1} = nchi/2, else it should be
  // 2^n2r = nchi.  NB If we don't compute many aP, we may not have
  // covered all the genus classes (either nchi or nchi/2) with at
  // least one P for which aP is nonzero. In that case we should
  // compute more aP.
  for (int i=0; i<n1ds; i++)
    {
      nflist[i].base_change_code();
      int cmd = nflist[i].is_CM();
      if (cmd)
        {
          int ngcl = nflist[i].genus_classes.size();
          INT D1(cmd);
          if (posmod(D1,4)!=1) D1*=4;
          if (div_disc(D1, Quad::disc)) // then we have an unramified self-twist
            {
              nflist[i].CMD = D1;
              if (2*ngcl<nchi)
                {
                  cout<<"Newform "<<i<<" appears to have self-twist by "<<D1
                      <<" but only "<<ngcl<<" genus classes have nonzero aP so far, out of "<<nchi/2<<" expected."
                      <<" Compute more aP!." << endl;
                }
              //assert (ngcl == nchi/2);
            }
          else // the usual case, not a self-twist
            {
              if (ngcl<nchi)
                {
                  cout<<"Newform "<<i<<" appears to not have self-twist, "
                      <<" and only "<<ngcl<<" genus classes have nonzero aP so far, out of "<<nchi<<" expected."
                      <<" Compute more aP!." << endl;
                }
              //assert (ngcl == nchi);
            }
        }
    }

  // Now we can split n1ds, and hence n2ds, by character:
  if (n2r==0) // trivial special case
    {
      new1dims.resize(1, n1ds);
      new2dims.resize(1, n2ds);
    }
  else
    {
      new1dims.resize(nchi, 0);
      for (int i=0; i<n1ds; i++)
        {
          // Test each newform for being unramified self-twist:
          INT cmd = nflist[i].CMD;
          // cout<<"form #"<<i<<": cmd="<<cmd<<endl;
          int j=0;
          if (cmd!=0) // then it is
            j = std::find(Quad::all_disc_factors.begin(), Quad::all_disc_factors.end(), cmd)
              - Quad::all_disc_factors.begin();
          new1dims[j] +=1;
        }
      // cout<<"new1dims = "<<new1dims<<endl;
      new2dims.resize(nchi, 0);
      for (int j=0; j<nchi; j++)
        {
          new2dims[j] = newdims[j] - new1dims[j];
          assert (new2dims[j]>=0 && "found more newforms in one component than the dimensions!");
        }
    }
}

void newforms::output_to_file(string eigfile) const
{
  int echo=0; // Set to 1 to echo what is written to the file for debugging
  ofstream out;
  out.open(eigfile.c_str());
  // Line 1
  out<<n1ds<<" "<<n2ds<<endl;
  if(echo) cout<<n1ds<<" "<<n2ds<<endl;
  if (n2r>0)
    {
      int nchi = 1<<n2r;
      for (int i=0; i<nchi; i++)
        out << new1dims[i] << " ";
      for (int i=0; i<nchi; i++)
        out << new2dims[i] << " ";
      out<<endl;
    }

  if(!n1ds)
    {
      out.close();
      return;
    }
  // Line 2
  out<<nap<<endl;  if(echo) cout<<nap<<endl;

  if (characteristic==0)
    {

  // Line 3: SFE

  for( const auto& f : nflist)
    {
      out<<setw(5)<<(f.sfe)<<" ";
      if(echo) cout<<setw(5)<<(f.sfe)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;

  // Line 4: pdot

  for( const auto& f : nflist)
    {
      out<<setw(5)<<(f.pdot)<<" ";
      if(echo) cout<<setw(5)<<(f.pdot)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;

  // Line 5: dp0

  for( const auto& f : nflist)
    {
      out<<setw(5)<<(f.dp0)<<" ";
      if(echo) cout<<setw(5)<<(f.dp0)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;

  // Line 6: cuspidal factor

  for( const auto& f : nflist)
    {
      out<<setw(5)<<(f.cuspidalfactor)<<" ";
      if(echo) cout<<setw(5)<<(f.cuspidalfactor)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;

  // Line 7: lambda

  for( const auto& f : nflist)
    {
      Quad lambda = f.lambda;
      out<<setw(5)<< lambda.re()<<" "<< lambda.im()<<" ";
      if(echo) cout<<setw(5)<< lambda.re()<<" "<< lambda.im()<<" ";
    }
  out<<endl;  if(echo) cout<<endl;

  // Line 8: lambdadot

  for( const auto& f : nflist)
    {
      out<<setw(5)<<(f.lambdadot)<<" ";
      if(echo) cout<<setw(5)<<(f.lambdadot)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;

  // Lines 9,10,11,12: a,b,c,d:

  Quad a;
  for( const auto& f : nflist)
    {
      a = f.a;
      out<<setw(5)<< a.re()<<" "<< a.im()<<" ";
      if(echo) cout<<setw(5)<< a.re()<<" "<< a.im()<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  for( const auto& f : nflist)
    {
      a = f.b;
      out<<setw(5)<< a.re()<<" "<< a.im()<<" ";
      if(echo) cout<<setw(5)<< a.re()<<" "<< a.im()<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  for( const auto& f : nflist)
    {
      a = f.c;
      out<<setw(5)<< a.re()<<" "<< a.im()<<" ";
      if(echo) cout<<setw(5)<< a.re()<<" "<< a.im()<<" ";
    }
  out<<endl;  if(echo) cout<<endl;
  for( const auto& f : nflist)
    {
      a = f.d;
      out<<setw(5)<< a.re()<<" "<< a.im()<<" ";
      if(echo) cout<<setw(5)<< a.re()<<" "<< a.im()<<" ";
    }
  out<<endl;  if(echo) cout<<endl;

  // Line 13: matdot

  for( const auto& f : nflist)
    {
      out<<setw(5)<<(f.matdot)<<" ";
      if(echo) cout<<setw(5)<<(f.matdot)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;

  // Line 14: bc

  for( const auto& f : nflist)
    {
      out<<setw(5)<<(f.bc)<<" ";
      if(echo) cout<<setw(5)<<(f.bc)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;

  // Line 15: cm code

  for( const auto& f : nflist)
    {
      out<<setw(5)<<(f.cm)<<" ";
      if(echo) cout<<setw(5)<<(f.cm)<<" ";
    }
  out<<endl;  if(echo) cout<<endl;

  // Line 16: self-twist code (only present when class number is even)

  if (n2r>0)
    {
      for( const auto& f : nflist)
        {
          out<<setw(5)<<(f.CMD)<<" ";
          if(echo) cout<<setw(5)<<(f.CMD)<<" ";
        }
      out<<endl;  if(echo) cout<<endl;
    }
  out<<endl;  if(echo) cout<<endl;

  // One line per bad prime: Atkin-Lehner eigenvalues

  for(int i=0; i<(int)badprimes.size(); i++)
    {
      for( const auto& f : nflist)
	{
	  out<<setw(5)<<(f.aqlist)[i]<<" ";
	  if(echo) cout<<setw(5)<<(f.aqlist)[i]<<" ";
	}
      out<<endl;
      if(echo) cout<<endl;
    }
  out<<endl;  if(echo) cout<<endl;
    } // end of if (characteristic==0) block

  // One line per prime: Fourier coefficients (=Hecke eigenvalues for good primes)

  for(int i=0; i<nap; i++)
    {
      for( const auto& f : nflist)
	{
	  out<<setw(5)<<(f.aplist)[i]<<" ";
	  if(echo) cout<<setw(5)<<(f.aplist)[i]<<" ";
	}
      out<<endl;      if(echo) cout<<endl;
    }
  out.close();
}

// For each newform we want a pivotal index j in [1,ngens] such that
// the j'th coordinate is nonzero, so that we can compute the Hecke
// eigenvalue a(p) by only computing the image of the j'th
// edge-generator.

void newforms::find_jlist()
{
  if(verbose>1)
    cout<<"Finding pivotal indices..."<<flush;

  int i, j, ok=0; j0=0;

  // First we see whether a single j works for all the newforms:

  for(j=1; (!ok)&&(j<=h1->ngens); j++)
    {
      ok=1;
      for (i=0; (i<n1ds)&&ok; i++)
        ok=(nflist[i].basis[j]!=0);
      if(ok) j0=j;
    }
  if(ok)
    {
      jlist.insert(j0);
      modsym m0 = h1->generator(j0);
      mjlist[j0] = m0;
      for (i=0; i<n1ds; i++)
	{
          newform& nfi = nflist[i];
	  nfi.j0 = j0;
	  nfi.m0 = m0;
          nfi.fac = I2long(nfi.basis[j0]);
          if(nfhmod!=0) nfi.facinv=invmod(nfi.fac, I2long(nfhmod));
	}
      if(verbose>1)
        cout<<"index j0="<<j0<<" works as a pivot for all newforms"<<endl;
      return;
    }

  if(verbose>1)
    cout<<"...failed to find a pivotal index which works for all newforms..." <<flush;

  // Instead, find a set of pivots to use:
  for (i=0; i<n1ds; i++)
    {
      newform& nfi = nflist[i];
      vec& bas = nfi.basis;
      j=1; while(bas[j]==0) j++;
      modsym m = h1->generator(j);
      jlist.insert(j); // jlist is a set so this will do nothing if it is already there
      nfi.j0 = j;
      nfi.m0 = m;
      mjlist[j] = m;
      nfi.fac = I2long(bas[j]);
      if(nfhmod!=0) nfi.facinv=invmod(nfi.fac, I2long(nfhmod));
    }
  if(verbose>1)
    cout<<"set of pivotal indices = "<<jlist<<endl;
}

//#define DEBUG_APVEC

// compute eigenvalues given the image images[j] for each j in jlist
vector<long> newforms::apvec_from_images(map<int,vec> images, pair<long,long> apbounds, const string& name)
{
  vector<long> apv(n1ds);

  for (int i=0; i<n1ds; i++)
    {
      int j0 = nflist[i].j0;
      int fac = nflist[i].fac;
      int top = I2long(images[j0][i+1]);
      long ap;
      // The eigenvalue is now top/fac (which should divide exactly)
      long nfhm = I2long(nfhmod);
      if(nfhm)
        ap=mod(xmodmul(top,nflist[i].facinv,nfhm), nfhm);
      else
        {
#ifdef DEBUG_APVEC
          cout << "ap   = " << top << "/" << fac << " = " <<top/fac<<endl;
#endif
          if (top%fac !=0)
            {
              cout<<"Problem in apvec: for newform #"<<(i+1)<<", with pivotal index "<<j0<<" and pivot "<<fac<<endl;
              cout<<"\timage list = "<<images[j0]<< " has "<<(i+1)<<" entry "<<top<<" which is not divisible by pivot "<<fac<<endl;
              cout<<flush;
            }
           ap = top/fac;
        }

      if (characteristic>0)
        apv[i] = posmod(ap, characteristic);
      else
        {
          apv[i] = ap;
          // check it is in range (in characteristic 0 only):
          if ((ap<apbounds.first)||(ap>apbounds.second))
            {
              cout<<"Error:  eigenvalue "<<ap<<" for operator "<<name
                  <<" for form # "<<(i+1)<<" is outside valid range "
                  << apbounds.first<<"..."<<apbounds.second<<endl;
              exit(1);
            }
        }
    }
  return apv;
}

// compute eigenvalue of op for each newform and check that it is within bounds
vector<long> newforms::apvec(const matop& op, pair<long,long> apbounds)
{
#ifdef DEBUG_APVEC
  cout<<"In apvec with operator "<<op.name()<<endl;
#endif
  // Compute the image images[j] of the j'th symbol under op, for all necessary j.
  map<int,vec> images;
  for( auto j : jlist) // between 1 and ngens inclusive
    images[j] = h1->applyop(op, mjlist[j], 1);

  vector<long> apv = apvec_from_images(images, apbounds, op.name());
#ifdef DEBUG_APVEC
  cout << "eigenvalue list = " << apv << endl;
#endif
  return apv;
}

// Special code for T(P) for Euclidean fields, for good P only

// The following utility does the following.  Given an integer ind:
// - if ind>0 it adds the ind'th row of pcd to imagej;
// - if ind<0 it subtracts the |ind|'th row of pcd from imagej;
// - if ind=0 it leaves imagej unchaged.
// if hmod is nonzero the vector addition is done modulo hmod.

void update(const mat& pcd, vec& imagej, long ind, scalar hmod)
{
  if (ind==0) return;
  vec part = (ind>0? pcd.row(ind): -pcd.row(-ind));
  imagej = reduce_mod_p(imagej + part, hmod);
}

// compute eigenvalue at P for each newform (good P, Euclidean) and check that it is in [minap,..,maxap]
vector<long> newforms::apvec_euclidean(Quadprime& P, pair<long,long> apbounds)
{
  assert (Quad::is_Euclidean && "field must be Euclidean in apvec_euclidean()");
  assert (val(P,N)==0 && "P must be good in apvec_euclidean()");
  Quad p = P.gen();

  //images[j] is the image of the j'th M-symbol
  map<int,vec> images;
  Quad a,c,q,u1,u2,u3;

  // Compute the image of the necessary M-symbols (hopefully only one)

  for( const auto& j : jlist)
    {
      vec imagej=vec(n1ds); // initialised to 0
      // Since this code is only used in the Euclidean case,
      // all symbols have type 0
      long s_number = h1->ER.gen(j);  // (c:d) symbol number
      Quad u, v;
      h1->P1.make_symb(s_number, u, v);

      // Now we compute the projected image of symbol s=(u:v)_t under
      // T(P).

      // This code is for Euclidean fields only, using Manin-Heilbronn
      // matrices: Loop over residues res mod P and for each res
      // compute several M-symbol image parts (u1:v1).  Accumulate the
      // associated vectors in vec imagej using the utility
      // update(projcoord, imagej, ind, nfhmod), where (u1:v1) is the
      // ind'th symbol.

      mat& pcd = h1->projcoord;

      // Matrix [1,0;0,p]
      long ind = h1->ER.coords(h1->index(u,p*v));
      update(pcd,imagej,ind,nfhmod);

      // Matrix [p,0;0,1]
      ind = h1->ER.coords(h1->index(p*u,v));
      update(pcd,imagej,ind,nfhmod);

      // Other matrices, several for each nonzero residue b mod p
      vector<Quad> resmodp = P.residues();
      for ( auto& b : resmodp)
        {
          if(b.is_zero()) continue; // handled above as special case
          a = -p;
          u1=u*p; u2=v-u*b;
          ind = h1->ER.coords(h1->index(u1,u2));
          update(pcd,imagej,ind,nfhmod);
          while(!b.is_zero())
            {
              q=a/b; c=a-b*q; u3=q*u2-u1;
              a=-b; b=c; u1=u2; u2=u3;
              ind = h1->ER.coords(h1->index(u1,u2));
              update(pcd,imagej,ind,nfhmod);
            }
        }
      images[j]=imagej;
    }

  vector<long> apv = apvec_from_images(images, apbounds, opname(P,N));
#ifdef DEBUG_APVEC
  cout << "eigenvalue list = " << apv << endl;
#endif
  return apv;
}

// apvec(P) returns a list of eigenvalues e for the prime P (one for each newform),
// specified as follows:
//
// - for good P with [P] square, e=a(P)
// - for good P with [P] non-square, e=|a(P)| (sign not yet determined), via a(P)^2, via a(P^2)
// - for bad P with Q=P^e||N and [Q] square, the A-L eigenvalue at Q
// - for bad P with Q=P^e||N and [Q] non-square, +1 (sign not yet determined)

// compute a[P] for each newform
vector<long> newforms::apvec(Quadprime& P)
{
#ifdef DEBUG_APVEC
  cout<<"In apvec with P = "<<P<<endl;
#endif
  long normp=I2long(P.norm());
  int i, vp = val(P,N);

  if ((characteristic>0) && ((vp>0) || (normp%characteristic ==0)))
    {
      vector<long> apv(n1ds, -999);
#ifdef DEBUG_APVEC
      cout << "ignored prime: ap list = " << apv << endl;
#endif
      return apv;
    }

  if (vp>0) // bad prime
    {
      long e = val(P,N);
      Qideal Pe = P;
      while (--e) Pe*=P;
      if (Pe.is_principal())
        {
          return apvec(AtkinLehnerQOp(P,N), {-1,1});
        }
      if (Pe.has_square_class())
        {
          Qideal A = Pe.sqrt_coprime_to(N);
          return apvec(AtkinLehnerQChiOp(P,A,N), {-1,1});
        }
      // Other values need to be found differently -- here we set them to 0
      vector<long> apv(n1ds, 0);
      return apv;
    }

  // now P is a good prime

  if (Quad::is_Euclidean)
    return apvec_euclidean(P, eigenvalue_range(P));

  if (P.is_principal())
    return apvec(HeckePOp(P,N), eigenvalue_range(P)); // T(P)

  if (P.has_square_class())
    {
      Qideal A = P.sqrt_coprime_to(N);
      return apvec(HeckePChiOp(P,A,N), eigenvalue_range(P)); // T(P)*T(A,A)
    }

  // Now the class number is even, and we deal with the newforms individually.
  // NB The call to eigenvalueHecke() may set the 'fake' flag on some newforms.
  vector<long> apv;
  for (i=0; i<n1ds; i++)
    apv.push_back(nflist[i].eigenvalueHecke(P, verbose));
  return apv;
}

// Strategy for operators used to automatically cut out 1-dimensional
// rational eigenspaces, in characteristic 0. NB We no longer include
// Atkin-Lehner operators W_Q in splitting, owing to the complications
// in computing oldform multiplicities.

// The first n2r (=Quad::class_group_2_rank) operators are T(A,A)
// where A runs over n2r ideals coprime to N whose classes generate
// the 2-torsion in the class group. These are involutions, but we
// only consider the +1-eigenspace since we only want to find
// eigenspaces with trivial unramified character.  When the class
// number is odd, then n2r=0 so there are none of these.  The list of
// ideals A used here is newforms::nulist.

// Then come nap operators for good primes P, where the constructor
// sets nap to the default MAXDEPTH (by default) and fills the array
// goodprimes with the first nap primes not dividing N.  The operator
// for P is *either* T(A,A)*T(P), when the class [P] is square and
// A^2*P is principal; *or* T(A,A)*T(P^2) when [P] is not square,
// and A*P is principal.
//
// In the first case the eigenvalues considered are integers a with
// |a|<=2*sqrt(N(P)).  In the second case we use the identity
//
// T(P^2) = T(P)^2 - N(P)T(P,P)
//
// to deduce that -- when the central character is trivial -- the
// eigenvalues satisfy
//
// a(P^2) = a(P)^2 - N(P)
//
// so the eigenvalues we consider for T(A,A)*T(P^2) are
// {a^2-N(P) : 0<=a<=2*sqrt(N(P))}.

matop newforms::h1matop(int i) // return the list of matrices defining the i'th operator
{
  assert (i>=0);
  if (h1matops[i].length()!=0)
    return h1matops[i];
  // else we have not already computed and stored it, so we do so now

  // cout<<"Computing matop "<<i<<"..."<<flush;
  if (i<n2r) // then we use T(A,A) where A is the i'th generator of the class group mod squares
    {
      h1matops[i] = CharOp(nulist[i], N);
      return h1matops[i];
    }
  // else we yield, for P the (i-n2r)'th good prime,
  // T(P)          if P is principal, or
  // T(A,A)*T(P)   if [P] is square with A^2*P principal, or
  // T(P^2)        if P^2 is principal, or
  // T(A,A)*T(P^2) where A*P is principal

  Quadprime P = goodprimes[i-n2r];
  if (P.is_principal())
    {
      h1matops[i] = HeckePOp(P, N);
      return h1matops[i];
    }
  if (P.has_square_class())
    {
      Qideal A = P.sqrt_coprime_to(N);
      h1matops[i] = HeckePChiOp(P, A, N);
      return h1matops[i];
    }
  Qideal P2 = P*P;
  if (P2.is_principal())
    {
      h1matops[i] = HeckeP2Op(P, N);
      return h1matops[i];
    }
  Qideal A = P.equivalent_mod_2_coprime_to(N,1);
  h1matops[i] = HeckeP2ChiOp(P, A, N);
  return h1matops[i];
}

// the list of possible (integer) eigenvalues for the i'th operator:
vector<long> newforms::eigrange(int i)
{
  assert (i>=0);
  if (eigranges[i].empty())
    {
      if (verbose>1)
        cout<<"Filling eigrange["<<i<<"]"<<endl;
      if (i<n2r)
        {
          eigranges[i] = {1};
        }
      else
        {
          Quadprime P = goodprimes[i-n2r];
          if (verbose>1)
            cout<<"Using goodprimes["<<i-n2r<<"] = "<<P<<endl;
          if (characteristic>0)
            eigranges[i] = range(0,characteristic-1);
          else
            eigranges[i] = good_eigrange(P);
          if (verbose>1)
            cout << "eigrange for P = " << P << " (norm "<<P.norm()<<"):\t" << eigranges[i] << endl;
        }
    }
  else
    {
      if (verbose>1)
        cout << "cached eigrange[" << i << "] is " << eigranges[i] << endl;
    }
  return eigranges[i];
}

// Conjugate newforms data: NB this *only* conjugates data needed for
// output_to_file(), not everything!
void newforms::conjugate(int debug)
{
  h1=0;
  of=0;
  // Conjugate level (and do nothing more if level is self-conjugate):
  Qideal Nbar = N.conj();
  if (N==Nbar)
    return;
  if (debug)
    cout<<"Conjugate level = "<<Nbar<<endl;
  N = Nbar;

  // Conjugate badprimes:

  if (debug)
    cout<<"Original bad primes = "<<badprimes<<endl;
  vector<Quadprime> conj_badprimes = Nbar.factorization().sorted_primes();
  for (int i=0; i<(int)badprimes.size(); i++)
    {
      bad_prime_conjugation_permutation.push_back(std::find(badprimes.begin(), badprimes.end(), conj_badprimes[i].conj()) - badprimes.begin());
    }
  if (debug)
    {
    cout<<"Conjugation permutation of bad primes = "<<bad_prime_conjugation_permutation<<endl;

    cout<<"Before conjugating individual newforms:"<<endl;
    display();
    }

  // Conjugate aq and ap for each newform:

  for (auto f = nflist.begin(); f!=nflist.end(); f++)
    {
      f->nf=this;
      f->conjugate(debug);
    }
  if (debug)
    {
      cout<<"After conjugating individual newforms:"<<endl;
      display();
    }
}

// Conjugate newform data: NB this *only* conjugates data needed for
// output_to_file(), not everything!
void newform::conjugate(int debug)
{
  if (debug)
    cout<<"Conjugating data for newform "<<index<<" at level "<<nf->N<<endl;

  // Conjugate AL eigs:

  vector<long> conj_aqlist; // temporary for permuted values
  for (int i=0; i<(int)aqlist.size(); i++)
    {
      conj_aqlist.push_back(aqlist[nf->bad_prime_conjugation_permutation[i]]);
    }
  if (debug)
    cout<<"Conjugate aqlist =  "<<conj_aqlist<<endl;
  aqlist = conj_aqlist;

  // Conjugate aP coeffs (for all P, not just good P):

  vector<long> conj_aplist; // temporary for permuted values
  for (int i=0; i<(int)aplist.size(); i++)
    {
      conj_aplist.push_back(aplist[Quadprimes::conjugate_index(i)]);
    }
  if (debug)
    cout<<"Conjugate aplist =  "<<conj_aplist<<endl;
  aplist = conj_aplist;

  // Conjugate integration matrix:

  a = a.conj();
  b = a.conj();
  c = a.conj();
  d = a.conj();
}

//end of newforms.cc
