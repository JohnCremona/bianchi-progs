// FILE HOMSPACE.CC: Implemention of class homspace

//#define USE_CRT // if using smats  mod MODULUS, try CRT-ing with another prime
                // NB this is experimental only

#include <eclib/msubspace.h>
#include <eclib/xmod.h>
#include "homspace.h"
const string W_opname("W");
const string T_opname("T");

int liftmats_chinese(const smat& m1, scalar pr1, const smat& m2, scalar pr2,
                     smat& m, scalar& dd);

int sign(int a)
{
  if (a) return (a>0? 1: -1);
  return 0;
}

// In add_rel(rel), rel is a list of (positive) (c:d)-symbol numbers i
// such that the corresponding symbols add to 0 in homology.  We use
// the map i -> j=coordindex[i] to convert this to a vector of
// dimension ngens, and append that as a new row to the relation
// matrix relmat.

void homspace::add_rel(const vector<int>& rel)
{
  vector<int>::const_iterator r;
  long c;
  if (verbose)
    {
      cout<<"Relation: ";
      for (r = rel.begin(); r!=rel.end(); r++)
        cout<<(*r)<<" "<<symbol(*r)<<" ";
      cout <<" --> ";
    }
#ifdef USE_SMATS
  svec relation(ngens);
#else
  vec relation(ngens);
#endif
  for (r = rel.begin(); r!=rel.end(); r++)
    {
      c = coordindex[*r];
      if(c)
#ifdef USE_SMATS
        relation.add(abs(c), sign(c));
#else
        relation[abs(c)] += sign(c);
#endif
    }
  if (verbose) cout<<relation<<endl;
#ifdef USE_SMATS
  if(relation.size()==0)
#else
  if(trivial(relation))
#endif
  return;
  numrel++;
  if(numrel<=maxnumrel)
    {
      make_primitive(relation);
      relmat.setrow(numrel,relation);
    }
  else
    cerr<<"Too many face relations (numrel="<<numrel
        <<", maxnumrel="<<maxnumrel<<")"<<endl;
}

homspace::homspace(const Quad& n, int hp, int cuspid, int verb) :symbdata(n)
{
  verbose=verb;
  cuspidal=cuspid;
  hmod = 0;
  if (verbose) symbdata::display();
  plusflag=hp;                  // Sets static level::plusflag = hp
  long i,j,k;
  ngens=0;
  coordindex.resize(nsymb);
  gens.resize(nsymb+1);
  //NB start of gens array is at 1 not 0

  edge_relations(); // sets coordindex

  if (verbose)
    {
      cout << "After 2-term relations, ngens = "<<ngens<<endl;
      cout << "gens = ";
      for (i=1; i<=ngens; i++) cout << gens[i] << " ";
      cout << endl;
      cout << "coordindex = \n";
      for (i=0; i<nsymb; i++)
	cout << i<<":\t"<<symbol(i)<<"\t"<<coordindex[i] << "\n";
      cout << endl;
    }

  face_relations(); //fills relmat with the relations

  if(verbose)
    {
      cout << "Finished face relations: ";
      cout << "number of relations = " << numrel;
      cout << " (bound was "<<maxnumrel<<")"<<endl;
    }

   vec pivs, npivs;
#if(USE_SMATS)
   smat_elim sme(relmat);
   int d1;
   smat ker = sme.kernel(npivs,pivs), sp;
   int ok = liftmat(ker,MODULUS,sp,d1);
   int mod2 = 1073741783; // 2^30-41
   bigint mmod = to_ZZ(MODULUS)*to_ZZ(mod2);
   if (!ok)
     {
       if(verbose)
         cout << "failed to lift modular kernel using modulus "
              << MODULUS << endl;
#ifdef USE_CRT
       if(verbose)
         cout << "repeating kernel computation, modulo " << mod2 << endl;
       smat_elim sme2(relmat,mod2);
       vec pivs2, npivs2;
       smat ker2 = sme2.kernel(npivs2,pivs2), sp2;
       ok = (pivs==pivs2);
       if (!ok)
         {
           cout<<"pivs do not agree:\npivs  = "<<pivs<<"\npivs2 = "<<pivs2<<endl;
         }
       else
         {
           if(verbose) cout << " pivs agree" << endl;
           ok = liftmats_chinese(ker,MODULUS,ker2,mod2,sp,d1);
         }
       if (ok)
         {
           if(verbose)
             cout << "success using CRT, combined modulus = "<<mmod
                  <<", denominator= " << d1 << "\n";
         }
       else
         {
           if(verbose)
             cout << "CRT combination with combined modulus "<<mmod<<" failed\n" << endl;
         }
#endif
     }
   if (ok)
     {
       denom1 = d1;
     }
   else
     {
       hmod = MODULUS;
       denom1 = 1;
     }
   relmat=smat(0,0);
   if(verbose>1)
     {
       cout << "kernel of relmat = " << sp.as_mat() << endl;
       cout << "pivots = "<<pivs << endl;
       cout << "denom = "<<d1 << endl;
     }
   rk = sp.ncols();
   coord.init(ngens+1,rk); // 0'th is unused
   for(i=1; i<=ngens; i++)
     coord.setrow(i,sp.row(i).as_vec());
   // if hmod>0, coord is only defined modulo hmod
   sp=smat(0,0); // clear space
#else
  relmat = relmat.slice(numrel,ngens);
  if(verbose) 
    {
      cout << "relmat = "; relmat.output_pari(cout); cout << endl;
    }
  msubspace sp = kernel(mat_m(relmat),0);
  rk = dim(sp);
  if(verbose>2) cout<<"coord = "<<basis(sp)<<endl;
  coord = basis(sp).shorten((int)1);
  pivs = pivots(sp);
  denom1 = I2int(denom(sp));
  relmat.init(); relation.init(); sp.clear(); 
#endif
  if (verbose)
    {
      cout << "rk = " << rk << endl;
      cout << "coord:" << coord;
      if (hmod)
	cout << "failed to lift, coord is only defined modulo "<<hmod<<endl;
      else
        cout << "lifted ok, denominator = " << denom1 << endl;
      cout << "pivots = " << pivs <<endl;
    }
  if (rk>0)
    {
      freegens.resize(rk);
      for (i=0; i<rk; i++) freegens[i] = gens[pivs[i+1]];
      if (verbose)
        { 
          cout << "freegens: ";
          for (i=0; i<rk; i++) cout << freegens[i] << " ";
          cout << endl;
        }
    }
  pivs.init();
  {
    cusplist cusps(2*rk,this);
     mat deltamat(2*rk,rk);
     for (i=0; i<rk; i++)
       {
         modsym m(symbol(freegens[i]));
         for (j=1; j>-3; j-=2)
           {
             RatQuad c = (j==1 ? m.beta() : m.alpha());
             k = cusps.index(c);   //adds automatically if new
             deltamat(k+1,i+1) += j;  // N.B. offset of 1
           }
       }
     ncusps=cusps.count();
     if(verbose)cout << "ncusps = " << ncusps << endl;
     kern = kernel(smat(deltamat));
     vec pivs, npivs;
     int d2;
     smat sk; 
    int ok = liftmat(smat_elim(deltamat).kernel(npivs,pivs),MODULUS,sk,d2);
     if (!ok)
       cout << "**!!!** failed to lift modular kernel\n" << endl;
     denom2=d2;
    tkernbas = transpose(kern.bas());         // dim(kern) x rank
    deltamat.init(0); // clear space.
    if(verbose>1)
      {
        cout<<"tkernbas = "<<tkernbas.as_mat()<<endl;
      }
  }
   if(cuspidal)
     dimension = dim(kern);
   else
     dimension = rk;
  denom3 = denom1 * denom2;
  freemods.resize(rk);
  needed.resize(rk);
 
  if (dimension>0)
    {
      const smat& basiskern = basis(kern);
      if (verbose)  cout << "Freemods:\n";
      for (i=0; i<rk; i++)
        {
          freemods[i] = modsym(symbol(freegens[i])) ;
          needed[i]   =  (cuspidal? ! trivial(basiskern.row(i+1).as_vec())
                          : 1);
          if (verbose) 
            {
              cout << i << ": " << freemods[i];
              if (!needed[i]) cout << " (not needed)";
              cout << endl;
            }
        }
      if (verbose)
        {
          cout << "Basis of ker(delta):\n";
          cout << basiskern;
          cout << "pivots: " << pivots(kern) << endl;
        }
    }
  if (verbose) cout << "Finished constructing homspace.\n";
}

// 2-term (edge) relations

void homspace::edge_relations()    // computes coordindex
{
  Quad unit = fundunit;
  long lenrel = Quad::nunits;
  if(!plusflag) {unit=fundunit*fundunit; lenrel/=2;}
  symbop eps(this,unit,0,0,1);
  symbop sof(this,0,-1,1,0);
  vector<int> a(lenrel), b(lenrel);
  vector<int> check(nsymb, 0);
  int j, k, triv;
  if(verbose) cout << "About to start on 2-term relations.\n";
  for (j=nsymb-1; j>=0; j--)
    {
      if (check[j]==0)
        { 
	  if(verbose>1) cout << "j = " << j << ":\t";
          a[0]=j; b[0]=sof(j); triv=(j==b[0]);
          for(k=1; k<lenrel; k++)
            {
              a[k]= eps(a[k-1]);
              b[k]= eps(b[k-1]);
              triv= triv | (j==b[k]);
            }
          for (k=0; k<lenrel; k++) check[a[k]]=check[b[k]]=1;
	  if(verbose>1)
	    {
	      cout<<"+:\t";
	      for (k=0; k<lenrel; k++) cout<<a[k]<<" ";
              cout<<endl;
	      cout<<"\t-:\t";
	      for (k=0; k<lenrel; k++) cout<<b[k]<<" ";
              cout<<endl;
	    }
          if (triv)
            for (k=0; k<lenrel; k++) coordindex[a[k]]=coordindex[b[k]]=0;
          else
            {   
              gens[++ngens] = j;
              for(k=0; k<lenrel; k++)
                {
                  coordindex[a[k]] =  ngens;
                  coordindex[b[k]] = -ngens;
                }
            }
        }
    }
}

//
// face relations
//

void homspace::face_relations()
{
  long field = Quad::d;

  maxnumrel=2*nsymb;
  numrel = 0;

#if(USE_SMATS)
  relmat = smat(maxnumrel,ngens);
#else
  relmat.init(maxnumrel,ngens);
#endif
  triangle_relation_0();

  switch (field)
    {
    case 1:
    case 3:
      {
        if(!plusflag)
          triangle_relation_1_3();
        break;
      }
    case 2:
      {
        square_relation_2();
        break;
      }
    case 7:
      {
        rectangle_relation_7();
        break;
      }
    case 11:
      {
        hexagon_relation_11();
        break;
      }
    default:
      {
        cerr<<"homspace not implemented for field "<<field<<endl;
        exit(1);
      }
    }
}

// triangle relation for all fields
void homspace::triangle_relation_0()
{
  if(verbose)
    {
      cout << "Face relation type 1 (triangles):\n";
    }
  vector<int> rel(3);
  long j, k;
  symbop tof(this,1,1,-1,0);
  symbop rof(this,0,1,1,0);
  vector<int> check(nsymb, 0);
  for (k=0; k<nsymb; k++)
    if (check[k]==0)
      {
        for(j=0; j<3; j++)
          {
            rel[j] = (j? tof(rel[j-1]): k);
            check[rel[j]] = 1;
            if (plusflag)
              check[rof(rel[j])] = 1;
          }
        add_rel(rel);
      }
  if(verbose)
    {
      cout << "After face relation type 1, number of relations = " << numrel <<"\n";
    }
}

// extra triangle relation for fields 1, 3
void homspace::triangle_relation_1_3()
{
  if(verbose)
    {
      cout << "Face relation type 2 (triangles):\n";
    }
  vector<int> rel(3);
  long j, k;

  Quad w(0,1);
  long field = Quad::d;
  symbop xof = (field==1? symbop(this,w,1,1,0): symbop(this,1,w,w-1,0));

  vector<int> check(nsymb, 0);
  for (k=0; k<nsymb; k++)
    if (check[k]==0)
      {
        for(j=0; j<3; j++)
          {
            rel[j] = (j? xof(rel[j-1]): k);
            check[rel[j]] = 1;
          }
        add_rel(rel);
      }
  if(verbose)
    {
      cout << "After face relation type 2, number of relations = " << numrel <<"\n";
    }
}

// extra square relation for field 2
void homspace::square_relation_2()
{
  if(verbose)
    {
      cout << "Face relation type 2 (squares):\n";
    }
  vector<int> rel(4);
  long j, k;

  Quad w(0,1);
  symbop uof(this,w,1,1,0);
  symbop sof(this,0,-1,1,0);
  symbop J(this,fundunit,0,0,1);

  vector<int> check(nsymb, 0);
  for (k=0; k<nsymb; k++)
    if (check[k]==0)
      {
        for(j=0; j<4; j++)
          {
            rel[j] = (j? uof(rel[j-1]): k);
            if(plusflag) // NB det(uof)=-1
              {
                check[rel[j]] = 1;
                check[sof(rel[j])] = 1;
              }
            else
              {
                if(j%2==0)
                  check[rel[j]] = 1;
              }
          }
        if (!plusflag) // since uof() has det -1
          {
            rel[1] = J(rel[1]);
            rel[3] = J(rel[3]);
          }
        add_rel(rel);
      }
  if(verbose)
    {
      cout << "After face relation type 2, number of relations = " << numrel <<"\n";
    }
}

// extra rectangle relation for field 7
void homspace::rectangle_relation_7()
{
  if(verbose)
    {
      cout << "Face relation type 2 (rectangles):\n";
    }
  vector<int> rel(4);
  long j, k;

  Quad w(0,1);
  symbop yof(this,1,-w,1-w,-1);
  symbop us7of(this,w,-1,1,0);
  symbop rof(this,0,1,1,0);

  vector<int> check(nsymb, 0);
  for (k=0; k<nsymb; k++)
    if (check[k]==0)
      {
        for(j=0; j<4; j+=2) // j=0,2; j+1=1,3
          {
            rel[j] = (j? yof(rel[j-2]): k);
            check[rel[j]] = 1;
            rel[j+1] = us7of(rel[j]);
            if (plusflag)
              check[rof(rel[j+1])]=1;
          }
        add_rel(rel);
      }
}

// extra hexagon relation for field 11
void homspace::hexagon_relation_11()
{
  if(verbose)
    {
      cout << "Face relation type 2 (hexagons):\n";
    }
  vector<int> rel(6);
  long j, k;

  Quad w(0,1);
  symbop xof(this,1,-w,1-w,-2);
  symbop us11of(this,w,-1,1,0);
  symbop sof(this,0,-1,1,0);
  symbop rof(this,0,1,1,0);

  vector<int> check(nsymb, 0);
  for (k=0; k<nsymb; k++)
    if (check[k]==0)
      {
        for(j=0; j<6; j+=2) // j=0,2,4; j+1=1,3,5
          {
            rel[j] = (j? xof(rel[j-2]): k);
            check[rel[j]] = 1;
            rel[j+1] = us11of(rel[j]);
            if(plusflag)
              check[rof(rel[j+1])]=1;
          }
        add_rel(rel);
      }
  if(verbose)
    {
      cout << "After face relation type 2, number of relations = " << numrel <<"\n";
    }
}

void homspace::kernel_delta()
{

}

vec homspace::chain(const symb& s) const  //=old getcoord
{
 long i= coordindex[index(s)];
 vec ans(coord.ncols());
 if (i) ans = sign(i)*(coord.row(abs(i)));
 return ans;
}

vec homspace::chaincd(const Quad& c, const Quad& d) const //=old getcoord2
{
 long i= coordindex[index2(c,d)];
 vec ans(coord.ncols());
 if (i) ans = sign(i)*(coord.row(abs(i)));
 return ans;
}

vec homspace::projchaincd(const Quad& c, const Quad& d) const 
{
 long i= coordindex[index2(c,d)];
 vec ans(projcoord.ncols());
 if (i) ans = sign(i)*(projcoord.row(abs(i)));
 return ans;
}

vec homspace::chain(const Quad& nn, const Quad& dd) const //=old qtovec2
{
   vec ans = chaincd(0,1), part;
   Quad c=0, d=1, e, a=nn, b=dd, q, f;
   while (b!=0)
   { q=a/b;
     f=a; a=-b; b= f-q*b;
     e=d; d= c; c= e+q*c; c=-c;
     part = chaincd(c,d);
     if(hmod)
       ans.addmodp(part,hmod);
     else
       ans += part;
   }
   if(hmod) ans=reduce_modp(ans,hmod);
   return ans;
}

vec homspace::projcycle(const Quad& nn, const Quad& dd) const  //
{
  vec ans = projchaincd(0,1), part;
   Quad c=0, d=1, e, a=nn, b=dd, q, f;
   while (b!=0)
   { q=a/b;
     f=a; a=-b; b= f-q*b;
     e=d; d= c; c= e+q*c; c=-c;
     part = projchaincd(c,d);
     if(hmod)
       ans.addmodp(part,hmod);
     else
       ans += part;
   }
   if(hmod) ans=reduce_modp(ans,hmod);
   return ans;
}

vec homspace::applyop(const matop& mlist, const RatQuad& q) const
{ vec ans(rk), part;
  long i=mlist.length();
  while (i--)
    {
      part = chain(mlist[i](q));
      if(hmod)
        ans.addmodp(part,hmod);
      else
        ans += part;
    }
  if(hmod) ans=reduce_modp(ans,hmod);
  return ans;
}

mat homspace::calcop(const string opname, const Quad& p, const matop& mlist, int dual, int display) const
{
  mat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     { vec colj = applyop(mlist,freemods[j]);
       if(hmod) colj=reduce_modp(colj,hmod);
       m.setcol(j+1,colj);
     }
  if(cuspidal) m = restrict_mat(smat(m),kern).as_mat();
  if(dual) m = transpose(m);
  if (display) cout << "Matrix of " << opname << "(" << p << ") = " << m;
  if (display && (dimension>1)) cout << endl;
  return m;
}

vec homspace::calcop_col(const string opname, const Quad& p, const matop& mlist, int j, int display) const
{
  vec colj = applyop(mlist,freemods[j-1]);
  if(hmod) colj=reduce_modp(colj,hmod);
  return colj;
}

mat homspace::calcop_cols(const string opname, const Quad& p, const matop& mlist, const vec& jlist, int display) const
{
  int i, j, d = dim(jlist);
  mat m(d,rk);
  for (i=1; i<=d; i++)
    {
      j = jlist[i];
      vec colj = applyop(mlist,freemods[j-1]);
      if(hmod) colj=reduce_modp(colj,hmod);
      m.setcol(i,colj);
     }
  return m;
}
 
smat homspace::s_calcop_cols(const string opname, const Quad& p, const matop& mlist, const vec& jlist, int display) const
{
  int i, j, d = dim(jlist);
  smat m(d,rk);
  for (i=1; i<=d; i++)
    {
      j = jlist[i];
      svec colj = applyop(mlist,freemods[j-1]);
      if(hmod) colj.reduce_mod_p(hmod);
      m.setrow(i,colj);
     }
  return m;
}
 
smat homspace::s_calcop(const string  opname, const Quad& p, const matop& mlist, 
			int dual, int display) const
{
  smat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     { svec colj = applyop(mlist,freemods[j]);
       if(hmod) colj.reduce_mod_p(hmod);
       m.setrow(j+1,colj);
     }
  if(cuspidal)
    {
      m = restrict_mat(transpose(m),kern);
      if(dual) m = transpose(m);
    }
  else if(!dual) {m=transpose(m);}
  if (display) 
    {
      cout << "Matrix of " << opname << "(" << p << ") = ";
      if (dimension>1) cout << "\n";
      cout<<m.as_mat();
    }
  return m;
}

mat homspace::calcop_restricted(const string opname, const Quad& p, const matop& mlist, const subspace& s, int dual, int display) const
{
  long d=dim(s);
  mat m(d,rk);
  for (long j=0; j<d; j++)
     {
       long jj = pivots(s)[j+1]-1;
       vec colj = applyop(mlist,freemods[jj]);
       if(hmod) colj=reduce_modp(colj,hmod);
       m.setrow(j+1,colj);
     }
  if(hmod)
    m = matmulmodp(m,basis(s),hmod);
  else
    m = m*basis(s);
  if(!dual) m=transpose(m); // dual is default for restricted ops
  if (display) cout << "Matrix of " << opname << "(" << p << ") = " << m;
  if (display && (dimension>1)) cout << endl;
  return m;
}
 
smat homspace::s_calcop_restricted(const string opname, const Quad& p, const matop& mlist, const ssubspace& s, int dual, int display) const
{
  long d=dim(s);
  smat m(d,rk);
  for (long j=1; j<=d; j++)
     { 
       long jj = pivots(s)[j];
       svec colj = applyop(mlist,freemods[jj-1]);
       if(hmod) colj.reduce_mod_p(hmod);
       m.setrow(j,colj);
     }
  m = mult_mod_p(m,basis(s),MODULUS);
  m.reduce_mod_p();
  if(!dual) m=transpose(m); // dual is default for restricted ops
  if (display)
    {
      cout << "Matrix of " << opname << "(" << p << ") = " << m.as_mat();
      if (dimension>1) cout << endl;
    }
  return m;
}
 
mat homspace::heckeop(const Quad& p, int dual, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return calcop(name,p,matlist,dual,display);
}
 
vec homspace::heckeop_col(const Quad& p, int j, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return calcop_col(name,p,matlist,j,display);
}
 
mat homspace::heckeop_cols(const Quad& p, const vec& jlist, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return calcop_cols(name,p,matlist,jlist,display);
}
 
smat homspace::s_heckeop_cols(const Quad& p, const vec& jlist, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return s_calcop_cols(name,p,matlist,jlist,display);
}
 
smat homspace::s_heckeop(const Quad& p, int dual, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return s_calcop(name,p,matlist,dual,display);
}
 
mat homspace::heckeop_restricted(const Quad& p, const subspace& s, int dual, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return calcop_restricted(name,p,matlist,s,dual,display);
}
 
smat homspace::s_heckeop_restricted(const Quad& p, const ssubspace& s, int dual, int display) const
{
 matop matlist(p,modulus);
 string name = (div(p,modulus)) ? W_opname : T_opname;
 return s_calcop_restricted(name,p,matlist,s,dual,display);
}
 
mat homspace::wop(const Quad& q, int dual, int display) const
{
 matop matlist(q,modulus);
 return calcop(W_opname,q,matlist,dual,display);
}
 
mat homspace::fricke(int dual, int display) const
{
 matop frickelist(modulus,modulus);
 return calcop(W_opname,modulus,frickelist,dual,display);
}

mat homspace::opmat(int i, int dual, int verb)
{
  if((i<0)||(i>=nap)) return mat(dimension);  // shouldn't happen
  Quad p = primelist[i];
  if(verbose) 
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p<<")...";
  return heckeop(p,dual,verb); // Automatically chooses W or T
}

vec homspace::opmat_col(int i, int j, int verb)
{
  if((i<0)||(i>=nap)) return vec(dimension);  // shouldn't happen
  Quad p = primelist[i];
  if(verbose) 
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p<<")...";
  return heckeop_col(p,j,verb); // Automatically chooses W or T
}

mat homspace::opmat_cols(int i, const vec& jlist, int verb)
{
  if((i<0)||(i>=nap)) return mat(dimension);  // shouldn't happen
  Quad p = primelist[i];
  if(verbose) 
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p<<")...";
  return heckeop_cols(p,jlist,verb); // Automatically chooses W or T
}

smat homspace::s_opmat_cols(int i, const vec& jlist, int verb)
{
  if((i<0)||(i>=nap)) return smat(dimension);  // shouldn't happen
  Quad p = primelist[i];
  if(verbose)
    cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p<<")..."<<flush;
  return s_heckeop_cols(p,jlist,verb); // Automatically chooses W or T
}

mat homspace::opmat_restricted(int i, const subspace& s, int dual, int verb)
{
  if((i<0)||(i>=nap)) return mat(dim(s));  // shouldn't happen
  Quad p = primelist[i];
  if(verbose) 
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p
	  <<") restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
  return heckeop_restricted(p,s,dual,verb); // Automatically chooses W or T
}

smat homspace::s_opmat(int i, int dual, int v)
{
  //  if(i==-1) return s_conj(dual,v);
  if((i<0)||(i>=nap)) 
    {
      return smat(dimension);  // shouldn't happen
    }
  Quad p = primelist[i];
  if(v) 
    {
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p<<")...";
      smat ans = s_heckeop(p,dual,0); // Automatically chooses W or T
      cout<<"done."<<endl;
      return ans;
    }
  else return s_heckeop(p,dual,0); // Automatically chooses W or T
}

smat homspace::s_opmat_restricted(int i, const ssubspace& s, int dual, int v)
{
  if((i<0)||(i>=nap)) 
    {
      return smat(dim(s));  // shouldn't happen
    }
  Quad p = primelist[i];
  if(v) 
    {
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p
	  <<") restricted to subspace of dimension "<<dim(s)<<" ..."<<flush;
      smat ans = s_heckeop_restricted(p,s,dual,0); // Automatically chooses W or T
      cout<<"done."<<endl;
      return ans;
    }
  else return s_heckeop_restricted(p,s,dual,0); // Automatically chooses W or T
}

vector<long> homspace::eigrange(long i)  // implementing virtal function in matmaker
{
  vector<long> ans;
  if((i<0)||(i>=nap)) return ans;  // shouldn't happen
  Quad p = primelist[i];
  long normp = quadnorm(p);
  if (verbose) 
    cout << "eigrange for p = " << p << ":\t";
  if(div(p,modulus))
    {
      vector<long> ans(2);
      ans[0]=-1;
      ans[1]=1;
      if (verbose) 
	cout << ans << endl;
      return ans;
    }
  else
    {
      long aplim=2;
      while (aplim*aplim<=4*normp) aplim++;
      aplim--;
      if(verbose)
	cout << "|ap| up to "<<aplim<<":\t";
      long ap, l = 2*aplim+1;
      vector<long> ans(l);
      ans[0]=0;
      for(ap=-aplim; ap<=aplim; ap++)
	ans[ap+aplim] = ap;
      if (verbose) 
	cout << ans << endl;
      return ans;
    }
}

vec homspace::maninvector(const Quad& p) const
{
  vector<Quad> resmodp=residues(p);
  vec ans = chain(0,p), part;             // =0, but sets the right length.
  vector<Quad>::const_iterator res=resmodp.begin();
  while(res!=resmodp.end())
    {
      part = chain(*res++,p);
      if(hmod)
        ans.addmodp(part,hmod);
      else
        ans += part;
    }
  if(hmod) ans=reduce_modp(ans,hmod);
  return ans;
}

vec homspace::manintwist(const Quad& lambda, const vector<Quad>& res, vector<int> chitable) const
{
  vec ans = chain(0,lambda), part;          // =0, but sets the right length.
  vector<int>::const_iterator chi=chitable.begin();
  vector<Quad>::const_iterator r=res.begin();
  while(r!=res.end())
   {
     part = (*chi++)*chain(*r++,lambda);
      if(hmod)
        ans.addmodp(part,hmod);
      else
        ans += part;
   }
  if(hmod) ans=reduce_modp(ans,hmod);
 return ans;
}

vec homspace::projmaninvector(const Quad& p) const    // Will only work after "proj"
{
  vector<Quad> resmodp=residues(p);
  vec ans = projcycle(0,p), part;         // =0, but sets the right length.
  vector<Quad>::const_iterator res=resmodp.begin();
  while(res!=resmodp.end())
    {
      part = projcycle(*res++,p);
      if(hmod)
        ans.addmodp(part,hmod);
      else
        ans += part;
    }
  if(hmod) ans=reduce_modp(ans,hmod);
  return ans;
}

vec homspace::newhecke(const Quad& p, const Quad& n, const Quad& d) const
                                     // Will only work after "proj"
{ 
  vec ans = projcycle(p*n,d), part;
  vector<Quad> resmodp=residues(p);  Quad dp = d*p;
  vector<Quad>::const_iterator res=resmodp.begin();
  while(res!=resmodp.end())
    {
      part = projcycle(n+d*(*res++),dp);
      if(hmod)
        ans.addmodp(part,hmod);
      else
        ans += part;
    }
  if(hmod) ans=reduce_modp(ans,hmod);
  return ans;
}

#if(USE_SMATS)
void mergeposval(long* pos, long* val, long& npos, long f)
{
  if(f==0) return;
  long af=abs(f), sf=sign(f), i=0,j;
  while((pos[i]<af)&&(i<npos)) i++;
  if(i==npos) //run off end
    {
      pos[npos]=af; val[npos]=sf; 
      npos++;
      return;
    }
  if(pos[i]==af) // change existing entry
    {
      val[i]+=sf;
      if(val[i]!=0) return;  // else close up
      for(j=i; j<npos; j++) {pos[j]=pos[j+1]; val[j]=val[j+1];}
      npos--;
      return;
    }
  // else pos[i-1] < af < pos[i], so insert
  for(j=npos; j>i; j--) {pos[j]=pos[j-1]; val[j]=val[j-1];}
  npos++;
  pos[i]=af; val[i]=sf;
  return;
}
#endif


#define DEBUG_CHINESE

int liftmats_chinese(const smat& m1, scalar pr1, const smat& m2, scalar pr2,
                     smat& m, scalar& dd)
{
  long modulus=(long)pr1*(long)pr2,n,d,mij;
  long nr,nc,u,v;
  float lim=floor(sqrt(modulus/2.0));

  dd = bezout(pr1,pr2,u,v); //==1
  if (dd!=1) return 0;

  // First time through: compute CRTs, common denominator and success flag
  m = m1; // NB We assume that m1 and m2 have nonzero entries in the same places
  for(nr=0; nr<m1.nro; nr++)
    for(nc=0; nc<m1.col[nr][0]; nc++)
      {
        mij = mod(v*m1.val[nr][nc],pr1)*pr2 + mod(u*m2.val[nr][nc],pr2)*pr1;
        mij = mod(mij,modulus);
#ifdef DEBUG_CHINESE
        if (((mij-m1.val[nr][nc])%pr1)||((mij-m2.val[nr][nc])%pr2))
          {
            cout<< "bad CRT(["<<m1.val[nr][nc]<<","<<m2.val[nr][nc]<<"],["<<pr1<<","<<pr2<<"]) = "<<mij<<endl;
          }
#endif
        m.val[nr][nc] = mij;
	if (modrat(mij,modulus,lim,n,d))
          dd=lcm(d,dd);
        else
          {
#ifdef DEBUG_CHINESE
            cout<<"CRT("<<m1.val[nr][nc]<<","<<m2.val[nr][nc]<<")="<<mij<<" (mod "<<modulus<<") fails to lift (lim="<<lim<<")\n";
            cout << "Problems encountered in chinese lifting of smat modulo "<<pr1<<" and "<<pr2<< endl;
#endif
            return 0;
          }
      }
  dd=abs(dd);
#ifdef DEBUG_CHINESE
  cout << "Common denominator = " << dd << "\n";
#endif
  // Second time through: rescale
  for(nr=0; nr<m.nro; nr++)
    for(nc=0; nc<m.col[nr][0]; nc++)
      {
        m.val[nr][nc] = mod(xmodmul((dd/d),(long)m.val[nr][nc],modulus),modulus);
      }
  return 1;
}

vec reduce_modp(const vec& v, const scalar& p)
{
  int i, d=dim(v);
  scalar ai, p2 = p>>1;
  vec ans(d);
  for(i=1; i<=dim(v); i++)
    {
      ai = v[i]%p;
      while( ai>p2) ai-=p;
      while(-ai>p2) ai+=p;
      ans[i] = ai;
    }
  return ans;
}

mat reduce_modp(const mat& m, const scalar& p)
{
  int i, j, nr=m.nrows(), nc=m.ncols();
  scalar aij, p2 = p>>1;
  mat ans(nr,nc);
  for(i=1; i<=nr; i++)
    for(j=1; j<=nc; j++)
      {
        aij = m(i,j)%p;
        while( aij>p2) aij-=p;
        while(-aij>p2) aij+=p;
        ans(i,j) = aij;
      }
  return ans;
}
