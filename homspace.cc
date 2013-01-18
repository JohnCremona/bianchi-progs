// FILE HOMSPACE.CC: Implemention of class homspace

//#define USE_SMATS   // define in Makefile to use sparse matrix methods

#include <eclib/msubspace.h>
#include "homspace.h"
const string W_opname("W");
const string T_opname("T");

void homspace::userel(vec& rel)
{
  long h = vecgcd(rel);
  if (h)
    {  
      if (h>1) rel/=h;
      if(numrel<maxnumrel)
        relmat.setrow(++numrel,rel);
      else 
        cerr<<"Too many relations! numrel="<<numrel
	    <<", maxnumrel="<<maxnumrel<<endl;
    }
}
 
homspace::homspace(const Quad& n, int hp, int verb) :symbdata(n)
{
  verbose=verb;
  cuspidal=0;
  if (verbose) symbdata::display();
  long field = Quad::d;
  plusflag=hp;                  // Sets static level::plusflag = hp
  long i,j,k,ngens=0;
  coordindex = new int[nsymb];
  int* check = new int[nsymb];
  int* gens = new int[nsymb+1];
  //NB start of gens array is at 1 not 0

// 2-term (edge) relations:

  Quad unit = fundunit;
  long lenrel = Quad::nunits;
  if(!plusflag) {unit=fundunit*fundunit; lenrel/=2;}
  symbop eps(this,unit,0,0,1);
  symbop sof(this,0,-1,1,0);
  int *a=new int[lenrel], *b=new int[lenrel],  triv;
  i=nsymb; while (i--) check[i]=0;
  if(verbose) cout << "About to start on 2-term relations.\n";
  for (j=0; j<nsymb; j++)  
    {
      if (check[j]==0)
        { 
	  if(verbose>1) cout << "j = " << j << ":\t";
          a[0]=j; b[0]=sof(j); triv=0;
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
	      for (k=0; k<lenrel; k++) cout<<a[k]<<" "; cout<<endl;
	      cout<<"\t-:\t";
	      for (k=0; k<lenrel; k++) cout<<b[k]<<" "; cout<<endl;
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

// end of 2-term relations

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
  delete[] a; delete[] b;
   
//
// face relations
//
  if(plusflag)  maxnumrel = 2*ngens+10;
  else          maxnumrel = 3*ngens+10;

  maxnumrel=2*nsymb;

#ifdef USE_SMATS
   smat relmat(maxnumrel,ngens);
   svec newrel(ngens);
#else
  relmat.init(maxnumrel,ngens);
  vec newrel(ngens);
#endif
  numrel = 0;
  int * rel = new int[6];  //max length
  long ij, fix;
//
// first the 3-term relation for all fields:
//
  if(verbose)
    {
      cout << "Face relation type 1 (triangles):\n";
    }
  symbop tof(this,1,1,-1,0);
  symbop rof(this,0,1,1,0);
  i=nsymb; while (i--) check[i]=0;
  for (k=0; k<nsymb; k++) if (check[k]==0)
    {
#ifdef USE_SMATS
      newrel.clear();
#else
      for (j=1; j<=ngens; j++) newrel[j]=0;
#endif
      rel[2]=tof(rel[1]=tof(rel[0]=k));
      if (verbose)   
//      cout << "Relation: " << rel[0]<<" "<<rel[1]<<" "<<rel[2]<<" --> ";
        cout<<"Relation: "<<rel[0]<<symbol(rel[0])<<" "<<rel[1]<<symbol(rel[1])<<" "<<rel[2]<<symbol(rel[2])<<" --> ";
      if(tof(rel[2])!=k)
	cout<<"Error in TS-relation!"<<endl;
      for (j=0; j<3; j++)
        {
          check[ij=rel[j]] = 1;
	  if (plusflag) check[rof(ij)] = 1;
          fix = coordindex[ij];
          if (verbose)  cout << fix << " ";
#ifdef USE_SMATS
	  if(fix) newrel.add(abs(fix),(fix>0?1:-1));
#else
          if (fix!=0) newrel[abs(fix)] += sign(fix);
#endif
        }
      if (verbose)  cout << endl;
#ifdef USE_SMATS
     if(newrel.size()!=0) 
       {
	 numrel++;
	 make_primitive(newrel);
	 if(numrel<=maxnumrel)
	   relmat.setrow(numrel,newrel);
	 else 
	   cout<<"Too many 3-term relations (numrel="<<numrel
	       <<", maxnumrel="<<maxnumrel<<")"<<endl;
       }
#else
      userel(newrel);
#endif
    }

if(verbose)
  {
    cout << "After face relation type 1, number of relations = " << numrel <<"\n";
    cout << "Face relation type 2:\n";
  }

//
// Now one extra depending on the field:
//
  Quad w(0,1);
  
  //Define extra ops needed:
  // field 1:
  symbop x1of(this,w,1,1,0);      
  // field 3:
  symbop x3of(this,1,w,w-1,0);
  // field 2:
  symbop uof(this,w,1,1,0);
  // field 7:
  symbop yof(this,1,-w,1-w,-1);
  symbop us7of(this,w,-1,1,0);
  // field 11:
  symbop xof(this,1,-w,1-w,-2);
  symbop us11of(this,w,-1,1,0);

  switch (field) 
    {
    case 1: case 3:
      if(plusflag) break;
      i=nsymb; while (i--) check[i]=0;
      for (k=0; k<nsymb; k++) if (check[k]==0)
        {
#ifdef USE_SMATS
	  newrel.clear();
#else
          for (j=1; j<=ngens; j++) newrel[j]=0;
#endif
          if(field==1) rel[2]=x1of(rel[1]=x1of(rel[0]=k));
          else         rel[2]=x3of(rel[1]=x3of(rel[0]=k));
          if (verbose)   
	    //  cout<<"Relation: "<<rel[0]<<" "<<rel[1]<<" "<<rel[2]<<" --> ";
            cout<<"Relation: "<<symbol(rel[0])<<" "<<symbol(rel[1])<<" "<<symbol(rel[2])<<" --> ";
          for (j=0; j<3; j++)
            {
              check[ij=rel[j]] = 1;
              fix = coordindex[ij];
              if (verbose)  cout << fix << " ";
#ifdef USE_SMATS
	      if(fix) newrel.add(abs(fix),(fix>0?1:-1));
#else
              if (fix!=0) newrel[abs(fix)] += sign(fix);
#endif
            }
          if (verbose)  cout << endl;
#ifdef USE_SMATS
	  if(newrel.size()!=0) 
	    {
	      numrel++;
	      make_primitive(newrel);
	      if(numrel<=maxnumrel)
		relmat.setrow(numrel,newrel);
	      else 
		cout<<"Too many n-term relations (numrel="<<numrel
		    <<", maxnumrel="<<maxnumrel<<")"<<endl;
	    }
#else
          userel(newrel);
#endif
        }
      break;
    case 2:
      i=nsymb; while (i--) check[i]=0;
      for (k=0; k<nsymb; k++) if (check[k]==0)
        {
#ifdef USE_SMATS
	  newrel.clear();
#else
          for (j=1; j<=ngens; j++) newrel[j]=0;
#endif
          rel[3]=uof(rel[2]=uof(rel[1]=uof(rel[0]=k)));
          if (verbose)   
            cout<<"Relation: "<<rel[0]<<" "<<rel[1]<<" "<<rel[2]<<" "<<rel[3]<<" --> ";
          for (j=0; j<4; j++)
            {
              check[ij=rel[j]] = 1;
              if (plusflag) check[sof(ij)] = 1;
              fix = coordindex[ij];
              if (verbose)  cout << fix << " ";
#ifdef USE_SMATS
	      if(fix) newrel.add(abs(fix),(fix>0?1:-1));
#else
              if (fix!=0) newrel[abs(fix)] += sign(fix);
#endif
            }
          if (verbose)  cout << endl;
#ifdef USE_SMATS
	  if(newrel.size()!=0) 
	    {
	      numrel++;
	      make_primitive(newrel);
	      if(numrel<=maxnumrel)
		relmat.setrow(numrel,newrel);
	      else 
		cout<<"Too many n-term relations (numrel="<<numrel
		    <<", maxnumrel="<<maxnumrel<<")"<<endl;
	    }
#else
          userel(newrel);
#endif
        }
      break;
    case 7:
      i=nsymb; while (i--) check[i]=0;
      for (k=0; k<nsymb; k++) if (check[k]==0)
        {
#ifdef USE_SMATS
	  newrel.clear();
#else
          for (j=1; j<=ngens; j++) newrel[j]=0;
#endif
          rel[2]=yof(rel[0]=k);
          rel[1]=us7of(rel[0]); 
          rel[3]=us7of(rel[2]);
          if (verbose)   
            cout<<"Relation: "<<rel[0]<<" "<<rel[1]<<" "<<rel[2]<<
                  " "<<rel[3]<<" --> ";
          check[k]=check[rel[2]]=1;
          check[rof(rel[1])]=check[rof(rel[3])]=1;
          for (j=0; j<4; j++)
            { 
              fix = coordindex[rel[j]];
              if (verbose)  cout << fix << " ";
#ifdef USE_SMATS
	      if(fix) newrel.add(abs(fix),(fix>0?1:-1));
#else
              if (fix!=0) newrel[abs(fix)] += sign(fix);
#endif
            }
          if (verbose)  cout << endl;
#ifdef USE_SMATS
	  if(newrel.size()!=0) 
	    {
	      numrel++;
	      make_primitive(newrel);
	      if(numrel<=maxnumrel)
		relmat.setrow(numrel,newrel);
	      else 
		cout<<"Too many 3-term relations (numrel="<<numrel
		    <<", maxnumrel="<<maxnumrel<<")"<<endl;
	    }
#else
          userel(newrel);
#endif
        }
      break;
    case 11:
      i=nsymb; while (i--) check[i]=0;
      for (k=0; k<nsymb; k++) if (check[k]==0)
        {
#ifdef USE_SMATS
	  newrel.clear();
#else
          for (j=1; j<=ngens; j++) newrel[j]=0;
#endif
          rel[4]=xof(rel[2]=xof(rel[0]=k));
          rel[1]=us11of(rel[0]); 
          rel[3]=us11of(rel[2]); 
          rel[5]=us11of(rel[4]);
          if (verbose)   
            cout<<"Relation: "<<rel[0]<<" "<<rel[1]<<" "<<rel[2]<<
                  " "<<rel[3]<<" "<<rel[4]<<" "<<rel[5]<<" --> ";
          check[k]=check[rel[2]]=check[rel[4]]=1;
          check[rof(rel[1])]=check[rof(rel[3])]=check[rof(rel[5])]=1;
          for (j=0; j<6; j++)
            { 
              fix = coordindex[rel[j]];
              if (verbose)  cout << fix << " ";
#ifdef USE_SMATS
	      if(fix) newrel.add(abs(fix),(fix>0?1:-1));
#else
              if (fix!=0) newrel[abs(fix)] += sign(fix);
#endif
            }
          if (verbose)  cout << endl;
#ifdef USE_SMATS
	  if(newrel.size()!=0) 
	    {
	      numrel++;
	      make_primitive(newrel);
	      if(numrel<=maxnumrel)
		relmat.setrow(numrel,newrel);
	      else 
		cout<<"Too many n-term relations (numrel="<<numrel
		    <<", maxnumrel="<<maxnumrel<<")"<<endl;
	    }
#else
          userel(newrel);
#endif
        }
      break;
    default: cerr<<"homspace not implemented for field "<<field<<endl;
      exit(1);
    }

  if(verbose) 
    {
      cout << "Finished face relations: ";
      cout << "number of relations = " << numrel;
      cout << " (bound was "<<maxnumrel<<")"<<endl;
#ifdef USE_SMATS
      cout << "relmat = \n" << relmat<<endl;
#endif
    }

   vec pivs, npivs;
#ifdef USE_SMATS
   smat_elim sme(relmat);
   relmat=smat(0,0);
   int d1;
   smat sp = liftmat(sme.kernel(npivs,pivs),MODULUS,d1);
   denom1=d1;
   if(verbose>1) 
     {
       cout << "kernel of relmat = " << sp.as_mat() << endl;
       cout << "pivots = "<<pivs << endl;
       cout << "denom = "<<d1 << endl;
     }
   rk = ncols(sp);
   coord.init(ngens+1,rk); // 0'th is unused
   for(i=1; i<=ngens; i++) 
     coord.setrow(i,sp.row(i).as_vec());
   sp=smat(0,0); // clear space
#else
  relmat = relmat.slice(numrel,ngens);
  if(verbose) 
    {
      cout << "relmat = "; relmat.output_pari(cout); cout << endl;
    }
  msubspace sp = kernel(mat_m(relmat),0);
  rk = dim(sp);
  coord = basis(sp).shorten((int)1);
  pivs = pivots(sp);
  denom1 = I2int(denom(sp));
  relmat.init(); newrel.init(); sp.clear(); 
#endif
  if (verbose) 
    {
      cout << "rk = " << rk << endl;
      cout << "coord:" << coord;
      if(denom1!=1) cout << "denominator = " << denom1 << endl;
      cout << "pivots = " << pivs <<endl;
    }
  if (rk>0)
    {
      freegens = new int[rk]; if (!freegens) exit(1);
      for (i=0; i<rk; i++) freegens[i] = gens[pivs[i+1]];
      if (verbose)
        { 
          cout << "freegens: ";
          for (i=0; i<rk; i++) cout << freegens[i] << " ";
          cout << endl;
        }
    }
  delete [] check; delete [] gens; delete [] rel; 
  pivs.init();
  {
     cusplist cusps(2*rk);
     mat deltamat(2*rk,rk);
     for (i=0; i<rk; i++)
       {
         modsym m(symbol(freegens[i]));
         for (j=1; j>-3; j-=2)
           {
             cusp c = (j==1 ? m.beta() : m.alpha());
             k = cusps.index(c);   //adds automatically if new
             deltamat(k+1,i+1) += j;  // N.B. offset of 1
           }
       }
     ncusps=cusps.count();
     if(verbose)cout << "ncusps = " << ncusps << endl;
     kern = kernel(smat(deltamat));
     vec pivs, npivs;
     int d2;
     smat sk = liftmat(smat_elim(deltamat).kernel(npivs,pivs),MODULUS,d2);
     denom2=d2;
  }
   if(cuspidal)
     dimension = dim(kern);
   else
     dimension = rk;
  denom3 = denom1 * denom2;
  freemods = new modsym[rk]; if (!freemods) exit(1);
  needed   = new int[rk];    if (!needed) exit(1);
  
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

vec homspace::chain(const symb& s) const  //=old getcoord
{
 long i= coordindex[index(s)];
 vec ans(ncols(coord));
 if (i) ans = sign(i)*(coord.row(abs(i)));
 return ans;
}

vec homspace::chaincd(const Quad& c, const Quad& d) const //=old getcoord2
{
 long i= coordindex[index2(c,d)];
 vec ans(ncols(coord));
 if (i) ans = sign(i)*(coord.row(abs(i)));
 return ans;
}

vec homspace::projchaincd(const Quad& c, const Quad& d) const 
{
 long i= coordindex[index2(c,d)];
 vec ans(ncols(projcoord));
 if (i) ans = sign(i)*(projcoord.row(abs(i)));
 return ans;
}

vec homspace::chain(const Quad& nn, const Quad& dd) const //=old qtovec2
{
   vec ans = chaincd(0,1);
   Quad c=0, d=1, e, a=nn, b=dd, q, f;
   while (b!=0)
   { q=a/b; 
     f=a; a=-b; b= f-q*b; 
     e=d; d= c; c= e+q*c; c=-c;
     ans += chaincd(c,d);
   }
   return ans;
}

vec homspace::projcycle(const Quad& nn, const Quad& dd) const  //
{
   vec ans = projchaincd(0,1);
   Quad c=0, d=1, e, a=nn, b=dd, q, f;
   while (b!=0)
   { q=a/b; 
     f=a; a=-b; b= f-q*b; 
     e=d; d= c; c= e+q*c; c=-c;
     ans += projchaincd(c,d);
   }
   return ans;
}

vec homspace::applyop(const matop& mlist, const RatQuad& q) const
{ vec ans(rk);  long i=mlist.length;
  while (i--) ans += chain(mlist[i](q));
  return ans;
}
 
mat homspace::calcop(const string opname, const Quad& p, const matop& mlist, int dual, int display) const
{
  mat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     { vec colj = applyop(mlist,freemods[j]);
       m.setcol(j+1,colj);
     }
  if(cuspidal) m = restrict_mat(smat(m),kern).as_mat();
  if(dual) {  m=transpose(m);}
  if (display) cout << "Matrix of " << opname << "(" << p << ") = " << m;
  if (display && (dimension>1)) cout << endl;
  return m;
}
 
smat homspace::s_calcop(const string  opname, const Quad& p, const matop& mlist, 
			int dual, int display) const
{
  smat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     { svec colj = applyop(mlist,freemods[j]);
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
       svec colj = applyop(mlist,freemods[jj]);
       m.setrow(j+1,colj.as_vec());
     }
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
       m.setrow(j,colj);
     }
  m = mult_mod_p(m,basis(s),MODULUS);
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

mat homspace::opmat(long i, int dual, int verb)
{
  if((i<0)||(i>=nap)) return mat(dimension);  // shouldn't happen
  Quad p = primelist[i];
  if(verbose) 
      cout<<"Computing " << ((div(p,modulus)) ? W_opname : T_opname) <<"("<<p<<")...";
  return heckeop(p,dual,verb); // Automatically chooses W or T
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
      ans[0]=1;
      ans[1]=-1;
      if (verbose) 
	cout << ans << endl;
      return ans;
    }
  else
    {
      long aplim=2;
      while (aplim*aplim<=4*normp) aplim++; aplim--;
      if(verbose)
	cout << "|ap| up to "<<aplim<<":\t";
      long ap, j, l = 2*aplim+1;
      vector<long> ans(l);
      ans[0]=0;
      for(ap=1, j=1; ap<=aplim; ap++)
	{
	  ans[j++] = ap;
	  ans[j++] = -ap;
	}
      if (verbose) 
	cout << ans << endl;
      return ans;
    }
}

vec homspace::maninvector(const Quad& p) const
{
  vector<Quad> resmodp=residues(p); 
  vec tvec = chain(0,p);             // =0, but sets the right length.
  vector<Quad>::const_iterator res=resmodp.begin();
  while(res!=resmodp.end())
    tvec += chain(*res++,p);
  //  return kernelpart(tvec);
  return tvec;
}

vec homspace::manintwist(const Quad& lambda, const vector<Quad>& res, int* chitable) const
{
 vec sum = chain(0,lambda);          // =0, but sets the right length.
 int *chi=chitable;
 vector<Quad>::const_iterator r=res.begin();
 while(r!=res.end())
   sum += (*chi++)*chain(*r++,lambda);
 // return kernelpart(sum);
 return sum;
}

vec homspace::projmaninvector(const Quad& p) const    // Will only work after "proj"
{
  vector<Quad> resmodp=residues(p); 
  vec tvec = projcycle(0,p);             // =0, but sets the right length.
  vector<Quad>::const_iterator res=resmodp.begin();
  while(res!=resmodp.end())
    tvec += projcycle(*res++,p);
  return tvec;
}

vec homspace::newhecke(const Quad& p, const Quad& n, const Quad& d) const
                                     // Will only work after "proj"
{ 
  vec tvec = projcycle(p*n,d);
  vector<Quad> resmodp=residues(p);  Quad dp = d*p;
  vector<Quad>::const_iterator res=resmodp.begin();
  while(res!=resmodp.end())
    tvec += projcycle(n+d*(*res++),dp);
  return tvec;
}

#ifdef USE_SMATS
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

#ifndef USE_XSPLIT

int startswith(const vector<long>& a, const vector<long>& b, long l)
{
  int ans=1;
  for (int i=0; (i<l)&&ans; i++) ans = (a(i)==b(i));
  return ans;
}

#define MODULUS 92681

splitter::splitter(homspace* h, long d, int dualflag, int meth, int v)
:h1(h), maxdepth(d), dual(dualflag), depth(0), method(meth), verbose(v)
{
  h1denom = h1->h1denom();
  nest = new subspace[d];
  nest[0] = subspace(h->h1dim());
  aplist = vector<long>(d);
  havemat = new int[d]; long i=d; while(i--) havemat[i]=0;
  opmats = new mat[d];
  vector<Quad>::const_iterator pr=quadprimes.begin();
  for(i=0; i<d; i++) 
    {
      if(i<level::npdivs) plist.push_back(level::plist(i));
      else 
        {
          while(div(*pr,level::modulus)) pr++; 
          plist.push_back(*pr++); 
        }
    }
}

void splitter::splitoff(const vector<long>& apl)
{
  long targetdim= 1;
  if(verbose)cout<<"Entering splitter, old depth = "<<depth<<endl;
  if(verbose&&(depth>0))cout<<"Old aplist = "<<aplist<<endl;
  if(verbose)cout<<"New aplist = "<<apl<<endl;
  while(!startswith(aplist,apl,depth)) nest[depth--]=subspace(1);
  if(verbose)cout<<"Starting depth = "<<depth<<endl;
  while((dim(nest[depth])>targetdim) && (depth<maxdepth))
    {
      if(verbose)cout<<"Increasing depth to "<<depth+1<<endl;
      aplist[depth]=apl(depth);
      if(!havemat[depth])
        {
          if(depth<level::npdivs)
            {
            if(verbose)cout<<"Computing W("<<plist[depth]<<")"<<endl;
            opmats[depth] = h1->wop(plist[depth]);
          }
          else
            {
            if(verbose)cout<<"Computing T("<<plist[depth]<<")"<<endl;
            opmats[depth] = h1->heckeop(plist[depth]);
          }
          if(dual)opmats[depth]=transpose(opmats[depth]);
          havemat[depth]=1;
        }
      if(method<3)
        nest[depth+1] = subeigenspace(opmats[depth],h1denom*aplist[depth],nest[depth],method);
      else
        nest[depth+1] = psubeigenspace(opmats[depth],h1denom*aplist[depth],nest[depth],MODULUS);
      depth++;
    }
  long finaldim=dim(nest[depth]);
  if(finaldim!=targetdim)
    {
      cerr<<"error in splitter::splitoff with aplist = "<<apl<<endl;
      cerr<<"final subspace has dimension "<<finaldim<<endl;
    }
  subspace s = nest[depth];
  basisvector = basis(s).col(1);
  if(method==3)
    basisvector = lift(basisvector,MODULUS);
  else 
    makeprimitive(basisvector);
}

#endif // ifndef USE_XSPLIT

