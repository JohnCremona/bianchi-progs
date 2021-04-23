// FILE homspace.cc: Implemention of class homspace

#include "homspace.h"

// define method for kernel procedures:
// 0 = exact, starts to report errors from about level (3[9,2+a]), norm 81
// 2 = modular, works;
#define KER_RELMAT_DEFAULT_METHOD 2
#undef KER_RELMAT_RUNTIME_PROMPT

#undef USE_OBS_ALGOL_RELNS
#undef USE_REDUNDANT_TRIANGLES

void rel_matrix::userel(vec& rel)
{
  long h = vecgcd(rel);
  if (h)
    {  
      if (h>1) rel/=h;
      if(numrel<maxnumrel)
        setrow(++numrel,rel);
      else 
        cerr<<"Too many relations!"<<endl;
    }
}

homspace::homspace(const Qideal& n, int hp, int verb) :symbdata(n)
{
  GL_vs_SL_flag = hp;                  // was:  static level::plusflag
  verbose=verb;
  if (verbose) symbdata::display(cout, verbose);

  switch (Field::d) {
  case 2: case 19: case 5:
    { if (verbose) cout << "Constructing homspace...\n";}
    break;
  default:
    { cerr<<"homspace not implemented for field "<<Field::d<<endl;
      exit(1);
    }
  }
  coordindex = new long[nsymb*geometry::nt()];
  // for debugging: to ensure that coordindex is assigned ok, init to daft val!
  // for (long pch=0; pch<nsymb*geometry::nt(); pch++) coordindex[pch]=-99;

  gens = new long[nsymb*geometry::nt()+1];
  //NB start of gens array is at 1 not 0  -- of course, usually way too long!

  relate_two_term();
  relate_faces();
  restrict_to_kerdelta();
  if (verbose) cout << "Finished constructing homspace.\n";
}

void homspace::relate_two_term()
{
  int* check = new int[nsymb];
  long i,j,k,l,m;
  ngens=0;

// 2-term (edge) relations:
  if(verbose) cout << endl << "About to start on 2-term relations.\n";

  symbop sof(this,0,-1,1,0);

  switch (Field::d)
    {
    case 2:
{
  symbop eps(this,-1,0,0,1);
  symbop eps2(this,1,0,0,1);
  long lenrel = (GL_vs_SL_flag) ? 2 : 1;
  long *a=new long[lenrel], *b=new long[lenrel]; int triv;
  i=nsymb; while (i--) check[i]=0;

  for (j=0; j<nsymb; j++)  
    {
      if (check[j]==0)
        { // if(verbose) cout << "j = " << j << "\n";
          a[0]=j; b[0]=sof(j); triv=0;
          for(k=1; k<lenrel; k++)
            {
              a[k]= (GL_vs_SL_flag) ? eps(a[k-1]) : eps2(a[k-1]); 
              b[k]= (GL_vs_SL_flag) ? eps(b[k-1]) : eps2(b[k-1]);
              triv= triv | (j==b[k]);
            }
          for (k=0; k<lenrel; k++) check[a[k]]=check[b[k]]=1;
          if (triv)
            for (k=0; k<lenrel; k++)
	      { coordindex[tysofs(a[k],0)]=0;
		coordindex[tysofs(b[k],0)]=0;
	      }
          else
            {   
              gens[++ngens] = tysofs(j,0);
              for(k=0; k<lenrel; k++)
                {
                  coordindex[tysofs(a[k],0)] =  ngens;
                  coordindex[tysofs(b[k],0)] = -ngens;
                }
            }
        }
    }
}

      break;   // end of edge-rels for d=2

    case 19:

// "type 0"
  {
    symbop jof(this,-1,0,0,1);
    //  symbop sof(this,0,-1,1,0);

    i=nsymb; while (i--) check[i]=0;
    for (j=0; j<nsymb; j++)  
      {
	if (check[j]==0)
	  { // if(verbose) cout << "j = " << j << "\n";
	    k=jof(j);
	    l=sof(j);
	    m=sof(k);
	    
	    // debug test
	    if (m != jof(l))
	      { cerr <<"ERROR: Unexpected V4 failure in type 0 2-term!"<<endl;}
	    
	    check[j]=check[k]=check[l]=check[m]=1;

	    // debug test
	    if ( ((j==l)&&(k!=m)) || ((j==m)&&(k!=l)) )
	      { cerr <<"Warning: unexpected result in relate_two_term!"<<endl;}

#ifndef USE_OBS_ALGOL_RELNS
	    if ( (j==l)||(j==m) )
#else
	      if ( (j==l)||(k==m) )     // obs as in algol prog
#endif
		{
		  coordindex[tysofs(j,0)]=0;
		  coordindex[tysofs(k,0)]=0;
		  coordindex[tysofs(l,0)]=0;
		  coordindex[tysofs(m,0)]=0;
		}
	      else
		{   
		  gens[++ngens] = tysofs(j,0);
		  coordindex[tysofs(j,0)]=ngens;
		  coordindex[tysofs(k,0)]=ngens;
		  coordindex[tysofs(l,0)]=-ngens;
		  coordindex[tysofs(m,0)]=-ngens;
		}
	  }
      }
  }

// "type 1"
  {
    symbop kof(this,-1,Quad(0,1),0,1);
    i=nsymb; while (i--) check[i]=0;
    for (j=0; j<nsymb; j++)  
      {
	if (check[j]==0)
	  { // if(verbose) cout << "j = " << j << "\n";
	    k=kof(j);
	    
	    check[j]=check[k]=1;
	    
	    gens[++ngens] = tysofs(j,1);
	    coordindex[tysofs(j,1)]=ngens;
	    coordindex[tysofs(k,1)]=ngens;
	  }
      }
  }

      break;   // end of edge-rels for d=19

    case 5:
{
//  symbop eps(this,fundunit,0,0,1);
//  symbop eps2(this,fundunit*fundunit,0,0,1);
  symbop jof(this,-1,0,0,1);
//  symbop sof(this,0,-1,1,0);
  symbop bmatrix(this, -5, Quad(2,2), Quad(-2,2), 5);
  symbop ematrix(this, -2, Quad(0,1), Quad(0,1), 2);
  symbop amatrix(this, -1, Quad(1,1), 0, 1);

// "type 0"
  i=nsymb; while (i--) check[i]=0;
  for (j=0; j<nsymb; j++)  
    {
      if (check[j]==0)
        { // if(verbose) cout << "j = " << j << "\n";
	  k=jof(j);
	  l=sof(j);
	  m=sof(k);

	  // debug test
	  if (m != jof(l))
	    { cerr <<"ERROR: Unexpected V4 failure in type 0 2-term!"<<endl;}
	  
          check[j]=check[k]=check[l]=check[m]=1;

	  // debug test
	  if ( ((j==l)&&(k!=m)) || ((j==m)&&(k!=l)) )
	    { cerr << "Warning: unexpected result in relate_two_term!"<<endl;}

#ifndef USE_OBS_ALGOL_RELNS
	  if ( (j==l)||(j==m) )
#else
	  if ( (j==l)||(k==m) )     // obs as in algol prog
#endif
	    {
	      coordindex[tysofs(j,0)]=0;
	      coordindex[tysofs(k,0)]=0;
	      coordindex[tysofs(l,0)]=0;
	      coordindex[tysofs(m,0)]=0;
	    }
          else
            {   
              gens[++ngens] = tysofs(j,0);
	      coordindex[tysofs(j,0)]=ngens;
	      coordindex[tysofs(k,0)]=ngens;
	      coordindex[tysofs(l,0)]=-ngens;
	      coordindex[tysofs(m,0)]=-ngens;
	    }
        }
    }

// "type 2"
  i=nsymb; while (i--) check[i]=0;
  for (j=0; j<nsymb; j++)  
    {
      if (check[j]==0)
        { // if(verbose) cout << "j = " << j << "\n";
	  k=bmatrix(j);
	  l=ematrix(j);
	  m=ematrix(k);

	  // debug test
	  if (m != bmatrix(l))
	    { cerr <<"ERROR: Unexpected V4 failure in type 2 2-term!"<<endl;}	  
          check[j]=check[k]=check[l]=check[m]=1;

	  // debug test
	  if ( ((j==l)&&(k!=m)) || ((j==m)&&(k!=l)) )
	    { cerr << "Warning: unexpected result in relate_two_term!"<<endl;}

#ifndef USE_OBS_ALGOL_RELNS
	  if ( (j==l)||(j==m) )
#else
	  if ( (j==l)||(k==m) )     // obs as in algol prog
#endif
	    {
	      coordindex[tysofs(j,2)]=0;
	      coordindex[tysofs(k,2)]=0;
	      coordindex[tysofs(l,2)]=0;
	      coordindex[tysofs(m,2)]=0;
	    }
          else
            {   
              gens[++ngens] = tysofs(j,2);
	      coordindex[tysofs(j,2)]=ngens;
	      coordindex[tysofs(k,2)]=ngens;
	      coordindex[tysofs(l,2)]=-ngens;
	      coordindex[tysofs(m,2)]=-ngens;
	    }
        }
    }

// "type 1"
  i=nsymb; while (i--) check[i]=0;
  for (j=0; j<nsymb; j++)  
    {
      if (check[j]==0)
        { // if(verbose) cout << "j = " << j << "\n";
	  k=amatrix(j);
	  
          check[j]=check[k]=1;

	  gens[++ngens] = tysofs(j,1);
	  coordindex[tysofs(j,1)]=ngens;
	  coordindex[tysofs(k,1)]=ngens;
        }
    }
}
      break;  // end of edge-rels for d=5

    default:
      cerr<<"Error: relate_two_term not implemented for field "<<Field::d<<" !";
      exit(1);
    }


  delete[] check;

// end of 2-term relations

  if (verbose)
    {
      cout << "After 2-term relations, ngens = "<<ngens<<endl;
      cout << "gens = ";
      for (i=1; i<=ngens; i++) cout << gens[i] << " ";
      cout << endl;
      cout << "coordindex = ";
      for (i=0; i<nsymb*geometry::nt(); i++) cout << coordindex[i] << " ";
      cout << endl;
    }
}


//
// face relations
//

void homspace::relate_faces()
{
  rel_matrix relmat(geometry::relmat_shape()*ngens,ngens);
                                          // should init to 0 -- enough rows??
  vec newrel(ngens);
  long* rel = new long[geometry::max_length_of_rel()];  //max length
  int* check = new int[nsymb];
  long i,j,k,ij, fix;

  switch (Field::d)
    {
    case 2:
{
//
// first the 3-term relation for all (Euclidean?) fields:
//
  if(verbose)
    {
      cout << "Face relation (I) + (TS) + ((TS)^2) :\n";
    }

  symbop tof(this,1,1,-1,0);
  symbop rof(this,0,1,1,0);
  i=nsymb; while (i--) check[i]=0;
  for (k=0; k<nsymb; k++) if (check[k]==0)
    {
      for (j=1; j<=ngens; j++) newrel[j]=0;
      rel[2]=tof(rel[1]=tof(rel[0]=k));
      if (verbose)   
        cout << "Relation: " << rel[0]<<" "<<rel[1]<<" "<<rel[2]<<" --> ";
      for (j=0; j<3; j++)
        {
          check[ij=rel[j]] = 1;  if (GL_vs_SL_flag) check[rof(ij)] = 1;
          fix = coordindex[tysofs(ij,0)];
          if (verbose)  cout << fix << " ";
          if (fix!=0) newrel[abs(fix)] += sign(fix);
        }
      if (verbose)  cout << endl;
      relmat.userel(newrel);
    }

if(verbose)
  {
    cout << "After face relation type 1, number of relations = " << relmat.get_numrel() <<"\n";
    cout << "Face relation type 2:\n";
  }

//
// Now one extra depending on the field:
//

// ****** NEEDS CORRECTING WHEN  GL_vs_SL_flag  IS NOT EQUAL TO  1  ******

      symbop sof(this,0,-1,1,0);
      symbop uof(this,Quad(0,1),1,1,0);
      i=nsymb; while (i--) check[i]=0;
      for (k=0; k<nsymb; k++) if (check[k]==0)
        {
          for (j=1; j<=ngens; j++) newrel[j]=0;
          rel[3]=uof(rel[2]=uof(rel[1]=uof(rel[0]=k)));
          if (verbose)   
            cout<<"Relation: "<<rel[0]<<" "<<rel[1]<<" "<<rel[2]<<" "<<rel[3]<<" --> ";
          for (j=0; j<4; j++)
            {
              check[ij=rel[j]] = 1;
              if (GL_vs_SL_flag) check[sof(ij)] = 1;
              fix = coordindex[tysofs(ij,0)];
              if (verbose)  cout << fix << " ";
              if (fix!=0) newrel[abs(fix)] += sign(fix);
            }
          if (verbose)  cout << endl;
          relmat.userel(newrel);
        }
}
      break;  // end of face-rels for d=2

    case 19:

// Face relation (I) + (TS) + ((TS)^2)
  {
    if(verbose) cout << "Face relation (I) + (TS) + ((TS)^2) :\n";

    symbop tsof(this,1,-1,1,0);
    
    symbop debug_ts2of(this,0,-1,1,-1);

    i=nsymb; while (i--) check[i]=0;
    for (k=0; k<nsymb; k++) if (check[k]==0)
      {
	for (j=1; j<=ngens; j++) newrel[j]=0;
	rel[2]=tsof(rel[1]=tsof(rel[0]=k));
	
	//debug
	if (rel[2] != debug_ts2of(k))
	  { cerr <<"ERROR: Group action failure in triangle reln!"<<endl; }
	if (rel[0] != tsof(rel[2]))
	  { cerr <<"ERROR: Triangle rotation failure in 3-term reln!"<<endl;}
	
	if (verbose)   
	  cout << "Relation: " << rel[0]<<" "<<rel[1]<<" "<<rel[2]<<" --> ";
	for (j=0; j<3; j++)
	  {
	    check[ij=rel[j]] = 1;
	    fix = coordindex[tysofs(ij,0)];
	    if (verbose)  cout << fix << " ";
	    if (fix!=0) newrel[abs(fix)] += sign(fix);
	  }
	if (verbose)  cout << endl;
	relmat.userel(newrel);
      }
    
    if(verbose)
      {
	cout<< "Current number of relations = " << relmat.get_numrel() <<endl;
      }
  }
  
// Face relation [I] + [N2] + [(N2)^2]
  {
    if(verbose) cout << "Face relation [I] + [N2] + [(N2)^2] :\n";

    symbop n_two_of(this,Quad(1,1),Quad(2,-1),2,Quad(0,-1));

    symbop debug_n_two_2of(this,Quad(0,1),Quad(2,-1),2,Quad(-1,-1));

    i=nsymb; while (i--) check[i]=0;
    for (k=0; k<nsymb; k++) if (check[k]==0)
      {
	for (j=1; j<=ngens; j++) newrel[j]=0;
	rel[2]=n_two_of(rel[1]=n_two_of(rel[0]=k));
      
	//debug
	if (rel[2] != debug_n_two_2of(k))
	  { cerr <<"ERROR: Group action failure in triangle reln!"<<endl; }
	if (rel[0] != n_two_of(rel[2]))
	  { cerr <<"ERROR: Triangle rotation failure in 3-term reln!"<<endl;}
      
	if (verbose)   
	  cout << "Relation: " << rel[0]<<" "<<rel[1]<<" "<<rel[2]<<" --> ";
	for (j=0; j<3; j++)
	  {
	    check[ij=rel[j]] = 1;
	    fix = coordindex[tysofs(ij,1)];
	    if (verbose)  cout << fix << " ";
	    if (fix!=0) newrel[abs(fix)] += sign(fix);
	  }
	if (verbose)  cout << endl;
	relmat.userel(newrel);
      }
  
    if(verbose)
      {
	cout<< "Current number of relations = " << relmat.get_numrel() <<endl;
      }
  }

// Face relation [I] + [N1] + [(N1)^2]
  {
    if(verbose) cout << "Face relation [I] + [N1] + [(N1)^2] :\n";

    symbop sof(this,0,-1,1,0);
    symbop n_one_of(this,Quad(-1,1),2,2,Quad(0,-1));
    
    symbop debug_n_one_2of(this,Quad(0,1),2,2,Quad(1,-1));

    i=nsymb; while (i--) check[i]=0;
    for (k=0; k<nsymb; k++) if (check[k]==0)
      {
	for (j=1; j<=ngens; j++) newrel[j]=0;
	rel[2]=n_one_of(rel[1]=n_one_of(rel[0]=k));
      
	//debug
	if (rel[2] != debug_n_one_2of(k))
	  { cerr <<"ERROR: Group action failure in triangle reln!"<<endl; }
	if (rel[0] != n_one_of(rel[2]))
	  { cerr <<"ERROR: Triangle rotation failure in 3-term reln!"<<endl;}
      
	if (verbose)   
	  cout << "Relation: " << rel[0]<<" "<<rel[1]<<" "<<rel[2]<<" --> ";
	for (j=0; j<3; j++)
	  {
	    check[ij=rel[j]] = 1;
	    check[sof(rel[j])] = 1;    // short cut - last triangle of p'ron
	    fix = coordindex[tysofs(ij,1)];
	    if (verbose)  cout << fix << " ";
	    if (fix!=0) newrel[abs(fix)] += sign(fix);
	  }
	if (verbose)  cout << endl;
	relmat.userel(newrel);
      }
  
    if(verbose)
      {
	cout<< "Current number of relations = " << relmat.get_numrel() <<endl;
      }
  }


// mixed 4-term reln  (I) + [N1] - (N1) + [S]
  {
    if (verbose) cout << "Mixed 4-term relation:\n";
    symbop sof(this,0,-1,1,0);
    symbop n_one_of(this,Quad(-1,1),2,2,Quad(0,-1));

    i=nsymb; while (i--) check[i]=0;
    for (k=0; k<nsymb; k++) if (check[k]==0)
      {
	for (j=1; j<=ngens; j++) newrel[j]=0;
	rel[0]=k;
	rel[1]=n_one_of(k);
	rel[2]=rel[1];
	rel[3]=sof(k);

	if (verbose)   
	  cout<<"Relation: "<<rel[0]<<" +"<<rel[1]<<" -"<<rel[2]<<
	    " +"<<rel[3]<<" --> ";
	check[k]=1;
	for (j=0; j<4; j++)
	  { 
	    long ty=(j%2);
	    fix = coordindex[tysofs(rel[j], ty )];
	    if (j==2) fix = -fix;
	    if (verbose)  cout << fix << " ";
	    if (fix!=0) newrel[abs(fix)] += sign(fix);
	  }
	if (verbose)  cout << endl;
	relmat.userel(newrel);
      }
  }

      break;  // end of face-rels for d=19

    case 5:
{
// Face relation (I) + (TS) + ((TS)^2)

  if(verbose) cout << "Face relation (I) + (TS) + ((TS)^2) :\n";

  symbop tsof(this,1,-1,1,0);

  symbop debug_ts2of(this,0,-1,1,-1);

#ifdef USE_REDUNDANT_TRIANGLES
  for (long casus=0; casus<=2; casus+=2)
#else
  for (long casus=0; casus<=0; casus+=2)   // JSB predicts casus 2 redundant!
#endif
    {
      i=nsymb; while (i--) check[i]=0;
      for (k=0; k<nsymb; k++) if (check[k]==0)
	{
	  for (j=1; j<=ngens; j++) newrel[j]=0;
	  rel[2]=tsof(rel[1]=tsof(rel[0]=k));
	  
	  //debug
	  if (rel[2] != debug_ts2of(k))
	    { cerr <<"ERROR: Group action failure in triangle reln!"<<endl; }
	  if (rel[0] != tsof(rel[2]))
	    { cerr <<"ERROR: Triangle rotation failure in 3-term reln!"<<endl;}

	  if (verbose)   
	    cout << "Relation: " << rel[0]<<" "<<rel[1]<<" "<<rel[2]<<" --> ";
	  for (j=0; j<3; j++)
	    {
	      check[ij=rel[j]] = 1;
	      fix = coordindex[tysofs(ij,casus)];
	      if (verbose)  cout << fix << " ";
	      if (fix!=0) newrel[abs(fix)] += sign(fix);
	    }
	  if (verbose)  cout << endl;
	  relmat.userel(newrel);
	}

      if(verbose)
	{
	  cout << "After type "<<casus<<" triangle relation,";
	  cout<< " number of relations = " << relmat.get_numrel() <<"\n";
	}
    }

// "type 1" 4-term reln ---  [I] + [B] + [B^2] + [B^3]
// that is,  (I)_1  -  (T^-1)_1  +  (B^2)_1  - (B^2T^-1)_1

  if (verbose) cout << "Type 1 4-term relation:\n";

  symbop dcmatrix(this,1,-1,0,1);
  symbop b2matrix(this,Quad(0,1),2,2,Quad(0,-1));
  
  i=nsymb; while (i--) check[i]=0;
  for (k=0; k<nsymb; k++) if (check[k]==0)
    {
      for (j=1; j<=ngens; j++) newrel[j]=0;
      rel[0]=k;
      rel[1]=dcmatrix(k);
      rel[2]=b2matrix(k);
      rel[3]=dcmatrix(rel[2]);
      if (verbose)   
	cout<<"Relation: "<<rel[0]<<" -"<<rel[1]<<" +"<<rel[2]<<" -"<<rel[3]<<" --> ";
      for (j=0; j<4; j++)
	{
	  check[ij=rel[j]] = 1;
	  fix = coordindex[tysofs(ij,1)];
	  if (j%2) fix = -fix;
	  if (verbose)  cout << fix << " ";
	  if (fix!=0) newrel[abs(fix)] += sign(fix);
	}
      if (verbose)  cout << endl;
      relmat.userel(newrel);
    }

  if(verbose)
    {
      cout<< "Current number of relations = " << relmat.get_numrel() <<endl;
    }

// mixed 4-term reln, using N-matrices ---   (I) -  [I]  +  (BS)  -  [BS]
// ... using only P-matrices                 e0  -   e1  -   e2   +  (TS)^2 e1
// that is,                               (I)_0  - (I)_1 - (I)_2  + ((TS)^2))_1

  if (verbose) cout << "Mixed 4-term relation:\n";

  symbop ts2of(this,0,-1,1,-1);

  i=nsymb; while (i--) check[i]=0;
  for (k=0; k<nsymb; k++) if (check[k]==0)
    {
      for (j=1; j<=ngens; j++) newrel[j]=0;
      rel[2]=rel[1]=rel[0]=k;
      rel[3]=ts2of(k);

      if (verbose)   
	cout<<"Relation: "<<rel[0]<<" -"<<rel[1]<<" -"<<rel[2]<<
	  " +"<<rel[3]<<" --> ";
      check[k]=1;
      for (j=0; j<4; j++)
	{ 
	  long ty;                                        //  j= 0  1  2  3
	  if (j==3) ty=1; else ty=j;                      // ty= 0  1  2  1
	  fix = coordindex[tysofs(rel[j], ty )];
	  if ((j==1)||(j==2)) fix = -fix;
	  if (verbose)  cout << fix << " ";
	  if (fix!=0) newrel[abs(fix)] += sign(fix);
	}
      if (verbose)  cout << endl;
      relmat.userel(newrel);
    }
}
      break;  // end of face-rels for d=5

    default:
      cerr<< "Error: relate_faces not implemented for field "<<Field::d<<" !";
      exit(1);
    }

  delete[] check;
  delete[] rel;
  relmat.trim();

  if(verbose) 
    {
      cout << "Finished face relations: ";
      cout << "number of relations = " << relmat.get_numrel() << endl;
      cout << "relmat = " << endl;
      relmat.output_pretty(cout);
    }


  int ker_method=KER_RELMAT_DEFAULT_METHOD;

#ifdef KER_RELMAT_RUNTIME_PROMPT
  if (verbose)
    { cout << "Input method for kernel(relmat): ";
      cin >> ker_method;
      cout << endl;
    }
#endif

  subspace sp = kernel(relmat, ker_method);
  rk = dim(sp);
  coord = basis(sp);
  vec pivs = pivots(sp);
  denom1 = denom(sp);

  if (verbose) 
    {
      cout << "rk = " << rk << endl;
      cout << "coord:" << endl;
      coord.output_pretty(cout);
      if(denom1!=1) cout << "denominator = " << denom1 << endl;
      cout << "pivots = " << pivs <<endl;
    }

  long *freegens;
  if (rk>0)
    {
      freemods = new n_modsym[rk];
      freegens = new long[rk];
      for (i=0; i<rk; i++)
	{
	  freegens[i] = gens[pivs[i+1]];
	  freemods[i] = modsymoftys(freegens[i]);
	}
      if (verbose)
        { 
          cout << "freegens: ";
          for (i=0; i<rk; i++) cout << freegens[i] << " ";
          cout << endl;
        }
    }
  delete[] freegens;
  delete[] gens;

//  relmat.init(); newrel.init();   // no need? destructors should suffice!
}

void homspace::restrict_to_kerdelta()
{
  needed = new int[rk];
  cusplist cusps(modulus, GL_vs_SL_flag, 2*rk);    // 2*rk is upper bd on no.
                                                   // of inequiv cusps arising
  mat deltamat(2*rk,rk);
  long i,j,k;
  for (i=0; i<rk; i++)
    {
      n_modsym m=freemods[i];
      for (j=1; j>-3; j-=2)
	{
	  n_cusp c = (j==1 ? m.beta() : m.alpha());
	  k = cusps.index(c);   //adds automatically if new
	  deltamat(k+1,i+1) += j;  // N.B. offset of 1
	}
    }
  ncusps=cusps.count();
  if(verbose)cout << "ncusps = " << ncusps << endl;

  deltamat=deltamat.slice(ncusps,rk);
  if (verbose)
    { cout << "sliced deltamat = " << endl;
      deltamat.output_pretty(cout);
    }

// the i'th col of deltamat denotes the image of the i'th freegen under the map
// delta: (chns/bdys) --> H_0(quotient space) = v.sp. with basis inequiv cusps

  int ker_method=KER_RELMAT_DEFAULT_METHOD;

#ifdef KER_RELMAT_RUNTIME_PROMPT
  if (verbose)
    { cout << "Input method for kernel(deltamat): ";
      cin >> ker_method;
      cout << endl;
    }
#endif

  kern = kernel(deltamat, ker_method);
  dimension = dim(kern);
  denom2 = denom(kern);
  denom3 = denom1 * denom2;

// the cols of basiskern now form a basis of the desired homology space H_1
// (a col expresses an elt of our basis for H_1 in terms of the freemods)

// a trivial row means that the corr freemod is not needed in expressing the
// basic H_1 classes in terms of freemods

  if (dimension>0)
    {
      const mat& basiskern = basis(kern);
      if (verbose)  cout << "Freemods:\n";
      for (i=0; i<rk; i++)
        {
          needed[i]   =  ! trivial(basiskern.row(i+1));
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
          basiskern.output_pretty(cout);
          cout << "pivots: " << pivots(kern) << endl;
        }
    }
}

vec homspace::chain_signed_symb(long k, long eps) const
{
 long i= coordindex[k];
 vec ans(rk);
 if (i) ans = (sign(i)*eps)*(coord.row(abs(i)));
 return ans;
}

vec homspace::reconvert(const n_cusp&z) const
// express {infty, z} in terms of the generators of homology
{
  vec ans(rk);
  pseudo_euclid ps(z);
// for debugging
//  if((z.a()==Quad(4,4))&&(z.c()==Quad(-22,11)))
//    {
//      cerr << "Debug event in homspace::reconvert has occurred!\n";
//    }
  while (ps.not_infty())         // process the modular symbol +/- (M)_t
    {
      long k = tysofs(numsymb2(ps.summand_m().c(), ps.summand_m().d()),
		      ps.summand_t() );
      ans += chain_signed_symb(k, ps.summand_eps());
    }
  return ans;
}

vec homspace::reconvert(const n_modsym&m) const
//  {a,b} = {infty, b} - {infty, a}
{
  vec ans = reconvert(m.beta()) - reconvert(m.alpha());
  return ans;
}

mat homspace::calc_conj() const
{
  mat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     {
       vec colj(rk);
       colj = reconvert(freemods[j].conj());
       m.setcol(j+1,colj);
     }
  mat ans = restrict(m,kern,1);
  return ans;
}

mat homspace::calcop(const matop& mlist) const
{
  mat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     {
       vec colj(rk);
       for (long i=0; i<mlist.getlength(); i++)
	 {
	   colj += reconvert(mlist(i)(freemods[j]));
	 }
       m.setcol(j+1,colj);
     }
  mat ans = restrict(m,kern,1);
  return ans;
}

mat homspace::chi_p(const Qideal&p) const
{
  matop mlist=Chi_p(p);
  return calcop(mlist);
}
  
mat homspace::hecke_p(const Qideal&p) const
{
  matop mlist=T_p(p);
  return calcop(mlist);
}

// Added 20-07-1998
mat homspace::hecke_u_princ(const Qideal&p) const
{
  matop mlist=U_p(p);
  return calcop(mlist);
}

mat homspace::hecke_npnp(const Qideal&p) const
{
  matop mlist=T_npnp(p);
  return calcop(mlist);
}

mat homspace::hecke_npnpchi(const Qideal&p, const Qideal&r) const
{
  matop mlist=T_npnp(p);
  n_mat chimat=Chi_p(r)[0];
  for (long i=0; i<mlist.getlength(); i++) {mlist[i]=mlist[i]*chimat;}
  return calcop(mlist);
}

mat homspace::hecke_npnq(const Qideal&p, const Qideal&q) const
{
  matop mlist=T_npnq(p,q);
/*
  { // experiment: try inverse matrices
    long imax=mlist.getlength();
    matop mlist2(imax);
    for (long i=0; i<imax; i++) { mlist2[i]=mlist[i].inv_mod_scalars(); }
    cerr << "mlist: " << endl << mlist << endl;
    cerr << "mlist2: " << endl << mlist2 << endl;
    mlist=mlist2;
    cerr << "mlist: " << endl << mlist << endl;
  }
*/
  return calcop(mlist);
}

mat homspace::hecke_npnqchi(const Qideal&p, const Qideal&q, const Qideal&r) const
{
  matop mlist=T_npnq(p,q);
  n_mat chimat=Chi_p(r)[0];
  for (long i=0; i<mlist.getlength(); i++) {mlist[i]=mlist[i]*chimat;}
  return calcop(mlist);
}

/*
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
     e=d; d= c; c=-q*c-e;
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
     e=d; d= c; c=-q*c-e;
     ans += projchaincd(c,d);
   }
   return ans;
}

vec homspace::applyop(const matop& mlist, const RatQuad& q) const
{ vec ans(rk);  long i=mlist.length;
  while (i--) ans += chain(mlist[i](q));
  return ans;
}
 
mat homspace::calcop(char* opname, const Quad& p, const matop& mlist, int display) const
{
  mat m(rk,rk);
  for (long j=0; j<rk; j++) if (needed[j])
     { vec colj = applyop(mlist,freemods[j]);
       m.setcol(j+1,colj);
     }
  mat ans = restrict(m,kern);
  if (display) cout << "Matrix of " << opname << "(" << p << ") = " << ans;
  if (display && (dimension>1)) cout << endl;
  return ans;
}
 
mat homspace::heckeop(const Quad& p,int display) const
{
 matop matlist(p,modulus);
 char* name = (div(p,modulus)) ? "W" : "T";
 return calcop(name,p,matlist,display);
}
 
mat homspace::wop(const Quad& q,int display) const
{
 matop matlist(q,modulus);
 return calcop("W",q,matlist,display);
}
 
mat homspace::fricke(int display) const
{
 matop frickelist(modulus,modulus);
 return calcop("W",modulus,frickelist,display);
}

mat homspace::opmat(long i, long v)
{
  if((i<0)||(i>=nap)) return mat(dimension);  // shouldn't happen
  Quad p = primelist(i);
  if(verbose) 
      cout<<"Computing " << ((div(p,modulus)) ? "W" : "T") <<"("<<p<<")...";
  return heckeop(p,v); // Automatically chooses W or T
}

longlist homspace::eigrange(long i)    // implementing virtual fn in matmaker
{
  if((i<0)||(i>=nap)) return longlist(0);  // shouldn't happen
  Quad p = primelist(i);
  long normp = p.norm();
  if (verbose) 
    cout << "eigrange for p = " << p << ":\t";
  if(div(p,modulus))
    {
      longlist ans(2);
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
      longlist ans(l);
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
  Quadlist resmodp=residues(p); 
  vec tvec = chain(0,p);             // =0, but sets the right length.
  for (Quadvar i(resmodp); i.ok(); i++) tvec += chain((Quad)i,p);
  return kernelpart(tvec);
}

vec homspace::manintwist(const Quad& lambda, const Quadlist& res, int* chitable) const
{
 vec sum = chain(0,lambda);          // =0, but sets the right length.
 long i=0;
 for (Quadvar r(res); r.ok(); r++, i++) 
    sum += chitable[i]*chain((Quad)r,lambda);
 return kernelpart(sum);
}

vec homspace::projmaninvector(const Quad& p) const    // Will only work after "proj"
{
  Quadlist resmodp=residues(p); 
  vec tvec = projcycle(0,p);             // =0, but sets the right length.
  for (Quadvar i(resmodp); i.ok(); i++) tvec += projcycle((Quad)i,p);
  return tvec;
}

vec homspace::newhecke(const Quad& p, const Quad& n, const Quad& d) const
                                     // Will only work after "proj"
{ 
  vec tvec = projcycle(p*n,d);
  Quadlist resmodp=residues(p);  Quad dp = d*p;
  for (Quadvar k(resmodp); k.ok(); k++) 
    tvec += projcycle(n+d*Quad(k),dp);
  return tvec;
}

*/

// END OF FILE homspace.cc
