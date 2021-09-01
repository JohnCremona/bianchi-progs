load "input.m";

function ToCode(x)   // Convert a numbered generator of the vector space
                     // into a standard code
   y:=x mod np;
   y:=y eq 0 select np else y; 
   z:=Z!(1+(x-y)/np);
   return <z,alphabet[y]>;   // the number tells you the basic path
                             // the letter which coset
end function;

// KI_g is OK/N (N being the level).
// PKI_g is P^1(OK/N), which is where the (c:d) symbols live.

function List(S,level_g,KI_g,PKI_g)  // Generates a list of cd-symbols
                                     // where S is the prime power 
                                     // factorisation of the level
   level_l:=&*[s:s in S];
   if level_l eq level_g then
      KI_l:=KI_g;
      PKI_l:=PKI_g;
   else
      KI_l:=quo<OK|level_l>;
      PKI_l:=ProjectiveSpace(KI_l,1);
   end if;

   cdlist:=[* *];
   for pi in S do
      KI_l2:=quo<OK|pi>;
      G:=FreeAbelianGroup(2);
      H:=sub<G|[Eltseq(g):g in Basis(pi)]>;
      cr1:=[OK!Eltseq(a):a in Transversal(G,H)];
      cr2:=[c:c in cr1 | CoprimeToIdeal(c,pi) eq false];

      crep1:=[Matrix([[1,0],[x,1]]):x in cr1] cat
             [Matrix([[0,-1],[1,x]]):x in cr2];
      PKI_l2:=ProjectiveSpace(KI_l2,1);
      cdlist cat:=[*[PKI_l2![KI_l2!m[2][1],KI_l2!m[2][2]]:m in crep1]*];
   end for;
     
   while #cdlist gt 1 do
      cdn:=#cdlist;
      B:=Bezout(S[cdn-1],S[cdn]);
      new:=[];
      level_l:=S[cdn-1]*S[cdn];
      KI_l:=quo<OK|level_l>;
      PKI_l:=ProjectiveSpace(KI_l,1);
      for X in cdlist[cdn-1] do
         for Y in cdlist[cdn] do      
            c:=OK!X[1]*B[2]+OK!Y[1]*B[1];
            d:=OK!X[2]*B[2]+OK!Y[2]*B[1];
//      c:=CRT(S[cdn-1],S[cdn],OK!X[1],OK!Y[1]);
//      d:=CRT(S[cdn-1],S[cdn],OK!X[2],OK!Y[2]);
// there doesn't seem to be any real difference between the CRT function 
// and my code
            new cat:=[PKI_l![KI_l!c,KI_l!d]];      
         end for;
      end for;
      Prune(~cdlist);
      cdlist[cdn-1]:=new;
      Prune(~S);
      S[cdn-1]:=level_l;
   end while;      

   assert #cdlist eq 1;
   cdlist:=cdlist[1]; 
   return cdlist,KI_l,PKI_l;  

end function;

function ToMatrix(cd,level_l,KI_l,PKI_l)  
// Construct a right coset representative matrix from a cd-symbol
// Based on Bygott pg 26, Prop 22

   gamma:=OK!cd[1];
   delta:=OK!cd[2]; 

   x1:=OK!(delta*w);
   x2:=OK!(delta);
   x3:=OK!(gamma*w);
   x4:=OK!(gamma);
     
   if CoprimeToElement(gamma,delta) eq true then
     
      N:=Matrix(Z,4,2,[x1[1],x1[2],
                       x2[1],x2[2],
                      -x3[1],-x3[2],
	              -x4[1],-x4[2]]);
   
      v:=Vector(Z,[1,0]);

      s:=Solution(N,v);  
      alpha:=s[1]*w+s[2];
      beta:=s[3]*w+s[4];

      M:=MK!Matrix([[alpha,beta],[gamma,delta]]);
      assert Determinant(M) eq 1;
      return M;

   else

      z1:=Basis(level_l)[1];
      z2:=Basis(level_l)[2];

      N:=Matrix(Z,6,2,[x1[1],x1[2],
                       x2[1],x2[2],
                      -x3[1],-x3[2],
	              -x4[1],-x4[2],
                       z1[1],z1[2],
                       z2[1],z2[2]]);
   
      v:=Vector(Z,[1,0]);

      s:=Solution(N,v);  
      alpha:=s[1]*w+s[2];
      beta:=s[3]*w+s[4];

      assert Determinant(Matrix([[alpha,beta],[gamma,delta]]))-1 in level_l;

      if CoprimeToIdeal(alpha,level_l) eq true then
         
         M:=RMatrixSpace(Z,4,2)!0;

         BasisA:=Basis(ideal<OK|alpha>);
         BasisL:=Basis(level_l);

         M[1][1]:=BasisA[1][1];
         M[1][2]:=BasisA[1][2];
         M[2][1]:=BasisA[2][1];
         M[2][2]:=BasisA[2][2];
         M[3][1]:=BasisL[1][1];
         M[3][2]:=BasisL[1][2];
         M[4][1]:=BasisL[2][1];
         M[4][2]:=BasisL[2][2];

         v:=Vector(Z,[1,0]);
         sol:=Solution(M,v);

         u1:=sol[1]*BasisA[1]+sol[2]*BasisA[2];
         u2:=sol[3]*BasisL[1]+sol[4]*BasisL[2];  
 
         a:=alpha;
         b:=u2+beta*u1;
         c:=-u2+gamma*u1;
         d:=(1+b*c)/a;
         assert d in OK;
         N:=MK!Matrix([[a,b],[c,d]]);

         assert Determinant(N) eq 1;
         assert PKI_l![KI_l!cd[1],KI_l!cd[2]] eq PKI_l![KI_l!c,KI_l!d];
         return N;

      else 
         if alpha ne 0 then
            S:=Support(level_l) diff Support(ideal<OK|alpha>);
            q:=ideal<OK|1>;
            for s in S do
               q *:=s;
            end for;
         else
	    q:=ideal<OK|1>;
         end if;
     
         q1:=Basis(q)[1];
         q2:=Basis(q)[2];
         tau:=q^-1;

         while tau+ideal<OK|alpha> ne ideal<OK|1> do
            tau:=(Random(-10,10)*q1+Random(-10,10)*q2)*q^(-1);
         end while;

         assert Denominator(tau) eq 1;
 
         tf,lambda:=IsPrincipal(q*tau);
         assert tf eq true;
         alpha:=alpha+lambda*gamma;
         beta:=beta+lambda*delta;

         assert CoprimeToIdeal(alpha,level_l) eq true;

         M:=RMatrixSpace(Z,4,2)!0;

         BasisA:=Basis(ideal<OK|alpha>);
         BasisL:=Basis(level_l);
  
         M[1][1]:=BasisA[1][1];
         M[1][2]:=BasisA[1][2];
         M[2][1]:=BasisA[2][1];
         M[2][2]:=BasisA[2][2];
         M[3][1]:=BasisL[1][1];
         M[3][2]:=BasisL[1][2];
         M[4][1]:=BasisL[2][1];
         M[4][2]:=BasisL[2][2];

         v:=Vector(Z,[1,0]);
         sol:=Solution(M,v);

         u1:=sol[1]*BasisA[1]+sol[2]*BasisA[2];
         u2:=sol[3]*BasisL[1]+sol[4]*BasisL[2];  

         a:=alpha;
         b:=u2+beta*u1;
         c:=-u2+gamma*u1; 
         d:=(1+b*c)/a;
         assert d in OK;
         N:=MK!Matrix([[a,b],[c,d]]);
         assert Determinant(N) eq 1;
         assert PKI_l![KI_l!cd[1],KI_l!cd[2]] eq PKI_l![KI_l!c,KI_l!d];
         return N;

      end if;
   end if;
end function;

function AreEquivalent(CuspA,CuspB,level_g)  // returns whether or not two
                                             // cusps are equivalent modulo
                                             // Gamma_0(level_g)
   if CuspA eq CuspB then
      return "true";
   end if;
   
   if CuspIdealClass(CuspA) ne CuspIdealClass(CuspB) then
      return "false";
   else
      d1:=Denominator(CuspA[1]);
      p1:=CuspA[1]*d1;
      q1:=CuspA[2]*d1;
  
      d2:=Denominator(CuspB[1]);
      p2:=CuspB[1]*d2;
      q2:=CuspB[2]*d2;

      A:=ideal<OK|p1,q1>;
      B:=ideal<OK|p2,q2>;

      tf,gen:=IsPrincipal(A/B);

      p1 /:=gen;             // write the first fraction so its numerator and
      q1 /:=gen;             // denominator generate the same ideal as the
                             // second fraction

      assert ideal<OK|p1,q1> eq ideal<OK|p2,q2>;

      z1:=Conj(p2);
      z2:=Conj(q2);

      C:=ideal<OK|z1,z2>;     // Conjugate ideal to B

      tf2,delta:=IsPrincipal(B*C);

      z1:=Basis(C)[1];
      z2:=Basis(C)[2];

      x1:=OK!(p1*z1);
      x2:=OK!(p1*z2);
      x3:=OK!(q1*z1);
      x4:=OK!(q1*z2);
      x5:=OK!(p2*z1);
      x6:=OK!(p2*z2);
      x7:=OK!(q2*z1);
      x8:=OK!(q2*z2);

      N1:=Matrix(Z,4,2,[x1[1],x1[2],
                        x2[1],x2[2],
                       -x3[1],-x3[2],
		       -x4[1],-x4[2]]);

      N2:=Matrix(Z,4,2,[x5[1],x5[2],
                        x6[1],x6[2],
                       -x7[1],-x7[2],
		       -x8[1],-x8[2]]);

      v1:=Vector(Z,[delta,0]);

      s1:=Solution(N1,v1);
      s2:=Solution(N2,v1);

      N3:=Matrix([[p1,s1[3]*z1+s1[4]*z2],
                  [q1,s1[1]*z1+s1[2]*z2]]);
      N4:=Matrix([[p2,s2[3]*z1+s2[4]*z2],
                  [q2,s2[1]*z1+s2[2]*z2]]);     

      X:=ideal<OK|q1*q2/delta>*ideal<OK|z1,z2>^2;
      Y:=ideal<OK|delta>*level_g; 
      J:=X+Y;

      M1:=N4*N3^-1;
      assert M1 in MK;
      assert Determinant(M1) eq 1;
      assert ActC(M1,CuspA) eq CuspB;

      llentry:=delta*M1[2,1];    

      if llentry in J then  
         return "true";
      else
         return "false";
      end if;      
   end if;
end function;

function DimH(level_g,verbose_flag)  // when verbose_flag is false,
    // just returns dimension of the  homology space. Otherwise the
    // following data is also outputed (where N=level_g):
    //
    // KI_g:  OK/N
    // PKI_g: P^1(OK/N)
    // V: ambient space of all (c:d)_i symbols
    // V2: Quotient of V modulo relations
    // V3: Cuspidal subspace of V2 (kernel of boundary map)
    // MSC: a modular symbol basis for V2 (in the form {alpha,beta})
    // CDlist: the list of (c:d) symbols
    // qmap: quotient map V -> V2
    
    if verbose_flag then
	print "Level = ",level_g;
    end if;  

   KI_g:=quo<OK|level_g>;
   PKI_g:=ProjectiveSpace(KI_g,1);
   MKI_g:=RMatrixSpace(KI_g,2,2);

   assert Norm(level_g) ne 1;  // check the level isn't trivial  
 
   CDlist:=List([s[1]^s[2]:s in Factorization(level_g)],level_g,KI_g,PKI_g);

   cdn:=#CDlist;
   if verbose_flag then
       print "Number of (c:d) symbols = ",cdn;
   end if;  
   
   crep:=[ToMatrix(cd,level_g,KI_g,PKI_g): cd in CDlist];

   CDlist:=[PKI_g![KI_g!m[2][1],KI_g!m[2][2]]:m in crep];

   function CDClass(tcd)     // given a cd-symbol returns its position
      for cdl in CDlist do   // in the above list
         if tcd eq cdl then
            return Position(CDlist,cdl);
         end if;
      end for;      
      error "help",<K!(OK!tcd[1]),K!(OK!tcd[2])>;
   end function;

   function ToSymbol(M)  // Convert the bottom row of a matrix to a cd symbol
      return PKI_g![KI_g!M[2][1],KI_g!M[2][2]];
   end function;

   Rels1:=[];                  // Create 1 relation for each coset rep
   for x in rels do            // using code
      for M in crep do
         dummy1:=[];
         for r in x do          
            tc:=ToSymbol(M*r[2]);       
            p:=CDClass(tc);
            dummy1 cat:=[<r[1],r[3],p>];  // sign,path,coset   
         end for;
         Rels1 cat:=[dummy1];
      end for;
   end for;

   Rels2:=[];                // Write relations in terms of simple 
   for x in Rels1 do         // numbered variables
      dummy2:=[];
      for r in x do
         dummy2 cat:=[<r[1],r[2]+(r[3]-1)*np>];   // sign,number         
      end for;
      Rels2 cat:=[dummy2];
   end for;

   V:=VectorSpace(Q,np*cdn);   // original vector space

   Rels3:=[&+[t[1]*V.t[2]:t in r]:r in Rels2];   // Encode the relations
   Rels4:=[r:r in Rels3 | r ne V!0];             // in the vector space
   M:=RMatrixSpace(Q,#Rels4,np*cdn)!0;           
   for i in [1..#Rels4] do
      M[i]:=Rels4[i];
   end for;
   EF:=EchelonForm(M);

   S:=sub<V| [EF[i]:i in [1..#Rels4] | EF[i] ne V!0]>;

   V2,qmap:=quo<V|S>;     // V modulo the relations and a map going from
                          // V to V2

   B:=Basis(S);
   pivots:=Set(Pivots(B));
   gens:=[i : i in [1..np*cdn] | i notin pivots];    // a basis for V2

   Gens:=[ToCode(i): i in gens];    // this basis converted into my code

   MS:=[];                              // Convert remaining symbols into
   for g in Gens do                     // standard (gamma)_alpha form
      pa:=Position(alphabet,g[2]);
      alpha:=EndPoints[pa];    
      gamma:=ToMatrix(CDlist[g[1]],level_g,KI_g,PKI_g);
      MS cat:=[<alpha,gamma>];
   end for;

   MSC:=[[ActC(m[2],m[1]),ActC(m[2],Cusps![1,0])]:m in MS]; 
   // convert symbols into modular symbols of the form {a,b}

   cc:=C;  // which is defined in core.m                             

   for c in MSC do                     // generates a list of inequivalent                                        
      for cp in c do                   // cusps mod Gamma_0(level_g)
         mcp:=Cusps![-cp[1],cp[2]];
         for class in cc do
            if AreEquivalent(cp,class,level_g) eq "true" then
               v:=0;
               break;
            elif AreEquivalent(mcp,class,level_g) eq "true" then
               v:=0;
               break;
            else
 	       v:=1;
            end if;
         end for;
         if v eq 1 then
     	    cc cat:=[Cusps![cp[1],cp[2]]];
         end if;
      end for;
   end for;

   function CuspClass(Cusp,level_g)   // given a cusp returns which cusp
                                      // in the list above it is
                                      // equivalent to
      for class in cc do
         mCusp:=Cusps![-Cusp[1],Cusp[2]];
         if AreEquivalent(Cusp,class,level_g) eq "true" then
            return Position(cc,class);
         elif AreEquivalent(mCusp,class,level_g) eq "true" then
            return Position(cc,class);
         end if;
      end for;
   end function;

   vec_ms:=V2;
   vec_c:=VectorSpace(Q,#cc);
   vb:=Basis(vec_c);
  
   image:=[vb[CuspClass(m[2],level_g)]-vb[CuspClass(m[1],level_g)]:m in MSC];
   
   bmap:=hom<vec_ms -> vec_c | image>;     // defines the boundary map

   V3:=Kernel(bmap);      //  the kernel of the boundary map on V2

   if verbose_flag then
      return <KI_g,PKI_g,V,V2,V3,MSC,CDlist,qmap>;
   else
      return Dimension(V3);
   end if;

end function;
