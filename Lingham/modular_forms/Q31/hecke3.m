load "homology.m";
load "region3.m";

function Fricke(N)

   tf,g:=IsPrincipal(N^3);

   F:=Factorization(N^2);

   z:=1*OK;

   for f in F do
      z:=z*f[1]^f[2];
   end for;
  
   c:=cm(hinv(z));

   if c eq 1 then        
   elif c eq 2 then
      z*:=ideal<OK|211,w+89>;
   elif c eq 3 then
      z*:=ideal<OK|211,w+121>;
   end if;

   tf,z:=IsPrincipal(z);

   x:=1*OK;

   for f in F do
      x:=x*f[1]^f[2];
   end for;
  
   c:=cm(hinv(x));

   if c eq 1 then   
   elif c eq 2 then
      x:=x*ideal<OK|233,w+96>;
   elif c eq 3 then
      x:=x*ideal<OK|233,w+136>;
   end if;

   tf,x:=IsPrincipal(x);

   x3:=K!Basis(N)[1];
   x4:=K!Basis(N)[2];

   X:=Matrix(Z,4,2,[(x3*x*g)[1],(x3*x*g)[2],               
                    (x4*x*g)[1],(x4*x*g)[2],  
                    (-x3*z)[1],(-x3*z)[2],
                    (-x4*z)[1],(-x4*z)[2]]); 
     
   v:=Vector(Z,[g[1],g[2]]);

   s:=Solution(X,v);  
   y:=s[3]*x3+s[4]*x4;
   W:=s[1]*x3+s[2]*x4;
  
   M:=Matrix([[x,y],[z,g*W]]);
  
   assert Determinant(M)*OK eq g*OK;   
   assert ideal<OK|M[1][1],M[2][1]> eq N^2;
   assert ideal<OK|M[1][2],M[2][2]> eq N;
   assert M[2][2] in N^2;
   assert M[2][1] in N;

   return M;

end function;

function AtkinLehner(N,pp)  // calculates a matrix representation of the
                            // Atkin-Lehner involution W_pp 

   tf,g:=IsPrincipal(pp^3);

   F1:=Factorization(pp^2);
   F2:=Factorization(N/pp);

   z:=1*OK;

   for f in F1 cat F2 do
      z:=z*f[1]^f[2];
   end for;
  
   c:=cm(hinv(z));

   if c eq 1 then        
   elif c eq 2 then
      z*:=ideal<OK|211,w+89>;
   elif c eq 3 then
      z*:=ideal<OK|211,w+121>;
   end if;

   tf,z:=IsPrincipal(z);

   x:=1*OK;

   for f in F1 do
      x:=x*f[1]^f[2];
   end for;
  
   c:=cm(hinv(x));

   if c eq 1 then   
   elif c eq 2 then
      x:=x*ideal<OK|233,w+96>;
   elif c eq 3 then
      x:=x*ideal<OK|233,w+136>;
   end if;

   tf,x:=IsPrincipal(x);

   x3:=K!Basis(pp)[1];
   x4:=K!Basis(pp)[2];

   X:=Matrix(Z,4,2,[(x3*x*g)[1],(x3*x*g)[2],               
                    (x4*x*g)[1],(x4*x*g)[2],  
                    (-x3*z)[1],(-x3*z)[2],
                    (-x4*z)[1],(-x4*z)[2]]); 
     
   v:=Vector(Z,[g[1],g[2]]);

   s:=Solution(X,v);  
   y:=s[3]*x3+s[4]*x4;
   W:=s[1]*x3+s[2]*x4;
  
   M:=Matrix([[x,y],[z,g*W]]);

   assert Determinant(M)*OK eq g*OK;   
   assert M[1][1] in pp^2;
   assert M[1][2] in pp;
   assert M[2][1] in N*pp;
   assert M[2][2] in pp^2;

   return [M];
        
end function;  

function PrincipalHecke(p)  // calculates a list of matrices whose formal sum
                            // defines T_p for p a principal ideal
   a,beta:=IsPrincipal(p);
   assert IsPrime(p) eq true;
   assert a eq true;
   F:=quo<OK|p>;
   G:=FreeAbelianGroup(2);
   H:=sub<G|[G!Eltseq(g):g in Basis(p)]>;
   cr:=[OK!Eltseq(a):a in Transversal(G,H)];
   Tp:=[Matrix([[beta,0],[0,1]])];
   for alpha in cr do
      Tp cat:=[Matrix([[1,alpha],[0,beta]])];
   end for;
   return Tp;
end function;

function NewBasis(b,level_g)
// Based on pg 192-193 of "A Course in Computational Algebraic Number Theory"
// Finds a basis of b whose first element is in level_g
// Needed for the function M_ab
   
   b1:=Basis(level_g meet b)[1];

   F:=Factorization(b1*OK);
   r:=#F;
   G:=[<F[i][1],Valuation(b,F[i][1])>:i in [1..r]];

   I:=&*[f[1]^(f[2]+1):f in G];

   A:=[];
   for i in [1..r] do
      A[i]:=I/(G[i][1]^(G[i][2]+1));
   end for;

   assert &+[a:a in A] meet OK eq ideal<OK|1>;
 
   M:=RMatrixSpace(Z,2*r,2)!0;
   for i in [1..r] do
      M[2*i-1][1]:=Basis(A[i])[1][1];
      M[2*i-1][2]:=Basis(A[i])[1][2];
      M[2*i][1]:=Basis(A[i])[2][1];
      M[2*i][2]:=Basis(A[i])[2][2];
   end for;

   v:=Vector(Z,[1,0]);
   sol:=Solution(M,v);

   U:=[];
   for i in [1..r] do
      u:=sol[2*i-1]*Basis(A[i])[1]+sol[2*i]*Basis(A[i])[2];
      U[i]:=K!u;
   end for;

   assert &+[U[i]:i in [1..r]] eq 1;

   B:=[];
   for i in [1..r] do
      t:=Basis(G[i][1]^G[i][2])[1];
      assert t notin G[i][1]^(G[i][2]+1);
      B[i]:=t;
   end for;

   b2:=&+[B[i]*U[i]:i in [1..r]];

   assert b1 in b meet level_g;
   assert ideal<OK|b1,b2> eq b;
   assert b1 in level_g; 

   return <K!b1,K!b2>;

end function;

function M_ab(A,B,level_g)  // finds a matrix which maps OK^2 to
                            // A+B and has lower left entry in level_g

   tf,gamma:=IsPrincipal(A*B);

   c,d:=Explode(NewBasis(A,level_g)); 

   x1:=Basis(B)[1];
   x2:=Basis(B)[2];

   N:=Matrix(Z,4,2,[(x1*d)[1],(x1*d)[2],
                    (x2*d)[1],(x2*d)[2],
                    (-x1*c)[1],(-x1*c)[2],
                    (-x2*c)[1],(-x2*c)[2]]);
   
   v:=Vector(Z,[gamma[1],gamma[2]]);

   s:=Solution(N,v);  
   a:=s[1]*x1+s[2]*x2;
   b:=s[3]*x1+s[4]*x2;

   M:=MK!Matrix([[d,b],[c,a]]);

   assert Determinant(M) eq gamma;
   assert M[2][1] in level_g;
   assert ideal<OK|M[1][1],M[2][1]> eq A;
   assert ideal<OK|M[1][2],M[2][2]> eq B;

   return M;

end function;

function HeckeNP(prime,level_g,KI_g,PKI_g) 
// calculates a list of matrices whose formal sum defines T_p for p a 
// non-principal ideal

   SMlist:=[];

   a:=prime^2;                  
   b:=prime;
   level_l:=a/b;
   Wlist,KI_l,PKI_l:=List([level_l],level_g,KI_g,PKI_g);    

   M:=RMatrixSpace(Z,4,2)!0;

   BasisL:=Basis(level_l);
   BasisG:=Basis(level_g);
      
   M[1][1]:=BasisL[1][1];
   M[1][2]:=BasisL[1][2];
   M[2][1]:=BasisL[2][1];
   M[2][2]:=BasisL[2][2];
   M[3][1]:=BasisG[1][1];
   M[3][2]:=BasisG[1][2];
   M[4][1]:=BasisG[2][1];
   M[4][2]:=BasisG[2][2];
      
   v:=Vector(Z,[1,0]);
   s:=Solution(M,v);

   x:=s[1]*BasisL[1]+s[2]*BasisL[2];
   y:=s[3]*BasisG[1]+s[4]*BasisG[2];
   
   assert x+y eq 1;           
      
   Wlist:=[[OK!w[1]*y,OK!w[2]*y+x]:w in Wlist]; 

   Wrep:=[ToMatrix(mark,level_l,KI_l,PKI_l):mark in Wlist];

   Mab:=M_ab(a,b,level_g);      
   SMlist cat:=[Mab*W:W in Wrep];       
 
   assert #SMlist eq Norm(prime)+1;

   for sm in SMlist do  
      assert sm[2][1] in level_g;
   end for;
 
   return SMlist;

end function;    

function ModToCD(Cusp,level_g,KI_g,PKI_g,CDlist,V,qmap,V2)  
// Converts modular symbol of the form {infinity,cusp} to an M-symbol 

   function ToSymbol(M)  // Convert the bottom row of a matrix to a cd symbol
      return PKI_g![KI_g!M[2][1],KI_g!M[2][2]];
   end function;

   function MP(a,b)   
      x:=Position(EndPoints,a);
      y:=Position(EndPoints,b);
      if y eq np+1 then
         return x;
      else
         return -y;
      end if;
   end function;

   function num(d)
      m1:=Position(CDlist,ToSymbol(d[1]));
      m2:=MP(d[2][1],d[2][2]);
      m3:=Sign(m2)*(Abs(m2)+(m1-1)*np);
      return m3;
   end function;

   B:=[];
   beta:=Cusp;
   if Cusp in [Cusps![w,2],Cusps![w-1,2]] then
      return qmap(-1*V.Position(EndPoints,Cusp));  
   elif Cusp eq Cusps![1,0] then
      return V2!0;
   else
      X:=Matrix([[1,0],[0,1]]);      
      i:=1;
         while beta notin C do            // C is defined in core.m 
            M:=Matrix([[1,0],[0,1]]);
            beta,N:=translation(beta);
            X:=N*X;
            if beta in C then 
            else
               M:=region(beta);
               beta:=ActC(M,beta);
               X:=M*X;
            end if;
            B[i]:=[X^(-1),M];
            i +:=1;
         end while;

      assert ActC(X,Cusp) in C;
      assert Abs(Z!Determinant(X)) eq 1;

      s:=ActC(X,Cusp);

      I:=Matrix([[1,0],[0,1]]);
      nB:=#B;

      if s eq Cusps![1,0] then

         D:=[<b[1],[ActC(b[2],s),s]>:b in B | b[2] ne I];
         F:=&+[V.num(d):d in D];  

         return qmap(F);

      else

         D:=[<b[1],[ActC(b[2],Cusps![1,0]),Cusps![1,0]]>:b in B | b[2] ne I];
         D cat:=[<B[nB][1],[Cusps![1,0],s]>];
         F:=&+[Sign(num(d))*V.Abs(num(d)):d in D];       
           
         return qmap(F);
      end if;
   end if;
end function;    

function Hecke(p,level_g,KI_g,PKI_g)   // given a prime ideal returns
                                       // the appropriate operator
   assert IsPrime(p) eq true;
   if Valuation(level_g,p) eq 0 then
      if IsPrincipal(p) eq true then
         return PrincipalHecke(p); 
      else
         return HeckeNP(p,level_g,KI_g,PKI_g);
      end if;
   else
      e:=Valuation(level_g,p);
      return AtkinLehner(level_g,p^e);   
   end if;
end function;

function Eigen(level_g,info_g,ns,primes,HM)   // calculates the eigenvalue of a
                                              // modular form for each operator
   P:=primes;

   KI_g:=info_g[1];
   PKI_g:=info_g[2];
   V:=info_g[3];
   V2:=info_g[4];  
   V3:=info_g[5];
   MSC:=info_g[6];
   CDlist:=info_g[7];
   qmap:=info_g[8];  

   eigen:=[];
   BKer:=Basis(ns);
   NKer:=#BKer;

   if NKer eq 1 then     // if the new space is 1 dimensional

// Now for each remaining modular symbol we check whether it occurs
// in the basis for the 1-D eiegenspace
// If it does we calculate its image under each of the operators

      for modsym in MSC do
         i:=Position(MSC,modsym);    
         if BKer[1][i] ne 0 then
            eigen[i]:=[&+[ModToCD(ActC(M,modsym[2]),level_g,
                          KI_g,PKI_g,CDlist,V,qmap,V2)-
                          ModToCD(ActC(M,modsym[1]),level_g,KI_g,PKI_g,
                          CDlist,V,qmap,V2):
                          M in HM[j]]:j in [1..#P]];
         end if;
      end for;

      function HeckeAction(x,p)   // returns the image of the generator of
                                  // a 1-D eigenspace under the operator 
                                  // corresponding to p
         j:=Position(primes,p);
         return &+[x[i]*eigen[i][j]:i in [1..#MSC] | x[i] ne 0];
      end function;
   
      values:=[* *];
        
      f2:=BKer[1];
      V4:=sub<V3|f2>;         // A 1-D subsapce of V3

      for p in P do;                     // Calculate the eigenvalues for
         eigen_v:=HeckeAction(f2,p);     // each operator
         if Valuation(level_g,p) gt 1 then
            if Coordinates(V4,eigen_v)[1] eq 1 then
               values[Position(P,p)]:="+0";
            else
               values[Position(P,p)]:="-0";
            end if;
         elif Valuation(level_g,p) eq 1 then
            values[Position(P,p)]:=-Coordinates(V4,eigen_v)[1];
         else        
            values[Position(P,p)]:=Coordinates(V4,eigen_v)[1];           
         end if;
      end for;     

      FM:=Fricke(level_g);

      feigen:=V2!0;
      for modsym in MSC do    
         i:=Position(MSC,modsym);      
         if f2[i] ne 0 then       
            feigen +:=f2[i]*(ModToCD(ActC(FM,modsym[2]),level_g,KI_g,PKI_g,CDlist,V,qmap,V2)-
                             ModToCD(ActC(FM,modsym[1]),level_g,KI_g,PKI_g,CDlist,V,qmap,V2));
         end if;
      end for;
 
      fi:=Coordinates(V4,feigen)[1]; 
     
      return values,fi;

   else // if the new space isn't 1 dimensional

      J:=1;
      while Valuation(level_g,primes[J]) ne 0 do
	 J+:=1;           
      end while;

      for modsym in MSC do
         i:=Position(MSC,modsym); 
         d:=0;
         for b in BKer do
            d:=d+Abs(b[i]);
         end for;
         if d ne 0 then        
            eigen[i]:=&+[ModToCD(ActC(M,modsym[2]),level_g,
                         KI_g,PKI_g,CDlist,V,qmap,V2)-
                         ModToCD(ActC(M,modsym[1]),level_g,KI_g,PKI_g,
                         CDlist,V,qmap,V2):
                         M in HM[J]];
         end if;
      end for;

      A:=ZeroMatrix(Q,NKer,NKer);   

      V4:=VectorSpace(Q,NKer);
      for b in BKer do     
         h:=Position(BKer,b);               
         eigen2:=&+[b[i]*eigen[i]:i in [1..#MSC] | b[i] ne 0];
         A[h]:=V4!Coordinates(ns,eigen2);    
      end for;

      E:=[e:e in Eigenvalues(A)];  // We identify the 1-D eigenspaces

      H:=[];
      if E eq [] then 
         return [];
      else                 
         H:=[Basis(Eigenspace(A,e[1])):e in E];     
      end if;        

      G:=[];
         
      for h in H do
         for hh in h do                  
            V6h:=&+[hh[i]*BKer[i]:i in [1..NKer]];
            G cat:=[V6h];
         end for;
      end for;
       
      ns_r:=sub<ns|G>;  // rational part of the newspace

      BKer_r:=Basis(ns_r); 
      NKer_r:=#BKer_r;  

      eigen:=[];

      for modsym in MSC do
         i:=Position(MSC,modsym); 
         d:=0;
         for b in BKer_r do
            d:=d+Abs(b[i]);
         end for;
         if d ne 0 then        
            eigen[i]:=[&+[ModToCD(ActC(M,modsym[2]),level_g,
                          KI_g,PKI_g,CDlist,V,qmap,V2)-
                          ModToCD(ActC(M,modsym[1]),level_g,KI_g,PKI_g,
                          CDlist,V,qmap,V2):
                          M in HM[j]]:j in [1..#primes]];
         end if;
      end for;

      if #H eq NKer_r then // if the space has already been split by first T_p

         r:=[**];  
         for f in H do 	                  
            V5g:=&+[f[1][i]*BKer_r[i]:i in [1..NKer_r]];
            V5:=sub<ns|V5g>;

            function HeckeAction(x,p)   // returns the image of the generator of
                                        // a 1-D eigenspace under the operator 
                                        // corresponding to p
               j:=Position(primes,p);
               return &+[x[i]*eigen[i][j]:i in [1..#MSC] | x[i] ne 0];
            end function;
   
            values:=[* *];

            for p in primes do                    
  	       eigen_v:=HeckeAction(V5g,p); 
               if Valuation(level_g,p) gt 1 then
                  if Coordinates(V5,eigen_v)[1] eq 1 then
                     values[Position(primes,p)]:="+0";
                  else
                     values[Position(primes,p)]:="-0";
                  end if;
               elif Valuation(level_g,p) eq 1 then
                  values[Position(primes,p)]:=-Coordinates(V5,eigen_v)[1];
               else        
                  values[Position(primes,p)]:=Coordinates(V5,eigen_v)[1];           
               end if;
            end for;     

            FM:=Fricke(level_g);

            feigen:=V2!0;
            for modsym in MSC do    
               i:=Position(MSC,modsym);      
               if V5g[i] ne 0 then       
                  feigen +:=V5g[i]*(ModToCD(ActC(FM,modsym[2]),level_g,KI_g,PKI_g,CDlist,V,qmap,V2)-
                                    ModToCD(ActC(FM,modsym[1]),level_g,KI_g,PKI_g,CDlist,V,qmap,V2));
               end if;
            end for;
 
            fi:=Coordinates(V5,feigen)[1]; 
                 
            r cat:=[*[*values,fi*]*];
           
         end for;

         return r;

      else // if we need to use more T_p to split the space

	 ns_us:=ns_r;
         J+:=1;
         r:=[**]; 

         while Dimension(ns_us) gt 0 do 

            BKer_s:=Basis(ns_us); 
            NKer_s:=#BKer_s;  
       
            HMJ:=HM[J]; 

            A:=ZeroMatrix(Q,NKer_s,NKer_s);   

            V4:=VectorSpace(Q,NKer_s);
            for b in BKer_s do     
               h:=Position(BKer_s,b);               
               eigen2:=&+[b[i]*eigen[i][J]:i in [1..#MSC] | b[i] ne 0];
               A[h]:=V4!Coordinates(ns_us,eigen2);    
            end for;

            E:=[e:e in Eigenvalues(A)];  // We identify the 1-D eigenspaces
            F:=[];
            H:=[];
            if E eq [] then 
            else
               F:=[Basis(Eigenspace(A,e[1]))[1]:e in E | e[2] eq 1];        
               H:=[Basis(Eigenspace(A,e[1])):e in E | e[2] gt 1];     
            end if;

         // and then proceed by working with these spaces in turn
                  
            for f in F do 
	                  
               V5g:=&+[f[i]*BKer_s[i]:i in [1..NKer_s]];
               V5:=sub<ns_us|V5g>;    

               function HeckeAction(x,p)   // returns the image of the generator of
                                           // a 1-D eigenspace under the operator 
                                           // corresponding to p
                  j:=Position(primes,p);
                  return &+[x[i]*eigen[i][j]:i in [1..#MSC] | x[i] ne 0];
               end function;
   
               values:=[* *];

               for p in primes do                    
  	          eigen_v:=HeckeAction(V5g,p);      
                  if Valuation(level_g,p) gt 1 then
                     if Coordinates(V5,eigen_v)[1] eq 1 then
                        values[Position(primes,p)]:="+0";
                     else
                        values[Position(primes,p)]:="-0";
                     end if;
                  elif Valuation(level_g,p) eq 1 then
                     values[Position(primes,p)]:=-Coordinates(V5,eigen_v)[1];
                  else        
                     values[Position(primes,p)]:=Coordinates(V5,eigen_v)[1];           
                  end if;
               end for;     

               FM:=Fricke(level_g);

               feigen:=V2!0;
               for modsym in MSC do    
                  i:=Position(MSC,modsym);      
                  if V5g[i] ne 0 then       
                     feigen +:=V5g[i]*(ModToCD(ActC(FM,modsym[2]),level_g,KI_g,PKI_g,CDlist,V,qmap,V2)-
                                       ModToCD(ActC(FM,modsym[1]),level_g,KI_g,PKI_g,CDlist,V,qmap,V2));
                  end if;
               end for;
 
               fi:=Coordinates(V5,feigen)[1]; 
                 
               r cat:=[*[*values,fi*]*];
         
            end for;

            G:=[];

            if H ne [] then
               for h in H do
                  for hh in h do                  
                     V6h:=&+[hh[i]*BKer_s[i]:i in [1..NKer_s]];
                     G cat:=[V6h];
                  end for;
               end for;
            end if;

            ns_us:=sub<ns_us|G>;
            J +:=1;

         end while;

         return r;

      end if;
   end if;

end function;



