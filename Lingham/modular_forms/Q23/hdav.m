info_g:=info;
level_g:=I;

   KI_g:=info_g[1];
   PKI_g:=info_g[2];
   V:=info_g[3];
   V2:=info_g[4];  
   V3:=info_g[5];
   MSC:=info_g[6];
   CDlist:=info_g[7];
   qmap:=info_g[8];  

   fd:=Dimension(V3);      

//   ns:=NewS(I,info);
   BKer:=Basis(ns); 
   snd:=Dimension(ns);       

   HM:=Hecke(ideal<OK|2,w>,level_g,KI_g,PKI_g); 

   eigen:=[];
                                           
   for modsym in MSC do
      i:=Position(MSC,modsym); 
      d:=0;
      for b in BKer do
         d:=d+Norm(b[i]);
      end for;
      if d ne 0 then
         eigen[i]:=V2!0;
         for M in HM do
            eigen[i] +:=ModToCD(ActC(M,modsym[2]),level_g,KI_g,PKI_g,CDlist,V,qmap,V2)-
            ModToCD(ActC(M,modsym[1]),level_g,KI_g,PKI_g,CDlist,V,qmap,V2);
         end for;
      end if;
   end for;

   L:=AlgebraicClosure();

   A:=ZeroMatrix(LL,snd,snd);   

   V4:=VectorSpace(LL,snd);
   for b in BKer do     
      h:=Position(BKer,b);               
      eigen2:=&+[b[i]*eigen[i]:i in [1..#MSC] | b[i] ne 0];
      A[h]:=V4!Coordinates(ns,eigen2);    
   end for;

   E:=[e:e in Eigenvalues(A)]; 
   H:=[];
   if E eq [] then 
      F:=[];
   else
      F:=[Basis(Eigenspace(A,e[1]))[1]:e in E | e[2] eq 1];        
      H:=[Basis(Eigenspace(A,e[1])):e in E | e[2] gt 1];     
   end if;

for f in F do
      V5g:=&+[f[i]*BKer[i]:i in [1..snd]];
      V5:=sub<ns|V5g>;

   P:=PrimeIdeals(50);

   HM:=[**];           // compile a list of all the operators in terms
                       // of matrices
   for p in P do
      HM cat:=[*Hecke(p,level_g,KI_g,PKI_g)*];
   end for;

// Now for each remaining modular symbol we check whether it occurs
// in the basis for the 1-D eiegenspace
// If it does we calculate its image under each of the operators

eigen:=[];

   for modsym in MSC do
      i:=Position(MSC,modsym);    
      eigen[i]:=[];
      if V5g[i] ne 0 then
         for j in [1..#P] do      
            eigen[i][j]:=V2!0;
            for M in HM[j] do
               eigen[i][j] +:=ModToCD(ActC(M,modsym[2]),level_g,
                              KI_g,PKI_g,CDlist,V,qmap,V2)-
                              ModToCD(ActC(M,modsym[1]),level_g,KI_g,PKI_g,
		              CDlist,V,qmap,V2);
            end for;
         end for;
      end if;
   end for;

   function HeckeAction(x,p)   // returns the image of the generator of
                               // a 1-D eigenspace under the operator 
                               // corresponding to p
      j:=Position(P,p);
      return &+[x[i]*eigen[i][j]:i in [1..#MSC] | x[i] ne 0];
   end function;
   
   values:=[* *];       

   for p in P do;                     // Calculate the eigenvalues for
      eigen_v:=HeckeAction(V5g,p);      // each operator
      if Valuation(level_g,p) gt 1 then
         if Coordinates(V5,eigen_v)[1] eq 1 then
            values[Position(P,p)]:="+0";
         else
            values[Position(P,p)]:="-0";
         end if;
      elif Valuation(level_g,p) eq 1 then
         values[Position(P,p)]:=-Coordinates(V5,eigen_v)[1];
      else        
         values[Position(P,p)]:=Coordinates(V5,eigen_v)[1];           
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
     
   values,fi;
   "";

end for;












