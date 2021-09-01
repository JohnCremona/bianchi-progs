function degen(N,d)          // calculates a matrix representation
                             // of the degeneracvy map from level
                             // N to level N/d

     if cm(hinv(d)) eq 1 then  // if d is principal we can use the
                               // classical definition

        tf,t:=IsPrincipal(d);
        return Matrix([[t,0],[0,1]]);
 
     else

        q1:=ideal<OK|193,w+55>;
        q2:=ideal<OK|193,w+137>;
 
        q:=hinv(d) eq hinv(q1) select q1 else q2;     

        tf,alpha:=IsPrincipal(q/d);    

        M1:=Matrix([[1,0],[0,K!alpha]]);

        a,c:=Explode(Bezout(q,N*d));     

        M2:=Matrix([[a,-1],[c,1]]);

        assert Determinant(M2) eq 1;    

        tf,beta:=IsPrincipal((q*d^2)^(-1));

        M3:=Matrix([[K!beta,0],[0,1]]);

        M:=M1*M2*M3;

        assert M1[1][1]*d+M1[2][1]*d eq d;
        assert M1[1][2]*d+M1[2][2]*d eq q;
        assert M2[1][1]*d+M2[2][1]*q eq d*q;
        assert M2[1][2]*d+M2[2][2]*q eq 1*OK;
        assert M3[1][1]*d*q+M3[2][1]*OK eq d^(-1);
        assert M3[1][2]*d*q+M3[2][2]*OK eq 1*OK;    
        assert M[1][1]*d+M[2][1]*d eq d^(-1);
        assert M[1][2]*d+M[2][2]*d eq 1*OK;     

        return M^(-1);

     end if;

end function;

function NewS(level_N,Ninfo)  // returns the newpart of V3

   KI_N:=Ninfo[1];
   PKI_N:=Ninfo[2];
   V_N:=Ninfo[3];
   V2_N:=Ninfo[4];  
   V3_N:=Ninfo[5];
   MSC_N:=Ninfo[6];
   CDlist_N:=Ninfo[7];
   qmap_N:=Ninfo[8];  

   NewSpace:=V3_N;

   for p in Support(level_N) do

      level_M:=level_N/p;     

      Minfo:=DimH(level_M,1);    

      KI_M:=Minfo[1];
      PKI_M:=Minfo[2];
      V_M:=Minfo[3];
      V2_M:=Minfo[4];  
      V3_M:=Minfo[5];
      MSC_M:=Minfo[6];
      CDlist_M:=Minfo[7];
      qmap_M:=Minfo[8];       

      B:=Basis(V3_N);
      image_p:=[];
      image_1:=[];

      A:=degen(level_N,p);

      for b in B do
         i:=Position(B,b); 
         image_p[i]:=V2_M!0;
         for modsym in MSC_N do
    	    cof:= b[Position(MSC_N,modsym)];
            if cof ne 0 then
               cd1:=ModToCD(ActC(A,modsym[1]),level_M,KI_M,PKI_M,CDlist_M,V_M,qmap_M,V2_M);
               cd2:=ModToCD(ActC(A,modsym[2]),level_M,KI_M,PKI_M,CDlist_M,V_M,qmap_M,V2_M);
               image_p[i] +:=cof*(cd2-cd1);
            end if;
         end for;
      end for;

      for b in B do
         i:=Position(B,b); 
         image_1[i]:=V2_M!0;
         for modsym in MSC_N do
            cof:= b[Position(MSC_N,modsym)];
            if cof ne 0 then
               cd1:=ModToCD(modsym[1],level_M,KI_M,PKI_M,CDlist_M,V_M,qmap_M,V2_M);
               cd2:=ModToCD(modsym[2],level_M,KI_M,PKI_M,CDlist_M,V_M,qmap_M,V2_M);
               image_1[i] +:=cof*(cd2-cd1);
            end if;
         end for;      
      end for;

      degenmap_p:=hom<V3_N -> V3_M | image_p>;
      degenmap_1:=hom<V3_N -> V3_M | image_1>;

      NewSpace:=NewSpace meet Kernel(degenmap_p) meet Kernel(degenmap_1);
   
   end for;

   return NewSpace;

end function;
