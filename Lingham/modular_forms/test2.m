load "ncore.m";

function ReJ(z)
   return Q!Trace(K!z)/2;
end function;

function w_part(z)
   d:=Denominator(z);
   return (OK!(d*z))[2]/d;
end function;

function NormList(bound)  // generates a list of elements of OK with norm leq 
                          // the bound
   return [K!a: a in 
   &join{&join{{K![(x-y)/2,y],K![(-x-y)/2,y],K![(x+y)/2,-y],K![(-x+y)/2,-y]}
	     : x in [0..Floor(Sqrt(4*bound-D*y^2))] | IsEven(x-y)}
             : y in [0..Floor(Sqrt(4*bound/D))]}];
end function;

function Pos(a) // for a in K determines whether a is "positive"
   return aa[1] gt 0 or aa[1] eq 0 and aa[2] ge 0 where aa:=Eltseq(a);
end function;

function NormList0(bound) // generates a list of "positive" elements of OK with
                          // norm leq the bound
   return [a : a in NormList(Q!bound) | Pos(a)];
end function;

function PrincipalCoveringLM(P) // Finds the centres of pricipal standard 
                                // hemispheres which cover the point P
   z:=K!P[1]; t:=Q!P[2];
   mulist:=[a : a in NormList0(1/t) | a ne 0];
   ans:=[];   
   for mu in mulist do
      bound:=Ceiling((Sqrt(Abs(1-t*Norm(mu)))+Sqrt(Norm(mu*z)))^2);         
      lambdalist:= NormList(bound);
      for lambda in lambdalist do
         if Norm(mu*z-lambda)+Norm(mu)*t le 1 and Norm(ideal<OK|lambda,mu>) eq 1
            then ans cat:=[[lambda,mu]];
         end if;
      end for;
   end for;
   return ans;
end function;

function Reduce(lambda,mu);      // Given a fraction lambda/mu it will
   I:=ideal<OK|lambda>;          // return a pair which give the same fraction
   J:=ideal<OK|mu>;              // and generate coprime ideals
   G:=GCD(I,J);
   IR:=I/G;
   JR:=J/G;
   flag1,g1:=IsPrincipal(IR);
   flag2,g2:=IsPrincipal(JR);
   L:=K!g1;
   M:=K!g2;
   if L/M eq lambda/mu then 
      return KK!<L,M>;
   else 
      return KK!<L,-M>;
   end if;
end function;

RF:=RealField(5);
CF:=ComplexField(5);
eval:=hom<L -> RF | RF!Sqrt(D)>;
evalC:=hom<L -> CF | CF.1*Sqrt(D)>;

function InitialList(bound)
   List:={};
   nl:=NormList(2*bound);
   n2:=NormList0(bound);
   n2:=[x:x in n2 | x ne 0];
   for x in nl do
      for y in n2 do
         if ReJ(x/y) ge -1 and ReJ(x/y) le 1 then
            if w_part(x/y) ge 0 and w_part(x/y) le 1 then
               if IsPrincipal(ideal<OK|x,y>) then
                  lm:=Reduce(x,y);
                  if ReJ(lm[1]/lm[2])+1/Sqrt(Norm(lm[2])) ge -1/2 and 
                     ReJ(lm[1]/lm[2])-1/Sqrt(Norm(lm[2])) le 1/2 then
                     if eval(iim(w_part(lm[1]/lm[2])))-1/Sqrt(Norm(lm[2])) 
                        le eval(v/b)                     
                        then List:=List join {lm};
                     end if;
                  end if;                
               end if;
            end if;
         end if;
      end for;
   end for;
   
   return List;

end function;

function CheckCovering(P,S)   // Checks whether a particular point is properly                              
   z:=K!P[1]; t:=Q!P[2];      // covered by one of the hemispheres in S
   ans:=[];
   for c in S do     
      lambda:=c[1];  
      mu:=c[2];
      if Norm(mu*z-lambda)+Norm(mu)*t lt 1
         then ans cat:=[[lambda,mu]];
      end if;
   end for;
   return ans;
   
end function;

function Check(S)      // Checks to see if the list contains any unecessary
   for c1 in S do      // hemispheres
      i:=0;
      I:=[];
      for c2 in (S diff {c1}) do
         for c3 in (S diff {c1,c2}) do

            n1:=Norm(c1[2]);
            n2:=Norm(c2[2]);
            n3:=Norm(c3[2]);

          eqn1:=(X-ReJ(c1[1]/c1[2]))^2+(Y-L!(w_part(c1[1]/c1[2])*v/(b/2)))^2+Z-1/n1;
          eqn2:=(X-ReJ(c2[1]/c2[2]))^2+(Y-L!(w_part(c2[1]/c2[2])*v/(b/2)))^2+Z-1/n2;
          eqn3:=(X-ReJ(c3[1]/c3[2]))^2+(Y-L!(w_part(c3[1]/c3[2])*v/(b/2)))^2+Z-1/n3;

            C:=Scheme(A,[eqn1,eqn2,eqn3]);

            if Dimension(C) eq 0 then
               P1:=Points(C);
               P2:=[p:p in P1 | Q!p[3] ge 0];
               if P2 ne [] then
                  p:=P2[1];
                  if p in I then
	          else                    
                     I cat:=[p];	    
                     if p[1] in Q and p[3] in Q and Q!p[3] gt 0 then
                        if Q!p[1] ge -1/2 and Q!p[1] le 1/2 then
                           if Q!(p[2]/v) ge 0 and Q!(p[2]/v) le 1/b then
  	        	      pc:=CheckCovering([im(L!(p[1]+p[2])),K!p[3]],S);
                              if pc eq [] then		 
  		                 i +:=1;
                              end if;
                           end if;
                        end if;
                     end if;
                  end if;
               end if;
            end if;
         end for;
      end for;
      if i lt 1 then
         Exclude (~S,c1);
      end if;
   end for;

   return S;

end function;

function Inter(S)     // Swan's algorithm                       
   I:=[];
   I2:=[];
   Si:=[];
   for c1 in S do
      for c2 in (S diff {c1}) do
         for c3 in (S diff {c1,c2}) do
            n1:=Norm(c1[2]);
            n2:=Norm(c2[2]);
            n3:=Norm(c3[2]);

          eqn1:=(X-ReJ(c1[1]/c1[2]))^2+(Y-L!(w_part(c1[1]/c1[2])*v/(b/2)))^2+Z-1/n1;
          eqn2:=(X-ReJ(c2[1]/c2[2]))^2+(Y-L!(w_part(c2[1]/c2[2])*v/(b/2)))^2+Z-1/n2;
          eqn3:=(X-ReJ(c3[1]/c3[2]))^2+(Y-L!(w_part(c3[1]/c3[2])*v/(b/2)))^2+Z-1/n3;

            C:=Scheme(A,[eqn1,eqn2,eqn3]);

            if Dimension(C) eq 0 then
               P1:=Points(C);
               P2:=[p:p in P1 | Q!p[3] ge 0];
               P3:=[p:p in P1 | Q!p[3] eq 0];

               if P3 ne [] then
                  p:=P3[1];                 
                  if Q!p[1] ge -1/2 and Q!p[1] le 1/2 then
                     if Q!(p[2]/v) ge 0 and Q!(p[2]/v) le 1/b then
		        z:=im(p[1]+p[2]);                                          
  		        Si cat:=[z];
                     end if;
                  end if;
               end if;

               if P2 ne [] then
                  p:=P2[1];
                  if p in I then
	          else                    
                     I cat:=[p];
                     if p[1] in Q and p[3] in Q and Q!p[3] gt 0 then
                        if Q!p[1] ge -1/2 and Q!p[1] le 1/2 then
                           if Q!(p[2]/v) ge 0 and Q!(p[2]/v) le 1/b then
			      x:=[im(L!(p[1]+p[2])),K!p[3]];
                              I2 cat:=[[K!x[1],K!x[2]]];
                              PC:=PrincipalCoveringLM(x);
                              S_p:=PC[1];
                              c:=S_p[1]/S_p[2];
                              n:=Norm(S_p[2]);
                              q:=im(p[1]+p[2]);  
                              z_1:=-(ReJ(q)-ReJ(c))^2-
                                   (L!(w_part(q)*v/(b/2)-L!(w_part(c)*v/(b/2))))^2+1/n;
                              for pc in PC do
                                 c:=pc[1]/pc[2];
                                 n:=Norm(pc[2]);
                                 q:=im(p[1]+p[2]);
                                 z:=-(ReJ(q)-ReJ(c))^2-
                                    (L!(w_part(q)*v/(b/2)-L!(w_part(c)*v/(b/2))))^2+1/n;
                                 if Q!z gt Q!z_1 then
                                    S_p:=pc;
                                    z_1:=z;
                                 end if;
                              end for;
                              if <S_p[1],S_p[2]> notin S then
"new point found";
			         S:=S join {<S_p[1],S_p[2]>};
                              end if;
                           end if;
                        end if;
                     end if;
                  end if;
               end if;
            end if;
         end for;
      end for;
   end for;

   return S,Set(I2),Set(Si);

end function;

function Boundary(I,S)   // Given a list of intersection points I this function                         
   I2:=[];               // throws out those that don't actually lie in the f.r. 
 
   for p in I do  
      if Q!Trace(K!p[1])/2 gt -1/2 then 
         pc:=CheckCovering(p,S);
         if pc eq [] then
            I2 cat:=[p];
         end if;
      end if;
   end for;

   return I2;

end function;

function SphereDenom(lm,P)    
   l:=K!lm[1]; m:=K!lm[2]; z:=K!P[1];
   return Norm(m*z-l)+P[2]*Norm(m);
end function;

conj:=hom<K -> K | 1-w>; 

function act(M,P)   // Returns the image of P under M
   z:=P[1]; t:=P[2];
   a:=K!M[1,1]; b:=K!M[1,2]; c:=K!M[2,1]; d:=K!M[2,2];
   det:=K!Determinant(M);
   denom:=SphereDenom([-d,c],P);
   if denom ne 0 
      then return <((a*z+b)*conj(c*z+d)+a*conj(c)*t)/denom,Norm(det)*t/denom^2>;
      else return Infinity();
   end if; 
end function;

function IsOnSphere(lm,P) // Returns whether P lies on the sphere centre lm
   return SphereDenom(lm,P) eq 1;
end function;

function AssociatedLM(P) // Finds the centres of standard hemispheres which go
                         // through the point P
   mulist:=[a : a in NormList0(1/P[2]) | a ne 0];
   ans:=[[1,0]];
   z:=P[1]; t:=Q!P[2];
   for mu in mulist do
      bound:=Ceiling((Sqrt(Abs(1-t*Norm(mu)))+Sqrt(Norm(mu*z)))^2); 
      lambdalist:= NormList(bound);
      for lambda in lambdalist do
         if Norm(mu*z-lambda)+Norm(mu)*t eq 1
            then ans cat:=[[lambda,mu]];
         end if;        
      end for;
   end for;
   return ans;
end function;

function Stab(P) // Finds the stabiliser of the point P
   LM:=AssociatedLM(P);
   mulist:=[lm[2] : lm in LM | lm[2] ne 0];    
   ans:=[GL(2,K)!Matrix([[1,0],[0,1]])];
   for mu in mulist do
      lambdalist:=[lm[1] : lm in LM | lm[2] eq mu];
      // look for pairs lambda1,lambda2 whose product is -1 or +1 mod mu
      for lambda1, lambda2 in lambdalist do
         beta:=-(lambda1*lambda2+1)/mu;
         if beta in OK then 
            ans cat:=[GL(2,K)!Matrix([[lambda1,beta],[mu,-lambda2]])]; 
         end if;
         beta:=-(lambda1*lambda2-1)/mu;
         if beta in OK 
            then ans cat:=[GL(2,K)!Matrix([[lambda1,beta],[mu,-lambda2]])]; 
         end if;
      end for;
   end for;
   cans:=[M : M in ans | act(M,P) eq P];
   sans:=Set(cans);
   return [M : M in sans];
end function;

function FromPtoQ(P,Q) // Finds the matrices that map P to Q 
   t:=P[2];
   if t ne Q[2] then return []; end if;
   if P eq Q then return Stab(P); end if;      
   LMP:=AssociatedLM(P);
   LMQ:=AssociatedLM(Q);
   mulist:={lm[2] : lm in LMP} meet {lm[2] : lm in LMQ};
   ans:=[];
   for mu in mulist do
      if mu eq 0 then // look for a translation
         beta:=Q[1]-P[1];
         if beta in OK 
            then ans cat:=[GL(2,K)!Matrix([[1,beta],[0,1]])];
         end if;
         else
         lambda1list:=[lm[1] : lm in LMP | lm[2] eq mu];
         lambda2list:=[lm[1] : lm in LMQ | lm[2] eq mu];
         // look for pairs lambda1,lambda2 whose product is -1 or +1 mod mu
         for lambda1 in lambda1list do
            for lambda2 in lambda2list do
               beta:=-(lambda1*lambda2+1)/mu;
               if beta in OK then 
                  ans cat:=[GL(2,K)!Matrix([[lambda2,beta],[mu,-lambda1]])]; 
               end if;
               beta:=-(lambda1*lambda2-1)/mu;
               if beta in OK then 
                  ans cat:=[GL(2,K)!Matrix([[lambda2,beta],[mu,-lambda1]])]; 
               end if;
            end for;
         end for;
      end if;
   end for;
return [M : M in ans | act(M,P) eq Q];
end function;

function ActC(M,C)        // Returns the image of a cusp under M
   a:=M[1,1]; b:=M[1,2]; 
   c:=M[2,1]; d:=M[2,2];
   if C[2] eq 0 then      
      return Cusps![a,c]; 
   else 
      z:=C[1];      
      return Cusps![(a*z+b),(c*z+d)];
   end if; 
end function;

function translation(c)
   assert c[2] eq 1;
   z:=c[1];
   x:=w_part(z)/w_part(w);
   if x ne 1/2 then
      x:=Round(x);
   else
      x:=0;
   end if;
   z -:=x*w;
   if w_part(z) lt 0 then
      z:=-z;
      u:=1;
   else
      u:=0;
   end if;

   y:=ReJ(z);
   y:=Floor(y+1/2);   
   z -:=y;
   if u eq 1 then
      M:=Matrix([[1,-y],[0,1]])*Matrix([[1,0],[0,-1]])*
         Matrix([[1,-x*w],[0,1]]);
      return <ActC(M,c),MK!M>;
   else
      M:=Matrix([[1,-x*w-y],[0,1]]);
      return <ActC(M,c),MK!M>;
   end if;
end function;

function InversionMatrix(c,NewCusps)  
// Construct the inversion matrix corresponding to the cusp c

   c:=Reduce(c[1],c[2]);

   gamma:=OK!c[1];
   delta:=OK!c[2]; 

   x1:=OK!(-gamma*w);
   x2:=OK!(-gamma);
   x3:=OK!(delta*w);
   x4:=OK!(delta);
     
   N:=Matrix(IR,4,2,[x1[1],x1[2],
                     x2[1],x2[2],
                    -x3[1],-x3[2],
             	    -x4[1],-x4[2]]);
   
   v:=Vector(IR,[1,0]);

   s:=Solution(N,v);  
   alpha:=s[1]*w+s[2];
   beta:=s[3]*w+s[4];

   M:=MK!Matrix([[alpha,beta],[delta,-gamma]]);
   assert Determinant(M) eq 1;

   t:=ActC(M,Cusps![c[1],c[2]]);
   assert t eq Cusps![1,0];

   nc:= ActC(M,[1,0]);
   N:=Matrix(OK,2,2,[[1,0],[0,1]]);
   if nc[2] eq 1 then
      p,N:=Explode(translation(nc));
      nc:=p;
   end if;

   assert nc in NewCusps;

   X:=MK!(N*M);

   assert ActC(X,Cusps![c[1],c[2]]) eq Cusps![1,0];
   assert ActC(X,Cusps![1,0]) in NewCusps;

   return X;

end function;

function actinfinity(M)   // Gives the image of infinity under the matrix M    
   a:=M[1,1]; c:=M[2,1]; 
   if c ne 0 
      then return Cusps![a,c];
      else return Cusps![1,0];
   end if; 
end function;

function Image(C1,C2,NewCusps)  // Finds a modular symbol of the form {c1,infinity} 
                                // equivalent to the symbol {C1,C2} and the matrix
                                // which maps one to the other 
                                // c1 lies in the fundamental region
     x:=C1[1];
     b:=Denominator(K!x);
     a:=x*b;
     y:=C2[1];
     d:=Denominator(K!y);
     c:=y*d;    

   if C2[2] eq 0 then
      t:=translation(C1);
      return <t[1],Cusps![1,0],t[2]^(-1)>;
   else

      A:=InversionMatrix([c,d],NewCusps);  
      z:=ActC(A,C1); 

      t:=translation(z);

      A:=t[2]*A; 

      c1:=ActC(A,C1);
      c2:=ActC(A,C2);

      assert ActC(A^(-1),c1) eq C1;
      assert actinfinity(A^(-1)) eq C2;

      return <c1,c2,A^(-1)>,1; 
  
   end if;
end function;

function ImageInf(C1,NewCusps) // Finds a modular symbol of the form 
                               // {e,infinity} equivalent to the symbol
                               // {infinity,C1} and the matrix
                               // which maps one to the other
                               // e lies in the fundamental region
   z:=C1[1];
   beta:=Denominator(z);
   alpha:=z*beta;
   A:=InversionMatrix([alpha,beta],NewCusps); 
   z:=actinfinity(A)[1];  
   d:=Denominator(z);
   t:=translation(Cusps![z*d,d]);
   A:=t[2]*A;
 
   return <actinfinity(A),ActC(A,C1),A^(-1)>,1; 
 
end function;

function Vorbits(V)      // divides the vertices of the f.r. into orbits
   II:=[];               // excluding the singular points
   for i in V do
      if i[2] ne 0 then
         II cat:=[KQ!<i[1],i[2]>];
      end if;
   end for;

   III:=[];
   for i in II do
      III cat:=[[i]];
   end for;

   for a in II do                 
      for b in II do              
         i:=Position(II,a);
         if a ne b then
            if FromPtoQ(a,b) ne [] then
               III[i] cat:=[b];         
            end if;
         end if;
      end for;
   end for;

   IV:=[Set(x):x in III];
   IV:=Set(IV);
   IV:=[x:x in IV];

   V:=[];

   for x in IV do
      V[Position(IV,x)]:=[y:y in x];
   end for;

   return V;
end function;

function RepPoly(VO,S)      // chooses representative polyhedra from each
   VII:=[];                 // orbit excluding those around the singular points
   for i in VO do           // and those which cross the boundary of the f.r. 
      for j in i do           
         VI:=[];               
         for n in S do
            if IsOnSphere(n,j) then	
               VI cat:=[Position(S,n)];
            end if;
         end for;
         if #VI ge 3 then
            VII cat:=[VI];
            break;
         end if;
      end for;
   end for;

   return VII;

end function;

function Vertices(P) // Given a point P this function calculates the centres of
                     // the principal hemispheres through it and then returns
                     // the points above these centres that lie on the covering
                     // surface 
   mulist:=[a : a in NormList0(1/P[2]) | a ne 0];
   ans:=[];
   z:=K!P[1]; t:=P[2];
   for mu in mulist do
      bound:=Ceiling((Sqrt(Abs(1-t*Norm(mu)))+Sqrt(Norm(mu*z)))^2);         
      lambdalist:= NormList(bound);
      for lambda in lambdalist do    
         if Norm(mu*z-lambda)+Norm(mu)*t le 1 and Norm(ideal<OK|lambda,mu>) eq 1
            then T:=(1/Norm(mu)); 
            ans cat:=[KQ!<lambda/mu,T>];
         end if;
      end for;
   end for;
   return ans;
end function;

time S1:=InitialList(16);
time S2:=Check(S1);

S2:={<-1,1>,<0,1>,<1,1>,<3,1-w>,<-3,w>,<w-2,3>,<w+1,3>,<2*w-3,4>,<2*w-1,4>,<2*w+1,4>,
     <2*w+1,4-w>,<2*w-3,w+3>,<w-6,w+3>,<w+5,4-w>,<5,1-w>,<-5,w>,<2*w-5,w+3>,<2*w+3,4-w>,
     <-3,w+1>,<w,3>,<3,2-w>,<w-1,3>}; 

S2c:=[<-1,1>,<0,1>,<1,1>,<-3,w>,<3,1-w>,<w-2,3>,<w+1,3>,<2*w-3,4>,<2*w-1,4>,<2*w+1,4>,
      <w-6,w+3>,<2*w+1,4-w>,<2*w-3,w+3>,<w+5,4-w>,<-5,w>,<5,1-w>,<2*w+3,4-w>,<2*w-5,w+3>,
      <w,2>,<w-1,2>];

S2a:=[Cusps![x[1],x[2]]:x in S2c];

S2b:=[x[1]/x[2]:x in S2a];

time S3,I1,SP1:=Inter(S2);
assert S2 eq S3;
SP2:=[x:x in SP1];
I2:=Boundary(I1,S3);               // Vertices of f.r.
S4:=[x:x in S3];                   // List of hemispheres forming the f.r.

S5:=[x:x in S4 | translation(Cusps![x[1],x[2]])[2] eq MI]; 

//S5b:=[x:x in S5a | x notin [<2*w-3,4>,<3,1-w>]];

NC:=[Cusps![s[1],s[2]]:s in S5];

// Inversion matrix corresponding to each hemisphere   
IM:=[<c,InversionMatrix(c,NC)>: c in S4];

S6:=[x[1]/x[2]:x in S5] cat SP2;         // List of generating edges
//S6:=[x[1]:x in EndPoints | x[2] ne 0];

VO:=Vorbits(I2);                         // Orbits of vertices under Gamma
RP:=RepPoly(VO,S2c);                     // Representative polyhedra


EndPoints:=[Cusps![2*w-3,4],
   	    Cusps![2*w-1,4],
            Cusps![9*w-13,20],     // (2*w+1)/(4-w)
            Cusps![9*w+4,20],      // (2*w-3)/(w+3)
            Cusps![w-2,3],
            Cusps![3*w-3,8],       // -3/w
            Cusps![3*w,8],         // 3/(1-w)
            Cusps![w,3],
  	    Cusps![0,1],
	    Cusps![w-1,3],
   	    Cusps![3,2-w],         // (3*w+3)/10
   	    Cusps![-3,w+1],        // (3*w-6)/10
            Cusps![w-1,2],
   	    Cusps![w,2],
	    Cusps![1,0]  	    
];

S6:=[x:x in EndPoints | x[2] ne 0];

function Rels(RepPoly)

for p in RepPoly do
<p[1],p[2],Infinity()>;
i,_,M1:=Explode(Image(S2a[p[1]],S2a[p[2]],NC));
a:=Position(S6,i);
i,_,M2:=Explode(Image(S2a[p[2]],Cusps![1,0],NC));
b:=Position(S6,i);
i,_,M3:=Explode(ImageInf(S2a[p[1]],NC));
c:=Position(S6,i);
printf "%o_%o\n",M1,a;
printf "%o_%o\n",M2,b;
printf "%o_%o\n",M3,c;
assert ActC(M1,S6[a]) eq Cusps![S2a[p[1]][1],S2a[p[1]][2]]; 
assert actinfinity(M1)[1] eq S2a[p[2]][1]/S2a[p[2]][2]; 
<ActC(M1,EndPoints[a]),actinfinity(M1)>;
<ActC(M2,EndPoints[b]),actinfinity(M2)>;
<ActC(M3,EndPoints[c]),actinfinity(M3)>;
"";
<p[2],p[3],Infinity()>;
i,_,M1:=Explode(Image(S2a[p[2]],S2a[p[3]],NC));
a:=Position(S6,i);
i,_,M2:=Explode(Image(S2a[p[3]],Cusps![1,0],NC));
b:=Position(S6,i);
i,_,M3:=Explode(ImageInf(S2a[p[2]],NC));
c:=Position(S6,i);
printf "%o_%o\n",M1,a;
printf "%o_%o\n",M2,b;
printf "%o_%o\n",M3,c;
<ActC(M1,EndPoints[a]),actinfinity(M1)>;
<ActC(M2,EndPoints[b]),actinfinity(M2)>;
<ActC(M3,EndPoints[c]),actinfinity(M3)>;
"";
<p[3],p[1],Infinity()>;
i,_,M1:=Explode(Image(S2a[p[3]],S2a[p[1]],NC));
a:=Position(S6,i);
i,_,M2:=Explode(Image(S2a[p[1]],Cusps![1,0],NC));
b:=Position(S6,i);
i,_,M3:=Explode(ImageInf(S2a[p[3]],NC));
c:=Position(S6,i);
printf "%o_%o\n",M1,a;
printf "%o_%o\n",M2,b;
printf "%o_%o\n",M3,c;
<ActC(M1,EndPoints[a]),actinfinity(M1)>;
<ActC(M2,EndPoints[b]),actinfinity(M2)>;
<ActC(M3,EndPoints[c]),actinfinity(M3)>;
"";
"";
end for;
return 0;
end function;


for v in VO do
   vv:=v[1];
   vv;
   "";
   for x in v do
      B:=Boundary(Vertices(x),S3);
      M:=FromPtoQ(x,vv);
      for m in M do
         for b in B do
            act(m,b);
         end for;
         actinfinity(m);
         "";
      end for;
   end for;
   "";
   "";
end for;

//Rels(RP);





