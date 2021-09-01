complex_w := Conjugate(K.1,1);
CF:=Parent(complex_w);
complex:=hom<K -> CF | complex_w>;  
// complex(x) should be the same as Conjugate(x,1)

// ReJ(a+b*w) = a+b/2 (in Q)

function ReJ(z)
   return Q!Trace(K!z)/2;
end function;

// ImQ(a+b*w) = b*Sqrt(23)/2 (in C)

function ImQ(z)
   return Im(complex(z));
end function;

// w_part(a+b*w) = b

function w_part(z)
    return  Eltseq(K!z)[2];
end function;

// ???

function line(P1,P2,z)
   d:=P2-P1;
   y:=ImQ(P1)+ImQ(d)*(ReJ(z)-ReJ(P1))/ReJ(d);
   return y;
end function;

// The set of inversion matrices

MSet:=[
       Matrix([[-2*w+1,-w-5],[-4,2*w-3]]),
       Matrix([[2*w-1,6],[4,1-2*w]]),
       Matrix([[-2*w+1,w-6],[-4,2*w+1]]),
       Matrix([[w+1,3],[2-w,-w-1]]),
       Matrix([[w-3,w+2],[w+2,3-w]]),
       Matrix([[w+2,3-w],[3-w,-w-2]]),
       Matrix([[w-2,3],[w+1,2-w]]),
       Matrix([[0,-1],[1,0]])
];

//  The set of cusps, centres of the inversion matrices

CSet:=[
       Cusps![2*w-3,4],
       Cusps![2*w-1,4],
       Cusps![2*w+1,4],
       Cusps![3*w-5,8],
       Cusps![5*w-3,12],
       Cusps![5*w-2,12],
       Cusps![3*w+2,8],
       Cusps![0,1]
];

// This tests that each inversion matrix takes the corresponding
// central cusp to infinity and takes infinity to another central cusp.

for m in MSet do
   assert ActC(m,[-m[2,2],m[2,1]]) eq Cusps![1,0];
   assert ActC(m,Cusps![1,0]) in CSet;
end for;

CSet_perm:= [ 2, 2, 2, 4, 5, 6, 7, 8 ];
assert forall{i : i in [1..#CSet] | ActC(MSet[i],CSet[i]) eq Cusps![1,0]};
assert forall{i : i in [1..#CSet] | ActC(MSet[i],Cusps![1,0]) eq CSet[CSet_perm[i]]};

function region1(z)
   if ImQ(z) ge line((11*w-17)/23,(w-1)/2,z) then 
      return MSet[1];
   elif ImQ(z) ge line((22*w-11)/46,(w-1)/2,z) then
      return MSet[2];
   elif ImQ(z) le line((73*w-71)/184,(w-1)/2,z) and 
        ImQ(z) gt line((w-2)/3,(73*w-71)/184,z) then
      return MSet[4];
   elif ImQ(z) gt line((73*w-71)/184,(2*w-1)/5,z) then
      return MSet[5];
   else return MSet[8];
   end if;
end function;

function region2(z)
   if ImQ(z) ge line((11*w+6)/23,w/2,z) then 
      return MSet[3];
   elif ImQ(z) ge line((22*w-11)/46,w/2,z) then
      return MSet[2];
   elif (ImQ(z) le line((73*w-2)/184,w/2,z) and 
        ImQ(z) gt line((w+1)/3,(73*w-2)/184,z)) then
      return MSet[7];
   elif ImQ(z) gt line((73*w-2)/184,(2*w-1)/5,z) then
      return MSet[6];
   else return MSet[8];
   end if;
end function;

// for a cusp c=[z,1], with z in the fundamental region of C/OK,
// returns the inversion matrix which optimally reduces z.

// The fundamental region is defined by the conditions
//
//  -1/2 <= ReJ(z)    <= 1/2
//    0  <= w_part(z) <= 1/2
//
//  i.e., for z=a+b*w in K,  -1/2 <= (a+b/2) <= 1/2 and 0<=b<=1/2.
//
//  i.e., for z=x+y*i in C,  -1/2 <= x <= 1/2 and 0<=y<=sqrt(23)/4.

function region(c)
   assert c[2] eq 1;
   z:=c[1];
   assert ReJ(z) ge -1/2 and ReJ(z) lt 1/2; 
   assert w_part(z) ge 0 and w_part(z) le 1/2; 

   if ReJ(z) lt 0 then
      return region1(z);
   else
      return region2(z);   
   end if;
end function;

// for a cusp c=[z,1], with z in K, 
// returns c'=M(c), M where M is a generalised translation matrix
// and c' is in the fundamental region of C/OK.
//
// "Generalized" means that det(M) may be -1.

function translation(c)
   assert c[2] eq 1;
   z:=c[1];
   x:=w_part(z);
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
      return ActC(M,c),M;
   else
      M:=Matrix([[1,-x*w-y],[0,1]]);
      return ActC(M,c),M;
   end if;
end function;

// gm(z), for z a cusp, returns a single pseudo-Euclidean transform of
// z, i.e. a translation to the fundamental region of C/OK followed by
// the approproate inversion matrix which reduces the (generalised)
// denominator of z.

function gm(z)
   return region(translation(z));
end function;

// Pseudo-Euclidean algorithm on a pair a,b in OK
//
// Returns matrix X,c with X in GL(2,OK) and c in C such that X(a:b) = c
// where C={[1:0],[w:2],[w-1:2]}.

function EuAl(a,b)
   Cusp:=Cusps![a,b];
   beta:=Cusp;
   if Cusp in C then
      return Matrix([[1,0],[0,1]]);  
   else
      X:=Matrix([[1,0],[0,1]]);
      while beta notin C do
         beta,N:=translation(beta);
         X:=N*X;
         if beta in C then
         else
            M:=region(beta);
            beta:=ActC(M,beta);
            X:=M*X;
         end if;
      end while;
      assert ActC(X,Cusp) in C;
      assert Determinant(X) in [-1,1];
      return X,ActC(X,Cusp);
   end if;
end function;

// returns the generalized denominator of a/b,
// i.e. the norm of the denominator ideal of <a>/<b>, i.e. <b>/<a,b>.

function psi(a,b)
   return Norm(ideal<OK|b>)/Norm(ideal<OK|a,b>);
end function;

// The key idea is that if
// c=[a:b] and gm(c)=c'=[a',b'] then
// psi(a',b')<psi(a,b)
