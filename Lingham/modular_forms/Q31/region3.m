CF:=ComplexField(30);
complex:=hom<K -> CF | (1+Sqrt(31)*CF.1)/2>;  

function ReJ(z)
   return Q!Trace(K!z)/2;
end function;

function ImQ(z)
   return Im(complex(z));
end function;

function w_part(z)
   d:=Denominator(z);
   return (OK!(d*z))[2]/d;
end function;

function line(P1,P2,z)
   d:=P2-P1;
   y:=ImQ(P1)+ImQ(d)*(ReJ(z)-ReJ(P1))/ReJ(d);
   return y;
end function;

MSet:=[
       Matrix([[2*w-1,w+7],[4,3-2*w]]),
       Matrix([[2*w-3,w+7],[4,1-2*w]]),
       Matrix([[2*w-1,8-w],[4,-2*w-1]]),
       Matrix([[-2*w-1,-8],[w-4,2*w+1]]),
       Matrix([[3-2*w,-8],[-w-3,2*w-3]]),
       Matrix([[w-2,w+1],[3,2-w]]),
       Matrix([[3,1-w],[-w,-3]]),
       Matrix([[-3,w],[w-1,3]]),
       Matrix([[w-2,3],[3,-w-1]]),
       Matrix([[0,-1],[1,0]])
];

CSet:=[
       Cusps![2*w-3,4],
       Cusps![2*w-1,4],
       Cusps![2*w+1,4],
       Cusps![9*w-13,20],     // (2*w+1)/(4-w)
       Cusps![9*w+4,20],      // (2*w-3)/(w+3)
       Cusps![w-2,3],
       Cusps![3*w-3,8],       // -3/w
       Cusps![3*w,8],         // 3/(1-w)
       Cusps![w+1,3],
       Cusps![0,1]
];

function region1(z)
   if ImQ(z) ge line((28*w-45)/62,(w-1)/2,z) then 
      return MSet[1];
   elif ImQ(z) ge line((14*w-7)/31,(w-1)/2,z) then
      return MSet[2];
   elif ImQ(z) ge line((3*w-5)/7,(105*w-161)/248,z) and 
        ImQ(z) ge line((105*w-161)/248,(w-1)/2,z) then
      return MSet[4];
   elif ImQ(z) gt line((9*w-20)/31,(30*w-46)/93,z) and 
        ImQ(z) le line((105*w-161)/248,(30*w-46)/93,z)  then
      return MSet[6];
   elif ImQ(z) gt line((30*w-46)/93,(32*w-16)/93,z) then
      return MSet[7];
   else return MSet[10];
   end if;
end function;

function region2(z)
   if ImQ(z) ge line((28*w+17)/62,w/2,z) then 
      return MSet[3];
   elif ImQ(z) ge line((14*w-7)/31,w/2,z) then
      return MSet[2];
   elif ImQ(z) ge line((3*w+2)/7,(105*w+56)/248,z) and 
        ImQ(z) ge line((105*w+56)/248,w/2,z) then
      return MSet[5];
   elif ImQ(z) gt line((9*w+11)/31,(30*w+16)/93,z) and 
        ImQ(z) le line((105*w+56)/248,(30*w+16)/93,z)  then
      return MSet[9];
   elif ImQ(z) gt line((30*w+16)/93,(32*w-16)/93,z) then
      return MSet[8];
   else return MSet[10];
   end if;
end function;

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
      return ActC(M,c),M;
   else
      M:=Matrix([[1,-x*w-y],[0,1]]);
      return ActC(M,c),M;
   end if;
end function;

function gm(z)
   return region(translation(z));
end function;

for m in MSet do
   assert ActC(m,[-m[2,2],m[2,1]]) eq Cusps![1,0];
   assert ActC(m,Cusps![1,0]) in CSet;
end for;

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

L<v>:=NumberField(x^2-31);
im:=hom<L -> K | 2*w-1>;
R:=RealField(5);
evalu:=hom<L -> R | R!Sqrt(31)>;

function test(x)
z:=im(x[1]+v*x[2]);
evalu(x[1]);
evalu(v*x[2]);
return CSet[Position(MSet,gm([z,1]))];
end function;
