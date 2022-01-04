// Define D to be a positive squarefree integer before loading this

//D:=31;
IR:=IntegerRing();
Q:=Rationals();
Zx<x>:=PolynomialRing(IR);
L<v>:=NumberField(x^2-D);

if D mod 4 eq 3 then
   L2:=NumberField(x^2+D);
   OL:=IntegerRing(L2);
   f:=MinimalPolynomial(OL.2,IR);
   K<w>:=NumberField(f);
   OK:=IntegerRing(K);  
   im:=hom<L -> K | 2*w-1>;
   b:=4;
   conj:=hom<K -> K | 1-w>; 
else
   K<w>:=NumberField(x^2+D);
   OK:=IntegerRing(K);
   im:=hom<L -> K | w>;
   b:=2;
   conj:=hom<K -> K | -w>; 
end if;

iim:=Inverse(im);
Cusps:=ProjectiveSpace(K,1);
A<X,Y,Z>:=AffineSpace(L,3);
MK:=RMatrixSpace(K,2,2);
KQ:=CartesianProduct(K,Q);
KK:=CartesianProduct(K,K);
MI:=Matrix([[1,0],[0,1]]);
MJ:=Matrix([[1,0],[0,-1]]);

