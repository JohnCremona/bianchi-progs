Z:=IntegerRing();
Q:=Rationals();
Zx<x>:=PolynomialRing(Z);
K<w>:=NumberField(x^2-x+8);
OK:=IntegerRing(K);
Cusps:=ProjectiveSpace(K,1);
MK:=RMatrixSpace(OK,2,2);

CG,h:=ClassGroup(OK);   // CG is the ideal class group as an abstract
                        // group 0,CG.1,2*CG.1
hinv:=Inverse(h);       // maps an ideal to an element of CG
ClassNumber:=#CG;
cm:=ClassMap(CG);       // maps an element of CG to {1,2,3}

function CuspIdealClass(cusp)
   d:=Denominator(cusp[1]);
   return hinv(ideal<OK|d*cusp[1],d*cusp[2]>);
end function;

Conj:=hom<OK -> OK|1,1-w>;  // the action of complex conjugation on
                            // elements of OK

C:=[Cusps![1,0],Cusps![w,2],Cusps![w-1,2]];  // Standard choice of cusp
                                             // class representatives
                                             // under the action of the
                                             // full group

function ActC(M,C)        // Returns the image of a cusp under M
   a:=M[1,1]; b:=M[1,2]; 
   c:=M[2,1]; d:=M[2,2];
   if C[2] eq 0 then      
      return Cusps![a,c]; 
   else 
      z:=C[1]/C[2];      
      return Cusps![(a*z+b),(c*z+d)];
   end if; 
end function;

function Bezout(a,b)        // Given two coprime ideals finds a pair
                            // of elements whose sum is equal to 1
  assert a+b eq ideal<OK|1>;

  M:=RMatrixSpace(Z,4,2)!0;

  BasisA:=Basis(a);
  BasisB:=Basis(b);

  M[1][1]:=BasisA[1][1];
  M[1][2]:=BasisA[1][2];
  M[2][1]:=BasisA[2][1];
  M[2][2]:=BasisA[2][2];

  M[3][1]:=BasisB[1][1];
  M[3][2]:=BasisB[1][2];
  M[4][1]:=BasisB[2][1];
  M[4][2]:=BasisB[2][2];

  v:=Vector(Z,[1,0]);
  sol:=Solution(M,v);

  u1:=sol[1]*BasisA[1]+sol[2]*BasisA[2];
  u2:=sol[3]*BasisB[1]+sol[4]*BasisB[2];

  assert u1+u2 eq 1;
  assert u1 in a;
  assert u1-1 in b;
  assert u2 in b;
  assert u2-1 in a;

  return [u1,u2];

end function;

function CoprimeToIdeal(element,ideal)
   return ideal<OK|element>+ideal eq 1*OK;
end function;

function CoprimeToElement(elementa,elementb)
   return ideal<OK|elementa>+ideal<OK|elementb> eq 1*OK;
end function;

function Pivots(B)   // find pivots of reduced basis.
   return [Min(Support(b)): b in B];
end function;

function oi(p);      // function for outputing ideals in a nice form
   tf,gen:=IsPrincipal(p);
   if tf eq true then
      return <gen[1]+gen[2]*w>;
   else
      B:=Basis(p);
      return <K!B[1],K!B[2]>;
   end if;
end function;

/* Attach("egrosNF.m"); */
/* Attach("ks4.m"); */
/* SE:=ECgrS(K,[ideal<OK|2,w>,ideal<OK|2,1-w>,ideal<OK|7,w+4>,ideal<OK|7,w+2>]:verb:=true); */
/* SE2:=[e:e in SE | Norm(Conductor(e)) le 200]; */
/* SE2:=Sort(SE2,func<Ea,Eb|Norm(Conductor(Ea))-Norm(Conductor(Eb))>); */

/* load "curves3.m"; */
/* load "mina4a6.m"; */

/* for e in SE2 do */
/* a4,a6:=Explode(mina4a6(e)); */
/* e1:=EllipticCurve([0,0,0,a4,a6]); */
/* e1; */
/* Norm(Conductor(e1)); */
/* ToF(e1,P); */
/* end for; */

