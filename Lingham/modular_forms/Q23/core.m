Z:=IntegerRing();
Q:=Rationals();
Zx<x>:=PolynomialRing(Z);
K<w>:=NumberField(x^2-x+6);
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
    return Automorphism(Cusps,Transpose(M))(C);
end function;

function Bezout(a,b)        // Given two coprime ideals finds a pair
                            // of elements whose sum is equal to 1
  assert a+b eq ideal<OK|1>;

  BasisA:=Basis(a);
  BasisB:=Basis(b);

  M := RMatrixSpace(Z,4,2) !
      &cat([Eltseq(x) : x in BasisA] cat [Eltseq(x) : x in BasisB]);

  v:=Vector(Z,[1,0]);
  sol:=Solution(M,v);

  u1:=sol[1]*BasisA[1]+sol[2]*BasisA[2];
  u2:=sol[3]*BasisB[1]+sol[4]*BasisB[2];

  assert u1+u2 eq 1;
  assert u1 in a;
  assert u2 in b;

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
