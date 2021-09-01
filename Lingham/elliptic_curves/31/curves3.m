Z:=IntegerRing();
Q:=Rationals();
Zx<x>:=PolynomialRing(Z);
K<w>:=NumberField(x^2-x+8);
OK:=IntegerRing(K);
CG,h:=ClassGroup(OK);   // CG is the ideal class group as an abstract
                        // group 0,CG.1,2*CG.1
hinv:=Inverse(h);       // maps an ideal to an element of CG
ClassNumber:=#CG;
cm:=ClassMap(CG);       // maps an element of CG to {1,2,3}

Conj:=hom<OK -> OK|1,1-w>;  // the action of complex conjugation on
                            // elements of OK

function oi(p);      // function for outputing ideals in a nice form
   tf,gen:=IsPrincipal(p);
   if tf eq true then
      return <gen[1]+gen[2]*w>;
   else
      B:=Basis(p);
      return <K!B[1],K!B[2]>;
   end if;
end function;

load "levels.m";

function Denom(x)
  den1 := LCM([Denominator(c) : c in Eltseq(x)]);
  x1 := OK!(den1*x);
  num := GCD(den1, GCD(ChangeUniverse(Eltseq(x1), Integers())));
  return ExactQuotient(den1, num), OK!(x1/num);
end function;

function MyResidueClassField(pid)
  F, m := ResidueClassField(OK, pid);
  p := Minimum(pid);
  e := ChineseRemainderTheorem(pid, ideal<OK|OK!p> / pid^RamificationIndex(pid),
                               OK!1, OK!0);
  f := function(x)
         
     if x eq 0 then 
        return F!0; 
     end if;

     v := Valuation(x, pid);
    
     if v gt 0 then 
        return F!0; 
     end if;

     den := Denom(x);
     v := Valuation(den, p);
     y := x * e^v; // y is now p-integral and = x mod pid
     den := Denom(y);

     return m(OK!(den*y))/F!den;

  end function;

  return F, map< K -> F | x :-> f(x), y :-> y @@ m >;

end function;


function ToF_s(E,p)
   
   rf,rmap:=MyResidueClassField(p);

   T,E_m:=LocalInformation(E,p);
  
   a1,a2,a3,a4,a6:=Explode(Coefficients(E_m));  
   a1:=rmap(a1);
   a2:=rmap(a2);
   a3:=rmap(a3);
   a4:=rmap(a4);
   a6:=rmap(a6);
  
   A<y,z>:=AffineSpace(rf,2);
   M<y,z>:=PolynomialRing(rf,2);
  
   f:=y^2+a1*z*y+a3*y-z^3-a2*z^2-a4*z-a6;
   SE:=Curve(A,f);   

   np:=1;
       
   for a in rf do
      for b in rf do
         if [a,b] in SE eq true then
            np +:=1;
         end if;
      end for;
   end for;

   ToF:=1+Norm(p)-np;

   return ToF;

end function;

function ToF(E,P)
   Cond:=Conductor(E);
   values:=[];  
  
   for p in P do
      if Valuation(Cond,p) eq 0 then
         T,E_m:=LocalInformation(E,p);
         values[Position(P,p)]:=TraceOfFrobenius(Reduction(E_m,p));
      else	 
	 values[Position(P,p)]:=ToF_s(E,p);
      end if;
   end for;
 
   return values;
 
end function;

function Econj(E)
   a1,a2,a3,a4,a6:=Explode(Coefficients(E));
   a1:=K!Conj(a1);
   a2:=K!Conj(a2);
   a3:=K!Conj(a3);
   a4:=K!Conj(a4);
   a6:=K!Conj(a6);
   E2:=EllipticCurve([a1,a2,a3,a4,a6]);

   return E2;

end function;


E1:=EllipticCurve([1,-w-1,1,w-3,1]);                          // 10
E2:=EllipticCurve([w+1,0,1,-w-1,0]);                          // 14
E3:=EllipticCurve([w,1-w,w,-w,0]);                            // 20
E4:=EllipticCurve([0,-1,0,3-w,-3]);                           // 32
E5:=EllipticCurve([0,0,0,621,-2916*w-10638]);                 // 32
E6:=EllipticCurve([w+1,1,1,-2*w-4,2-w]);                      // 49
E7:=EllipticCurve([0,0,0,-25920*w-10395,-3981312*w+3297078]); // 50
E8:=EllipticCurve([0,w+1,0,w-2,-1]);                          // 64
E9:=EllipticCurve([1-w,1-w,1,w,1-w]);                         // 70
E10:=EllipticCurve([w+1,0,1,0,0]);                            // 70
E11:=EllipticCurve([0,0,0,-324*w+3429,-25812*w-17226]);       // 80
E12:=EllipticCurve([0,0,0,1917,-11664*w+5886]);               // 80
E13:=EllipticCurve([w,w,w,-w-3,2-2*w]);                       // 80
E14:=EllipticCurve([0,0,0,1728*w-18603,-95364*w+495882]);     // 80
E15:=EllipticCurve([w+1,w-1,0,-2*w-3,3]);                     // 82
E16:=EllipticCurve([1,w+1,1,w-2,-1]);                         // 98
E17:=EllipticCurve([0,0,0,-675*w-4968,28782*w+123984]);       // 100
E18:=EllipticCurve([0,w-1,w,-5,1-w]);                         // 100
E19:=EllipticCurve([0,0,0,648*w+216,972*w-24732]);            // 100
E20:=EllipticCurve([0,0,0,-1485*w-29835,-77166*w-1798146]);   // 100
E21:=EllipticCurve([0,0,2-w,3*w-8,-w-1]);                     // 112
E22:=EllipticCurve([2-w,1,0,-4,w]);                           // 112
E23:=EllipticCurve([K|0,-1,1,0,0]);                           // 121
E24:=EllipticCurve([K|1,-1,1,-1,1]);                          // 124
E25:=EllipticCurve([1-w,-1,0,-1,3-2*w]);                      // 128
E26:=EllipticCurve([0,w+1,0,w,0]);                            // 128
E27:=EllipticCurve([1,1,1,w-5,2*w-5]);                        // 140
E28:=EllipticCurve([1,1,1,-w,1-w]);                           // 142
E29:=EllipticCurve([0,0,0,-1971*w-3672,-78354*w-22896]);      // 144
E30:=EllipticCurve([w+1,1,0,w,0]);                            // 160
E31:=EllipticCurve([w+1,1,w+1,-2*w+2,3-w]);                   // 160
E32:=EllipticCurve([0,0,0,-1323*w+5400,-347382*w-148176]);    // 160
E33:=EllipticCurve([1,w,1,1,2*w+2]);                          // 190
E34:=EllipticCurve([K|1,0,1,-1,0]);                           // 196
E35:=EllipticCurve([0,0,0,17739*w+35613,-715878*w+6711606]);  // 196
E36:=EllipticCurve([0,0,0,-5427*w-19224,-2419794*w+3211920]); // 200
E37:=EllipticCurve([0,0,0,108*w+756,29052*w+44388]);          // 200

E:=[E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11,E12,E13,E14,E15,E16,E17,E18,E19,E20,E21,E22,E23,E24,
    E25,E26,E27,E28,E29,E30,E31,E32,E33,E34,E35,E36,E37];

P:=PrimeIdeals(50);

for e in E do
   e;
   printf "Conductor: %o\n",oi(Conductor(e));
   printf "Rank: %o\n",RankBound(e);
   ToF(e,P);
   ""; 
end for;
