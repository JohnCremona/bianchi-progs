load "core.m";
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

E1:=EllipticCurve([w,-1,0,-w-6,0]);                             //6
E2:=EllipticCurve([w,-w+1,1,-1,0]);                             //26
E3:=EllipticCurve([1,-w-1,1,1,w-2]);                            //27
E4:=EllipticCurve([0,-w,0,4*w-1,-5]);                           //32
E5:=EllipticCurve([w+1,-w,w+1,-w-1,0]);                         //36
E6:=EllipticCurve([0,w-1,0,5*w+6,0]);                           //48
E7:=EllipticCurve([w+1,1,0,-1,0]);                              //52
E8:=EllipticCurve([w+1,w-1,w,-2*w,0]);                          //54
E9:=EllipticCurve([0,w,w+1,2*w-6,-4]);                          //54
E10:=EllipticCurve([w+1,-1,0,-2-2*w,3-w]);                      //62
E11:=EllipticCurve([0,-w-1,0,3,-2*w-5]);                        //64
E12:=EllipticCurve([w+1,1,w+1,0,0]);                            //72
E13:=EllipticCurve([w,w+1,0,-4-2*w,3+w]);                       //72
E14:=EllipticCurve([w+1,w+1,1,-w-4,w+5]);                       //78
E15:=EllipticCurve([w,0,1,6,-2*w-4]);                           //82
E16:=EllipticCurve([K!1,-1,0,-10,-12]);                         //92
E17:=EllipticCurve([0,-w-1,0,5*w-22,-6*w+48]);                  //96
E18:=EllipticCurve([w,w-1,0,-w-1,0]);                           //96
E19:=EllipticCurve([w,w,w,-w-5,-2*w-4]);                        //96
E20:=EllipticCurve([w+1,-w-1,0,-w,3*w+1]);                      //108
E21:=EllipticCurve([0,0,w,2*w-3,w+3]);                          //108
E22:=EllipticCurve([0,0,0,552*w+1221,4888*w-34762]);            //108
E23:=EllipticCurve([1,w-1,0,2*w-3,5]);                          //108
E24:=EllipticCurve([0,0,0,-53160*w-43995,-5067640*w+19402006]); //108
E25:=EllipticCurve([w+1,0,w+1,2*w-1,-6*w-8]);                   //116
E26:=EllipticCurve([K|0,-1,1,0,0]);                             //121
E27:=EllipticCurve([1,-1+w,1-w,3+w,-9+7*w]);                    //124
E28:=EllipticCurve([0,w+1,0,w-2,-1]);                           //128
E29:=EllipticCurve([1,w-1,w+1,-4,-w]);                          //141
E30:=EllipticCurve([-w+2,w,-w+2,3*w+3,2*w+4]);                  //144
E31:=EllipticCurve([w,w-1,0,1,0]);                              //144
E32:=EllipticCurve([0,0,0,-18927*w-14202,1857222*w-1211004]);   //144
E33:=EllipticCurve([0,w,0,-w-3,-2*w+3]);                        //144
E34:=EllipticCurve([w+1,-1,w+1,-5*w+3,2*w+2]);                  //144
E35:=EllipticCurve([0,0,0,864*w-26811,-95472*w+1553094]);       //144
E36:=EllipticCurve([1,0,w,-2,-2*w+4]);                          //156
E37:=EllipticCurve([w+1,w-1,0,-2*w-6,-2*w+6]);                  //162
E38:=EllipticCurve([w+1,-1,1,-4*w-2,w+5]);                      //162
E39:=EllipticCurve([w+1,w+1,1,-7,-3*w+5]);                      //162
E40:=EllipticCurve([1,-1,w,-2*w,4]);                            //162
E41:=EllipticCurve([0,-w-1,1,w-2,1]);                           //177
E42:=EllipticCurve([w,1-w,0,-4*w-4,0]);                         //192
E43:=EllipticCurve([w,0,w,-1-w,-w]);                            //192
E44:=EllipticCurve([w,-1,w,2,1-w]);                             //192
E45:=EllipticCurve([0,w-2,0,-8,-w-5]);                          //192
E46:=EllipticCurve([K|1,0,1,-1,0]);                             //196

E:=[E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11,E12,E13,E14,E15,E16,E17,E18,E19,E20,
    E21,E22,E23,E24,E25,E26,E27,E28,E29,E30,E31,E32,E33,E34,
    E35,E36,E37,E38,E39,E40,E41,E42,E43,E44,E45,E46];

P:=PrimeIdeals(50);

for e in E do
   e;
   printf "Conductor: %o\n",oi(Conductor(e));
   printf "Rank: %o\n",RankBound(e);
   ToF(e,P);
   ""; 
end for;



