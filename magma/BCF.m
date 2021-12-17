// before loading this define K and OK

load "sorting.m";

Qx<x> := PolynomialRing(Rationals());

function make_k(d)
  if d mod 4 eq 0 then d:= d/4; end if;
  if d lt 0 then
  if d mod 4 eq 1 then
    return NumberField(x^2-x+(1-d)/4);
  else
    return NumberField(x^2-d);
  end if;
  else
   if d mod 4 eq 3 then
     return NumberField(x^2-x+(1+d)/4);
   else
     return NumberField(x^2+d);
   end if;
 end if;
end function;

K := make_k(d);
OK := Integers(K);
print "K = ", K;

Vdata := VoronoiData(BianchiCuspForms(K,1*OK));

function S2(N)
  return BianchiCuspForms(K,N : VorData:=Vdata);
end function;

function dimS2(N)
  return Dimension(S2(N));
end function;

procedure dimtable(minnorm, maxnorm, include_zero)
  for n in [minnorm..maxnorm] do
     S := IdealsOfNorm(K, n);
     for i in [1..#S] do
           N := S[i];
           label:=IdealLabel(N);
           d := dimS2(N);
           if d gt 0 or include_zero then
              print label, d;
           end if;
     end for;
  end for;
end procedure;


