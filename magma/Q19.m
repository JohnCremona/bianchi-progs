Qx<x> := PolynomialRing(Rationals());
K<w> := NumberField(x^2-x+5);
OK := Integers(K);
w:= OK.2;

function S2(N)
 return BianchiCuspForms(K,N*OK);
end function;

function dimS2(N)
  return Dimension(S2(N));
end function;

procedure dimtable(minnorm, maxnorm)
  for n in [minnorm..maxnorm] do
     fl, res := NormEquation(OK, n);
     if fl then
        for i in [1..#res] do
           N := res[i];
           d := dimS2(N);
           if d gt 0 then
              print n, K!N, d;
           end if;
        end for;
     end if;
  end for;
end procedure;
