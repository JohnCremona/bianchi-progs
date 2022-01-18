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

procedure dimtable(d, minnorm, maxnorm, include_zero)
  K := make_k(d);
  OK := Integers(K);
  print "K = ", K;
  Vdata := VoronoiData(BianchiCuspForms(K,1*OK));
  outputfilename := "dims" cat Sprintf("%o",d) cat "DY.txt";
  print "Output file is ", outputfilename;
  SetColumns(0);

  for n in [minnorm..maxnorm] do
     S := IdealsOfNorm(K, n);
     for i in [1..#S] do
           N := S[i];
           label:=IdealLabel(N);
           d := Dimension(BianchiCuspForms(K,N : VorData:=Vdata));
           if d gt 0 or include_zero then
              print label, d;
              PrintFile(outputfilename, Sprintf("%o %o",label, d));
           end if;
     end for;
  end for;
end procedure;


