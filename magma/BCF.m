// before loading this define K and OK

load "sorting.m";

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


