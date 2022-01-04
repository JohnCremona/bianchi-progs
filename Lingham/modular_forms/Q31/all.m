load "core3.m";
load "input3.m";
//load "levels.m";
load "homology_new.m"; // loads "input.m"
//load "hecke3.m";    // loads "homology.m", "region.m"
//load "bianchi3.m";

load "sorting.m";

procedure dimtable(minnorm, maxnorm, include_zero)
  for n in [minnorm..maxnorm] do
     S := IdealsOfNorm(K, n);
     for i in [1..#S] do
           N := S[i];
           label:=IdealLabel(N);
           d := DimH(N, false);
           if d gt 0 or include_zero then
              print label, d;
           end if;
     end for;
  end for;
end procedure;
