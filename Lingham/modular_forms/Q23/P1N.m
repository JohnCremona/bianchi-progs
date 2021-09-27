// Given an integral ideal N, returns a list of elements of P^1(N) as
// pairs of elements of O (the order of N)

function P1Nlist(N)
  O := Order(N);
  one := O!1;
  if Norm(N) eq 1 then return [[O!0, one]]; end if;
  FN := Factorization(N);
  if #FN eq 1 then
    OmodN := quo<O|N>;
    res := [O!x: x in OmodN];
    nires := [O!x: x in OmodN | not IsUnit(x)];
    return [[x,one]: x in res] cat [[1,x]: x in nires];
  else
    Q := FN[1][1]^FN[1][2];
    M := N/Q;
    P1Q := P1Nlist(Q);
    P1M := P1Nlist(M);
    return &cat[[[CRT(Q,M,y[1],x[1]), CRT(Q,M,y[2],x[2])]: x in P1M]: y in P1Q];
  end if;
end function;

// Given an integral ideal N, returns a list of elements of P^1(N)

function P1N(N)
  L := P1Nlist(N);
  O := Order(N);
  OmodN := quo<O|N>;
  P1 := ProjectiveSpace(OmodN, 1);
for x in L do
print x;
print [OmodN!xi: xi in x];
print P1![OmodN!xi: xi in x];
end for;
return [P1![OmodN!xi:xi in x]: x in L];
end function;

function P1(N)
  O := Order(N);
  OmodN := quo<O|N>;
  res := [O!x: x in OmodN];
  P := ProjectiveSpace(OmodN, 1);
return [P![x,y]: x in res, y in res | x*O+y*O+N eq 1*O];
end function;
