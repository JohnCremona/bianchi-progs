Qx<x> := PolynomialRing(Rationals());
K<w> := NumberField(x^2-x+12);
assert Discriminant(K) eq -47;
OK := Integers(K);
w:= OK.2;
