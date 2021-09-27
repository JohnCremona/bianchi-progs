Qx<x> := PolynomialRing(Rationals());
K<w> := NumberField(x^2-x+6);
assert Discriminant(K) eq -23;
OK := Integers(K);
w:= OK.2;
