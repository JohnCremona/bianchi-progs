Qx<x> := PolynomialRing(Rationals());
K<w> := NumberField(x^2+5);
assert Discriminant(K) eq -20;
OK := Integers(K);
w:= OK.2;
