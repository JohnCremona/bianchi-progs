from sage.all import QQ, NumberField, polygen, NFCusp, infinity, CC, RR, Matrix

dk = -43
x = polygen(QQ)
k = NumberField(x**2-x+(1-dk)//4, 'w')
assert k.discriminant() == dk
w = k.gen()
wbar = 1-w
zero = k(0)
one = k(1)
Ok = k.ring_of_integers()
emb = next(e for e in k.embeddings(CC) if e(w).imag()>0)
rootd=RR(-dk).sqrt()

J = Matrix(2,2,[-k(1), k(0), k(0), k(1)])
Smat = Matrix(2,2,[k(0), k(-1), k(1), k(0)])
inf = NFCusp(k,infinity)
tri0 = [NFCusp(k,0), inf, NFCusp(k,1)]


M_alphas = []
all_alphas = []
n_alphas = 0

add_alpha(0,-1,1,0)
add_alpha(w-1,5,2,-w) # alpha[1] = w/2
add_alpha(w,5,2,1-w)  # alpha[2] = (w-1)/2
add_four_alphas(3, w, 1-w)
add_two_alphas(3, 1+w)

all_alphas = [NFCusp(k,-M[1,1]/M[1,0]) for M in M_alphas]
n_alphas = len(all_alphas)
alpha_denoms = list(Set(M[1,0] for M in M_alphas))
alphas = dict([(d,[alpha for M,alpha in zip(M_alphas, all_alphas) if M[1,0]==d]) for d in alpha_denoms])
alpha_inv = [all_alphas.index(apply(M,cusp(oo))) for M in M_alphas]
assert all([alpha_inv[alpha_inv[i]]==i for i in range(n_alphas)])

for i in range(n_alphas):
     r1,s = frac(all_alphas[i])
     r2,s2 = frac(all_alphas[alpha_inv[i]])
     assert s==s2
     assert (r1*r2+1)/s in Ok

# triangles from EW:

EWT = []
# from P_v_5 (3.3.2):
EWT.append([cusp(oo), cusp(w/3), cusp((w+1)/3)])
EWT.append([cusp(0), cusp(3*w/11), cusp(3*(w+1)/13)])
EWT.append([cusp(3*w/11), cusp(w/3), cusp((5*w+2)/17)])
EWT.append([cusp((w+1)/4), cusp( (5*w+2)/17), cusp( 4*(w+1)/13)])
# from P_v_4 (3.3.3):
EWT.append([cusp(0), cusp(3*w/11), cusp(3*(w-1)/11)])
EWT.append([cusp(infinity), cusp( (w-1)/3), cusp( w/3)])
# from P_v_3 (3.3.4):
EWT.append([cusp((w+1)/3), cusp( (5*w+3)/13), cusp( 5*(w+1)/13)])
EWT.append([cusp(infinity), cusp( w/2), cusp( (w+1)/2)])
# from P_v_6 (3.3.5):
EWT.append([cusp(0), cusp(1), cusp(infinity)])
EWT.append([cusp((w+1)/4), cusp( 3*(w+1)/13), cusp( (4*w+5)/17)])
EWT.append([cusp((w+2)/4), cusp( (3*w+7)/13), cusp( (4*w+8)/17)])
EWT.append([cusp((w+1)/3), cusp( 4*(w+1)/13), cusp( (4*w+5)/13)])

Tlist, triangles = make_triangles()
print("{} triangles (before eliminating congruences)".format(len(triangles)))
assert len(triangles)==5
T0,t0=reduce_triangles(Tlist, triangles)
print("{} triangles (after eliminating congruences)".format(len(t0)))
assert len(t0)==2
assert [t for t in EWT if not check_poly_in_list(t,[tri0]+t0)] == []

# cyclic triangles for alpha with M_alpha of order 3, i.e. alpha=r/s with r^3=+-1 (mod s) so alpha' = (r+1)/s or (r-1)/s:
# i i'  r s r'
# 0 0   0 1 0    special case
# 1 2   w 2 w-1  special case
# 9 11  w 4 w-1  use
# 10 12 -w 4 1-w congruent to previous


# cyclic triangles for alpha with M_alpha of order 3, i.e. alpha=r/s with r^3=+-1 (mod s) so alpha' = (r+1)/s or (r-1)/s.
# Equivalently, trace(M_alpha) = +/-1
cyclic_triangles = [i for i,M in enumerate(M_alphas) if i>2 and M.trace() in [1,-1]]
print("cyclic triangle indices: {}".format(cyclic_triangles))
# []
# Hence there are no more cyclic triangles

# squares from EW:
EWS = []
# from P_v_5 (3.3.2):
EWS.append([cusp(oo), cusp(0), cusp(3*w/11), cusp(w/3)])
EWS.append([cusp(w/3), cusp((w+1)/3), cusp(4*(w+1)/13), cusp((5*w+2)/17)])
EWS.append([cusp(3*w/11), cusp((5*w+2)/17), cusp((w+1)/4), cusp(3*(w+1)/13)])
# from P_v_4 (3.3.3):
#EWS.append([cusp(oo), cusp(0), cusp(3*w/11), cusp(w/3)])
EWS.append([cusp(oo), cusp(0), cusp(3*(w-1)/11), cusp((w-1)/3)])
EWS.append([cusp(w/3), cusp(3*w/11), cusp(3*(w-1)/11), cusp((w-1)/3)])
# from P_v_3 (3.3.4):
EWS.append([cusp(oo), cusp((w+1)/3), cusp((5*w+3)/13), cusp(w/2)])
EWS.append([cusp(oo), cusp((w+1)/3), cusp((5*w+5)/13), cusp((w+1)/2)])
EWS.append([cusp(w/2), cusp((w+1)/2), cusp((5*w+5)/13), cusp((5*w+3)/13)])

print("{} squares (before eliminating congruences)".format(len(EWS)))
assert len(EWS)==8

# Of these 8 squares, the SL2-equivalence classes are {0,1,2,3}, {4}, {5,6,7}
assert all(poly_equiv(EWS[i], EWS[0], sign=1) for i in [1,2,3])
assert all(poly_equiv(EWS[i], EWS[5], sign=1) for i in [6,7])
assert not poly_equiv(EWS[0], EWS[4], sign=1)
assert not poly_equiv(EWS[5], EWS[4], sign=1)
assert not poly_equiv(EWS[0], EWS[5], sign=1)

# Under GL2, 4 is equivalent to {5,6,7} but not {0,1,2,3}
assert not poly_equiv(EWS[0], EWS[5], sign=0)
assert not poly_equiv(EWS[0], EWS[4], sign=0)
assert poly_equiv(EWS[4], EWS[5], sign=0)

squares = []
for i,S in enumerate(EWS):
    if not any(poly_equiv(S,s,sign=0) for s in squares):
        print("Square {} is new".format(i))
        squares.append(S)

print("{} squares (after eliminating congruences)".format(len(squares)))
assert len(squares)==2

# After the above we only have 2 squares (there would be 3 with sign=1).

# First is [oo,0,3/(1-w),w/3], map by S to get
S0 = apply(Smat, squares[0])
#  = [0,oo,(w-1)/3, -3/w] = [alpha_0, oo, alpha_6, S(alpha_3)]

E0=[std_edge(e) for e in poly_edges(S0)]
S0types = [e[0] for e in E0] # = [0, 4, 2, 3]
S0mats = [e[1] for e in E0]
assert [apply(M,[all_alphas[i],cusp(oo)]) for M,i in zip(S0mats, S0types)] == poly_edges(S0)
assert S0mats[0] == Tmat(zero) # identity
assert S0mats[1] == J*M_alphas[3]*J
assert S0mats[2] == Matrix(2,2,[-3, 2*w - 2, w, 7])
assert S0mats[3] == Smat

assert [all_alphas[0], cusp(oo)] == S0[:2]
assert apply(M_alphas[4], [all_alphas[4], cusp(oo)]) == S0[1:3]
U1 = Matrix(2,2,[-3, 2*w - 2, w, 7])
assert U1.det()==1
assert apply(U1,[all_alphas[2], cusp(oo)]) == S0[2:4]
assert apply(M_alphas[0],[all_alphas[3], cusp(oo)]) == S0[3:]+S0[:1]

# NB The image of this under J is not +1-equivalent to either square.

# Second is [w/3, 3/(1-w), -3/w, (w-1)/3] = [alpha[3], S(alpha[6]), S(alpha[3]), alpha[6]]
# (so S gives a self-congruence, cycling 2 steps).
S1 = squares[1]

assert apply(Smat,S1) == cycle_poly(S1,2)
# Apply a transform to get
U, S1 = std_poly(squares[1])
# S1 = [w/2, oo, (2*w-1)/3, (w-7)/(w+1)]

E1=[std_edge(e) for e in poly_edges(S1)]
S1types = [e[0] for e in E1] # = [1,8,1,8]
S1mats = [e[1] for e in E1]
assert S1mats[0] == Tmat(zero) # identity
assert S1mats[1] == -Tmat(w) * M_alphas[8]
assert S1mats[2] == -U.inverse()*Smat*U
assert S1mats[3] == M_alphas[2]*Tmat(w)
assert [apply(M,[all_alphas[i],cusp(oo)]) for M,i in zip(S1mats, S1types)] == poly_edges(S1)
