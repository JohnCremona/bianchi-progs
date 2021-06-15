from sage.all import QQ, NumberField, polygen

dk = -67
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

n_alphas = 0
M_alphas = []
all_alphas = []

add_alpha(0,-1,1,0)
add_alpha(w-1,8,2,-w) # alpha[1] = w/2
add_alpha(w,8,2,1-w)  # alpha[2] = (w-1)/2
add_four_alphas(3, w, 1-w)
add_two_alphas(3, 1+w)
add_four_alphas(4, w, w-1);
add_four_alphas(4, w+1, 2-w);
add_four_alphas(3-w, w+6, 2+w);
add_four_alphas(2+w, 7-w, 3-w);

alpha_denoms = list(Set(M[1,0] for M in M_alphas))
alphas = dict([(d,[alpha for M,alpha in zip(M_alphas, all_alphas) if M[1,0]==d]) for d in alpha_denoms])

alpha_inv = [all_alphas.index(apply(M,cusp(oo))) for M in M_alphas]
assert all([alpha_inv[alpha_inv[i]]==i for i in range(25)])
for i in range(n_alphas):
     r1,s = frac(all_alphas[i])
     r2,s2 = frac(all_alphas[alpha_inv[i]])
     assert s==s2
     assert (r1*r2+1)/s in Ok

def make_edge_matrix(a,b):
    return Matrix(2,2,[a.numerator(),b.numerator(),a.denominator(),b.denominator()])

pairs = [(i,j) for i,a in enumerate(all_alphas) for j,b in enumerate(all_alphas) if j>i and posdet(make_edge_matrix(a,b)) in alpha_denoms]

# EW's triangles:
EWT = []
# from P_v_3 (fig 3.4.2): all equi to tri2
EWT.append([cusp(oo), cusp(w/2),cusp((w+1)/2)])
EWT.append([cusp((8*w+3)/19), cusp((12*w+7)/29),cusp((2*w+1)/5)])
EWT.append([cusp((8*w+8)/19), cusp((12*w+10)/29),cusp((2*w+2)/5)])
EWT.append([cusp((7*w+5)/19), cusp((7*w+7)/19),cusp((w+1)/3)])
# from P_v_2 (fig 3.4.3):
EWT.append([cusp(oo), cusp(w/3),cusp((w+1)/3)])
EWT.append([cusp(w/2), cusp(7*w/17),cusp((8*w+3)/19)])
EWT.append([cusp((2*w+1)/5), cusp((7*w+5)/19),cusp((11*w+4)/29)])
EWT.append([cusp(w/3), cusp(7*w/17),cusp((11*w+4)/29)])
# from P_v_5 (fig 3.4.4):
EWT.append([cusp(oo), cusp(w/3),cusp((w+1)/4)])
EWT.append([cusp(oo), cusp(w/4),cusp((w+1)/4)])
EWT.append([cusp(5*w/17), cusp(w/3),cusp((w+1)/4)])
EWT.append([cusp(5*w/17), cusp(w/4),cusp((w+1)/4)])
# from P_v_6 (fig 3.4.5):
EWT.append([cusp(5*(w+1)/19), cusp(w/3),cusp((w+1)/3)])
EWT.append([cusp(5*(w+1)/19), cusp(w/3),cusp((w+1)/4)])
# from P_v_9 (fig 3.4.6):
EWT.append([cusp(0), cusp(4*w/17),cusp(4*(w+1)/19)])
# from P_v_8 (fig 3.4.7):
EWT.append([cusp(infinity), cusp((w-1)/4),cusp(w/4)])
EWT.append([cusp(0), cusp(4*(w-1)/17),cusp(4*w/17)])
# from P_v_11 (fig 3.4.8):
EWT.append([cusp(infinity), cusp((w-1)/2), cusp(w/2)])
EWT.append([cusp(infinity), cusp(w/2), cusp((w-7)/(w+2))])
EWT.append([cusp(infinity), cusp((w-7)/(w+2)), cusp((w+6)/(3-w))])
EWT.append([cusp(infinity), cusp((w+6)/(3-w)), cusp((w-1)/2)])
EWT.append([cusp((2*w-1)/5), cusp((w-1)/2), cusp(w/2)])
EWT.append([cusp((2*w-1)/5), cusp(w/2), cusp((w-7)/(w+2))])
EWT.append([cusp((2*w-1)/5), cusp((w-7)/(w+2)), cusp((w+6)/(3-w))])
EWT.append([cusp((2*w-1)/5), cusp((w+6)/(3-w)), cusp((w-1)/2)])

Tlist, triangles = make_triangles()
print("{} triangles (before eliminating congruences)".format(len(triangles)))
assert len(triangles)==31
T0,t0=reduce_triangles(Tlist, triangles)
print("{} triangles (after eliminating congruences)".format(len(t0)))
assert len(t0)==8
assert [t for t in EWT if not check_poly_in_list(t,[tri0]+t0)] == []

# cyclic triangles for alpha with M_alpha of order 3, i.e. alpha=r/s with r^3=+-1 (mod s) so alpha' = (r+1)/s or (r-1)/s.
# Equivalently, trace(M_alpha) = +/-1
cyclic_triangles = [i for i,M in enumerate(M_alphas) if i>2 and M.trace() in [1,-1]]
print("cyclic triangle indices: {}".format(cyclic_triangles))
# [9,10,11,12]
# Hence triangle (9,11,9) is cyclic via M_alphas[9]

#EW's squares:
EWS = []
# from P_v_2 (fig 3.4.3):
EWS.append([cusp(oo), cusp(w/2), cusp(7*w/17), cusp(w/3)])
EWS.append([cusp((2*w+1)/5), cusp((11*w+4)/29), cusp(7*w/17), cusp((8*w+3)/19)])
EWS.append([cusp(w/3), cusp((11*w+4)/29), cusp((7*w+5)/19), cusp((w+1)/3)])
# from P_v_9 (fig 3.4.6):
EWS.append([cusp(infinity), cusp(0), cusp(4*w/17), cusp(w/4)])
EWS.append([cusp(infinity), cusp(0), cusp(4*(w+1)/19), cusp((w+1)/4)])
EWS.append([cusp(w/4), cusp((w+1)/4), cusp(4*(w+1)/19), cusp(4*w/17)])
# from P_v_8 (fig 3.4.7):
EWS.append([cusp(infinity), cusp(0), cusp(4*(w-1)/17), cusp((w-1)/4)])
EWS.append([cusp((w-1)/4), cusp(w/4), cusp(4*w/17), cusp(4*(w-1)/17)])

squares = []
for i,S in enumerate(EWS):
    if not any(poly_equiv(S,s) for s in squares):
        print("Square {} is new".format(i))
        squares.append(S)

# After the above we have 3 squares:

S0 = squares[0] # [oo,w/2,7/(1-w),w/3]
S1 = squares[1] # [oo,0,4/(1-w),w/4]
S2 = squares[2] # [oo,0,4/(2-w),(w+1)/4]

assert [all_alphas[3], cusp(oo)] == squares[0][3:] + squares[0][:1]
assert apply(M_alphas[2], [all_alphas[2], cusp(oo)]) == squares[0][0:2]
U1 = Matrix(2,2,[7,2*w,1-w,5])
assert U1.det()==1
assert apply(U1,[all_alphas[4], cusp(oo)]) == squares[0][1:3]
assert apply(M_alphas[5],[all_alphas[12], cusp(oo)]) == squares[0][2:4]

# Yasaki data

# vertices
DYV = [ None,
        cusp((-w + 2)/(-4)),
        cusp((-w + 3)/(-4)),
        cusp((w + 3)/(-w + 2)),
        cusp((-2*w + 5)/(-8)),
        cusp((5)/(-w - 1)),
        cusp((w + 8)/(-2*w + 1)),
        cusp(oo),
        cusp((-w + 2)/(-3)),
        cusp((w + 2)/(-w + 3)),
        cusp((5)/(-w - 2)),
        cusp((2*w + 5)/(-2*w + 5)),
        cusp((-w + 8)/(-w - 6)),
        cusp((10)/(-2*w - 3)),
        cusp((2*w)/(-w + 7)),
        cusp((-4*w + 2)/(2*w - 15)),
        cusp((-w + 3)/(-5)),
        cusp((-5*w + 5)/(2*w - 20)),
        cusp((2*w + 13)/(-4*w + 2)),
        cusp((-2*w + 15)/(-2*w - 13)),
        cusp((3*w + 15)/(-5*w + 5)),
        cusp((-3*w + 18)/(-2*w - 18)),
        cusp((2*w + 18)/(-5*w)),
        cusp((-2*w + 20)/(-3*w - 15)),
        cusp((5*w)/(-3*w + 18)),
        cusp((1)/(-1)),
        cusp((w + 2)/(-w + 2)),
        cusp((0)/(1))
]

# triangles
tv = [ [1,2,4], [3,5,6], [1,2,7], [3,5,8], [2,10,12], [4,11,13],
       [1,9,14], [9,10,16], [17,18,21], [15,22,23], [19,20,24],
       [10,25,26], [7,10,25], [2,10,26], [2,7,10], [1,2,7], [1,7,9],
       [7,9,10], [2,7,10], [7,9,10], [7,10,25], [9,10,16],
       [10,16,25], [7,9,27], [16,25,27], [9,16,27], [7,25,27] ]

DYT = [[DYV[i] for i in t] for t in tv]

# Check that all DY's triangles are congruent to one of ours:
assert [t for t in DYT if not check_poly_in_list(t,[tri0]+t0)] == []
print("All DY's triangles are in our list")

# squares
sv = [ [2,3,6,4], [1,2,3,5], [1,4,6,5], [2,3,8,7], [1,5,8,7],
       [1,2,3,5], [1,2,10,9], [2,4,11,12], [1,4,13,14], [2,7,25,26],
       [1,2,10,9] ]

DYS = [[DYV[i] for i in s] for s in sv]

# Check that all DY's squares are congruent to one of ours:
assert [s for s in DYS if not check_poly_in_list(s,squares)] == []
print("All DY's squares are in our list")

