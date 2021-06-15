from sage.all import QQ, NumberField, polygen, NFCusp

dk = -163
x = polygen(QQ)
k = NumberField(x**2-x+(1-dk)//4, 'w')
w = k.gen()
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
add_alpha(w-1,20,2,-w) # alpha[1] = w/2
add_alpha(w,20,2,1-w)  # alpha[2] = (w-1)/2
add_four_alphas(3, w, 1-w)
add_two_alphas(3, 1+w)
add_four_alphas(4, w, w-1)
add_four_alphas(4, w+1, 2-w)
assert n_alphas==17

add_four_alphas(5, w, w-1)
add_four_alphas(5, 2*w, 2-2*w)
add_four_alphas(5, w+2, 1-2*w)
add_four_alphas(5, w-2, 2+2*w)
add_four_alphas(5, w+1, 1+2*w)
assert n_alphas==37

# 12 alphas with s=6 (norm 36)

add_four_alphas(6, w, 1-w)
add_four_alphas(6, 1+w, w-2)
add_four_alphas(6, 2+w, 3-w)
assert n_alphas==49

# 4 alphas with s=w (norm 41), and 4 conjugates of these

add_four_alphas(w, 12, 17)
add_four_alphas(1-w, 12, 17)
assert n_alphas==57

# 4 alphas with s=1+w (norm 43), and 4 conjugates of these

add_four_alphas(1+w, 12, w-17)
add_four_alphas(2-w, 12, -w-16)
assert n_alphas==65

# 8 alphas with s=2+w (norm 47), and 8 conjugates of these

add_four_alphas(2+w, 7, 18-w)
add_four_alphas(2+w, w-11, w-16)

add_four_alphas(3-w, 7, 17+w)
add_four_alphas(3-w, -w-10,-w-15)
assert n_alphas==81

# 10 alphas with s=7 (norm 49)

add_four_alphas(7, w+3, 2*w-1)
add_four_alphas(7, 2*w+1, 3-2*w)
add_two_alphas(7, 2+3*w)
assert n_alphas==91

# 4 alphas with s=3+w (norm 53), and 4 conjugates of these

add_four_alphas(3+w, 5-w, w-17)
add_four_alphas(4-w, w+4, -w-16)
assert n_alphas==99

alpha_denoms = list(Set(k(M[1,0]) for M in M_alphas))
alphas = dict([(d,[alpha for M,alpha in zip(M_alphas, all_alphas) if M[1,0]==d]) for d in alpha_denoms])

alpha_inv = [all_alphas.index(apply(M,cusp(oo))) for M in M_alphas]
assert all([alpha_inv[alpha_inv[i]]==i for i in range(25)])

def make_edge_matrix(a,b):
    return Matrix(2,2,[a.numerator(),b.numerator(),a.denominator(),b.denominator()])

pairs = [(i,j) for i,a in enumerate(all_alphas) for j,b in enumerate(all_alphas) if j>i and posdet(make_edge_matrix(a,b)) in alpha_denoms]

Tlist, triangles = make_triangles()
print("{} triangles (before eliminating congruences)".format(len(triangles)))
assert len(triangles) == 173
T0,t0=reduce_triangles(Tlist, triangles)
print("{} triangles (after eliminating congruences)".format(len(t0)))
assert len(t0) == 36

# cyclic triangles for alpha with M_alpha of order 3, i.e. alpha=r/s with r^3=+-1 (mod s) so alpha' = (r+1)/s or (r-1)/s.
# Equivalently, trace(M_alpha) = +/-1
cyclic_triangles = [i for i,M in enumerate(M_alphas) if i>2 and M.trace() in [1,-1]]
print("cyclic triangle indices: {}".format(cyclic_triangles))
# [9, 10, 11, 12, 17, 18, 19, 20]
# Hence triangle (i,i',i) is cyclic via M_alphas[i] for i = 9,17.
# (9, 11, 9), (17, 19, 17)

# Yasaki data

# vertices

DYV = [ None, cusp((1)/(-1)), cusp((w+2)/(-w+4)), cusp((w+3)/(-w +3)),
        cusp((-w+5)/(-8)), cusp((8)/(-w-4)), cusp((-2*w+2)/(w - 13)),
        cusp((w+2)/(-w+5)), cusp((-w+4)/(-8)), cusp((14)/(-2*w - 7)),
        cusp((w+3)/(-w+4)), cusp((7)/(-w-3)), cusp((-2*w+9)/(-16)),
        cusp((w+10)/(-2*w+1)), cusp((2*w + 5)/(-2*w+9)),
        cusp((7)/(-w-4)), cusp((-w+12)/(-w-12)),
        cusp((5*w+1)/(-3*w+28)), cusp((3*w+2)/(-2*w+16)),
        cusp((5*w+8)/(-4*w+25)), cusp((-w+33)/(-4*w-20)),
        cusp((2*w+35)/(-6*w-4)), cusp((3*w+9)/(-3*w+13)),
        cusp((-3*w+27)/(-2*w-29)), cusp((22)/(-3*w-10)),
        cusp((6)/(-w-3)), cusp((0)/(1)), cusp((-w+5)/(-9)),
        cusp((-7*w+7)/(4*w-46)), cusp((4*w+35)/(-8*w+4)),
        cusp((-3*w+42)/(-4*w-42)), cusp((w+1)/(-w+6)),
        cusp((7)/(-w-5)), cusp((w+8)/(-2*w+1)),
        cusp((3*w+7)/(-3*w+14)), cusp((21)/(-3*w-11)),
        cusp((-3*w+5)/(w-19)), cusp((-5*w + 6)/(2*w-31)),
        cusp((4*w+27)/(-6*w+11)), cusp((-4*w + 31)/(-2*w-36)),
        cusp((w+32)/(-5*w-8)), cusp((-w+19)/(-2*w-14)),
        cusp((w+1)/(-w+5)), cusp((6)/(-w-4)), cusp((-4*w+2)/(2*w-23)),
        cusp((-w+4)/(-7)), cusp((2*w+14)/(-3*w+5)),
        cusp((-2*w+16)/(-w-18)), cusp((w+18)/(-3*w-2)),
        cusp((2*w+21)/(-4*w+2)), cusp((-2*w+23)/(-2*w-21)), cusp(oo) ]

# triangles
tv = [ [1,2,4], [1,2,3], [1,3,5], [1,4,5], [7,8,9], [6,8,9], [6,7,9],
       [6,7,8], [11,12,13], [10,12,13], [10,11,13], [10,11,12],
       [14,15,16], [4,15,16], [4,14,16], [4,14,15], [11,17,19],
       [18,20,21], [11,18,22], [19,21,23], [11,13,22], [10,11,13],
       [10,13,24], [13,22,24], [8,25,26], [7,25,26], [7,8,26],
       [7,8,25], [7,15,27], [28,29,30], [27,32,33], [31,32,33],
       [27,31,33], [27,31,32], [4,5,10], [1,4,5], [1,5,10], [1,4,10],
       [2,4,15], [1,4,15], [1,2,15], [1,2,4], [8,9,11], [7,9,11],
       [7,8,9], [7,8,11], [10,11,12], [10,12,14], [9,11,12],
       [9,12,14], [7,9,34], [14,15,35], [4,14,15], [10,14,15],
       [4,10,15], [4,10,14], [11,18,36], [11,17,37], [36,38,40],
       [18,20,39], [11,18,22], [10,24,41], [8,11,26], [7,11,26],
       [7,8,26], [7,8,11], [27,31,32], [7,15,27], [1,15,42],
       [7,26,43], [4,10,15], [1,4,15], [1,10,15], [1,4,10],
       [10,14,15], [7,9,11], [10,11,45], [11,18,36], [44,46,47],
       [45,46,48], [10,41,47], [18,41,49], [44,49,50], [36,48,50],
       [1,10,15], [7,11,26], [10,11,45], [1,26,51] ]

DYT = [[DYV[i] for i in t] for t in tv]

# Check that all DY's triangles are congruent to one of ours:
assert [t for t in DYT if not check_poly_in_list(t,[tri0]+t0)] == []
print("All DY's triangles are in our list")

# squares
sv = [ [2,3,5,4], [11,18,21,19], [11,17,20,18], [17,19,21,20],
       [11,19,23,22], [11,18,21,19], [18,21,23,22], [10,11,22,24],
       [15,27,29,28], [7,15,28,30], [7,27,29,30], [9,11,10,14],
       [7,15,35,34], [7,9,14,15], [9,14,35,34], [11,17,20,18],
       [11,36,40,37], [18,36,38,39], [10,11,22,24], [10,11,18,41],
       [18,22,24,41], [7,27,31,43], [15,27,32,42], [1,15,7,26],
       [7,15,14,9], [7,11,10,15], [9,14,10,11], [10,11,18,41],
       [11,36,48,45], [10,45,46,47], [44,46,48,50], [41,47,44,49],
       [18,36,50,49], [7,11,10,15], [1,15,7,26], [1,10,11,26],
       [1,10,11,26], [1,10,45,51], [11,26,51,45]]


DYS = [[DYV[i] for i in s] for s in sv]

print("DY has {} squares".format(len(DYS)))
assert len(DYS) == 39

squares = []
for i,S in enumerate(DYS):
    SS = std_poly(S)[1]
    if not any(poly_equiv(SS,s) for s in squares):
        print("Square {} is new".format(i))
        squares.append(SS)

print("Number of incongruent squares is {}".format(len(squares)))
assert len(squares) == 13
for S in squares:
    ((i,j,kk,l),(x,y,z)) = square_parameters(S)
    print("add_square({},{},{},{}, {},{},{});".format(i,j,kk,l,x,y,z))

# add_square(42,36,8,30, 0,0,0);
# add_square(17,43,17,43, 0,0,0);
# add_square(38,0,37,19, 0,0,0);
# add_square(41,35,7,29, 0,0,0);
# add_square(11,31,89,33, 0,0,0);
# add_square(1,90,1,90, -w,w,w);
# add_square(88,9,87,8, 0,0,0);
# add_square(1,89,1,89, 1,-1);
# add_square(62,8,59,2, -1,0,0);
# add_square(11,17,11,17, 0,0,0);
# add_square(2,66,0,75, 0,0,0);
# add_square(50,9,55,2, 0,0,0);
# add_square(83,11,81,0, 0,0,0);
