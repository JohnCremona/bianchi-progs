from sage.all import QQ, NumberField, polygen, NFCusp

dk = -19
x = polygen(QQ)
k = NumberField(x**2-x+(1-dk)//4, 'w')
assert k.discriminant() == dk
w = k.gen()
Ok = k.ring_of_integers()
emb = next(e for e in k.embeddings(CC) if e(w).imag()>0)
rootd=RR(-dk).sqrt()
Ireps=[c.ideal() for c in k.class_group()]

J = Matrix(2,2,[-k(1), k(0), k(0), k(1)])
Smat = Matrix(2,2,[k(0), k(-1), k(1), k(0)])
inf = NFCusp(k,infinity)
tri0 = [NFCusp(k,0), inf, NFCusp(k,1)]


M_alphas = []
all_alphas = []
n_alphas = 0
add_alpha(0,-1,1,0)   # alpha[0] = 0
add_alpha(w-1,2,2,-w) # alpha[1] = w/2
add_alpha(w,2,2,1-w)  # alpha[2] = (w-1)/2

all_alphas = [NFCusp(k,-M[1,1]/M[1,0]) for M in M_alphas]
alpha_denoms = list(Set(M[1,0] for M in M_alphas))
alphas = dict([(d,[alpha for M,alpha in zip(M_alphas, all_alphas) if M[1,0]==d]) for d in alpha_denoms])
alpha_inv = [all_alphas.index(apply(M,cusp(oo))) for M in M_alphas]
assert all([alpha_inv[alpha_inv[i]]==i for i in range(n_alphas)])

for i in range(n_alphas):
     r1,s = frac(all_alphas[i])
     r2,s2 = frac(all_alphas[alpha_inv[i]])
     assert s==s2
     assert (r1*r2+1)/s in Ok


all_alphas = sum(alphas.values(), [])

# Whitley data

# triangles
EWT = []
# from triangular prism (3.2.2)
EWT.append([cusp(oo), cusp(w/2), cusp((w-1)/2)])
EWT.append([cusp(0), cusp(-2/w), cusp(-2/(w-1))])
# from cuboctahedron (3.2.3)
EWT.append([cusp(oo), cusp(0), cusp(1)])
EWT.append([cusp(oo), cusp(w/2), cusp((w+1)/2)])
EWT.append([cusp(w/2), cusp(2/(1-w)), cusp((w-2)/(w+1))])
EWT.append([cusp((w+1)/2), cusp(3/(2-w)), cusp((w-2)/w)])
EWT.append([cusp((w-2)/(w+1)), cusp(3/(2-w)), cusp((w+1)/3)])
EWT.append([cusp((w+1)/3), cusp(2/(2-w)), cusp((w-1)/(w+1))])
EWT.append([cusp(2/(1-w)), cusp(2/(2-w)), cusp(0)])
EWT.append([cusp((w-2)/w), cusp((w-1)/(w+1)), cusp(1)])

# Yasaki data (his indices start at 1 so the 0'th is a dummy)
# [
# <1, -1>,
# <1, 0>,
# <w + 1, -w + 1>,
# <-w + 2, -2>,
# <w + 2, -w>,
# <-w + 3, -2>,
# <-w + 1, -2>,
# <w + 1, -w + 2>,
# <-w + 2, -3>,
# <2, -w - 1>,
# <2, -w>,
# <3, -w - 1>,
# <0, 1>,
# <w, -w + 2>
# ]

DYV = [None, cusp(-1), cusp(oo), cusp((w + 1)/(-w + 1)),
       cusp(-(-w + 2)/2), cusp(-(w + 2)/w), cusp(-(-w + 3)/2),
       cusp(-(-w + 1)/2), cusp((w + 1)/(-w + 2)), cusp(-(-w +2)/3),
       cusp(2/(-w - 1)), cusp(2/(-w)), cusp(3/(-w - 1)),
       cusp(0), cusp(w/(-w + 2))]

tv = [[1,3,5],[2,4,6],[2,4,7],[7,8,11],[8,9,12],[3,4,12],[10,11,13],[1,2,13],[1,3,14],[9,10,14]]
DYT = [[DYV[i] for i in t] for t in tv]

sv = [ [1,2,6,5], [1,2,4,3], [3,4,6,5], [1,2,4,3], [4,7,8,12],
       [8,9,10,11], [2,7,11,13], [3,14,9,12], [1,14,10,13]]
DYS = [[DYV[i] for i in s] for s in sv]

Tlist, triangles = make_triangles()
print("{} triangles (before eliminating congruences)".format(len(triangles)))
assert len(triangles)==1
T0,t0=reduce_triangles(Tlist, triangles)
print("{} triangles (after eliminating congruences)".format(len(t0)))
assert len(t0)==1
# Check that all DY's triangles are congruent to one of ours:
assert [t for t in DYT if not check_poly_in_list(t,[tri0]+t0)] == []
print("All DY's triangles are in our list")
# Check that all EW's triangles are congruent to one of ours:
assert [t for t in EWT if not check_poly_in_list(t,[tri0]+t0)] == []
print("All EW's triangles are in our list")

# cyclic triangles for alpha with M_alpha of order 3, i.e. alpha=r/s with r^3=+-1 (mod s) so alpha' = (r+1)/s or (r-1)/s.
# Equivalently, trace(M_alpha) = +/-1
cyclic_triangles = [i for i,M in enumerate(M_alphas) if i>2 and M.trace() in [1,-1]]
print("cyclic triangle indices: {}".format(cyclic_triangles))
# []
# Hence there are no more cyclic triangles

print("{} squares (before eliminating congruences)".format(len(DYS)))
assert len(DYS)==9

squares = []
for i,S in enumerate(DYS):
    if not any(poly_equiv(S,s,sign=0) for s in squares):
        print("Square {} is new".format(i))
        squares.append(S)

print("{} squares (after eliminating congruences)".format(len(squares)))
assert len(squares)==1
