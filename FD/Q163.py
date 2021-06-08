from sage.all import QQ, NumberField, polygen, NFCusp

dk = -163
x = polygen(QQ)
k = NumberField(x**2-x+(1-dk)//4, 'w')
w = k.gen()
zero = k(0)
one = k(1)
Ok = k.ring_of_integers()

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

