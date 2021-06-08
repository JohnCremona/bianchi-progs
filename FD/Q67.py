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

def make_edge_matrix(a,b):
    return Matrix(2,2,[a.numerator(),b.numerator(),a.denominator(),b.denominator()])

pairs = [(i,j) for i,a in enumerate(all_alphas) for j,b in enumerate(all_alphas) if j>i and posdet(make_edge_matrix(a,b)) in alpha_denoms]

tri0 = [cusp(0),cusp(infinity),cusp(1)]

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
