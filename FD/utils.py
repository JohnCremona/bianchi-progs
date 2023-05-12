from sage.all import (oo, NFCusp, Matrix, cached_function)

def nf(x):
    """
    Returns the number field of a cusp or field element
    """
    try:
        return x.number_field()
    except AttributeError:
        return x.parent()

Imat = Matrix(2,2,[1,0,0,1])
Jmat = Matrix(2,2,[-1,0,0,1])
Smat = Matrix(2,2,[0,-1,1,0])

def Tmat(x):
    return Matrix(2,2, [1,x,0,1])

def to_k(c, k=None):
    if k is None:
        k = nf(c)
    return k(c.numerator()/c.denominator()) if c.denominator() else oo

def round(z): # z in k
    from sage.all import RR, CC
    k = nf(z)
    emb = next(e for e in k.embeddings(CC) if e(k.gen()).imag()>0)
    x, y = emb(z)
    rootd=RR(k.discriminant().abs()).sqrt()
    b = (2*y/rootd).round()
    a = ((2*x).round() - b)//2
    return k([a,b])

def ispos(a):
    return a[1]>0 or (a[1]==0 and a[0]>0)

def makepos(a):
    return a if ispos(a) else -a

def posgen(I):
    a = I.number_field()(I.gens_reduced()[0])
    return makepos(a)

def posdet(M):
    return makepos(M.det())

def kgcd(a,b, k=None):
    if k is None:
        k = nf(a)
    I = k.ideal(a,b)
    return posgen(I) if I.is_principal() else k(1)

def frac(x):
    """
    Return coprime [a,b] (if they exist) such that x=a/b, otherwise any such [a,b].
    x can be an element of k or a cusp.
    """
    if x == oo:
        return [1, 0]
    a = x.numerator()
    b = x.denominator()
    c = kgcd(a,b)
    k = nf(a)
    return [k(a/c),k(b/c)]

#We need to label cusps with something hashable, to use graph
#functionality, which can be converted back into cusps later.
#Unfortunately the NFCusp class is not hashable.

def cusp_label(a):
    n,d = a.numerator(), a.denominator()
    return (-n,-d) if list(d)<[0,0] else (n,d)

NFCusp._repr_ = lambda c: "[{}:{}]".format(*cusp_label(c)) if c.denominator() else "oo"

@cached_function
def smallest_ideal_class_representatives(k):
    """
    Return a list of ideals representing the ideal classes, each minimal in its class.
    """
    C = k.class_group()
    IC = [c.ideal() for c in C]
    maxnorm = max(I.norm() for I in IC)
    Ilist = k.ideals_of_bdd_norm(maxnorm)
    # make sure that each ideal is the one of smallest norm in its class
    for (c,(i,I)) in zip(C, enumerate(IC)):
        Inorm = I.norm()
        for n,Inlist in Ilist.items():
            for In in Inlist:
                if n<Inorm and C(In)==c:
                    print("Replacing {} of norm {} with equivalent {} of smaller norm {}".format(I,Inorm,In,n))
                    I = In
                    Inorm = n
                    IC[i] = I
    return IC

def cusp(x, k=None, Ireps=None):
    """
    Return x as a reduced cusp [a:b], or [1:0] if x is oo.
    """
    if x==oo:
        if k is None:
            raise RuntimeError("cusp(oo) must specify the number field")
        return NFCusp(k, 1, 0)

    if k is None:
        k = nf(x)

    if Ireps is None:
        Ireps = smallest_ideal_class_representatives(k)

    a, b = frac(k(x))
    if not ispos(b):
        a=-a;
        b=-b
    return NFCusp(k, a, b, lreps=Ireps)

def cusp2(a,b, k=None, Ireps=None):
    """
    Return a/b as a reduced cusp [a:b], or [1:0] if b=0.
    """
    return cusp(a/b if b else oo, k, Ireps)

def tri0(k):  # universal
    return [cusp(a, k) for a in [0, oo, 1]]
def tri1(k):  # universal for non-Euclidean class#1
    w = k.gen()
    return [cusp(a, k) for a in [oo, w/2, (w-1)/2]]
def tri2(k):  # universal for non-Euclidean class#1
    w = k.gen()
    return [cusp(a, k) for a in [oo, w/2, (w+1)/2]]

def translate_cusp(c, x):
    """
    translate the cusp c by x
    """
    return apply(Tmat(x),c)

def negate_cusp(c):
    """
    negate the cusp c to -c
    """
    return apply(Jmat,c)

def polygon(xlist):
    """
    Convert a list of k union {oo} to a list of cusps.
    """
    return [cusp(x) for x in xlist]

def cusp_from_string(s, k):
    if s == 'oo':
        return cusp(oo, k)
    a,b = [k(c) for c in s[1:-1].split(":")]
    return cusp2(a,b, k)

def make_poly_from_edges(F, k):
    """Use for F a face as returned by the graph function G.faces(), so is
    a list of ordered pairs such that the second element of each is
    the first of the next, cyclically, and each element of the pair is
    a pair (a,b) for the cusp [a:b].
    """
    return [cusp_from_string(e[0], k) for e in F]

def edge(x,y):
    """
    Convert a pair of k union {oo} elements to a list of 2 cusps.
    """
    return polygon([x,y])

def is_integral(M):
    """
    Return True iff all entries of M are integral.
    """
    return all(c.is_integral() for c in M.list())

def edge_to_mat(e):
    """Return M such that M(0)=e[0], M(oo)=e[1], with det(M) "positive".
    """
    a,c = frac(e[1])
    b,d = frac(e[0])
    M = Matrix(2,2,[a,b,c,d])
    return M if ispos(M.det()) else M * Jmat

def edet(e, pos=False):
    """
    Return det(M) where M(0)=e[0], M(oo)=e[1], with det(M) "positive" if pos.
    """
    a,c = frac(e[1])
    b,d = frac(e[0])
    det = a*d-b*c
    return makepos(det) if pos else det

def apply(U, P):
    try:
        return P.apply(U.list())
    except (AttributeError, TypeError):
        return [c.apply(U.list()) for c in P]

def alpha_index(a, alphas):
    """
    Return i if a==alphas[i], else -1.
    """
    try:
        return alphas.index(a)
    except ValueError:
        return -1

def alpha_index_with_translation(a, alphas):
    """
    Return (i,x) if a==alphas[i]+x with x integral, else (-1,0).
    """
    ka = to_k(a)
    for i, alpha in enumerate(alphas):
        x = ka - to_k(alpha)
        if x.is_integral():
            return (i, x)
    return (-1, 0)

def sigma_index(s, sigmas):
    """
    Return i if s==sigmas[i], else -1.
    """
    try:
        return sigmas.index(s)
    except ValueError:
        return -1

def sigma_index_with_translation(s, sigmas):
    """
    Return (i,x) if a==sigmas[i]+x with x integral, else (-1,0).
    """
    if s==sigmas[0]:
        return (0, 0)
    ks = to_k(s)
    for i, sigma in enumerate(sigmas[1:]):
        x = ks - to_k(sigma)
        if x.is_integral():
            return (i+1, x)
    return (-1, 0)

def make_M_alpha(a, alphas):
    inf = cusp(oo, nf(a))
    M = Matrix(2,2,a.ABmatrix()).inverse()
    j,x = alpha_index_with_translation(apply(M,inf), alphas)
    assert j>=0
    M = Tmat(-x)*M
    return M, j

def make_M_alphas(alphas):
    t = [make_M_alpha(a, alphas) for a in alphas]
    return [t[0] for t in t], [t[1] for t in t]

def std_edge(e, alphas):
    """Return (i, U) where e=U(e0), U in SL(2,Ok) and
    e0={alpha[i],oo}, else return (-1, None). Note that U,u0 exist
    and are unique for edges of the tessellation.
    """
    #print("standardising edge {}".format(e))
    k = nf(e[1])
    inf = cusp(oo,k)
    U = Matrix(2,2, e[1].ABmatrix())
    e1 = apply(U.inverse(), e)
    #print(" - e1 = {}".format(e1))
    assert e1[1]==inf
    i, x = alpha_index_with_translation(e1[0], alphas)
    #print(" - (i,x) = {}".format((i,x)))
    if i<0:
        return (-1, None)
    assert e1[0]==translate_cusp(alphas[i], x)
    U = U*Tmat(x)
    assert e == apply(U,[alphas[i],inf])
    return (i, U)

def edge_equiv(e1,e2):
    """Return U in SL(2,O) such that U maps e1 to e2, otherwise False.
    Note that such a U exists if and only if there exists V with
    det(V)=-1 mapping e1 to e2.
    """
    M1 = edge_to_mat(e1) # M1{0,oo}=e1
    M2 = edge_to_mat(e2) # M2{0,oo}=e2
    U = M2*M1.inverse()  # U(e1)=e2
    return U if U.det()==1 and is_integral(U) else False

def std_cusp(a):
    return cusp(oo) if a.is_infinity() else cusp(to_k(a))

def std_poly(P, alphas):
    """Return (U,P') such that U(P')=P, with the first edge of P' a basic
    edge and det(U)=1.  If the first edge of P is not the transform of
    a basic edge, return (False, None)
    """
    i, U = std_edge(P[:2], alphas)
    if i<0:
        return (False, None)
    return U, apply(U.inverse(),P)

def std_pol(poly, alphas):
    return [std_cusp(a) for a in std_poly(poly, alphas)[1]]

# The next two are only for polygons with all vertices principal

def poly_gl2_orbit_reps(polys, alphas):
    reps = []
    for poly in polys:
        poly = std_poly(poly, alphas)[1]
        if not check_poly_in_list(poly, reps):
            reps.append(poly)
    return reps

def poly_sl2_orbit_reps(polys, alphas):
    reps = []
    for poly in polys:
        poly = std_poly(poly, alphas)[1]
        if not check_poly_in_list(poly, reps, 1):
            reps.append(poly)
    return reps

def poly_edges(P):
    """
    Return the list of edges of P, including the wrap-around edge from last vertex to first.
    """
    npts = len(P)
    ee = [[P[i],P[i+1]] for i in range(npts-1)] + [[P[npts-1],P[0]]]
    return ee

def cycle_poly(P, n=1):
    """Return a new polygon with the vertices of P cycled round n steps
    (first moved to the end).  Negative n works fine!
    """
    return P[n:] + P[:n]

def reverse_poly(P):
    """
    Return a new polygon with the vertices of P in reverse order.
    """
    return P[::-1]

def conj_cusp(c):
    """
    Return the conjugate cusp
    """
    return NFCusp(nf(c), [x.conjugate() for x in frac(c)])

def conj_poly(P):
    """
    Return a new polygon with the vertices the conjugates of those of P (in the same order)
    """
    return [conj_cusp(c) for c in P]

def poly_equiv_exact(P1,P2, sign=1):
    """Test whether the polygons P1, P2 are GL(2,O)-equivalent (or
    SL(2,O)-equivalent if sign) in the exact sense that there exists a
    matrix in G/SL(2,O) mapping the ordered list of vertices of P1 to
    those of P2.

    Returns either a transforming matrix U, or False
    """
    npts = len(P1)
    if len(P2) != npts:
        return False

    e1 = P1[:2] # first edge of P1
    e2 = P2[:2] # first edge of P2
    U = edge_equiv(e1,e2) # det(U)=1
    if U and apply(U,P1)==P2:
        return U
    if sign != +1:
        P3 = apply(Jmat,P2)
        e3 = P3[:2]
        U = edge_equiv(e1,e3) # det(U)=1
        if U and apply(U,P1)==P3: # then JU(P1)=Jmat(P3)=P2
            return Jmat*U
    return False

def poly_equiv(P1,P2, cycle=True, reverse=True, sign=1):
    """Test whether the polygons P1, P2 are GL(2,O)-equivalent (or
    SL(2,O)-equivalent if sign).  If cycle is True (default), allow P2
    to be cycled, i.e. test if P1 is equivalent to a cycled version of
    P2.  If reverse is True (default), allow P2 to be reversed.

    So for triangles, with the defaults, up to 6 exact transforms are
    tested.

    Return either U, a matrix such that U(P1) = P2 (or a cycled or
    reversed copy of P2), or False.

    If sign==+1 unly U with det(U)=1 will be returned, otherwise
    det(U)=-1 is also possible.

    """
    npts = len(P1)
    if len(P2) != npts:
        return False

    if reverse:
        U = poly_equiv(P1, P2, cycle, False, sign)
        return U if U else poly_equiv(P1, reverse_poly(P2), cycle, False, sign)

    # Now reverse is False

    for n in range(npts if cycle else 1):
        U = poly_equiv_exact(P1, cycle_poly(P2, n), sign)
        if U:
            return U
    return False

def is_cyclic(T):
    return poly_equiv(T, cycle_poly(T), False)

def add_alpha(a,b,c,d, alphas, M_alphas):
    k = nf(a)
    assert a*d-b*c==1
    M_alphas.append(Matrix(2,2,[a,b,c,d]))
    alphas.append(cusp(-d/c, k))

def add_two_alphas(s, r, eps, alphas, M_alphas): # r^2=eps (mod s)
    t = (r*r*eps-1)/s
    add_alpha(-eps*r, t, s, -r, alphas, M_alphas) # alpha =  r/s
    add_alpha( eps*r, t, s,  r, alphas, M_alphas) # alpha = -r/s

def add_four_alphas(s, r1, r2, alphas, M_alphas): # r1*r2=-1 (mod s)
    t = -(r1*r2+1)/s
    add_alpha( r2, t, s, -r1, alphas, M_alphas) # alpha =  r1/s
    add_alpha(-r2, t, s,  r1, alphas, M_alphas) # alpha = -r1/s
    add_alpha( r1, t, s, -r2, alphas, M_alphas) # alpha =  r2/s
    add_alpha(-r1, t, s,  r2, alphas, M_alphas) # alpha = -r2/s

def make_triangles(alphas, M_alphas):
    """Return two lists, Tlist and triangles where each entry in Tlist is
    a triple (i,j,n) of alpha-indices such that M_alphas[i](alphas[j])
    = alphas[n] + x with x in Ok, and the corresponding entry in
    triangles is (alphas[i], oo, alphas[j]).  This is the triangle
    whose edges are

    {alphas[i], oo}
    {oo, alphas[j]} = image of {alphas[j'], oo} under M_alphas[j']
    {alphas[j], alphas[i]} = image of {alphas[n], oo} under M_alphas[i']*T^x

    """
    n_alphas = len(alphas)
    Tlist = []
    for i, Ma in enumerate(M_alphas):
        for j in range(i+1, n_alphas):
            n = alpha_index_with_translation(apply(Ma, alphas[j]))[0]
            if n>=0:
                print("({},{},{}) : <{},{},{}>".format(i,j,n,to_k(alphas[i]), oo, to_k(alphas[j])))
                Tlist.append((i,j,n))
    k = nf(alphas[0])
    return Tlist, [[alphas[i],NFCusp(k,oo),alphas[j]] for i,j,_ in Tlist]

def check_poly_in_list(T, polys, sign=0):
    """Return True if T is GL(2,Ok)-congruent to a polygon in polys or
    the reverse of one.  For SL(2,Ok)-equivalence, set sign=1.
    """
    #print("Testing whether {} is in list {}".format(T, polys))
    return any(poly_equiv(T,t, sign=sign) for t in polys)

def reduce_triangles(Tlist, triangles):
    Tlist0 = []
    triangles0 = []
    for trip, tri in zip(Tlist, triangles):
        if not check_poly_in_list(tri, triangles0):
            Tlist0.append(trip)
            triangles0.append(tri)
    return Tlist0, triangles0

def square_parameters(S, alphas, M_alphas, alpha_inv):
    """
    For S a square, returns [[i,j,k,l],[x,y,z]] where
    """
    k = nf(alphas[0])
    d = -k.disc().squarefree_part()
    inf = cusp(oo, k)

    S = std_poly(S, alphas)[1]
    assert S[1] == inf
    i = alpha_index(S[0], alphas)
    assert i!=-1
    jd, z = alpha_index_with_translation(S[2], alphas)
    assert jd!=-1
    j = alpha_inv[jd]
    beta = S[3]
    Mi=M_alphas[i];
    l, y = alpha_index_with_translation(apply(Mi, beta), alphas)
    assert l!=-1
    Tz = Matrix(2,2,[1,z,0,1])
    kd, x = alpha_index_with_translation(apply(M_alphas[jd]*Tz.inverse(), beta), alphas)
    assert kd!=-1
    kk = alpha_inv[kd]
    Mj=M_alphas[j]; Mk=M_alphas[kk]; Ml=M_alphas[l]
    alpha1 = x + to_k(apply(Mk, inf))
    alpha2 = y + to_k(apply(Ml.inverse(), inf))
    assert apply(Mi*Tz*Mj,cusp(alpha1)) == cusp(alpha2)
    xr,xi = x
    yr,yi = y
    zr,zi = z
    print("{} Q {} {} {} {} {} {} {} {} {} {}".format(d, i, j, kk, l, xr,xi, yr,yi, zr,zi))
    return [[i,j,kk,l],[x,y,z]]

def aaa_triangle_parameters(T, alphas, M_alphas, strict=True):
    """For T a triangle with all vertices principal cusps, returns
    [[i,j,k],u] where M_alphas[i](alpha[j]+u) = alpha[k] + x with x
    integral, where T has vertices [alpha_i, oo, alpha_j].

    We consider rotations and reflections.
    """
    k = nf(T[0])
    d = -k.disc().squarefree_part()
    if k.class_number()==1:
        for T0 in [tri0(k), tri1(k), tri2(k)]:
            if poly_equiv(T, T0, sign=0):
                print("universal triangle")
                return T0
    inf = cusp(oo, k)
    T = std_poly(T, alphas)[1]
    assert T[1] == inf
    for Jloop in range(2):
        if Jloop:
            T = std_poly(apply(Jmat,T), alphas)[1]
        for Rloop in range(3):
            if Rloop:
                T = std_poly(cycle_poly(T), alphas)[1]
            i = alpha_index(T[0], alphas)
            j, u = alpha_index_with_translation(T[2], alphas)
            if i>=0 and j>=0:
                k, x = alpha_index_with_translation(apply(M_alphas[i], T[2]), alphas)
                assert k!=-1
                assert translate_cusp(alphas[k],x) == apply(M_alphas[i], translate_cusp(alphas[j],u))
                if u==0 or not strict:
                    ur, ui = u
                    print("{} T {} {} {} {} {}".format(d, i, j, k, ur, ui))
                    return [[i,j,k], u]
    return aaa_triangle_parameters(T, alphas, M_alphas, False)

def symmetries(S):
    n = len(S)
    return [U for U in [poly_equiv_exact(S, cycle_poly(S,i), -1) for i in range(n)]
            +  [poly_equiv_exact(S, reverse_poly(cycle_poly(S,i)), -1) for i in range(n)] if U]

def std_aas_triangle(T, alphas):
    """Standard form for an aas-triangle: cycle so that the singular
    vertex is the last, then transform so that the first edge is a
    basic edge.  The result will have the form [alpha,oo,sigma+u]
    """
    while T[2].ideal()==1:
        T = cycle_poly(T)
    return std_poly(T, alphas)[1]

def aas_triangle_gl2_orbit_reps(Tlist, alphas):
    reps = []
    for T in Tlist:
        T = std_aas_triangle(T, alphas)
        if T not in reps and std_aas_triangle(apply(Jmat,T), alphas) not in reps:
            T = std_aas_triangle(reverse_poly(T), alphas)
            if T not in reps and std_aas_triangle(apply(Jmat,T), alphas) not in reps:
                reps.append(T)
    return reps

def aas_triangle_sl2_orbit_reps(Tlist, alphas):
    reps = []
    for T in Tlist:
        T = std_aas_triangle(T, alphas)
        if T not in reps:
            T = std_aas_triangle(reverse_poly(T), alphas)
            if T not in reps:
                reps.append(T)
    return reps

def aas_symmetries(T):
    """
    T should be a triangle with 3rd vertex only non-principal
    """
    return [Imat] + [U for U in [poly_equiv_exact(T, [T[1],T[0],T[2]], -1) ] if U]

def aas_triangle_parameters(T, alphas, M_alphas, sigmas):
    """
    For T a triangle with two vertices principal cusps and one singular,
    returns [[i,j,k],u] where M_alphas[i](sigmas[j]+u) = sigmas[k] + x
    with x integral, where T has vertices [alpha_i, oo, sigma_j+u].

    """
    k = nf(T[0])
    d = -k.disc().squarefree_part()
    T = std_aas_triangle(T, alphas)
    assert T[1] == cusp(oo, k)
    i = alpha_index(T[0], alphas)
    assert i!=-1
    j, u = sigma_index_with_translation(T[2], sigmas)
    assert j!=-1
    k, x = sigma_index_with_translation(apply(M_alphas[i], T[2]), sigmas)
    assert k!=-1
    ur, ui = u
    print("{} U {} {} {} {} {}".format(d, i, j, k, ur, ui))
    return [[i,j,k], u]

def hexagon_parameters(H, alphas, M_alphas):
    """
    For a principal hexagon [a_i, oo, a_j, b_2, gamma, b_1]
    """
    k = nf(H[0])
    d = -k.disc().squarefree_part()
    H = std_poly(H, alphas)[1]
    assert H[1]==cusp(oo, k)
    i = alpha_index(H[0], alphas)
    assert i>=0
    kk, x1 = alpha_index_with_translation(apply(M_alphas[i], H[5]), alphas)
    m, y1 = alpha_index_with_translation(apply(M_alphas[kk]*Tmat(-x1)*M_alphas[i], H[4]), alphas)
    j, u = alpha_index_with_translation(H[2], alphas)
    assert j>=0
    l, x2 = alpha_index_with_translation(apply(M_alphas[j]*Tmat(-u), H[3]), alphas)
    n, y2 = alpha_index_with_translation(apply(M_alphas[l]*Tmat(-x2)*M_alphas[j]*Tmat(-u), H[4]), alphas)
    gamma1 = apply(M_alphas[i].inverse()*Tmat(x1)*M_alphas[kk].inverse()*Tmat(y1), alphas[m])
    gamma2 = apply(Tmat(u)*M_alphas[j].inverse()*Tmat(x2)*M_alphas[l].inverse()*Tmat(y2), alphas[n])
    #print(gamma1, gamma2)
    assert gamma1==gamma2
    ur,ui = u
    x1r,x1i = x1
    y1r,y1i = y1
    x2r,x2i = x2
    y2r,y2i = y2
    print("{} H {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(d, i, j, kk, l, m, n, ur,ui, x1r,x1i, y1r,y1i, x2r,x2i, y2r,y2i))
    return [[i,j,kk,l,m,n],[u,x1,y1,x2,y2]]
