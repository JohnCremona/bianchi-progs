from sage.all import infinity, oo, NFCusp, Matrix
try:
    assert k
except NameError:
    from Q43 import k, Ok, w, zero, one, emb, rootd, alphas, M_alphas, all_alphas, n_alphas

J = Matrix(2,2,[-one, zero, zero, one])
Smat = Matrix(2,2,[zero, -one, one, zero])
inf = cusp(oo)
def Tmat(x):
    return Matrix(2,2, [one,x,zero,one])

def to_k(c):
    return k(c.numerator()/c.denominator()) if c.denominator() else oo

def round(z): # z in k
    x, y = emb(z)
    b = (2*y/rootd).round()
    a = ((2*x).round() - b)//2
    return a+b*w

def ispos(a):
    return a[1]>0 or (a[1]==0 and a[0]>0)

def makepos(a):
    return a if ispos(a) else -a

def posgen(I):
    a = I.number_field()(I.gens_reduced()[0])
    return makepos(a)

def posdet(M):
    return makepos(M.det())

def kgcd(a,b):
    return posgen(k.ideal(a,b))

def frac(x):
    """
    Return coprime a,b such that x=a/b.  x can be an element of k or a cusp.
    """
    if x == infinity:
        return [one, zero]
    a = x.numerator()
    b = x.denominator()
    c = kgcd(a,b)
    return [k(a/c),k(b/c)]

def cusp(x):
    """
    Return x as a reduced cusp [a:b], or [1:0] if x is oo.
    """
    return NFCusp(k,infinity) if x==infinity else NFCusp(k,frac(k(x)))

def polygon(xlist):
    """
    Convert a list of k union {oo} to a list of cusps.
    """
    return [cusp(x) for x in xlist]

def edge(x,y):
    """
    Convert a pair of k union {oo} elements to a list of 2 cusps.
    """
    return polygon([x,y])

def is_integral(M):
    """
    Return True iff all entries of M are integral.
    """
    return all(c in Ok for c in M.list())

def edge_to_mat(e):
    """Return M such that M(0)=e[0], M(oo)=e[1], with det(M) "positive".
    """
    a,c = frac(e[1])
    b,d = frac(e[0])
    M = Matrix(2,2,[a,b,c,d])
    return M if ispos(M.det()) else M * J

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

def alpha_index(a):
    """
    Return i if a==all_alphas[i], else -1.
    """
    try:
        return all_alphas.index(a)
    except ValueError:
        return -1

def alpha_index_with_translation(a):
    """
    Return (i,x) if a==all_alphas[i]+x with x integral, else (-1,0).
    """
    ka = to_k(a)
    for i, alpha in enumerate(all_alphas):
        x = ka - to_k(alpha)
        if x in Ok:
            return (i, x)
    return (-1, 0)

def std_edge(e, sign=+1):
    """When sign=+1, return (i, U) where e=U(e0), U in SL(2,Ok) and
    e0={alpha[i],oo}, else return (-1, None). Note that U,u0 are exist
    and are unique for edges of the tessellation.

    When sign=-1, return (j,V) with det(V)=-1 and e=V{alpha[j],oo}.
    To do this we apply the sign=+1 case to J(e) and adjust the
    resulting U.

    Usually V=U*J and alpha[j]=-alpha[i], except in the cases where
    alpha[i]=w/2, alpha[j]=(w-1)/2 (or vice versa) when
    V=U*[-1,w;0,1].  But we do not need these.

    """
    d = edet(e, True)
    if d not in alphas:
        d = -d
        if d not in alphas:
            return (-1, None)
    Me = edge_to_mat(e if sign==+1 else apply(J,e))
    for alpha in alphas[d]:
        #print(alpha_index(alpha), alpha)
        Ma = edge_to_mat([alpha, inf])
        U = Me * Ma.inverse()
        if is_integral(U):
            return (alpha_index(alpha), U if sign==+1 else J*U)
    return (-1, None)

def edge_equiv(e1,e2):
    """Return U in SL(2,O) such that U maps e1 to e2, otherwise False.
    Note that such a U exists if and only if there exists V with
    det(V)=-1 mapping e1 to e2.
    """
    i1, U1 = std_edge(e1) # e1 = U1{alpha[i1],oo}
    assert i1>=0
    i2, U2 = std_edge(e2) # e2 = U2{alpha[i2],oo}
    assert i2>=0
    if i1 != i2:
        return False
    # Now e2 = U(e1) with U=U2*U1^{-1} and det(U)=1, but U may not be integral
    U = U2*U1.inverse()
    return U if is_integral(U) else False

def std_poly(P):
    """Return (U,P') such that U(P')=P, with the first edge of P' a basic
    edge and det(U)=1.  If the first edge of P is not the transform of
    a basic edge, return (False, None)
    """
    i, U = std_edge(P[:2])
    if i<0:
        return (False, None)
    return U, apply(U.inverse(),P)

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
        P3 = apply(J,P2)
        e3 = P3[:2]
        U = edge_equiv(e1,e3) # det(U)=1
        if U and apply(U,P1)==P3: # then JU(P1)=J(P3)=P2
            return J*U
    return False

def poly_equiv(P1,P2, cycle=True, reverse=True, sign=1):
    """Test whether the polygons P1, P2 are GL(2,O)-equivalent (of
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

def add_alpha(a,b,c,d):
    global M_alphas, all_alphas, n_alphas
    assert a*d-b*c==1
    M_alphas.append(Matrix(2,2,[a,b,c,d]))
    all_alphas.append(cusp(-d/c))
    n_alphas += 1

def add_two_alphas(s, r): # r^2=-1 (mod s)
    t = -(r*r+1)/s
    add_alpha( r, t, s, -r) # alpha =  r/s
    add_alpha(-r, t, s,  r) # alpha = -r/s

def add_four_alphas(s, r1, r2): # r1*r2=-1 (mod s)
    t = -(r1*r2+1)/s
    add_alpha( r2, t, s, -r1) # alpha =  r1/s
    add_alpha(-r2, t, s,  r1) # alpha = -r1/s
    add_alpha( r1, t, s, -r2) # alpha =  r2/s
    add_alpha(-r1, t, s,  r2) # alpha = -r2/s

tri0 = [cusp(0), inf, cusp(1)]

def make_triangles():
    """Return two lists, Tlist and triangles where each entry in Tlist is
    a triple (i,j,n) of alpha-indices such that M_alphas[i](alphas[j])
    = alphas[n] + x with x in Ok, and the corresponding entry in
    triangles is (alphas[i], oo, alphas[j]).  This is the triangle
    whose edges are

    {alphas[i], oo}
    {oo, alphas[j]} = image of {alphas[j'], oo} under M_alphas[j']
    {alphas[j], alphas[i]} = image of {alphas[n], oo} under M_alphas[i']*T^x

    """
    Tlist = []
    for i, Ma in enumerate(M_alphas):
        for j in range(i+1, n_alphas):
            n = alpha_index_with_translation(apply(Ma, all_alphas[j]))[0]
            if n>=0:
                print("({},{},{}) : <{},{},{}>".format(i,j,n,to_k(all_alphas[i]), oo, to_k(all_alphas[j])))
                Tlist.append((i,j,n))
    return Tlist, [[all_alphas[i],NFCusp(k,oo),all_alphas[j]] for i,j,_ in Tlist]

def check_poly_in_list(T, polys, sign=0):
    """Return True if T is GL(2,Ok)-congruent to a polygon in polys or
    the reverse of one.  For SL(2,Ok)-equivalence, set sign=1.
    """
    return any(poly_equiv(T,t, sign=sign) for t in polys)

def reduce_triangles(Tlist, triangles):
    Tlist0 = []
    triangles0 = []
    for trip, tri in zip(Tlist, triangles):
        if not check_poly_in_list(tri, triangles0):
            Tlist0.append(trip)
            triangles0.append(tri)
    return Tlist0, triangles0
