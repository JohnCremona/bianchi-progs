# Functions for working with H3 and hemispheres etc.

from itertools import chain #, combinations

from sage.all import (Infinity, Matrix, ZZ, QQ, RR, CC, NumberField,
                      Graph, srange, Set, sign, var, implicit_plot3d, NFCusp, Integer, oo,
                      infinity, polygen, point, line, circle)

from utils import (nf, to_k, cusp, cusp_label, Imat, apply,
                   translate_cusp, negate_cusp, conj_cusp,
                   smallest_ideal_class_representatives)

from alphas import precomputed_alphas

def make_k(dk):
    """
    Given a negative fundamental discriminant, constructs the associated imaginary quadratic field
    and returns a dict containing this and useful other data
    """
    x = polygen(QQ)
    if dk%4==1:
        k = NumberField(x**2-x+(1-dk)//4, 'w')
    else:
        k = NumberField(x**2-dk//4, 'w')
    assert k.discriminant() == dk
    w = k.gen()
    emb = next(e for e in k.embeddings(CC) if e(w).imag()>0)
    return {'k': k, 'dk': dk, 'w': w, 'wbar': w.trace()-w, 'Ok': k.ring_of_integers(),
            'emb': emb, 'Ymax': emb(w).imag()/2,
            'Ireps': [c.ideal() for c in k.class_group()]}

# Points of H_3 are represented as pairs [z,t2] where z is in k and t2
# in QQ is the square of the height (so the actual point coordinates
# are (z,sqrt(t2))).

# Each principal cusp alpha=r/s with (r,s)=(1) determines the
# hemisphere S_alpha with equation |z-alpha|^2+t^2=1/|s|^2, or
# N(s*z-r)+N(s)*t^2=1.

def radius_squared(alpha):
    """
    For a principal cusp alpha, return the square radius of S_alpha.
    """
    return 1/alpha.denominator().norm()

def cusp_to_point(alpha):
    """
    For a principal cusp alpha = a in k, return the point [a,
    radius_squared(alpha)].
    """
    return [to_k(alpha), radius_squared(alpha)]

def tri_inter(a0, a1, a2):
    """Returns the triple intersection point of the hemispheres S_a_i,
    where a0, a1, a2 are principal cusps, if there is one, as a pair
    [z,t2] where z is in k and t2 in QQ is the square of the vertical
    coordinate.
    """
    alist = [a0,a1,a2]
    # Check the cusps are principal, not infinity, and with unit ideal
    assert all((not a.is_infinity()) and (a.ideal()==1) for a in alist)
    # Define the square radii and centres
    rho0, rho1, rho2 = [radius_squared(a) for a in alist]
    al0, al1, al2 = [to_k(a) for a in alist]
    n0, n1, n2 = [a.norm() for a in [al0, al1, al2]]
    #
    delta = al1*(al0-al2).conjugate() + al2*(al1-al0).conjugate() + al0*(al2-al1).conjugate()
    if delta==0:
        return None
    z = (al1*(n0-n2+rho2-rho0) + al2*(n1-n0+rho0-rho1) + al0*(n2-n1+rho1-rho2)) / delta
    t2 = rho0 - n0 - z.norm() + 2*(al0*z.conjugate()).real()
    assert t2 == rho1 - n1 - z.norm() + 2*(al1*z.conjugate()).real()
    assert t2 == rho2 - n2 - z.norm() + 2*(al2*z.conjugate()).real()
    return None if t2<0 else [z,t2]

def bi_inter(a1, a2):
    """Returns the point on the intersection of the hemispheres S_a_i
    (where a1, a2 are principal cusps) which is on the line from a1 to
    a2, as a pair [z,t2] where z is in k and t2 in QQ is the square of
    the vertical coordinate.

    Use: when both S_a_i pass through a singular point.
    """
    alist = [a1,a2]
    # Check the cusps are principal, not infinity, and with unit ideal
    assert all((not a.is_infinity()) and (a.ideal()==1) for a in alist)
    # Define the square radii and centres
    rho1, rho2 = [radius_squared(a) for a in alist]
    al1, al2 = [to_k(a) for a in alist]
    n1, n2 = [a.norm() for a in [al1, al2]]
    #
    delta = al2-al1
    z = ((al1+al2) + (rho1-rho2)*delta/delta.norm())/2
    t2 = rho1 - n1 - z.norm() + 2*(al1*z.conjugate()).real()
    assert t2 == rho2 - n2 - z.norm() + 2*(al2*z.conjugate()).real()
    return None if t2<0 else [z,t2]

def is_under(P, a):
    """
    Returns -1,0,+1 according as P is over, on, under S_a (a principal)
    """
    z, t2 = P
    ad = a.denominator()
    al = a.numerator()/ad
    return sign(1/ad.norm() - (z-al).norm() - t2)

def is_inside(a, b, strict=False):
    """Returns True iff a is inside (or strictly inside) the circle
    centred on b, where a,b are cusps with b principal.
    """
    k = nf(a)
    d2 = (to_k(a,k)-to_k(b,k)).norm()
    r2 = radius_squared(b)
    if strict:
        return d2 < r2
    else:
        return d2 <= r2

def covering_hemispheres1(P, option=None):
    """For P=[z,t2] in H_3, returns a list of cusps alpha such that P lies
    on or under S_alpha.

    If option is 'exact' only returns alpha for which P is on S_alpha exactly.
    If option is 'strict' only returns alpha for which P is strictly under S_alpha.
    Otherwise (default), returns alpha for which P is under or on S_alpha.

    """
    alphas = []
    z, t2 = P
    k = z.parent()
    a = z.numerator()   # in O_K
    b = z.denominator() # in Z
    sbound = (1/t2).floor()
    for snorm in range(1,1+sbound):
        umax = b*b*(1-snorm*t2)
        for s in k.elements_of_norm(snorm):
            #print("s = {}".format(s))
            if option=='exact':
                urange = [umax] if umax in ZZ else []
            else:
                urange = srange(umax.floor()+1)
            sa = s*a
            #print("umax={}, urange={}".format(umax,list(urange)))
            for unorm in urange:
                if unorm<umax or option != 'strict':
                    for u in k.elements_of_norm(unorm):
                        #print("  u = {}".format(u))
                        for rb in [sa+u, sa-u] if u else [sa]:
                            r = rb/b
                            #print("    r = {}".format(r))
                            if r.is_integral() and k.ideal(r,s)==1:
                                alphas.append(cusp(r/s, k))
    return alphas

def covering_hemispheres2(P, option=None, debug=False):
    """For P=[z,t2] in H_3, returns a list of cusps alpha such that P lies
    on or under S_alpha.

    If option is 'exact' only returns alpha for which P is on S_alpha exactly.
    If option is 'strict' only returns alpha for which P is strictly under S_alpha.
    Otherwise (default), returns alpha for which P is under or on S_alpha.

    """
    alphas = []
    z, t2 = P
    k = z.parent()
    a = z.numerator()   # in O_K
    b = z.denominator() # in Z
    sbound = (1/t2).floor()
    if debug:
        print("t2={} so bound on N(s) = {}".format(t2, sbound))
    for snorm in srange(1,1+sbound):
        for s in k.elements_of_norm(snorm):
            sz = s*z
            d1 = 1/snorm - t2
            assert d1>=0
            if debug:
                print("s = {}, norm {}: d1 = {}".format(s, snorm, d1))
            rbound = ((RR(sz.norm()).sqrt()+1)**2).floor()
            if debug:
                print("Bound on N(r) = {}".format(rbound))
            for rnorm in srange(1+rbound):
                for r in k.elements_of_norm(rnorm):
                    if k.ideal(r,s)!=1:
                        continue
                    for pm in [-1,1]:
                        a = pm*r/s
                        d = d1 - (a-z).norm()
                        if debug and d>=0:
                            print("a = {}, d = {}".format(a, d))
                        # we need d==0 for exact, d>0 for strict, else d>=0
                        ok = (d>0) if option=='strict' else (d==0) if option=='exact' else (d>=0)
                        if ok:
                            a = cusp(a,k)
                            if debug:
                                print(" OK {}".format(a))
                            alphas.append(a)
    return alphas

def covering_hemispheres_test(P, option=None):
    res1 = covering_hemispheres1(P, option)
    res2 = covering_hemispheres2(P, option)
    if sorted(res1) != sorted(res2):
        print("old and new disagree for P={}".format(P))
    return res1

covering_hemispheres = covering_hemispheres2

def hemispheres_through(P):
    return covering_hemispheres(P, 'exact')

def properly_covering_hemispheres(P):
    return covering_hemispheres(P, 'strict')

def is_maximal(P):
    return len(properly_covering_hemispheres(P))==0

def apply3d(M, P):
    """
    Return M(P) where M is in SL(2,O_K) and P=[z,t2] in H3.
    """
    z, t2 = P
    k = z.parent()
    try:
        a, b, c, d = [k(r) for r in M.list()]
    except AttributeError:
        a, b, c, d = [k(r) for r in M]
    n = (c*z+d).norm() + c.norm()*t2
    new_z = ((a*z+b)*(c*z+d).conjugate() + a*c.conjugate()*t2) / n
    new_t2 = t2 / n**2
    return [new_z, new_t2]

def infinity_matrix(a, P=None, Plist=None):
    """For a principal cusp a, returns M_a in GL(2,O_K) with M_a(a)=oo.
    If P is given in H3 it should be an interior point on S_a and then
    the matrix will be adjusted by premultiplying by a translation so
    that M_a(P) is in Plist.
    """
    M0 = Matrix(2,2,a.ABmatrix()).inverse()
    if P is None:
        return M0
    else:
        Q = apply3d(M0,P)
        if Q in Plist:
            return M0
        for R in Plist:
            z = R[0]-Q[0]
            if z.is_integral():
                M = Matrix([[1,z],[0,1]])*M0
                if apply3d(M,P) not in Plist:
                    print("P = {}".format(P))
                    print("M = {}".format(M))
                    print("M(P) = {}".format(apply3d(M,P)))
                assert apply3d(M,P) in Plist
                return M
        raise RuntimeError("infinity_matrix failed")


def singular_points_in_class(I, IC=None, verbose=False):
    """Given an ideal I, return a list of singular points of class [I]
    (one representative for each orbit under integral translations).

    Uses the new characterization of singular points as a/b for b one
    nonzero element of minimal norm in one non-principal ideal I in
    each ideal class, where I=(a,b).

    IC can be set to a list of ideal class representatives.

    """
    k = I.number_field()
    if I.is_principal():
        return [NFCusp(k, infinity)]
    if IC is None:
        IC = smallest_ideal_class_representatives(k)
    sigmas = []
    Inorm = I.norm()
    Ibar = k.ideal(Inorm)/I
    s = k(I.norm())
    slist = [s]
    if I!=Ibar:
        I2 = I*I
        if I2.is_principal():
            s2 = I2.gens_reduced()[0]
            assert s.norm()==s2.norm()
            slist.append(s2)
    if verbose:
        print("Ideal class #{}: denominators {}".format(IC.index(I), slist))
    for s in slist:
        rlist = [r for r in k.ideal(s).residues() if k.ideal(r,s) == I]
        ss = [cusp(reduce_mod_Ok(r/s), k, IC) for r in rlist]
        if verbose:
            print(" - denominator s = {}, numerators {}, sigmas {}".format(s, rlist, ss))
        sigmas += ss
    return sigmas

def singular_points_by_class(IC, verbose=False):
    """Return a list of lists of singular points, one sublist for each
    nontrivial ideal class, representative for each orbit under
    integral translations.

    Uses the new characterization of singular points as a/b for b one
    nonzero element of minimal norm in one non-principal ideal I in
    each ideal class, where I=(a,b).

    """
    return [singular_points_in_class(I, IC=IC, verbose=verbose) for I in IC]

def singular_points_new(k, verbose=False):
    """Return a list of singular points, one representative for each
    orbit under integral translations.

    Uses the new characterization of singular points as a/b for b one
    nonzero element of minimal norm in one non-principal ideal I in
    each ideal class, where I=(a,b).
    """
    return sum(singular_points_by_class(smallest_ideal_class_representatives(k), verbose), [])

def ab_to_k(k,ab):
    """MA's code returns each singular point in the form (a,b) with a,b
    rational, representing a+b*sqrt(-d) with d squarefree.  We convert
    to an element of k, assuming that k's defining polynomial is
    either X^2+d or X^2-X+(d+1)/4.
    """
    w = k.gen()
    rootd = 2*w-1 if k.discriminant()%4 else w
    a,b = ab
    return a+b*rootd

def singular_points_MA(k):
    """
    Singular points from MA's code
    """
    if k.class_number()==1:
        return []
    from FundDomands import singular_points, reduce_ab_mod_ok
    S = singular_points(k)
    # include negatives:
    S = S + [[-ab[0],-ab[1]] for ab in S]
    # reduce mod O_k
    S = [reduce_ab_mod_ok(k, ab) for ab in S]
    # convert to field elements
    S = [ab_to_k(k,ab) for ab in S]
    # remove repeats
    S = list(set(S))
    # convert into cusps whose ideals are standardised, and prepend oo
    IC = smallest_ideal_class_representatives(k)
    return [cusp(oo,k)] + [cusp(s,k,IC) for s in S]

def differ_by_integer(s,t):
    """
    If s,t are cusps, return True iff s-t is integral
    """
    if s.is_infinity():
        return t.is_infinity()
    if t.is_infinity():
        return False
    ks = s.numerator()/s.denominator()
    kt = t.numerator()/t.denominator()
    return (ks-kt).is_integral()

def test_singular_points(dmin, dmax, verbose=False):
    x = polygen(QQ)
    for d in srange(dmin,dmax+1):
        if not d.is_squarefree():
            continue
        k = NumberField(x**2-x+(d+1)//4 if d%4==3 else x**2+d, 'w')
        h = k.class_number()
        if h==1:
            continue
        if verbose:
            print("d={}, {} has class number {}".format(d, k, h))
        sigmas = singular_points_new(k)
        if verbose:
            print("New sigmas: {}".format(sigmas))
        old_sigmas = singular_points_MA(k)
        if verbose:
            print("Old sigmas: {}".format(old_sigmas))
        diff1 = [s for s in sigmas if not any(differ_by_integer(s,t) for t in old_sigmas)]
        diff2 = [s for s in old_sigmas if not any(differ_by_integer(s,t) for t in sigmas)]
        ok = True
        if diff1:
            ok = False
            print("d={}: sigmas from new code not in old: {}".format(d,diff1))
        if diff2:
            ok = False
            print("d={}: sigmas from old code not in new: {}".format(d,diff2))
        if ok:
            print("Old and new agree for d={}".format(d))


def tau(P1, P2):
    """Given P_i=[alpha_i,rho_i^2] for i=1,2, where alpha_i=r_i/s_i are
    principal cusps defining hemispheres (or circles) with square
    radius rho_i^2=1/N(s_i), return +2, +1, 0, -1, -2 as follows:

    +2 if they do not intersect and are external to each other
    +1 if they are externally tangent
    0  if they intersect in two distinct points
    -1 if they are internally tangent (or equal)
    -2 if they do not intersect and one is inside the other
    """
    a1, r1sq = P1
    a2, r2sq = P2
    d1 = (a1-a2).norm() - (r1sq + r2sq)
    d2 = d1**2 - 4*r1sq*r2sq
    if d2 < 0:
        return 0
    return sign(d1) * (1 if d2==0 else 2)

def circles_intersect(P1,P2):
    return tau(P1,P2)==0

def circles_tangent(P1,P2, exterior=True):
    return tau(P1,P2) == (+1 if exterior else -1)

def circle_inside_circle(P1,P2, strict=True):
    # t = P1[1]<P2[1] and tau(P1,P2)==-2
    # if strict or t:
    #     return t
    # return (P1[1]<P2[1] and tau(P1,P2)==-1)

    t1 = (P1[1]<P2[1]) if strict else (P1[1]<=P2[1])
    t2 = tau(P1,P2) in ([-2] if strict else [-2,-1])
    return t1 and t2

def xy_coords(alpha):
    """
    alpha = x+y*sqrt(-d) in k = Q(w) with either w=sqrt(-d) or w=(1+sqrt(-d))/2
    """
    x, y = list(alpha)
    if alpha.parent().gen().trace():
        y /=2
        x +=y
    return (x,y)

def reduce_mod_Ok(alpha):
    """
    Return in integer translate of alpha whose xy-coords satisfy
    -1/2 < x <= 1/2 and
    -1/2 < y <= 1/2 (even discriminant, w=sqrt(-d))
    -1/4 < y <= 1/4 (odd discriminant, w=(1+sqrt(-d))/2)
    """
    k = alpha.parent()
    w = k.gen()
    y = xy_coords(alpha)[1]
    r = 2 if w.trace() else 1
    alpha -= (r*y).round('down')*w
    x = xy_coords(alpha)[0]
    alpha -= x.round('down')
    assert in_rectangle(alpha)
    return alpha

def slope2(x,y):
    """
    Function used to order nonzero (x,y) in R^2 via their argument, going clockwise around the origin:

    (+,-) < (0,-) < (-,-) < (-,0) < (-,+) < (0,+) < (+,+) < (+,0)
    """
    return (sign(y), x/y) if y else (sign(x), Infinity)

def slope(alpha, centre=0):
    """
    As above for elements of an imaginary quadratic field k, assuming
    k=Q(w) with either w=sqrt(-d) or w=(1+sqrt(-d))/2.
    """
    return slope2(*xy_coords(alpha-centre))

def slope_before(es1, es2):
    e1, s1 = es1
    e2, s2 = es2
    return (e1==e2 and s1<=s2) or (e1!=e2 and s1>s2)

def in_first_half(alpha1, alpha2, centre=0):
    """
    Return True if the clockwise angle from alpha1 round to alpha2 is < pi.
    """
    return slope_before(slope(alpha1, centre), slope(alpha2, centre))

# plotting functions taken essentiall from MA

def plot1hemi(kdata, H):
    """
    kdata is a dict with keys 'k' (the field), 'emb' (embedding of k into CC)
    H = [z, rsq] with z in k defines a hemisphere
    """
    X, Y, Z = var('X, Y, Z')
    Ymax = kdata['Ymax']
    Xmax = 0.5
    x0, y0 = kdata['emb'](H[0])
    eq = (X - x0)**2 + (Y - y0)**2 + Z**2 - H[1]
    return implicit_plot3d(eq, (Y, -Ymax, Ymax ),  (X, -Xmax, Xmax), (Z, 0, 1), plot_points=60, aspect_ratio=1, color='lightgreen')

def plot_Bianchi_diagram(k, Hlist):
    """
    Hlist is a list of hemispheres H = [z,rsq] with z in k and square radius rsq
    """
    kdata = make_k(k.discriminant())
    return sum([plot1hemi(kdata, H) for H in Hlist])

def circ(c,r, fill):
    return circle(c, r,
                  aspect_ratio=1,
                  edgecolor='blue' if fill else 'black',
                  thickness=1 if fill else 2,
                  alpha = 0.2,
                  fill=fill)

def disc(c,r):
    return circle(c, r,
                  aspect_ratio=1, fill=True, rgbcolor='blue', alpha = 0.2)

def plot_circles_and_points(cc, pp, fill=False):
    circles = [circ(c, r, fill) for c, r in cc]
    points = [point(P, rgbcolor='red', pointsize=50) for P in pp]
    return sum(circles) + sum(points)

def plot_circles(alist, fill=False):
    k = nf(alist[0])
    emb = next(e for e in k.embeddings(CC) if e(k.gen()).imag()>0)
    A = [list(emb(to_k(a, k))) for a in alist]
    R = [RR(radius_squared(a)).sqrt() for a in alist]
    circles = [(c,r) for c,r in zip(A,R)]
    return plot_circles_and_points(circles, [], fill)

def plot_FunDomain_projection(k, alphas, sigmas, fill=False):
    w = k.gen()
    D = k.discriminant().abs()
    emb = next(e for e in k.embeddings(CC) if e(w).imag()>0)
    rootd = emb(w).imag()
    Ymax = 3*rootd/(4 if ZZ(D).mod(4) == 3  else 2)
    Xmax = 3*0.5

    triplets, extra_alphas = alpha_triples(alphas)

    A = [list(emb(to_k(a))) for a in alphas+extra_alphas]
    R = [RR(radius_squared(a)).sqrt() for a in alphas+extra_alphas]
    #print("circle centres: {}".format(A))

    S = [list(emb(to_k(s))) for s in sigmas if not s.is_infinity()]
    #print("singular points: {}".format(S))

    C = [list(emb(P[2][0])) for P in triplets]
    #print(triplets)
    #print("corners: {}".format(C))

    circles  = [(c, r) for c, r in zip(A,R)]
    proj = plot_circles_and_points(circles, S, fill=fill)

    z = w-ZZ(1)/2 if ZZ(D).mod(4)==3 else w
    TL=list(emb((-1+z)/2))
    TR=list(emb((1+z)/2))
    BR=list(emb((1-z)/2))
    BL=list(emb((-1-z)/2))
    lines = [line([TL,TR], rgbcolor='black'),
             line([TR,BR], rgbcolor='black'),
             line([BR,BL], rgbcolor='black'),
             line([BL,TL], rgbcolor='black')]

    proj += sum(lines)
    proj.set_axes_range(-Xmax, Xmax, -Ymax, Ymax)
    return proj

def is_redundant(P, alphas):
    """Return True iff P is strictly covered by any of the hemispheres
    S_a for a in the list alphas.
    """
    return any(is_under(P,a)==1 for a in alphas)

def triple_intersections(alphas):
    """Given a list of principal cusps alpha (all reduced mod O_k) return
    a list of "corners" P = [z,tsq] each the intersection of an S_a
    with at least two other S_{b+t} with z in the fundamental
    rectangle and tsq>0.

    Let u = (w-wbar)/2.  The fundamental rectangle F has TR corner at
    (u+1)/2 and BL corner minus this.  Using symmetries (negation and
    conjugation) we can work with the quarter-rectangle F4 with the
    same TR and BL=0.  To recover F from F4 take the union of
    z,-z,zbar,-zbar for z in F4.

    The 9 quarter-rectangles adjacent to F4 consist of

    -z, zbar, 1-zbar; -zbar, z, 1-zbar; u-z, u_zbar, u+1-zbar

    for z in F4.
    """
    w = alphas[0].number_field().gen()
    u = (w-w.conjugate())/2

    # Extract the alphas in F4:
    alphas4 = [a for a in alphas if cusp_in_quarter_rectangle(a)]

    # Extend these by 8 translations:
    def nbrs(a):
        k = nf(a)
        w = k.gen()
        z = to_k(a, k)
        cz = z.conjugate()
        zlist = [-z, cz, 1-z, -cz, 1-cz, w-z, w+cz,
                 cz+w-1 if w.trace() else 1+w-z]
        alist = [cusp(z2, k) for z2 in zlist]
        for b in alist:
            if not b.ideal()==1:
                print("cusp {} is a neighbour of principal cusp {} but is not principal".format(b,a))
        return alist

    xalphas4 = sum([nbrs(a) for a in alphas4], alphas4)
    n = len(xalphas4)

    # convert each cusp to a point P = [z,tsq] with tsq the square
    # radius of S_a:

    Alist = [cusp_to_point(a) for a in xalphas4]

    corners4 = []
    for a, A in zip(alphas4, Alist):
        bb = [(b,B) for b,B in zip(xalphas4, Alist) if circles_intersect(A,B)]
        for b,B in bb:
            cc = [c for c,C in bb if circles_intersect(B,C)]
            # now each pair of {a,b,c} intersect
            for c in cc:
                P = tri_inter(a, b, c)
                if P and P[1] and in_quarter_rectangle(P[0]) and not is_redundant(P, xalphas4) and P not in corners4:
                    corners4.append(P)

    # These corners are in F4, so we apply symmetries to get all those in F:
    corners = []
    for P in corners4:
        z = P[0]
        zbar = z.conjugate()
        for z2 in [z, -z, zbar, -zbar]:
            if in_rectangle(z2):
                P2 = [z2, P[1]]
                if P2 not in corners:
                    corners.append(P2)
    return corners

def alpha_triples(alphas):
    """Given a list of principal cusps
    alpha (all reduced mod O_k)
    return (1) a list of
    [tsq,(a1,a2,a3),P] where each ai is
    the translate of an alpha, P =
    [z,tsq] is a "corner", the triple
    intersection of the S_ai with P
    in the fundamental rectangle and
    tsq>0; (2) a list of the extra translates required

    """
    w = alphas[0].number_field().gen()
    corners = []
    triples = []
    alpha_translates = []
    # Extend the alphas by 8 translations:
    xalphas = alphas + sum([[translate_cusp(a,t) for t in [-w-1,-w,1-w,-1,1,-1+w,w,1+w]] for a in alphas], [])
    n = len(xalphas)
    # convert each cusp to a point
    # [a,tsq] with tsq the square
    # radius of S_a:
    Alist = [cusp_to_point(a) for a in xalphas]

    # get a list of pairs {i,j} with
    # i<j such that S_ai and S_aj
    # intersect properly:

    ij_list = [{i,j} for i,ai in
               enumerate(Alist) for j,aj in
               enumerate(Alist) if i<j and
               circles_intersect(ai, aj)]

    for i,j in ij_list:
        ai = xalphas[i]
        aj = xalphas[j]
        for k, ak in enumerate(alphas):# in range(max(i,j)+1, n):
            if {i,k} in ij_list and {j,k} in ij_list:
                #ak = xalphas[k]
                P = tri_inter(ai, aj, ak)
                if P and P[1] and in_rectangle(P[0]) and not is_redundant(P, xalphas):
                    if P not in corners:
                        trip = [P[1],(ai,aj,ak),P]
                        triples.append(trip)
                        corners.append(P)
                        for a in trip[1]:
                            if a not in alpha_translates and a not in alphas:
                                alpha_translates.append(a)
    triples.sort(key = lambda t:t[0])
    return triples, alpha_translates

def alpha_doubles(alphas, sigmas):
    """Given alphas, a list of principal cusps alpha (all reduced mod O_k)
    and sigmas, a list of all singular points sigma in one nontrivial
    ideal class, returns

    (1) a list of [rsq,(a1,a2),P] where each ai is the translate of an
    alpha, P = [z,rsq] is the bi-intersection of S_a1, S_a2 with z in
    the fundamental rectangle and rsq>0;

    (2) a list of the extra alpha-translates required.

    """
    w = alphas[0].number_field().gen()
    SP = [[to_k(s),0] for s in sigmas if not s.is_infinity()]
    # get list of alpha translates a such that at least one sigma is on S_a:
    xalphas = sum([[translate_cusp(a,t) for t in [a+b*w for a in [-1,0,1] for b in [-1,0,1]]] for a in alphas], [])
    xalphas = [a for a in xalphas if any(is_under(P,a)==0 for P in SP)]
    # convert each principal cusp to
    # a point [a,tsq] with tsq the
    # square radius of S_a:
    Alist = [cusp_to_point(a) for a in xalphas]
    #print("Alist: {}".format(Alist))

    # Now we consider each pair (ai,aj) with i<j which intersect properly:
    ij_list = [{i,j} for i,ai in
               enumerate(Alist) for j,aj in
               enumerate(Alist) if i<j and
               circles_intersect(ai, aj)]
    #print("ij_list: {}".format(ij_list))

    doublets = []
    alpha_translates = []

    for i,j in ij_list:
        ai = xalphas[i]
        aj = xalphas[j]
        P = bi_inter(ai, aj)
        if P and P[1] and not is_redundant(P, xalphas):
            sig = [s for s in SP if is_under(s,ai)==0 and is_under(s,aj)==0]
            assert len(sig)==1
            dub = [P[1],(cusp(sig[0][0]),ai,aj),P]
            if dub not in doublets:
                doublets.append(dub)
                for a in dub[1]:
                    if a not in alpha_translates and a not in alphas:
                        alpha_translates.append(a)
    return doublets, alpha_translates

def orbit_polyhedron(orb, Plist, Pverts, Pmats):
    i = orb[0]
    P = Plist[i]
    E = []
    for j in orb:
        Q = Plist[j]
        QV  = Pverts[j]
        MQP = [M for M in Pmats[j] if apply3d(M,Q)==P]
        V = [[apply(M, a) for a in QV] for M in MQP]
        E += sum([[[t[0],t[i]] for i in range(1,len(QV))] for t in V], [])
    G = Graph([[cusp_label(a),cusp_label(b)] for a,b in E])
    return G

def principal_polyhedra(alphas, debug=False):
    k = alphas[0].number_field()
    triplets, extra_alphas = alpha_triples(alphas)
    # only used for the 3d plot:
    hemispheres = [cusp_to_point(a) for a in alphas+extra_alphas]

    if debug:
        print("{} triplets:".format(len(triplets)))
        for t in triplets:
            print(t)

    Plist = [t[2] for t in triplets]
    if debug:
        print("Plist: {}".format(Plist))
    Psupps = [hemispheres_through(P) for P in Plist]
    if debug:
        print("Psupps: {}".format(Psupps))
    Pmats = [[Imat] + [infinity_matrix(a, P, Plist) for a in Psupp] for P, Psupp in zip(Plist, Psupps)]
    if debug:
        print("Pmats: {}".format(Pmats))
    Pverts = [[cusp(oo,k)] + list(t[1]) for t in triplets]
    orbits = set()
    for i,P in enumerate(Plist):
        Qlist = [apply3d(M,P) for M in Pmats[i]]
        orb = Set([Plist.index(Q) for Q in Qlist])
        if orb not in orbits:
            print("New orbit from P_{}={}: {}".format(i,P,orb))
            orbits.add(orb)
    orbits = [list(orb) for orb in orbits]
    print("Found {} orbits:".format(orbits))
    polyhedra = [orbit_polyhedron(orb, Plist, Pverts, Pmats) for orb in orbits]
    print("Constructed {} polyhedra".format(len(polyhedra)))
    print("Faces: {}".format([[len(F) for F in G.faces()] for G in polyhedra]))
    return polyhedra, hemispheres

def singular_polyhedra(alphas, sigmas, debug=False):
    k = alphas[0].number_field()
    doublets, extra_alphas = alpha_doubles(alphas, sigmas)

    if debug:
        print("{} doublets:".format(len(doublets)))
        for t in doublets:
            print(t)

    Rlist = [d[2] for d in doublets]
    if debug:
        print("Rlist: {}".format(Rlist))
    Rsupps = [hemispheres_through(R) for R in Rlist]
    Rmats = [[Imat] + [infinity_matrix(a, R, Rlist) for a in Rsupp] for R,Rsupp in zip(Rlist, Rsupps)]
    Rverts = [[cusp(oo,k)] + list(t[1]) for t in doublets]
    orbits = set()
    for i,P in enumerate(Rlist):
        #print("i={}, R={}, Rmats={}".format(i,P,Rmats))
        Qlist = [apply3d(M,P) for M in Rmats[i]]
        orb = Set([Rlist.index(Q) for Q in Qlist])
        if orb not in orbits:
            print("New orbit from R_{}={}: {}".format(i,P,orb))
        orbits.add(orb)
    orbit_reps = [list(orb) for orb in orbits]
    print("Found {} orbits with representative points {}:".format(len(orbits), orbit_reps))
    polyhedra = [orbit_polyhedron(orb, Rlist, Rverts, Rmats) for orb in orbit_reps]
    print("Constructed {} polyhedra".format(len(polyhedra)))
    print("Faces: {}".format([[len(F) for F in G.faces()] for G in polyhedra]))
    #faces = sum([G.faces() for G in polyhedra],[])
    return polyhedra

def all_polyhedra(k, debug=False):
    alphas = precomputed_alphas(k)
    sigmas = singular_points_by_class(smallest_ideal_class_representatives(k))[1:]
    polys, hemis = principal_polyhedra(alphas, debug)
    polys += sum([singular_polyhedra(alphas, sigs, debug) for sigs in sigmas], [])
    return polys, hemis

def is_poly_principal(T):
    return all(a.ideal()==1 for a in T)

half = Integer(1)/2

def xy_in_rectangle(xy, f):
    """
    f = 1 or 2
    """
    x,y = xy
    fy = f*y
    return -half<x and x<= half and -half<fy and fy<=half

def xy_in_quarter_rectangle(xy, f):
    """
    f = 1 or 2
    """
    x,y = xy
    fy = f*y
    return 0<=x and x<= half and 0<=fy and fy<=half

def in_rectangle(a):
    f = 1 + a.parent().disc()%2
    return xy_in_rectangle(xy_coords(a), f)

def in_quarter_rectangle(a):
    f = 1 + nf(a).disc()%2
    return xy_in_quarter_rectangle(xy_coords(a), f)

def cusp_in_rectangle(a):
    return in_rectangle(to_k(a))

def cusp_in_quarter_rectangle(a):
    return in_quarter_rectangle(to_k(a))

def is_sigma_surrounded(sigma, alist, debug=False):
    """Given a singular point s and a candidate list of principal cusps
    alist, tests whether the discs S_{a+t} for a in alist and t in Ok
    completely surround sigma.

    Returns either (True, xlist) with xlist a list of all a+t needed,
    or (False, [])

    """

    k = nf(alist[0])
    w = k.gen()

    # convert s to a point
    s = to_k(sigma, k)
    if debug:
        print("s = {}".format(s))

    # extend the candidate list by including offsets:

    offsets = [-1-w,-w,1-w,-1,1,-1+w,w,1+w]
    alist = sum([[translate_cusp(b,t) for t in offsets] for b in alist], alist)

    # extract the relevant alphas, if any:

    alist = [a for a in alist if is_under([s,0], a)==0]
    alist = [a for a in alist if not any(circle_inside_circle(cusp_to_point(a), cusp_to_point(b), False)
                                         for b in alist if b!=a)]

    if debug:
        print(" relevant alphas: {}".format(alist))

    Alist = [cusp_to_point(a) for a in alist]
    # sort these by slope:
    Alist.sort(key=lambda a: slope(a[0], s))
    Aslopes = [slope(a[0], s) for a in Alist]
    if debug:
        print(" Alist (sorted) = {}".format(Alist))
        print(" Relative slopes: {}".format(Aslopes))

    for i, t2 in enumerate(Aslopes):
        t1 = Aslopes[i-1]
        if not slope_before(t1, t2):
            if debug:
                print(" !Failure around {} between {} and {}".format(s, alist[i-1], alist[i]))
            return False, []
    return True, alist

def are_sigmas_surrounded(sigmas, alist, debug=False):
    """Given a list of singular points s and a candidate list of principal
    cusps alist, tests whether the discs S_{a+t} for a in alist and t
    in Ok completely surround all sigmas.

    Returns either (True, xlist) with xlist a list of all a+t needed,
    or (False, [])

    """
    xlist = []
    for s in sigmas:
        if s.is_infinity():
            continue
        ok, xlist1 = is_sigma_surrounded(s, alist, debug)
        if not ok:
            if debug:
                print("{} is not surrounded".format(s))
            return False, s
        if debug:
            print("{} is surrounded by {}".format(s, xlist1))
        for a in xlist1:
            if a not in xlist:
                xlist.append(a)

    if debug:
        print("All sigmas are surrounded, by {}".format(xlist))
    return True, xlist

def tri_det(a1, a2, a3):
    return Matrix(3,3,[a1,a2,a3, a1.conjugate(), a2.conjugate(), a3.conjugate(), 1, 1, 1]).det()

def intersection_points_in_k(a1,a2):
    """Given principal cusps a1,a2 returns a list of 0, 1 or 2 points (in
    k) where the circles S_a1, S_a2 intersect.
    """
    k = nf(a1)
    alist = [a1,a2]
    # Check the cusps are principal, not infinity, and with unit ideal
    assert all((not a.is_infinity()) and (a.ideal()==1) for a in alist)
    # Define the square radii and centres
    r1sq, r2sq = [radius_squared(a) for a in alist]
    al1, al2 = [to_k(a, k) for a in alist]
    delta = al2-al1
    n = delta.norm()
    d1 = n - (r1sq + r2sq)
    d2 = d1**2 - 4*r1sq*r2sq
    if d2 > 0:
        return []
    z = ((al1+al2) + (r1sq-r2sq)/delta.conjugate())/2
    return [z + r/(2*delta.conjugate()) for r in k(d2).sqrt(all=True, extend=False)]

def intersection_points_in_CC(a1,a2):
    """Given principal cusps a1,a2 returns a list of 0, 1 or 2 points (in
    CC) where the circles S_a1, S_a2 intersect.
    """
    k = nf(a1)
    emb = next(e for e in k.embeddings(CC) if e(k.gen()).imag()>0)
    alist = [a1,a2]
    # Check the cusps are principal, not infinity, and with unit ideal
    assert all((not a.is_infinity()) and (a.ideal()==1) for a in alist)
    # Define the square radii and centres
    r1sq, r2sq = [radius_squared(a) for a in alist]
    al1, al2 = [to_k(a, k) for a in alist]
    delta = al2-al1
    n = delta.norm()
    d1 = n - (r1sq + r2sq)
    d2 = d1**2 - 4*r1sq*r2sq
    if d2 > 0:
        return []
    z = emb(((al1+al2) + (r1sq-r2sq)/delta.conjugate())/2)
    if d2 == 0:
        return [z]
    rd2 = CC(d2).sqrt() # pure imaginary
    z1 = z + rd2/(2*emb(delta.conjugate()))
    z2 = 2*z-z1 # = z - rd2/(2*emb(delta.conjugate()))
    return [z1,z2]

def show_intersection(a1,a2):
    zz = intersection_points_in_CC(a1,a2)
    if len(zz)==2:
        zz.append((zz[0]+zz[1])/2)
    points = [list(z) for z in zz]
    k = nf(a1)
    emb = next(e for e in k.embeddings(CC) if e(k.gen()).imag()>0)
    A = [list(emb(to_k(a, k))) for a in [a1,a2]]
    R = [RR(radius_squared(a)).sqrt() for a in [a1,a2]]
    circles = [(c,r) for c,r in zip(A,R)]
    return plot_circles_and_points(circles, points, True)

def are_intersection_points_covered_by_one(a1, a2, a, plot=False):
    """Given principal cusps a1, a2, a such that the circles S_a1 and
    S_a2 intersect in distinct points, test whether S_a covers either
    or both.

    Returns 0 if neither, 2 if both, +1 or -1 if just one.  The signs
    are consistent so that if a returns +1 and a' returns -1 then each
    intersection point is covered by either S_a or S_a'.
    """
    k = nf(a1)
    emb = next(e for e in k.embeddings(CC) if e(k.gen()).imag()>0)
    alist = [a1,a2,a]
    # Check the cusps are principal, not infinity, and with unit ideal
    assert all((not a.is_infinity()) and (a.ideal()==1) for a in alist)

    # Define the square radii and centres
    r1sq, r2sq, rsq = [radius_squared(a) for a in alist]
    al1, al2, al = [to_k(a, k) for a in alist]
    n1, n2 = [a.norm() for a in [al1, al2]]
    #
    delta = al2-al1
    n = delta.norm()
    z0 = ((al1+al2) + (r1sq-r2sq)/delta.conjugate())/2
    d1 = n - (r1sq + r2sq)
    d2 = d1**2 - 4*r1sq*r2sq
    if d2 >= 0:
        raise RuntimeError("cusps {} and {} have non-intersecting circles")

    if plot:
        points = [list(z) for z in intersection_points_in_CC(a1,a2)]
        circle = (list(emb(to_k(a, k))), RR(radius_squared(a)).sqrt())
        pic = plot_circles([a1,a2], False) + plot_circles_and_points([circle], points, True)
        pic.show()
        input("press Enter...")

    T = 2 * n * (rsq - (z0-al).norm()) + d2/2 # rational
    T2 = T**2
    D = tri_det(al, al2, al1) # pure imaginary
    D2 = QQ(D**2)             # negative rational
    d2D2 = d2*D2              # positive rational

    # the covering condition is \pm sqrt(d2)*D < T
    #print("T = {}, D = {}, d2 = {}".format(T,D,d2))

    if d2D2 < T2:
        return 2 if T>0 else 0 if T<0 else '?'
    if d2D2 > T2:
        u = QQ(D/(w-w.conjugate()))
        return -1 if u>0 else +1 if u<0 else 0
    return 0

def is_singular(s, sigmas):
    from utils import sigma_index_with_translation
    return sigma_index_with_translation(s, sigmas)[0]!=-1

def translates(a):
    return [translate_cusp(a,t) for t in [-w-1,-w,-w+1,-1,0,1,w-1,w,w+1]]

def is_inside_one(z, alist):
    """Test whether the cusp z is strictly inside at least one S_a for a
    in alist.  If so return True, a; otherwise return False, None.
    """
    try:
        a = next(a for a in alist if is_inside(z, a, strict=True))
        return True, a
    except StopIteration:
        return False, None

def are_intersection_points_covered(a0, a1, alist, sigmas, debug=False):
    """Given principal cusps a0, a1 whose circles S_a0, S_a1 intersect,
    and a list of principal cusps alist each of whose circles S_a also
    intersects S_a0, test whether each of the two intersection points
    of S_a0 and S_a1 is either singular or strictly inside one of the
    S_a.

    We treat as a special case when the two intersection points are in
    k.  If not, the code still uses exact arithmetic.

    """
    k = nf(a0)
    z_in_k = intersection_points_in_k(a0,a1)
    if z_in_k:
        zz = [cusp(z, k) for z in z_in_k]
        if debug:
            print("intersection points in k: {}".format(z_in_k))
        # check that each is *either* singular *or* contained in some S_a2
        for z in zz:
            if is_singular(z, sigmas):
                if debug:
                    print("{} is ok: singular".format(z))
            else:
                ok, a1 = is_inside_one(z, alist)
                if ok:
                    if debug:
                        print("{} is ok: inside S_{}".format(z, a1))
                else:
                    return False
        return True

    # Now the intersection points are not in k. Check that either one
    # S_a covers both, or two cover one each:
    t = 0 # will hold +1 or -1 if we have covered only one of the two
    for a2 in alist:
        if a2 == a1:
            continue
        t2 = are_intersection_points_covered_by_one(a0, a1, a2, plot=False)
        if debug:
            print("a0={}, a1={}, a2={}:  t2={}, t={}".format(a0, a1,a2,t2,t))
        if t2: # it is 2, +1 or -1
            assert t2 in [-1,1,2]
            if debug:
                are_intersection_points_covered_by_one(a0, a1, a2, plot=True)
            if t2==2 or ([t,t2] in [[1,-1],[-1,1]]):
                if debug:
                    print("t={}, t2={}, about to return True".format(t,t2))
                return True
            assert t2 in [-1,1] and t in [0,t2]
            if debug:
                print("t={}, t2={}, setting t to {}".format(t,t2,t2))
            t = t2
    return False

def is_alpha_surrounded(a0, alist, sigmas, debug=False, plot=False):
    """Given a principal cusp a0, a candidate list of principal cusps
    alist, tests whether the boundary of the disc S_a0 is contained in
    the union of the translates S_{b+t} for b in alist, apart from any
    singular points on the boundary.  It suffices to consider all b+t
    such that S_{b+t} intersects S_a in two points and check that each
    of the points is either singular or contained in some other
    S_{b+t}.  This is simplest when the intersection points are in k;
    if not then the method still uses exact arithmetic in k
    throughout.

    Returns either (True, xlist) with xlist a list of all b+t
    needed, or (False, None)

    """
    k = nf(alist[0])
    w = k.gen()
    emb = next(e for e in k.embeddings(CC) if e(w).imag()>0)

    # convert a0 to a point with radius
    A0 = cusp_to_point(a0)
    if debug:
        print("A0 = {}".format(A0))

    # extend the candidate list by including offsets:

    alist = sum([translates(b) for b in alist], [])

    # check if S_a0 is strictly entirely contained in one S_alpha:
    if any(circle_inside_circle(A0, cusp_to_point(b), True) for b in alist):
        if debug:
            a1 = next(b for b in alist if circle_inside_circle(A0, cusp_to_point(b), True))
            print(" ok: circle {} is entirely inside circle {}".format(A0, cusp_to_point(a1)))
        return True

    # extract the relevant alphas, if any, namely those for which
    # S_alpha and S_a0 properly intersect:

    alist = [a for a in alist if circles_intersect(A0, cusp_to_point(a))]
    # alist = [a for a in alist if not any(circle_inside_circle(cusp_to_point(a), cusp_to_point(b), False)
    #                                      for b in alist if b!=a)]

    if debug:
        print(" relevant alphas: {}".format(alist))

    if debug and plot:
        pic = plot_circles([a0], False) + plot_circles(alist, True)
        pic.show(figsize=[30,30])
        input("press Enter...")

    for i, a1 in enumerate(alist):
        if debug:
            print("\nTesting intersection points of {} and {}".format(a0,a1))
        ok = are_intersection_points_covered(a0, a1, alist, sigmas, debug)
        if ok:
            if debug:
                print(" - ok: intersection points of {} and {} are covered".format(a0,a1))
        else:
            if debug:
                print(" - not ok: intersection points of {} and {} are not covered".format(a0,a1))
            return False
    if debug:
        print("OK: all intersection points of {} and {} are covered".format(a0, a1))

    return True

def are_alphas_surrounded(alist, slist, verbose=False, debug=False):
    """Given alist, a candidate list of principal cusps, and slist, a
    complete list of singular points, tests whether the boundary of
    every disc S_a is contained in the union of the translates of the
    S_b for b in alist, apart from any singular points on the
    boundary.

    Returns either (True, xlist) with xlist a list of any translates
    needed, or (False, [])

    """
    for i, a in enumerate(alist):
        if verbose or debug:
            print("Testing alpha #{}/{} = {}".format(i+1, len(alist), a))
        ok = is_alpha_surrounded(a, alist, slist, debug)
        if not ok:
            if verbose or debug:
                print(" no, {} is not surrounded".format(a))
            return False
        if verbose or debug:
            print(" ok, {} is surrounded".format(a))
    return True

def next_norm(k, n):
    """
    Returns the smallest integer m>=n which is a norm from O_k
    """
    while not k.elements_of_norm(n):
        n+=1
    return n

def elements_of_norm(k, n):
    return iter(k.elements_of_norm(n))

def elements_of_norm_upto(k, n, start=1):
    return chain(*(iter(k.elements_of_norm(n)) for n in range(start, n+1)))

def reduced_numerators(s):
    k = s.parent()
    one = ZZ(1)
    for r in k.ideal(s).residues():
        if k.ideal(r,s).norm() == one:
            yield r

def principal_cusps_iter(k, maxnorm_s):
    """
    Iterator yielding all principal r/s with N(s)<=maxnorm_s
    """
    for s in elements_of_norm_upto(k, maxnorm_s):
        for r in reduced_numerators(s):
            a = reduce_mod_Ok(r/s)
            yield cusp(a, k)

def principal_cusps_up_to(k, maxn, fussy=True):
    """List of all principal r/s with N(s)<=maxnorm_s, omitting any whose
    circles are contained in an earlier circle.  Since we loop through
    circles in decreasing order of radius, no circle can be contained
    in a later one.

    If fussy, automatically increment maxn to the next integer which is a norm.
    """
    alist = []
    Alist = []
    if fussy:
        maxn0 = maxn
        maxn = next_norm(k, maxn0)
        if maxn != maxn0:
            print(" increasing maxn to {} since there are no elements of norm {}".format(maxn, list(range(maxn0,maxn))))
    for a in principal_cusps_iter(k, maxn):
        A = cusp_to_point(a)
        #print("Testing {} = {} against {}".format(a, A, alist))
        if not any(circle_inside_circle(A, B, False) for B in Alist):
            #print("appending {} = {} to {}".format(a, A, alist))
            alist.append(a)
            Alist.append(A)
    return alist

def find_covering_alphas(k, sigmas=None, verbose=False):
    """Returns a finite list of principal cusps a such that the S_{a+t}
    for all integral t cover CC apart from singular points.

    For n>=1 successively, we test as a candidate set all a=r/s with
    r,s coprime, r reduced mod s, N(s)<=n (omitting any for which S_a
    is contained in any earlier S_a') until we succeed.

    sigmas can be set to a list of singular points (up to
    translation), otherwise these will be computed.

    Returns maxn, alphas, sigmas

    Other functions will then (1) saturate the set, (2) discard
    redundancies.

    """
    if sigmas is None:
        sigmas = singular_points_new(k)
    ok = False
    maxn = 0
    while not ok:
        maxn = next_norm(k, maxn+1)
        if verbose:
            print("Testing max norm {}".format(maxn))
        alphas = principal_cusps_up_to(k, maxn)
        ok = are_alphas_surrounded(alphas, sigmas, debug=False)
        if verbose and not ok:
            print("{} fails, continuing...".format(maxn))
        if ok:
            if verbose:
                print("Success with max norm {}!".format(maxn))
            return maxn, alphas, sigmas

def point_translates(P):
    return [[P[0]+t,P[1]] for t in [-1-w,-1,-1+w,-1,0,1,-1+w,w,1+w]]

def nverts(a, plist):
    return len([P for P in plist if is_under(P,a)==0])

def saturate_covering_alphas(k, alphas, debug=False):
    """Given a covering set of alphas as produced by
    find_covering_alphas(), add extras if necessary so that they are
    "saturated", i.e. define the extended fundamental domain.

    By Swan, we need to find the points P in H^3 with positive height
    where at least 3 hemispheres S_a intersect, and for each P check
    whether P is properly covered by an S_a for a not in the set of
    alphas (up to translation).  If so, we need to add a to the set of
    alphas.  If none, then we have the fundamental region (and can go
    on to discard any redundant alphas).

    """
    sat = False
    checked_points = []
    alphas1 = [a for a in alphas] # copy so original list unchanged
    while not sat:
        n = max(a.denominator().norm() for a in alphas1)
        m = next_norm(k, n+1)
        all_points = triple_intersections(alphas1)
        if debug:
            print("Found {} potential vertices".format(len(all_points)))
        points = [P for P in all_points if P[1]<=1/m]
        if debug:
            print(" -- of which {} are low enough to be properly covered by a new alpha".format(len(points)))
        points = [P for P in points if in_quarter_rectangle(P[0])]
        if debug:
            print(" -- of which {} lie in the first quadrant".format(len(points)))
        points = [P for P in points if P not in checked_points]
        if debug:
            print(" -- of which {} have not already been checked".format(len(points)))
        sat = True        # will be set to False if we find out that the alphas are not already saturated
        extra_alphas = [] # will be filled with any extra alphas needed
        for P in points:
            if debug:
                print(" - checking P = {}".format(P))
            extras = properly_covering_hemispheres(P)
            if extras:
                sat = False
                hts = [radius_squared(a) - (P[0]-to_k(a)).norm() for a in extras]
                m = max(hts)
                extras0 = [a for a,h in zip(extras, hts) if h==m]
                if debug:
                    print("   - found properly covering {} with norms {}".format(extras, [a.denominator().norm() for a in extras]))
                    print("     max height above P (height {}) is {}, for {}".format(P[1], m, extras0))
                for a in extras0:
                    ca = conj_cusp(a)
                    for b in [a, negate_cusp(a), ca, negate_cusp(ca)]:
                        if cusp_in_rectangle(b) and b not in extra_alphas:
                            extra_alphas.append(b)
            else:
                print("   - OK, no properly covering alphas found")
                checked_points.append(P)
        if sat:
            print(" alphas are saturated")
        else:
            print(" alphas not saturated, {} extras needed: {}".format(len(extra_alphas), extra_alphas))
        alphas1 += extra_alphas

    # Now delete any alphas with <3 vertices, allowing for translates
    pointsx = []
    for P in all_points:
        for Q in point_translates(P):
            if Q not in pointsx:
                pointsx.append(Q)
    alphas1 = [a for a in xalphas if nverts(a, pointsx)>=3]
    points1 = triple_intersections(alphas1)
    return alphas1, points1

def reduce_alphas_mod_Ok(alist):
    """Rahm's list of alpha = lambda/mu in k includes repeats (up to
    translation by Ok).

    This function returns a list with no repeats, as cusps.
    """
    a0 = next(a for a in alist if not a in QQ)
    k = a0.parent()
    Ireps = smallest_ideal_class_representatives(k)
    alist = [k(a) for a in alist]
    alphas = []
    for a in alist:
        if not any((a-b).is_integral() for b in alphas):
            alphas.append(reduce_mod_Ok(a))
        else:
            print("omitting a={} which is a translate of {}".format(a, next(b for b in alist if (a-b).is_integral())))
    return [cusp(a, k, Ireps) for a in alphas]

def cong_mod(r1, r2, s):
    return ((r1-r2)/s).is_integral()

def find_edge_pairs(alphas, debug=False):
    from utils import add_two_alphas, add_four_alphas, nf, ispos

    k = nf(alphas[0])
    A1 = [a for a in alphas if a.denominator()==1]
    A2 = [a for a in alphas if a.denominator()==2]
    A = [a for a in alphas if a.denominator() not in [1,2,3]]
    S = list(set(k(a.denominator()) for a in A))
    S.sort(key = lambda z: z.norm())
    if debug:
        print("Denominator 1: {}".format(A1))
        print("Denominator 2: {}".format(A2))
        print("Other denominators: {}".format(S))
        for s in S:
            print("s = {}: numerators {}".format(s, [a for a in A if a.denominator()==s]))

    new_alphas = []
    M_alphas = []
    pluspairs = []
    minuspairs = []
    fours = []

    for s in S:
        if debug:
            print("s = {}".format(s))
        As = [a for a in A if a.denominator()==s]
        for a in As:
            if debug:
                print("  a = {}".format(a))
            r = k(a.numerator())
            if debug:
                print("  r = {}".format(r))
            rs = (r,s)
            mrs = (-r,s)
            rsq = r*r
            if cong_mod(rsq, +1, s):
                if not any(pair in pluspairs for pair in (rs, mrs)):
                    if ispos(r):
                        if debug:
                            print("  - adding plus pair {}".format(rs))
                        pluspairs.append(rs)
                        add_two_alphas(s, r, +1, new_alphas, M_alphas)
                    else:
                        if debug:
                            print("  - adding plus pair {}".format(mrs))
                        pluspairs.append(mrs)
                        add_two_alphas(s, -r, +1, new_alphas, M_alphas)
                continue
            if cong_mod(rsq, -1, s):
                if not any(pair in minuspairs for pair in (rs, mrs)):
                    if ispos(r):
                        if debug:
                            print("  - adding minus pair {}".format(rs))
                        minuspairs.append(rs)
                        add_two_alphas(s, r, -1, new_alphas, M_alphas)
                    else:
                        if debug:
                            print("  - adding minus pair {}".format(mrs))
                        minuspairs.append(mrs)
                        add_two_alphas(s, -r, -1, new_alphas, M_alphas)
                continue
            if debug:
                print("  - looking for a foursome")
            try:
                adash = next(ad for ad in As if cong_mod(r*ad.numerator(), -1, s))
                rdash = k(adash.numerator())
                rds = (rdash,s)
                mrds = (-rdash,s)
                if not any(pair in fours for pair in (rs, mrs, rds, mrds)):
                    if ispos(r):
                        if debug:
                            print("  - adding foursome {}".format(rs))
                        fours.append(rs)
                        add_four_alphas(s, r, rdash, new_alphas, M_alphas)
                    else:
                        if debug:
                            print("  - adding foursome {}".format(mrs))
                        fours.append(mrs)
                        add_four_alphas(s, -r, -rdash, new_alphas, M_alphas)
            except StopIteration:
                print("no negative inverse found for {} mod {}".format(r, s))
    print("alphas with denominator 1: {}".format(A1))
    print("alphas with denominator 2: {}".format(A2))
    print("plus pairs: {}".format(pluspairs))
    print("minus pairs: {}".format(minuspairs))
    print("fours: {}".format(fours))
    return A1, A2, new_alphas, M_alphas, pluspairs, minuspairs, fours
