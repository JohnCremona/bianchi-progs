# Functions for working with H3 and hemispheres etc.

from itertools import groupby

from sage.all import (Infinity, Matrix, vector, ZZ, QQ, RR, CC,
                      Graph, Set, sign, var, implicit_plot3d, oo,
                      point, line, circle, show,
                      Factorization)

from utils import (make_k, nf, to_k, cusp, Imat, apply,
                   translate_cusp, negate_cusp, conj_cusp,
                   smallest_ideal_class_representatives, cusp_class_index,
                   alpha_index_with_translation, sigma_index_with_translation,
                   elements_of_norm, elements_of_norm_upto, next_norm,
                   xy_coords, in_rectangle, cusp_in_rectangle,
                   in_quarter_rectangle, cusp_in_quarter_rectangle,
                   reduce_mod_Ok, singular_points, singular_points_by_class
)

from alphas import precomputed_alphas

from polyhedra import poly_type, poly_types

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

def tri_inter_points(a0, a1, a2):
    """Returns the triple intersection point of the hemispheres S_a_i,
    where a0, a1, a2 are principal cusps, if there is one, as a pair
    [z,t2] where z is in k and t2 in QQ is the square of the vertical
    coordinate.

    In this version, each ai = [z,rsq] where z=r/s is principal and rsq=N(s)
    """
    alist = [a0,a1,a2]
    # Define the square radii and centres
    rho0, rho1, rho2 = [a[1] for a in alist]
    al0, al1, al2 = [a[0] for a in alist]
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

def tri_inter_cusps(a0, a1, a2):
    """Returns the triple intersection point of the hemispheres S_a_i,
    where a0, a1, a2 are principal cusps, if there is one, as a pair
    [z,t2] where z is in k and t2 in QQ is the square of the vertical
    coordinate.
    """
    t = [[to_k(a), radius_squared(a)] for a in [a0,a1,a2]]
    return tri_inter_points(*t)

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

def covering_hemispheres1(P, option=None, norm_s_lb=1, debug=False):
    """For P=[z,t2] in H_3, returns a list of principal cusps alpha=r/s
    such that P lies on or under S_alpha, and N(s)>=norm_s_lb.

    If option is 'exact' only returns alpha for which P is on S_alpha
    exactly, i.e. is_under(P,alpha)==0.

    If option is 'strict' only returns alpha for which P is strictly
    under S_alpha, i.e. is_under(P,alpha)==1.

    Otherwise (default), returns alpha for which P is under or on
    S_alpha, i.e. is_under(P,alpha)>=0.

    Method: alpha=r/s with r,s coprime and |sz-r|^2 + |s|^2t^2 <=1 (or
    <1 or =1).  Write z=a/b with a in O_K, b in Z and u=sa-rb.  Then
    N(u) <= b^2 (1-N(s)T^2), so (1) N(s)<=1/t^2 and (2) given s, N(u)
    <= b^2 (1-N(s)t^2), (3) r = (sa-u)/b if integral.  We loop over s
    and then u (up to sign), testing both signs in (3).

    If 'strict', all these inequalities are strict; if 'exact', they
    are all equality, so in particular N(s) must satisfy b^2*t2*N(s)
    integral. (Here t2 is rational.)

    """
    if debug:
        print(f"Finding covering hemispheres for {P = } with {option = }")
    alphas = []
    z, t2 = P
    k = z.parent()
    a = z.numerator()   # in O_K
    b = z.denominator() # in Z
    sbound = (1/t2).floor()
    if debug:
        print(f"{t2 = } so bound on N(s) is {sbound}")
    for s in elements_of_norm_upto(k, sbound, norm_s_lb):
        snorm = s.norm()
        if debug:
            print(f"{s = } with norm {snorm}")
        umax = b*b*(1-snorm*t2)
        sa = s*a
        if option=='exact':
            ulist = elements_of_norm(k, umax) # will give [] if umax not integral
        else:
            umax = umax.floor()
            if option=='strict':
                if debug:
                    print(f"{umax = }")
                ulist = elements_of_norm_upto(k, umax-1, 0)
            else:
                ulist = elements_of_norm_upto(k, umax, 0)

        for u in ulist:
            if debug:
                print(f"  {u = }")
            for rb in [sa+u, sa-u] if u else [sa]:
                r = rb/b
                if debug:
                    print(f"    {r = }")
                if r.is_integral() and k.ideal(r,s)==1:
                    alphas.append(cusp(r/s, k))
    if debug:
        print(f"Covering hemispheres are S_alpha for alpha = {alphas}")
    covering_codes = [0] if option=='exact' else [1] if option=='strict' else [0,1]
    if not all(is_under(P,a)in covering_codes for a in alphas):
        print(f"Problem in covering_hemispheres1({P}, {option})")
        print(f" - returned cusp list includes {[a for a in alphas if is_under(P,a) not in covering_codes]}")
        print(f"   for which is_under(P,a) is {[is_under(P,a) for a in alphas if is_under(P,a) not in covering_codes]}")
    return alphas

# The following function does the same as the preceding one but is slower for option 'exact'.

def covering_hemispheres2(P, option=None, norm_s_lb=1, debug=False):
    """For P=[z,t2] in H_3, returns a list of principal cusps alpha =r/s
    such that P lies on or under S_alpha, and N(s)>=norm_s_lb.

    If option is 'exact' only returns alpha for which P is on S_alpha exactly.
    If option is 'strict' only returns alpha for which P is strictly under S_alpha.
    Otherwise (default), returns alpha for which P is under or on S_alpha.

    """
    alphas = []
    z, t2 = P
    k = z.parent()
    a = z.numerator()   # in O_K
    sbound = (1/t2).floor()
    if debug:
        print(f"{t2 = } so bound on N(s) is {sbound}")
    for s in elements_of_norm_upto(k, sbound, norm_s_lb):
        snorm = s.norm()
        sz = s*z
        d1 = 1/snorm - t2
        assert d1>=0
        rbound = ((RR(sz.norm()).sqrt()+1)**2).floor()
        if debug:
            print(f"{s = }, {snorm = }: {d1 = }, bound on N(r) is {rbound}")
        for r in elements_of_norm_upto(k, rbound, 0):
            rnorm = r.norm()
            if snorm.gcd(rnorm)>1 and k.ideal(r,s)!=1:
                continue
            for pm in [-1,1] if r else [1]:
                a = pm*r/s
                d = d1 - (a-z).norm()
                if debug and d>=0:
                    print(f"{a = }, {d = }")
                # we need d==0 for exact, d>0 for strict, else d>=0
                ok = (d>0) if option=='strict' else (d==0) if option=='exact' else (d>=0)
                if ok:
                    a = cusp(a,k)
                    if debug:
                        print(f" OK {a}")
                    alphas.append(a)
    if debug:
        print(f"Covering hemispheres are S_alpha for alpha = {alphas}")
    covering_codes = [0] if option=='exact' else [1] if option=='strict' else [0,1]
    assert all(is_under(P,a) in covering_codes for a in alphas)
    return alphas

def covering_hemispheres_test(P, option=None, norm_s_lb=1):
    res1 = covering_hemispheres1(P, option, norm_s_lb)
    res2 = covering_hemispheres2(P, option, norm_s_lb)
    if sorted(res1) != sorted(res2):
        print("old and new disagree for P={}".format(P))
    return res1

def covering_hemispheres(P, option, norm_s_lb=1, debug=False):
    L1 = sorted(covering_hemispheres1(P,option, norm_s_lb))
    if not debug:
        return L1
    L2 = sorted(covering_hemispheres2(P,option, norm_s_lb))
    if L1!=L2:
        print(f"Inconsistent results for covering_hemispheres({P},{option})!")
        print(f"Method 1 gives {L1}")
        print(f"Method 2 gives {L2}")
        assert L1==L2
    return L2

def hemispheres_through(P):
    return covering_hemispheres1(P, 'exact')

def properly_covering_hemispheres(P, norm_s_lb=1):
    return covering_hemispheres2(P, 'strict', norm_s_lb)

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
        assert P[1]==Q[1]
        if Q in Plist:
            return M0
        for R in Plist:
            if R[1]!=Q[1]:
                continue
            z = R[0]-Q[0]
            if z.is_integral():
                M = Matrix([[1,z],[0,1]])*M0
                if apply3d(M,P) not in Plist:
                    print("P = {}".format(P))
                    print("M = {}".format(M))
                    print("M(P) = {}".format(apply3d(M,P)))
                assert apply3d(M,P) in Plist
                return M
        raise RuntimeError("infinity_matrix failed with a={}, P={}, Plist={}".format(a,P,Plist))

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

# plotting functions taken essentially from MA's code

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
    return implicit_plot3d(eq, (Y, -Ymax, Ymax ),  (X, -Xmax, Xmax), (Z, 0, 1), plot_points=60, aspect_ratio=1, color='lightgreen', name=str(kdata['k'].discriminant()))

def plot_Bianchi_diagram(kdata, Hlist):
    """
    Hlist is a list of hemispheres H = [z,rsq] with z in k and square radius rsq
    """
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

def plot_circles_and_points(cc, pp1=[], pp2=[], pp3=[], lines=[], fill=False):
    circles = [circ(c, r, fill) for c, r in cc]
    points1 = [point(P, rgbcolor='black', pointsize=30) for P in pp1]
    points2 = [point(P, rgbcolor='blue', pointsize=30) for P in pp2]
    points3 = [point(P, rgbcolor='green', pointsize=30) for P in pp3]
    ll = [line(l, rgbcolor='black') for l in lines]
    return sum(circles) + sum(points1) + sum(points2) + sum(points3) + sum(ll)

def plot_circles(alist, fill=False, lines=False):
    k = nf(alist[0])
    emb = next(e for e in k.embeddings(CC) if e(k.gen()).imag()>0)
    A = [list(emb(to_k(a, k))) for a in alist]
    R = [RR(radius_squared(a)).sqrt() for a in alist]
    circles = [(c,r) for c,r in zip(A,R)]
    if lines:
        from itertools import combinations
        #ll1 = list(combinations(A,2))
        ll2 = [[list(z) for z in intersection_points_in_CC(*pq)] for pq in combinations(alist,2)]
        ll = ll2
    else:
        ll = []
    return plot_circles_and_points(circles, lines=ll, fill=fill)

def plot_FunDomain_projection(k, alphas, sigmas, fill=False):
    w = k.gen()
    D = k.discriminant().abs()
    emb = next(e for e in k.embeddings(CC) if e(w).imag()>0)
    rootd = emb(w).imag()
    Ymax = (2.5*rootd)/(4 if ZZ(D).mod(4) == 3  else 2)
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
    proj = plot_circles_and_points(circles, S, A, C, fill=fill)

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

def triple_intersections(alphas, debug=False):
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

    -z, zbar, 1-zbar; -zbar, z, 1-zbar; w-z, w+zbar, and either w+1-z or w+zbar-1

    for z in F4.
    """
    K = nf(alphas[0])
    w = K.gen()
    t = w.trace()

    # Extract the alphas in F4:
    A = [to_k(a, K) for a in alphas]
    AlA4 = [(al,a) for al,a in zip(alphas,A) if in_quarter_rectangle(a)]
    #Al4 = [a[0] for a in AlA4] # as cusps
    A4 = [a[1] for a in AlA4]  # as elements of K
    if debug:
        print(f"{len(A4) = }")

    # Extend these by 8 translations:
    def nbrs(z):
        cz = z.conjugate()
        return [z,-z, cz, 1-z, -cz, 1-cz, w-z, w+cz, w+cz-1 if t else w+1-z]

    XA4 = sum([nbrs(a) for a in A4], [])
    XAl4 = [cusp(a,K) for a in XA4]
    if debug:
        print(f"{len(XA4) = }")

    # convert each cusp z to a point P = [z,tsq] with tsq the square
    # radius of S_z:

    from utils import frac
    Alist = [[a,1/frac(a)[1].norm()] for a in XA4]

    # For each i get a list of j>i for which S_ai and S_aj intersect properly
    i2j = {}
    for i, Ai in enumerate(Alist):
        i2j[i] = [i+1+j for j,Aj in enumerate(Alist[i+1:]) if  circles_intersect(Ai, Aj)]
    if debug:
        print(" finished making i2j")

    # Hence make a list of triples (i,j,k) with i<j<k with pairwise proper intersections
    ijk_list = []
    for i, j_list in i2j.items():
        ijk_list += [(i,j,k) for j in j_list for k in j_list if k in i2j[j]]
    if debug:
        print(f" finished making ijk_list: {len(ijk_list)} triples")

    corners = []

    for (i,j,k) in ijk_list:
        P = tri_inter_points(Alist[i], Alist[j], Alist[k])
        if not P:
            continue
        z, t2 = P
        if not t2:
            continue
        # if debug:
        #     print(f" found {P = }")
        if not in_quarter_rectangle(z):
            continue
        if P in corners:
            continue
        if is_redundant(P, XAl4):
            continue
        # These corners are in F4, so we apply symmetries to get all those in F:
        zbar = z.conjugate()
        for z2 in [z, -z, zbar, -zbar]:
            if in_rectangle(z2):
                P2 = [z2, t2]
                if P2 not in corners:
                    if debug:
                        print(f" adding {P2 = } from (i,j,k)={(i,j,k)}")
                    corners.append(P2)

    if debug:
        print(f" returning {len(corners)} corners")
    return corners

def alpha_triples(alphas, debug=False):
    """
    Given a list of principal cusps alpha (all reduced mod O_k) return
    (1) a list of [tsq,(a1,a2,a3),P] where each ai is the translate of
    an alpha, P = [z,tsq] is a "corner", the triple intersection of
    the S_ai, with P in the fundamental rectangle and tsq>0; (2) a list
    of the extra translates required.
    """
    n0 = len(alphas)
    if debug:
        print(f"Finding triple intersection points of {n0} alphas")
    k = nf(alphas[0])
    corners = []
    triples = []
    alpha_translates = []
    # Extend the alphas by 8 translations:
    #w = k.gen()
    # xalphas = sum([[translate_cusp(a,t) for t in
    #                 [a+b*w for a in [-1,0,1] for b in [-1,0,1]]] for a in alphas],
    # [])
    #xalphas = alphas + sum([[translate_cusp(a,t) for t in [-w-1,-w,1-w,-1,1,-1+w,w,1+w]] for a in alphas], [])
    xalphas = sum([translates(a) for a in alphas], [])
    n = len(xalphas)
    if debug:
        print(f"After adding translates number of alphas is {n}")

    # convert each cusp to a point
    # [a,tsq] with tsq the square
    # radius of S_a:
    Alist = [cusp_to_point(a) for a in xalphas]

    # For each i get a list of j>i for which S_ai and S_aj intersect properly
    i2j = {}
    n2 = 0
    for i, ai in enumerate(Alist):
        i2j[i] = j_list = [i+1+j for j,aj in enumerate(Alist[i+1:]) if  circles_intersect(ai, aj)]
        n2 += len(j_list)
    if debug:
        print(f"{n2} ordered pairs intersect properly")
        #print(f"{i2j}")

    # Hence make a list of triples (i,j,k) with i<j<k with pairwise proper intersections
    ijk_list = []
    for i in range(n):
        j_list = i2j[i]
        ijk_list += [(i,j,k) for j in j_list for k in j_list if k in i2j[j]]
    n3 = len(ijk_list)
    if debug:
        print(f"{n3} ordered triples have mutual 2-way intersections")
        #print(f"{ijk_list}")

    for (i,j,k) in ijk_list:
        ai = xalphas[i]
        aj = xalphas[j]
        ak = xalphas[k]
        P = tri_inter_cusps(ai, aj, ak)
        if P and P[1] and in_rectangle(P[0]) and not is_redundant(P, xalphas):
            if P not in corners:
                if debug:
                    print(f"triple ({i},{j},{k}) gives corner {P}")
                trip = [P[1],(ai,aj,ak),P]
                triples.append(trip)
                corners.append(P)
                for a in trip[1]:
                    if a not in alpha_translates and a not in alphas:
                        alpha_translates.append(a)
            else:
                if debug:
                    print(f"triple ({i},{j},{k}) gives repeat of corner {P}")
    triples.sort(key = lambda t:t[0])
    return triples, alpha_translates

def sigma_triples(alphas, sigmas):
    """Given alphas, and sigmas (which can be a complete set or just those
    in one ideal class), returns a list of [rsq, (s, a1,a2), R] where
    each ai is the translate of an alpha, s is a sigma on both S_a1
    and S_a2, and R = [z,rsq] is the bi-intersection (with rsq>0) of
    S_a1, S_a2.

    """
    k = nf(alphas[0])
    # get list of finite sigmas as full points:
    xsigmas = [s for s in sigmas if not s.is_infinity()]
    SP = [[to_k(s),0] for s in xsigmas]

    # get list of alpha translates a such that at least one sigma is on S_a:
    w = alphas[0].number_field().gen()
    xalphas = sum([[translate_cusp(al, a+b*w) for a in [-1,0,1] for b in [-1,0,1]] for al in alphas], [])
    xalphas = [a for a in xalphas if any(is_under(P,a)==0 for P in SP)]

    triples = []

    for s,S in zip(xsigmas, SP):
        # find the alphas a such that S_a passes through s:
        alist = [a for a in xalphas if is_under(S,a)==0]
        # sort these by slope:
        alist.sort(key=lambda a: slope(to_k(a,k), S[0]))

        for i, ai in enumerate(alist):
            aj = alist[i-1]
            R = bi_inter(ai, aj)
            assert R and R[1]
            # test whether this corner R is a translate of one we already have
            old = False
            for t in triples:
                x = t[1][0] - R[0]
                if x.is_integral(): # we have a repeat corner, up to translation
                    old = True
                    R[0] += x
                    trip = [[translate_cusp(c, x) for c in (s,ai,aj)], R]
                    break
            if not old:
                trip = [[s,ai,aj], R]
            triples.append(trip)

    return triples

def orbit_polyhedron(orb, Plist, Pverts, Pmats,  debug=False):
    if debug:
        print(f"\nConstructing orbit polyhedron from {orb=}")
    i = orb[0]
    P = Plist[i]
    if debug:
        print(f"Base point P{i} = {P}")
    E = []
    for j in orb:
        Q = Plist[j]
        QV  = Pverts[j]
        if debug:
            print(f"\n P{j}={Q} with {len(QV)} vertices {QV} and {len(Pmats[j])} matrices ")
            print(f" which map P{j} to {[apply3d(M,Q) for M in Pmats[j]]}")
        # matrices which map Q to P:
        MQP = [M for M in Pmats[j] if apply3d(M,Q)==P]
        if debug:
            print(f" of which {len(MQP)} map P{j} to P{i}")
        # images of Q's vertices under these:
        V = [[apply(M, a) for a in QV] for M in MQP]
        edges = sum([[[str(t[0]), str(t[i])] for i in range(1,len(QV))] for t in V], [])
        if debug:
            print(f" and map P{j}'s vertices to")
            for v in V:
                print(v)
            print(f" -- so adding {len(edges)} edges: {edges}")
        E += edges
    # check edge relation is symmetric:
    # if not all(list(reversed(e)) in E for e in E):
    #     print("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    #     print(f"some edges not symmetric: {[e for e in E if list(reversed(e)) not in E]}")
    #     print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    G = Graph(E)
    Ptype = poly_type(G)
    if Ptype == 'unknown':
        print("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print(f"unrecognised polyhedron from {orb=}")
        print(f" - base point {P} with vertices {Pverts[i]}")
        print(f"{G = }")
        verts = G.vertices(sort=False)
        print(f" - {len(verts)} vertices: {verts}")
        edges = [(verts.index(v), verts.index(w)) for (v,w) in E]
        print(f" - {len(E)} edges: {E}")
        print(f"{edges = }")
        try:
            faces = G.faces()
            print(f" - faces {faces}")
            print(f" - face sizes {[len(F) for F in faces]}")
        except ValueError:
            print(" - graph is not planar")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
        if debug:
            show(G)
    else:
        if debug:
            print(f"Constructed a {Ptype}")
    return G

def principal_polyhedra(alphas, debug=False):
    print("Constructing principal polyhedra")
    k = nf(alphas[0])
    triplets, extra_alphas = alpha_triples(alphas, debug)

    if debug:
        print(f"Extra alphas: {[(a,alpha_index_with_translation(a, alphas)) for a in extra_alphas]}")
    # only used for the 3d plot:
    hemispheres = [cusp_to_point(a) for a in alphas+extra_alphas]

    if debug:
        print(f"{len(triplets)} triplets:")
        for t in triplets:
            print(t)

    Plist = [t[2] for t in triplets]
    if debug:
        print(f"{Plist = }")
        print(" - now finding Psupps: hemispheres through each P...")

    Psupps = [hemispheres_through(P) for P in Plist]
    if debug:
        print(f"Hemispheres through each P: {Psupps}")
        print(" - now finding matrices for each P...")

    # For each a with S_a through P find an inversion matrix M_a which
    # maps P to another Q in Plist:
    Pmats = [[Imat] +
             [infinity_matrix(a, P, Plist) for a in Psupp]
             for P, Psupp in zip(Plist, Psupps)]

    # Each corner P defines a subset of a polyhedron with vertices oo
    # and those translated alphas a which are in P's support, i.e. for
    # which S_a pass through P:
    if debug:
        print(" - now finding vertex set for each P...")
    Pverts = [[cusp(oo,k)] +
              [v for v in Psupp if alpha_index_with_translation(v, alphas)[0] >=0]
              for P,Psupp in zip(Plist, Psupps)]
    if debug:
        print(f"{Pverts = }")
        for P, Pv in zip(Plist, Pverts):
            print(f"{P=}: {[alpha_index_with_translation(v, alphas) for v in Pv[1:]]}")
        print(" - now finding the orbits...")

    # partition the P into orbits:
    orbits = set()
    used_Pi = Set()
    for i,P in enumerate(Plist):
        if i in used_Pi:
            continue
        Qlist = [apply3d(M,P) for M in Pmats[i]]
        orb = Set([Plist.index(Q) for Q in Qlist])
        if debug:
            print(f" + new orbit from P_{i}={P}: {orb}")
        used_Pi = used_Pi.union(orb)
        orbits.add(orb)
    orbits = [list(orb) for orb in orbits]
    if debug:
        print(f"Found {len(orbits)} orbits")
        print(f"{orbits = }")
        print(" - now finding the polyhedron from each orbit...")

    polyhedra = [orbit_polyhedron(orb, Plist, Pverts, Pmats, debug=debug) for orb in orbits]
    if debug:
        npoly = len(polyhedra)
        poly = "polyhedra" if npoly>1 else "polyhedron"
        print(f"Constructed {npoly} {poly}")
    for G in polyhedra:
        try:
            G.faces()
            #print(f"  {[len(F) for F in faces]}")
        except ValueError:
            print("   (not planar)")
    return polyhedra, hemispheres

def singular_polyhedra(alphas, sigmas, debug=False):
    print("Constructing polyhedra from one ideal class, sigmas {}".format(sigmas))
    k = nf(alphas[0])
    triplets = sigma_triples(alphas, sigmas)

    if debug:
        print("{} triplets from sigmas {}:".format(len(triplets), sigmas))
        for t in triplets:
            print(t)

    # we cannot just do Rlist = [t[1] for t in triplets] since there are repeats
    Rlist = list(k for k,_ in groupby(sorted([t[1] for t in triplets])))
    if debug:
        print(f"{Rlist = }")

    Rsupps = [hemispheres_through(R) for R in Rlist]
    if debug:
        print(f"Hemispheres through each R: {Rsupps}")
        print(" - now finding matrices for each R...")

    Rmats = [[Imat] +
             [infinity_matrix(a, R, Rlist) for a in Rsupp]
             for R, Rsupp in zip(Rlist, Rsupps)]

    if debug:
        print(" - now finding vertex set for each R...")
    Rverts = [[cusp(oo,k)] +
              [v for v in Rsupp if alpha_index_with_translation(v, alphas)[0] >=0] +
              sum([[v for v in t[0]] for t in triplets if t[1]==R], [])
              for R,Rsupp in zip(Rlist, Rsupps)]
    if debug:
        print(f"{Rverts = }")
        print(" - now finding the orbits...")

    # partition the R into orbits:
    orbits = set()
    used_Ri = Set()
    for i,R in enumerate(Rlist):
        if i in used_Ri:
            continue
        Qlist = [apply3d(M,R) for M in Rmats[i]]
        orb = Set([Rlist.index(Q) for Q in Qlist])
        if debug:
            print(f" + new orbit from R_{i}={R}: {orb}")
        used_Ri = used_Ri.union(orb)
        orbits.add(orb)
    orbits = [list(orb) for orb in orbits]
    if debug:
        print(f"Found {len(orbits)} orbits")
        print(f"{orbits = }")
        print(" - now finding the polyhedron from each orbit...")

    polyhedra = [orbit_polyhedron(orb, Rlist, Rverts, Rmats, debug=debug) for orb in orbits]
    if debug:
        npoly = len(polyhedra)
        poly = "polyhedra" if npoly>1 else "polyhedron"
        print(f"Constructed {npoly} {poly}")
    for G in polyhedra:
        try:
            G.faces()
            #print(f"  {[len(F) for F in faces]}")
        except ValueError:
            print("   (not planar)")
    return polyhedra

def all_polyhedra(k, alphas=None, debug=False):
    if alphas is None:
        alphas = precomputed_alphas(k)
    sigmas = singular_points_by_class(smallest_ideal_class_representatives(k))[1:]
    polys, hemis = principal_polyhedra(alphas, debug)
    polys += sum([singular_polyhedra(alphas, sigs, debug) for sigs in sigmas], [])
    return polys, hemis

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
    return plot_circles_and_points(circles, points, fill=True)

def are_intersection_points_covered_by_one(a1, a2, a, plot=False):
    """Given principal cusps a1, a2, a such that the circles S_a1 and
    S_a2 intersect in distinct points, test whether S_a covers either
    or both.

    Returns 0 if neither, 2 if both, +1 or -1 if just one.  The signs
    are consistent so that if a returns +1 and a' returns -1 then each
    intersection point is covered by either S_a or S_a'.
    """
    k = nf(a1)
    w = k.gen()
    emb = next(e for e in k.embeddings(CC) if e(w).imag()>0)
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
        pic = plot_circles([a1,a2], fill=False) + plot_circles_and_points([circle], points, fill=True)
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
    return sigma_index_with_translation(s, sigmas)[0]!=-1

def translates(a):
    w = nf(a).gen()
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

def is_alpha_surrounded(a0, alist, sigmas, pairs_ok=[], debug=False, plot=False):
    """Given a principal cusp a0, a candidate list of principal cusps
    alist, tests whether the boundary of the disc S_a0 is contained in
    the union of the S_{b} for b in alist, apart from any
    singular points on the boundary.  It suffices to consider all b
    such that S_{b} intersects S_a in two points and check that each
    of the points is either singular or contained in some other
    S_{b}.  This is simplest when the intersection points are in k;
    if not then the method still uses exact arithmetic in k
    throughout.

    pairs_ok is a list of pairs (a1,a2) whose intersection points are
    known to be covered, so can be skipped.

    Returns (True/False, new_pairs_ok) where new_pairs_ok is an
    updated list of pairs whose intersections have been shown to be
    covered.

    """
    # convert a0 to a point with radius
    A0 = cusp_to_point(a0)
    if debug:
        print("A0 = {}".format(A0))
    Alist = [cusp_to_point(b) for b in alist]

    # check if S_a0 is strictly entirely contained in one S_alpha:
    if any(circle_inside_circle(A0, A, True) for A in Alist):
        if debug:
            A1 = next(A for A in Alist if circle_inside_circle(A0, A, True))
            print(" ok: circle {} is entirely inside circle {}".format(A0, A1))
        return True, pairs_ok

    # extract the relevant alphas, if any, namely those for which
    # S_alpha and S_a0 properly intersect:

    alist = [a for a,A in zip(alist,Alist) if circles_intersect(A0, A)]
    if debug:
        print(" relevant alphas: {}".format(alist))

        if plot:
            pic = plot_circles([a0], fill=False) + plot_circles(alist, fill=True)
            pic.show(figsize=[30,30])
            input("press Enter...")

    a0_pairs_ok = [pr for pr in pairs_ok if a0 in pr]
    new_pairs_ok = pairs_ok.copy()
    all_ok = True
    for i, a1 in enumerate(alist):
        pair = [a0,a1]
        pair.sort()
        if pair in a0_pairs_ok:
            if debug:
                print("\nSkipping pair {}".format(pair))
            continue
        if debug:
            print("\nTesting intersection points of {}".format(pair))
        ok = are_intersection_points_covered(a0, a1, alist, sigmas, debug)
        if ok:
            new_pairs_ok.append(pair)
            if debug:
                print(" - ok: intersection points of {} and {} are covered".format(a0,a1))
        else:
            if debug:
                print(" - not ok: intersection points of {} and {} are not covered".format(a0,a1))
            all_ok = False
    if debug:
        if all_ok:
            print("OK: all intersection points of {} are covered".format(a0))
        else:
            print("No: not all intersection points of {} are covered".format(a0))
    return all_ok, new_pairs_ok

def are_alphas_surrounded(alist_ok, alist_open, slist, pairs_ok=[],
                          verbose=False, debug=False):
    """Given alist_ok and alist_open, lists of principal cusps, and slist,
    a complete list of singular points, tests whether the boundary of
    every disc S_a for a in alist_open is contained in the union of
    the translates of the S_b for b in alist_ok+alist_open, apart from
    any singular points on the boundary.

    Any a which pass are added to a new copy of alist_ok, while any
    which fail are added to a new alist_open, so success means that
    the latter is empty.  This allows for incremental testing by
    adding more a to alist_open.

    pairs_ok is list of pairs (a1,a2) whose intersection points are
    known to be covered.

    Returns (True/False, new_alist_ok, new_alist_open, new_pairs_ok).

    NB All a in alist_open will be tested, i.e. we carry on after a
    failure.

    """
    alist = alist_ok + alist_open
    alist = sum([translates(b) for b in alist], [])
    new_alist_ok = alist_ok.copy()
    new_alist_open = []
    all_ok = True
    for i, a in enumerate(alist_open):
        if not all_ok: # then return False so don't check remaining alphas
            new_alist_open.append(a)
            continue
        if cusp_in_quarter_rectangle(a):
            if verbose or debug:
                print(f"Testing alpha #{i+1}/{len(alist_open)} = {a}", end="...")
            ok, new_pairs_ok = is_alpha_surrounded(a, alist, slist, pairs_ok, debug)
            pairs_ok = new_pairs_ok
            if ok:
                if verbose or debug:
                    print(" ok! surrounded")
            else:
                all_ok = False
                if verbose or debug:
                    print(" no, not surrounded")
        else:
            ok = True
        if ok:
            new_alist_ok.append(a)
        else:
            new_alist_open.append(a)
    return all_ok, new_alist_ok, new_alist_open, pairs_ok

def reduced_numerators(s):
    k = s.parent()
    one = ZZ(1)
    for r in k.ideal(s).residues():
        if k.ideal(r,s).norm() == one:
            yield r

def principal_cusps_iter(k, maxnorm_s):
    """
    Return iterator yielding all principal r/s with N(s)<=maxnorm_s
    """
    for s in elements_of_norm_upto(k, maxnorm_s):
        for r in reduced_numerators(s):
            yield cusp(reduce_mod_Ok(r/s), k)

def principal_cusps_of_norm(k, norm_s):
    """
    Return iterator yielding all principal r/s with N(s)=maxnorm_s
    """
    for s in elements_of_norm(k, norm_s):
        for r in reduced_numerators(s):
            yield cusp(reduce_mod_Ok(r/s), k)

def principal_cusps_up_to(k, maxn, fussy=True):
    """List of all principal r/s with N(s)<=maxn, omitting any whose
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
        sigmas = singular_points(k)
    ok = False
    maxn = 0
    alphas_ok = []
    alphas_open = []
    pairs_ok = []
    Alist = []
    first = True #False
    while not ok:
        nc = 0 # number of new alphas added to list
        if first:
            maxn = -k.discriminant()//4
            new_cusps = principal_cusps_up_to(k, maxn)
            first = False
        else:
            maxn = next_norm(k, maxn+1)
            new_cusps = principal_cusps_of_norm(k, maxn)
        for a in new_cusps: #principal_cusps_of_norm(k, maxn):
            A = cusp_to_point(a)
            if not any(circle_inside_circle(A, B, False) for B in Alist):
                if cusp_in_quarter_rectangle(a):
                    alphas_open.append(a)
                    nc += 1
                else:
                    alphas_ok.append(a)
                Alist.append(A)
        if verbose and nc:
            s = "up to " if first else ""
            print(f"Adding {nc} alphas of norm {s}{maxn} (plus symmetrics); #alphas={len(alphas_ok)+len(alphas_open)} of which {len(alphas_ok)} are proved surrounded so far")
        if nc==0:
            continue
        ok, new_alphas_ok, new_alphas_open, new_pairs_ok = are_alphas_surrounded(alphas_ok, alphas_open, sigmas, pairs_ok, verbose=verbose, debug=False)
        if ok:
            if verbose:
                print("Success using {} alphas of with max norm {}!".format(len(new_alphas_ok), maxn))
            return maxn, new_alphas_ok, sigmas
        else:
            alphas_ok = new_alphas_ok
            alphas_open = new_alphas_open
            pairs_ok = new_pairs_ok
            if verbose:
                print("Some alphas are not surrounded, continuing...")

def point_translates(P):
    w = P[0].parent().gen()
    return [[P[0]+a+b*w, P[1]] for a in [-1,0,1] for b in [-1,0,1]]

def nverts(a, plist):
    return sum([1 for P in plist if is_under(P,a)==0])

def saturate_covering_alphas(k, alphas, sigmas, maxn=1, debug=False, verbose=False):
    """Given a covering set of alphas as produced by
    find_covering_alphas(), add extras if necessary so that they are
    "saturated", i.e. define the extended fundamental domain.

    By Swan, we need to find the points P in H^3 with positive height
    where at least 3 hemispheres S_a intersect, and for each P check
    whether P is properly covered by an S_a for a not in the set of
    alphas (up to translation).  If so, we need to add a to the set of
    alphas.  If none, then we have the fundamental region (and can go
    on to discard any redundant alphas).

    We assume that we have already considered all alpha=r/s with N(s)<=maxn.

    At the end we discard any alphas with <3 vertices (including translates and singular points).

    NB We assume that the initial list of alphas contains all r/s with
    N(s) up to some bound:
    """
    sat = False
    checked_points = []
    # copy so original list unchanged
    alphas1 = alphas.copy()
    while not sat:
        all_points = triple_intersections(alphas1, debug=debug)
        if debug:
            print(f"Found {len(all_points)} potential vertices")
        points = [P for P in all_points if maxn*P[1]<1]
        if debug:
            print(f" -- of which {len(all_points)} are low enough to be properly covered by a new alpha")
        points = [P for P in all_points if in_quarter_rectangle(P[0])]
        if debug:
            print(f" -- of which {len(points)} lie in the first quadrant")
        points = [P for P in points if P not in checked_points]
        if debug:
            print(f" -- of which {len(points)} have not already been checked")
        sat = True        # will be set to False if we find out that the alphas are not already saturated
        extra_alphas = [] # will be filled with any extra alphas needed
        points.sort(key=lambda P: -P[1]) # highest first
        for P in points:
            # if not sat:
            #     if debug:
            #         print("Restarting...")
            #     break
            if debug:
                print(f" - checking {P = }", end="...")
            extras = properly_covering_hemispheres(P)
            if debug:
                print(" done", end="...")
            #extras = properly_covering_hemispheres(P, maxn+1)
            #extras = covering_hemispheres2(P, 'strict', norm_s_lb=maxn+1)
            if extras:
                sat = False
                hts = [radius_squared(a) - (P[0]-to_k(a)).norm() for a in extras]
                m = max(hts)
                extras0 = [a for a,h in zip(extras, hts) if h==m]
                norms0 = [a.denominator().norm() for a in extras0]
                if debug:
                    print(f"   - found {len(extras)} properly covering hemispheres, with norms {[a.denominator().norm() for a in extras]}")
                    print(f"     max height above P (height {P[1]}) is {m}, for {extras0} with norms {norms0}")
                for a in extras0:
                    ca = conj_cusp(a)
                    for b in [a, negate_cusp(a), ca, negate_cusp(ca)]:
                        if cusp_in_rectangle(b) and b not in extra_alphas:
                            if debug:
                                print(f" - adding alpha {b}")
                            extra_alphas.append(b)
            else:
                if debug:
                    print("   - OK, no properly covering alphas found")
                checked_points.append(P)
        if verbose:
            if sat:
                m = max([a.denominator().norm() for a in alphas1])
                print(f" alphas are saturated! {len(alphas1)} alphas with max norm {m}")
            else:
                m = max([a.denominator().norm() for a in extra_alphas])
                print(f" alphas not saturated, {len(extra_alphas)} extras needed: {extra_alphas} (with norms at most {m})")
        alphas1 += extra_alphas

    m = max([a.denominator().norm() for a in alphas1])
    if verbose:
        print(f"After saturation we now have {len(alphas1)} alphas with max norm {m}")

    # Now delete any alphas with <3 vertices, allowing for translates
    pointsx = []
    for P in all_points+[[to_k(s,k),0] for s in sigmas if not s.is_infinity()]:
        for Q in point_translates(P):
            if Q not in pointsx:
                pointsx.append(Q)
    nv = [nverts(a, pointsx) for a in alphas1]
    # if verbose:
    #     print(f"# vertices for these alphas: {nv}")
    alphas1 = [a for a,n in zip(alphas1,nv) if n>=3]
    m = max([a.denominator().norm() for a in alphas1])
    if verbose:
        print(f"After removing alphas which go through <3 vertices, we now have {len(alphas1)} alphas with max norm {m}")
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

def denom_2_alphas(k):
    d = -k.discriminant().squarefree_part() # for compatibility with C++'s d
    if d in [1, 2, 3, 7, 11]:
        return []
    w = k.gen()
    d8 = d%8
    alist = []
    if d8 in [1,3,5]:
        alist.append(cusp(w/2,k))
    if d8 in [2,6]:
        alist.append(cusp((1+w)/2,k))
    if d8 == 3:
        alist.append(cusp((w-1)/2,k))
    return alist

def denom_2_sigmas(k):
    d = -k.discriminant().squarefree_part() # for compatibility with C++'s d
    if d in [1, 2, 3, 7, 11, 19, 43, 67, 163]:
        return []
    w = k.gen()
    d8 = d%8
    slist = []
    if d8 in [1,5]:
        slist.append(cusp((1+w)/2,k))
    if d8 in [2,6]:
        slist.append(cusp(w/2,k))
    if d8 == 7:
        slist.append(cusp(w/2,k))
        slist.append(cusp((1-w)/2,k))
    return slist

def denom_3_alphas(k):
    d = -k.discriminant().squarefree_part() # for compatibility with C++'s d
    if d in [1, 2, 3, 7, 11, 5, 6, 15, 19, 23]:
        return []

    w = k.gen()
    d12 = d%12

    if d12==3:
        alist = [w, w-1]
    if d12==7:
        alist = [w, 1-w, 1+w] if d>31 else [1+w]
    if d12==11:
        alist = [1+w]
    if d12 in [1,10]:
        alist = [w, 1+w, 1-w]
    if d12 in [2,5]:
        alist = [w]
    if d12 in [6,9]:
        alist = [1+w, w-1]

    return sum([[cusp(a/3,k), cusp(-a/3,k)] for a in alist], [])

def denom_3_sigmas(k):
    d = -k.discriminant().squarefree_part() # for compatibility with C++'s d
    if d in [1, 2, 3, 7, 11, 19, 43, 67, 163, 5]:
        return []
    w = k.gen()
    d12 = d%12
    slist = []
    if d12 in [2, 5]:
        slist.append(cusp((1+w)/3,k))
        slist.append(cusp((-1-w)/3,k))
        slist.append(cusp((1-w)/3,k))
        slist.append(cusp((w-1)/3,k))
    if d12 == 11:
        # if d==35:
        #     slist.append(cusp(w/3,k))
        #     slist.append(cusp(-w/3,k))
        if d>=35:
            slist.append(cusp(w/3, k))
            slist.append(cusp(-w/3, k))
            slist.append(cusp((w-1)/3,k))
            slist.append(cusp((1-w)/3,k))
    if d12 == 3 and d>15:
        slist.append(cusp((1+w)/3,k))
        slist.append(cusp((-1-w)/3,k))
    if d12 in [6,9] and d>6:
        slist.append(cusp(w/3,k))
        slist.append(cusp(-w/3,k))
    return slist

def alpha_in_list(a, alist, up_to_translation=True):
    if up_to_translation:
        return alpha_index_with_translation(a, alist)[0]>=0
    else:
        return a in alist

def compare_alpha_lists(alist1, alist2):
    return len(alist1)==len(alist2) and all(alpha_in_list(a,alist2) for a in alist1) and all(alpha_in_list(a,alist1) for a in alist2)

def find_edge_pairs(kdata, alphas, sigmas, debug=False, geout=None):
    from utils import ispos, add_two_alphas, add_four_alphas, half

    k = kdata['k']
    w = kdata['w']
    d = kdata['d'] # for compatibility with C++'s d

    # Extract the a for which 2*a or 3*a is integral, which we treat
    # separately:
    A1 = [a for a in alphas if to_k(a,k).is_integral()]
    assert A1 == [cusp(0,k)]
    A2 = [a for a in alphas if a not in A1 and (2*to_k(a,k)).is_integral()]
    A2exp = denom_2_alphas(k)
    if not compare_alpha_lists(A2, A2exp):
        print("*******************denom 2 alphas are {}, expected {}".format(A2, A2exp))
        raise ValueError
    A2 = A2exp # use the expected list for consistent normalisation and order
    A12 = A1 + A2
    A3 = [a for a in alphas if not alpha_in_list(a, A12) and (3*to_k(a,k)).is_integral()]
    A3exp = denom_3_alphas(k)
    if not compare_alpha_lists(A3, A3exp):
        print("*******************denom 3 alphas are {}, expected {}".format(A3, A3exp))
        raise ValueError
    A3 = A3exp # use the expected list for consistent normalisation and order
    A123 = A12 + A3

    # For a such that 2*a, 3*a are not integral and we make sure that we
    # have complete sets of {a,-a} pairs, not just up to translation:

    A = []
    for a in alphas:
        da = a.denominator()
        if not ispos(da):
            a = cusp(to_k(a))
        ma = negate_cusp(a)
        if not alpha_in_list(a, A123) and not alpha_in_list(ma, A):
            r,i = to_k(a,k)
            if w.trace()==0:
                if r<0 and i==half:
                    a = cusp(k([r,-half]), k)
                elif i<0 and r==half:
                    a = cusp(k([-half,i]), k)
            else:
                if i==half:
                    if r>0:
                        a = cusp(k([r,i-1]), k)
                    elif r<-half:
                        a = cusp(k([r+1,i-1]), k)
                elif 2*r+i==1 and i<0:
                    a = cusp(k([r-1,i]), k)
            r,i = to_k(a,k)
            if i>0:
                A.append(a)
                A.append(negate_cusp(a))
            else:
                A.append(negate_cusp(a))
                A.append(a)

    S = list(set(k(a.denominator()) for a in A))
    S.sort(key = lambda z: z.norm())
    if debug:
        print("Denominator 1,2,3: {}".format(A123))
        print("Other denominators: {}".format(S))
        for s in S:
            print("s = {}: numerators {}".format(s, [a for a in A if a.denominator()==s]))

    new_alphas = []
    M_alphas = []
    pluspairs = []
    minuspairs = []
    fours = []
    long_fours = []

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
                    else:
                        if debug:
                            print("  - adding plus pair {}".format(mrs))
                        pluspairs.append(mrs)
                continue
            if cong_mod(rsq, -1, s):
                if not any(pair in minuspairs for pair in (rs, mrs)):
                    if ispos(r):
                        if debug:
                            print("  - adding minus pair {}".format(rs))
                        minuspairs.append(rs)
                    else:
                        if debug:
                            print("  - adding minus pair {}".format(mrs))
                        minuspairs.append(mrs)
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
                            print("  - adding foursome {}".format((r,s,rdash)))
                        fours.append(rs)
                        long_fours.append((s,r,rdash))
                    else:
                        if debug:
                            print("  - adding foursome {}".format((-r,s,rdash)))
                        fours.append(mrs)
                        long_fours.append((s,-r,-rdash))
            except StopIteration:
                print("no negative inverse found for {} mod {}".format(r, s))
                raise ValueError

    for r,s in pluspairs:
        add_two_alphas(s, r, +1, new_alphas, M_alphas)
    for r,s in minuspairs:
        add_two_alphas(s, r, -1, new_alphas, M_alphas)
    for s, r1, r2 in long_fours:
        add_four_alphas(s, r1, r2, new_alphas, M_alphas)

    # Process the sigmas, standardising those with denominator 2 and 3 and putting the rest into +/- pairs

    # Extract the s with denominator 2 or 3, which we treat
    # separately:
    S2 = [s for s in sigmas if (not s.is_infinity()) and (2*to_k(s,k)).is_integral()]
    S2exp = denom_2_sigmas(k)
    if not compare_alpha_lists(S2, S2exp):
        print("*******************denom 2 sigmas are {}, expected {}".format(S2, S2exp))
    S2 = S2exp # use the expected list for consistent normalisation and order
    S3 = [s for s in sigmas if (not s.is_infinity()) and (not alpha_in_list(s, S2)) and (3*to_k(s,k)).is_integral()]
    S3exp = denom_3_sigmas(k)
    if not compare_alpha_lists(S3, S3exp):
        print("*******************denom 3 sigmas are {}, expected {}".format(S3, S3exp))
    S3 = S3exp # use the expected list for consistent normalisation and order
    S23 = S2 + S3

    S = []
    S_mod_neg = []
    for s in sigmas:
        ms = negate_cusp(s)
        if not s.is_infinity() and not alpha_in_list(s, S23) and not alpha_in_list(ms, S):
            r,i = to_k(s,k)
            if w.trace()==0:
                if r<0 and i==half:
                    s = cusp(k([r,-half]), k)
                elif i<0 and r==half:
                    s = cusp(k([-half,i]), k)
            else:
                if i==half:
                    if r>0:
                        s = cusp(k([r,i-1]), k)
                    elif r<-half:
                        s = cusp(k([r+1,i-1]), k)
                elif 2*r+i==1 and i<0:
                    s = cusp(k([r-1,i]), k)
            r,i = to_k(s,k)
            neg_s = negate_cusp(s)
            if i>0:
                S_mod_neg.append(s)
                S.append(s)
                S.append(neg_s)
            else:
                S_mod_neg.append(neg_s)
                S.append(neg_s)
                S.append(s)

    if debug:
        print(f"alphas with denominator | 2: {A2}")
        print(f"alphas with denominator | 3: {A3}")
        print(f"plus pairs: {pluspairs}")
        print(f"minus pairs: {minuspairs}")
        print(f"fours: {fours}")
        print(f"sigmas with denominator 2: {S2}")
        print(f"sigmas with denominator 3: {S3}")
        print(f"other (finite) sigmas (up to sign): {S_mod_neg}")

    new_sigmas = [cusp(oo,k)] + S2 + S3 + S

    geodata_file = f"geodata_{d}.dat"
    # print("//////////////////////////////")
    if geout:
        print(f"// tessellation edge data will be output to {geodata_file}")
    st = f"0\n0 {d=}\n0"
    if geout:
        geout.write(st+"\n")
    # else:
    #     print(st)
    for r,s in pluspairs:
        sr, si = s
        r1r, r1i = r
        r2r, r2i = -r
        st = f"{d} A {sr} {si} {r1r} {r1i} {r2r} {r2i}"
        if geout:
            geout.write(st+"\n")
        # else:
        #     print(st)
    for r,s in minuspairs:
        sr, si = s
        r1r, r1i = r
        r2r, r2i = r
        st = f"{d} A {sr} {si} {r1r} {r1i} {r2r} {r2i}"
        if geout:
            geout.write(st+"\n")
        # else:
        #     print(st)
    for s, r1, r2 in long_fours:
        sr, si = s
        r1r, r1i = r1
        r2r, r2i = r2
        st = f"{d} A {sr} {si} {r1r} {r1i} {r2r} {r2i}"
        if geout:
            geout.write(st+"\n")
        # else:
        #     print(st)
    for s in S_mod_neg:
        sr, si = s.denominator()
        rr, ri = s.numerator()
        st = f"{d} S {rr} {ri} {sr} {si}"
        if geout:
            geout.write(st+"\n")
    #     else:
    #         print(st)
    # print("//////////////////////////////")

    # for homology edge relation computation include alphas with denom 1,2,3:
    zero = k(0)
    one = k(1)
    two = k(2)
    three = k(3)
    minuspairs.append((zero,one))

    if d%4 == 1:
        minuspairs.append((w,two))     # w^2=-1 (mod 2)
        # could also go into pluspairs
    if d%4 == 2:
        minuspairs.append((w+1,two))   # (w+1)^2=-1 (mod 2)
        # could also go into pluspairs
    if d%8 == 3:
        long_fours.append((two,w,w+1)) # w(w+1)=-1 (mod 2)

    d12 = d%12
    if d12 in [1, 10]:
        minuspairs.append((w, three))          # w^2=-1 (mod 3)
        long_fours.append((three, 1+w, 1-w))   # (1+w)(1-w)=-1 (mod 3)
    if d12 == 7 and d>19:
        minuspairs.append((1+w, three))        # (1+w)^2=-1 (mod 3)
        if d>31:
            long_fours.append((three, w, 1-w)) # w(1-w)=-1 (mod 3)
    if d12 in [2, 5] and d>5:
        pluspairs.append((w, three))           # w^2=+1 (mod 3)
    if d12 == 11 and d>23:
        pluspairs.append((1+w, three))         # (1+w)^2=+1 (mod 3)
    if d12 == 3 and d>15:
        long_fours.append((three, w, w-1))     # w(w-1)=-1 (mod 3)
    if d12 in [6, 9] and d>6:
        long_fours.append((three, w+1, w-1))   # (w+1)(w-1)=-1 (mod 3)

    assert all(s.divides(r*r-1) for r,s in pluspairs)
    assert all(s.divides(r*r+1) for r,s in minuspairs)
    assert all(s.divides(r1*r2+1) for s,r1,r2 in long_fours)

    return A123, new_alphas, new_sigmas, pluspairs, minuspairs, long_fours

#
# From scratch:
#
def alpha_sigma_data(kdata, verbose=False, geout=None):
    k = kdata['k']
    d = kdata['d']
    if verbose:
        print(f"{k = }, class number {k.class_number()}")
    sigmas = singular_points(k)
    if verbose:
        print(f"{len(sigmas)} singular points: {sigmas}")
    maxn, alphas0, sigmas = find_covering_alphas(k, sigmas, verbose=verbose)
    if verbose:
        print(f"{len(alphas0)} covering alphas, max denom norm {maxn}")
    alphas1, points = saturate_covering_alphas(k, alphas0, sigmas, maxn, debug=verbose, verbose=verbose)
    maxn = max(a.denominator().norm() for a in alphas1)
    if verbose:
        print(f"{len(alphas1)} fundamental domain alphas, max denom norm {maxn}")
        print(f"{len(points)} fundamental vertices, min square height = {min(P[1] for P in points)}")
    # A2, new_alphas, M_alphas, pluspairs, minuspairs, long_fours
    data = find_edge_pairs(kdata, alphas1, sigmas, geout=geout)
    alphas2 = data[0] + data[1]
    new_sigmas = data[2]
    # for adding to precomputed alphas in alphas.py:
    alpha_string = "alphalist[{}] = [".format(d) + ", ".join([f"({a.numerator()})/({a.denominator()})" for a in alphas2]) + "]\n"
    alpha_string = alpha_string.replace(" ", "").replace('w','t').replace(",(",", (").replace("="," = ")
    alpha_string = alpha_string.replace("(0)/(1)", "0")
    alpha_file = f"alphas_{d}.py"
    print(f"Writing to {alpha_file} for inserting into alphas.py")
    with open(alpha_file, 'w') as aout:
        aout.write(alpha_string+"\n")
    if False: #verbose:
        print(alpha_string)
        sigma_string = "sigmas: [" + ", ".join([f"({s.numerator()})/({s.denominator()})" for s in new_sigmas]) + "]\n"
        print(sigma_string)
    plus_pairs = data[3]
    minus_pairs = data[4]
    fours = data[5]
    return alphas2, new_sigmas, plus_pairs, minus_pairs, fours

def edge_boundary_vector(e):
    v = vector(ZZ, e[0].number_field().class_number())
    v[cusp_class_index(e[1])] +=1
    v[cusp_class_index(e[0])] -=1
    return v

def edge_boundary_matrix(kdata, alphas, sigmas):
    """Return matrix of the boundary map from edges to cusps (for H^3
    modulo GL(2,O_k) or SL(2,O_k)).  The edge basis consists of first
    the [a,oo] for a in alphas (whose boundary is trivial), then the
    [s,oo] for s in sigmas[1:].  The cusp basis is indexed by ideal
    classes.
    """
    n = len(alphas)
    M = Matrix(ZZ, kdata['h'], n+len(sigmas)-1)
    for j,s in enumerate(sigmas):
        if j:
            M[0,n+j-1] = -1
            M[cusp_class_index(s), n+j-1] = 1
    return M

def edge_index(e, alphas, sigmas):
    """For e=[a1,a2] with a1,a2 principal return i where [oo,alphas[i]]
    is SL(2,Ok)-equivalent to e.

    For e=[a,s] with a principal and s singular, return i+len(alphas)-1,
    where [oo,sigmas[i]] is SL(2,Ok)-equivalent to e. Here i>=1.

    For e=[s,a] with a principal and s singular, return
    -(i+len(alphas)-1)<0, where [sigmas[i],oo] is SL(2,Ok)-equivalent to
    e.

    Raise an error if both e[0] amd e[1] are singular.

    Note that we assume sigmas[0]=oo which is not a singular point;
    the number of singular points is len(sigmas)-1.

    """
    a1, a2 = e
    if a1.ideal().is_principal(): # cases {a1,a2}, {a,s}
        U = Matrix(2,2, a1.ABmatrix())
        a1, a2 = apply(U.inverse(), e)
        assert a1.is_infinity()
        if a2.ideal().is_principal():
            i, x = alpha_index_with_translation(a2, alphas)
            assert i>=0
            return i
        else:
            i, x = sigma_index_with_translation(a2, sigmas)
            assert i>0
            return i+len(alphas)-1
    if a2.ideal().is_principal(): # case {s,a} = -{a,s}
        return - edge_index([e[1],e[0]], alphas, sigmas)
    # case {s1,s2} not required
    raise RuntimeError(f"edge {e} has neither end principal")

def face_boundary_vector(F, alphas, sigmas):
    v = vector(ZZ, len(alphas)+len(sigmas)-1)
    for i in range(len(F)):
        j = edge_index((F[i-1],F[i]), alphas, sigmas)
        #print(f" index of edge {(F[i-1],F[i])} is {j}")
        if j>=0:
            v[j] +=1
        else:
            v[-j] -=1
    return v

def face_boundary_matrix(faces, alphas, sigmas, plus_pairs, minus_pairs, fours, group, debug=False):
    nrows = len(alphas)+len(sigmas)-1
    ncols = len(faces) + len(plus_pairs) + len(minus_pairs) + len(fours)
    if group=="GL2":
        ncols += nrows # {a,oo}={-a,oo} and {s,oo}={-s,oo}
    else:
        ncols += len(faces) # negating all faces
        ncols += len(fours) # 2nd rel per four
    if debug:
        print(f"Face boundary matrix has {nrows=} and {ncols=}")
    M = Matrix(ZZ, nrows, ncols)
    n = 0 # = number of cols filled so far

    if debug:
        print("Edge identifications:")
    # edge relations
    # (0) relations {a,oo}={-a,oo} and {s,oo}={-s,oo} if GL2
    if group=="GL2":
        if debug:
            print(" alpha negations:")
        for i,a in enumerate(alphas):
            ma = negate_cusp(a)
            j = alpha_index_with_translation(ma, alphas)[0]
            M[i, n+i] +=1
            M[j, n+i] -=1
            if debug:
                print(f"alpha={a} --> column {n+i}: {M.column(n+i)}")
        n += len(alphas)
        if debug:
            print(" sigma negations:")
        m = len(alphas) ## offset into rows
        for i,s in enumerate(sigmas[1:]):
            ms = negate_cusp(s)
            j = sigma_index_with_translation(ms, sigmas)[0] -1 # omit sigmas[0]=oo
            M[m+i, n+i] +=1
            M[m+j, n+i] -=1
            if debug:
                print(f"sigma={s} --> column {n+i}: {M.column(n+i)}")
        n += len(sigmas)-1
    if debug:
        print(f" - after edge identifications (0), filled {n} columns out of {ncols} and rank is {M.rank()}")

    k = alphas[0].number_field()

    # (1) relations {a,oo}+{-a,oo}=0 for a=r/s, (r,s) in plus_pairs
    if debug:
        print(" + pairs:")
    for i, (r1,s) in enumerate(plus_pairs):
        for r in (r1,-r1):
            a = cusp(r/s, k)
            j = alpha_index_with_translation(a, alphas)[0]
            M[j, n+i] +=1
        if debug:
            print(f"Column {n+i}: {M.column(n+i)}")
    n += len(plus_pairs)
    if debug:
        print(f" - after edge identifications (1), filled {n} columns out of {ncols} and rank is {M.rank()}")

    # (2) relations 2{a,oo}=0 for a=r/s, (r,s) in minus_pairs
    if debug:
        print(" - pairs:")
    for i, (r,s) in enumerate(minus_pairs):
        a = cusp(r/s, k)
        j = alpha_index_with_translation(a, alphas)[0]
        M[j, n+i] +=2
        if debug:
            print(f"Column {n+i}: {M.column(n+i)}")
    n += len(minus_pairs)
    if debug:
        print(f" - after edge identifications (2), filled {n} columns out of {ncols} and rank is {M.rank()}")

    # (3) relations {a1,oo}+{a2,oo}=0 for a1=r1/s, a2=r2/s, (s,r1,r2) in fours
    if debug:
        print(" fours (1):")
    for i, (s,r1,r2) in enumerate(fours):
        for r in (r1,r2):
            a = cusp(r/s, k)
            j = alpha_index_with_translation(a, alphas)[0]
            M[j, n+i] +=1
        if debug:
            print(f"Column {n+i}: {M.column(n+i)}")
    n += len(fours)
    if group!="GL2":
        if debug:
            print(" fours (2):")
        for i, (s,r1,r2) in enumerate(fours):
            for r in (r1,r2):
                a = cusp(-r/s, k)
                j = alpha_index_with_translation(a, alphas)[0]
                M[j, n+i] +=1
            if debug:
                print(f"Column {n+i}: {M.column(n+i)}")
        n += len(fours)
    if debug:
        print(f" - after edge identifications (3), filled {n} columns out of {ncols} and rank is {M.rank()}")

    # face relations
    for i,F in enumerate(faces):
        M.set_column(n+i,face_boundary_vector(F, alphas, sigmas))
        if debug:
            #print(f"Face {i}: {F}")
            print(f"Column {n+i}: {M.column(n+i)}")
    n += len(faces) # = number of cols filled so far

    # the faces are up to GL2-equivalence, so if SL2 we need to also negate them
    if group!="GL2":
        for i,F in enumerate(faces):
            mF = [negate_cusp(a) for a in F]
            M.set_column(n+i,face_boundary_vector(mF, alphas, sigmas))
            if debug:
                #print(f"Negated face {i}: {mF}")
                print(f"Column {n+i}: {M.column(n+i)}")
        n += len(faces) # update number of cols filled so far
    if debug:
        print(f" - after face relations, filled {n} columns out of {ncols} and rank is {M.rank()}")

    assert n==M.ncols()
    return M

def compute_homology(kdata, alphas, sigmas, plus_pairs, minus_pairs, fours, faces, debug=False):
    M10 = edge_boundary_matrix(kdata, alphas, sigmas)
    #print(f"edge boundary matrix:\n{M10}")
    D,U,V = M10.smith_form(transformation=True)
    #print(f" - smith form:\n{D}")
    assert U*M10*V == D
    r = D.rank()
    #print(f" - rank: {r}")
    Vinv = V.inverse().change_ring(ZZ)
    hom = {}
    for group in ["GL2", "SL2"]:
        #print(f"{group} integral homology data for D={kdata['dk']}:")
        M21 = face_boundary_matrix(faces, alphas, sigmas, plus_pairs, minus_pairs, fours, group=group, debug=debug)
        #print(f"face boundary matrix:\n{M21}")
        assert M10*M21 == 0
        M = Vinv*M21
        assert D*M == 0
        #print(f" - after edge basis change:\n{M}")
        M = M.submatrix(r,0)
        #print(f" - after trimming:\n{M}")
        invariants = [d for d in M.elementary_divisors() if d!=1]
        #print(f" - invariants: {invariants}")
        rank = invariants.count(0)
        #print(f" - rank: {rank}")
        torsion_invariants = [d for d in invariants if d>1]
        invs = Factorization((p,1) for p in torsion_invariants)
        #print(f" - torsion invariants: {invs}")
        hom[group] = {'invariants': invariants,
                      'rank': rank,
                      'torsion_invariants': torsion_invariants,
                      'torsion_factors': invs}
    return hom

def tessellation(d, verbose=0, plot2D=False, plot3D=False, browser="/usr/bin/firefox"):
    from utils import (make_M_alphas,
                       make_poly_from_edges,
                       poly_gl2_orbit_reps,
                       aas_triangle_gl2_orbit_reps,
                       oriented_faces, face_index,
                       polyhedron_relation,
                       polygon_parameters,
                       is_poly_principal)

    kdata = make_k(d)
    k = kdata['k']
    if verbose:
        print("Field: {}".format(k))
        print("Discriminant: {}".format(k.discriminant()))
        print("Class number: {}".format(k.class_number()))

    geodata_file = f"geodata_{d}.dat"
    alphas = precomputed_alphas(d)
    if alphas:
        if verbose:
            print("using precomputed alphas")
        sigmas = singular_points(k)
        with open(geodata_file, 'w') as geout:
            data = find_edge_pairs(kdata, alphas, sigmas, geout=geout)
        alphas = data[0] + data[1]
        sigmas = data[2]
        plus_pairs = data[3]
        minus_pairs = data[4]
        fours = data[5]
    else:
        if verbose:
            print("computing alphas from scratch")
        with open(geodata_file, 'w') as geout:
            alphas, sigmas, plus_pairs, minus_pairs, fours = alpha_sigma_data(kdata, verbose=verbose, geout=geout)

    M_alphas, alpha_inv = make_M_alphas(alphas)
    if verbose:
        print("{} alphas".format(len(alphas)))
        print("{} sigmas".format(len(sigmas)))

    if plot2D:
        print("plotting projection of fundamental domain")
        show(plot_FunDomain_projection(k, alphas, sigmas))


    polyhedra, hemis = all_polyhedra(k, alphas, verbose>1)
    npoly = len(polyhedra)
    poly = "polyhedra" if npoly>1 else "polyhedron"
    print(f"{npoly} {poly} constructed")
    if plot3D:
        print("plotting fundamental domain")
        from sage.misc.viewer import viewer
        viewer.browser(browser)
        show(plot_Bianchi_diagram(kdata,hemis))

    pt = poly_types(polyhedra)
    nunk = pt['unknown']
    if nunk:
        poly = "polyhedra have" if nunk>1 else "polyhedron has"
        print(f"{nunk} {poly} unknown type!")
        if verbose:
            for G in polyhedra:
                if poly_type(G) == 'unknown':
                    show(G)

    for pol,num in pt.items():
        if num:
            print(f" {pol}: {num}")

    all_faces = sum([[make_poly_from_edges(f,k) for f in oriented_faces(G)] for G in polyhedra], [])
    triangles = [f for f in all_faces if len(f)==3]
    squares = [f for f in all_faces if len(f)==4]
    hexagons = [f for f in all_faces if len(f)==6]
    aaa_triangles = [T for T in triangles if is_poly_principal(T)]
    aas_triangles = [T for T in triangles if not is_poly_principal(T)]
    faces = aaa_triangles + aas_triangles + squares + hexagons

    print(f"All {len(faces)} polyhedron faces:")
    print(f" {len(aaa_triangles)} aaa-triangles,")
    print(f" {len(aas_triangles)} aas-triangles,")
    print(f" {len(squares)} squares,")
    print(f" {len(hexagons)} hexagons")

    if verbose:
        print()
        print("Finding GL2-orbits of faces...")

    aaa_triangles0 = poly_gl2_orbit_reps(aaa_triangles)
    aas_triangles0 = aas_triangle_gl2_orbit_reps(aas_triangles, alphas)
    squares0 = poly_gl2_orbit_reps(squares)
    hexagons0 = poly_gl2_orbit_reps(hexagons)
    faces = aaa_triangles0 + aas_triangles0 + squares0 + hexagons0

    print(f"{len(faces)} GL2-orbits of faces:")
    print(f" {len(aaa_triangles0)} aaa-triangles,")
    print(f" {len(aas_triangles0)} aas-triangles,")
    print(f" {len(squares0)} squares,")
    print(f" {len(hexagons0)} hexagons")

    if verbose:
        print()
        print("Finding redundant faces modulo polyhedron relations...")
    all_faces = hexagons0 + squares0 + aaa_triangles0 + aas_triangles0
    if False:
        print("List of GL2-inequivalent faces:")
        for F in all_faces:
            print(F)
        print(f"Face types: {[len(F) for F in all_faces]}")
        for P in polyhedra:
            print(f"Polyhedron {poly_type(P)} with oriented faces")
            for F in oriented_faces(P):
                print(f"{F = } with index {face_index(make_poly_from_edges(F,k), all_faces)}")
            print("and relation")
            print(polyhedron_relation(P, all_faces, k))
        M = Matrix([polyhedron_relation(P, all_faces, k) for P in polyhedra])
    redundant_faces = [] # [r.trailing_support() for r in H.rows() if r.trailing_coefficient()==1]

    # delete square faces of square prisms
    # delete hexagonal faces of hexagonal caps
    # delete aaa-traingles from aaas tetrahedra
    for poly in polyhedra:
        if poly_type(poly) == "square pyramid":
            squ = next(make_poly_from_edges(f,k) for f in oriented_faces(poly) if len(f)==4)
            ind, _ = face_index(squ, all_faces)
            if ind not in redundant_faces:
                if verbose:
                    print("omitting redundant square face from a square pyramid")
                redundant_faces.append(ind)
            continue
        if poly_type(poly) == "hexagonal cap":
            hexa = next(make_poly_from_edges(f,k) for f in oriented_faces(poly) if len(f)==6)
            ind, _ = face_index(hexa, all_faces)
            if ind not in redundant_faces:
                if verbose:
                    print("omitting redundant hexagonal face from a hexagonal cap")
                redundant_faces.append(ind)
            continue
        if poly_type(poly) == "tetrahedron":
            faces = [make_poly_from_edges(f,k) for f in oriented_faces(poly)]
            if not all(is_poly_principal(t) for t in faces):
                tri = next(t for t in faces if is_poly_principal(t))
                ind, _ = face_index(tri, all_faces)
                if ind not in redundant_faces:
                    if verbose:
                        print("omitting redundant aaa-triangle face from an aaas-tetrahedron")
                    redundant_faces.append(ind)
    if verbose:
        print(f"Redundant faces: {redundant_faces} ({len(redundant_faces)} out of {len(all_faces)})")
    i0 = 0 # index offset into list of all faces
    hexagons1 = [P for i,P in enumerate(hexagons0) if i+i0 not in redundant_faces]
    i0 += len(hexagons0)
    squares1 = [P for i,P in enumerate(squares0) if i+i0 not in redundant_faces]
    i0 += len(squares0)
    aaa_triangles1 = [P for i,P in enumerate(aaa_triangles0) if i+i0 not in redundant_faces]
    i0 += len(aaa_triangles0)
    aas_triangles1 = [P for i,P in enumerate(aas_triangles0) if i+i0 not in redundant_faces]
    faces = aaa_triangles1 + aas_triangles1 + squares1 + hexagons1

    print(f"{len(faces)} independent faces modulo polyhedron relations:")
    print(f" {len(aaa_triangles1)} aaa-triangles,")
    print(f" {len(aas_triangles1)} aas-triangles,")
    print(f" {len(squares1)} squares,")
    print(f" {len(hexagons1)} hexagons")

    print()

    geodata_file = f"geodata_{d}.dat"
    with open(geodata_file, 'a') as geout:
        for P in faces:
            polygon_parameters(P, alphas, M_alphas, alpha_inv, sigmas, geout=geout)

    with open(f"tessellation_{d}.txt", 'w') as tess_out:
        for pol,num in pt.items():
            if num:
                tess_out.write(f"{pol}: {num}\n")
        for P in polyhedra:
            from utils import cusp_from_string
            tess_out.write(poly_type(P)+"\n")
            V = [cusp_from_string(v, k) for v in P.vertices()]
            tess_out.write(f"  vertices: {V}\n")
            E = [[cusp_from_string(v, k) for v in e[:2]] for e in P.edges()]
            tess_out.write(f"  edges: {E}\n")
            F = [[[cusp_from_string(v, k) for v in e] for e in f] for f in P.faces()]
            tess_out.write(f"  faces: {F}\n")

    hom = compute_homology(kdata, alphas, sigmas, plus_pairs, minus_pairs, fours, faces)
    for group in ["GL2", "SL2"]:
        print(f"{group} integral homology data for D={kdata['dk']}:")
        print(f" - rank: {hom[group]['rank']}")
        print(f" - torsion invariants: {hom[group]['torsion_factors']}")

    return alphas, sigmas, faces, polyhedra, hom

def faces_from_file(kdata):
    from alphas import precomputed_alphas
    from utils import make_M_alphas, geodat_decoders, tri0, tri1, tri2

    d = kdata['d']
    k = kdata['k']
    w = kdata['w']
    alphas = precomputed_alphas(d)
    assert alphas
    sigmas = singular_points(k)
    alphas0, alphas1, sigmas, plus_pairs, minus_pairs, fours = find_edge_pairs(kdata, alphas, sigmas, geout=None)
    alphas = alphas0 + alphas1
    M_alphas, alpha_inv = make_M_alphas(alphas)

    # Faces for class number 1 fields treated separately in C++ code so not output to geodata:

    faces = [tri0(k)]          # triangle  [0,oo,1]
    if d in [19, 43, 67, 163]: # triangles [oo, w/2, (w-1)/2], [oo, w/2, (w+1)/2]
        faces.append(tri1(k))
        faces.append(tri2(k))
    if d in [2, 7]:
        faces.append([cusp(c,k) for c in [w/2,0,oo,w]])
    if d==11:
        faces.append([cusp(c,k) for c in [2*w/3,w/2,w/3,0,oo,w]])

    with open(f"../geodata/geodata_{d}.dat") as gin:
        for L in gin:
            cols = L.split()
            if int(cols[0])!=d:
                continue
            face_type = cols[1]
            params = [int(c) for c in cols[2:]]
            if face_type not in ['T', 'U', 'Q', 'H']:
                continue
            face = geodat_decoders[face_type](kdata, params, alphas, sigmas, M_alphas, alpha_inv)
            try:
                face_boundary_vector(face, alphas, sigmas)
                faces.append(face)
            except:
                print(f"Invalid face {face} from {L}")
    return alphas, sigmas, plus_pairs, minus_pairs, fours, faces

def integral_homology(d):
    kdata = make_k(d)
    alphas, sigmas, plus_pairs, minus_pairs, fours, faces = faces_from_file(kdata)
    hom = compute_homology(kdata, alphas, sigmas, plus_pairs, minus_pairs, fours, faces)
    for group in ["GL2", "SL2"]:
        print(f"{group} integral homology data for D={kdata['dk']}:")
        print(f" - invariants: {hom[group]['invariants']}")
        print(f" - rank: {hom[group]['rank']}")
        print(f" - torsion invariants: {hom[group]['torsion_factors']}")
