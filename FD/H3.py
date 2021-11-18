# Functions for working with H3 and hemispheres etc.

from itertools import chain #, combinations

from sage.all import (Infinity, Matrix, ZZ, QQ, CC, NumberField,
                      Graph, srange, Set, sign, var, implicit_plot3d, NFCusp, Integer, oo,
                      infinity, polygen)

from utils import (to_k, cusp, cusp2, cusp_label, Imat, apply,
                   translate_cusp, smallest_ideal_class_representatives)

# Precomputed alphas for some fields

def precomputed_alphas(k_or_d):

    if k_or_d in ZZ:
        # then d=k>0 should be square-free
        x = polygen(QQ)
        d = k_or_d
        if d%4==3:
            k = NumberField(x**2-x+(1+d)//4, 'w')
        else:
            k = NumberField(x**2+d, 'w')
    else:
        k = k_or_d
        d = -k.discriminant().squarefree_part()

    w = k.gen()
    IC = smallest_ideal_class_representatives(k)
    cusp = lambda a: NFCusp(k, k(a), lreps=IC)

    alphas_dict = {}
    alphas_dict[1] = [cusp(0)]
    alphas_dict[2] = [cusp(0)]
    alphas_dict[3] = [cusp(0)]
    alphas_dict[7] = [cusp(0)]
    alphas_dict[11] = [cusp(0)]
    alphas_dict[19] = [cusp(a) for a in [0, w/2,(w-1)/2]]

    alphas_dict[43] = [cusp(a) for a in [0, w/2, (-1+w)/2, w/3, -w/3,
                                         (1-w)/3, (-1+w)/3, (1+w)/3, (-1-w)/3]]

    alphas_dict[5] = [cusp(a) for a in [0, w/2, (w-4)/(2*w), (4-w)/(2*w)]]

    alphas_dict[6] = [cusp(a) for a in [0, (w+1)/2, 5/(2*w), -5/(2*w)]]

    alphas_dict[10] = [cusp(a) for a in [0, (w+1)/2, w/3, -w/3, (1+w)/3, (-1-w)/3, (1-w)/3,
                                         (-1+w)/3, 9/(2*w), -9/(2*w)]]

    alphas_dict[23] = [cusp(a) for a in [0, (1+2*w)/4, (-1-2*w)/4,
                                         (1+w)/3, (-1-w)/3, (2-w)/(1+w),
                                         (-2+w)/(1+w), (1+w)/(2-w),
                                         (-1-w)/(2-w), (-3+w)/(2+w),
                                         (3-w)/(2+w), (-2-w)/(3-w),
                                         (2+w)/(3-w)]]

    alphas_dict[31] = [cusp(a) for a in [0, (1+2*w)/4, (-1-2*w)/4, w/3,
                                         (-w)/3, (1-w)/3, (-1+w)/3,
                                         (1+w)/3, (-1-w)/3, 3/w, (-3)/w,
                                         3/(1-w), (-3)/(1-w), 3/(1+w),
                                         (-3)/(1+w), 3/(2-w), -3/(2-w),
                                         (-6+w)/(3+w), (6-w)/(3+w),
                                         (5+w)/(4-w), (-5-w)/(4-w)]]

    alphas_dict[67] = [cusp(a) for a in [0, w/2, (-1+w)/2, w/3, (-w)/3,
                                         (1-w)/3, (-1+w)/3, (1+w)/3, (-1-w)/3, w/4, (-w)/4, (-1+w)/4, (1-w)/4,
                                         (1+w)/4, (-1-w)/4, (2-w)/4, (-2+w)/4, (6+w)/(3-w), (-6-w)/(3-w),
                                         (2+w)/(3-w), (-2-w)/(3-w), (7-w)/(2+w), (-7+w)/(2+w), (3-w)/(2+w),
                                         (-3+w)/(2+w)]]

    alphas_dict[163] = [cusp(a) for a in [0, (w)/2, (-1+w)/2, (w)/3,
                                          (-w)/3, (1-w)/3, (-1+w)/3, (1+w)/3, (-1-w)/3, (w)/4, (-w)/4, (-1+w)/4,
                                          (1-w)/4, (1+w)/4, (-1-w)/4, (2-w)/4, (-2+w)/4, (w)/5, (-w)/5,
                                          (-1+w)/5, (1-w)/5, (2*w)/5, (-2*w)/5, (2-2*w)/5, (-2+2*w)/5, (2+w)/5,
                                          (-2-w)/5, (1-2*w)/5, (-1+2*w)/5, (-2+w)/5, (2-w)/5, (2+2*w)/5, (-2-2*w)/5,
                                          (1+w)/5, (-1-w)/5, (1+2*w)/5, (-1-2*w)/5, (w)/6, (-w)/6, (1-w)/6,
                                          (-1+w)/6, (1+w)/6, (-1-w)/6, (-2+w)/6, (2-w)/6, (2+w)/6, (-2-w)/6,
                                          (3-w)/6, (-3+w)/6, (12)/(w), (-12)/(w), (17)/(w), (-17)/(w),
                                          (12)/(1-w), (-12)/(1-w), (17)/(1-w), (-17)/(1-w), (12)/(1+w),
                                          (-12)/(1+w), (-17+w)/(1+w), (17-w)/(1+w), (12)/(2-w), (-12)/(2-w),
                                          (-16-w)/(2-w), (16+w)/(2-w), 7/(2+w), (-7)/(2+w), (18-w)/(2+w),
                                          (-18+w)/(2+w), (-11+w)/(2+w), (11-w)/(2+w), (-16+w)/(2+w),
                                          (16-w)/(2+w), 7/(3-w), (-7)/(3-w), (17+w)/(3-w), (-17-w)/(3-w),
                                          (-10-w)/(3-w), (10+w)/(3-w), (-15-w)/(3-w), (15+w)/(3-w), (3+w)/7,
                                          (-3-w)/7, (-1+2*w)/7, (1-2*w)/7, (1+2*w)/7, (-1-2*w)/7, (3-2*w)/7,
                                          (-3+2*w)/7, (2+3*w)/7, (-2-3*w)/7, (5-w)/(3+w), (-5+w)/(3+w),
                                          (-17+w)/(3+w), (17-w)/(3+w), (4+w)/(4-w), (-4-w)/(4-w), (-16-w)/(4-w),
                                          (16+w)/(4-w)]]

    if d in alphas_dict:
        return alphas_dict[d]

    print("No precomputed alphas for field {}".format(k))
    return None

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

def covering_hemispheres(P, option=None):
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
    -1 if they are internally tangent
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

def circle_inside_circle(P1,P2):
    return P1[1]<P2[1] and tau(P1,P2)==-2

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

def in_first_half(alpha1, alpha2, centre=0):
    """
    Return True if the clockwise angle from alpha1 round to alpha2 is < pi.
    """
    e1, s1 = slope(alpha1, centre)
    e2, s2 = slope(alpha2, centre)
    if s1==Infinity:
        return e1*e2<0 and s2!=Infinity
    if s2==Infinity:
        return e1*e2>0
    return e1*e2*(s1-s2)<0

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

def circ(c,r):
    return circle(c, r,
                  aspect_ratio=1, edgecolor='blue', alpha = 0.2)

def disc(c,r):
    return circle(c, r,
                  aspect_ratio=1, feill=True, rgbcolor='blue', alpha = 0.2)

def plot_FunDomain_projection(k, alphas, sigmas):
    w = k.gen()
    D = k.discriminant().abs()
    emb = next(e for e in k.embeddings(CC) if e(w).imag()>0)
    rootd = emb(w).imag()
    Ymax = 3*rootd/(4 if ZZ(D).mod(4) == 3  else 2)
    Xmax = 3*0.5

    triplets, extra_alphas = alpha_triples(alphas)

#    A = [xy_coords(to_k(a,k)) for a in alphas]
#    A = [(xy[0], xy[1]*rootd) for xy in A]
    A = [list(emb(to_k(a))) for a in alphas+extra_alphas]
    R = [RR(radius_squared(a)).sqrt() for a in alphas+extra_alphas]
    #print("circle centres: {}".format(A))

#    S = [xy_coords(to_k(s,k)) for s in sigmas]
#    S = [(xy[0], xy[1]*rootd) for xy in S]
    S = [list(emb(to_k(s))) for s in sigmas if not s.is_infinity()]
    #print("singular points: {}".format(S))

    C = [list(emb(P[2][0])) for P in triplets]
    #print(triplets)
    #print("corners: {}".format(C))

    circles  = [circ(centre, radius) for centre, radius in zip(A,R)]
    centres  = [point(P, rgbcolor = 'blue', pointsize=22) for P in A]
    sing_pts = [point(P, rgbcolor = 'red', pointsize=50) for P in S]
    corners  = [point(P, rgbcolor = 'yellow', pointsize=22) for P in C]

    z = w-ZZ(1)/2 if ZZ(D).mod(4)==3 else w
    TL=list(emb((-1+z)/2))
    TR=list(emb((1+z)/2))
    BR=list(emb((1-z)/2))
    BL=list(emb((-1-z)/2))
    lines = [line([TL,TR], rgbcolor='black'),
             line([TR,BR], rgbcolor='black'),
             line([BR,BL], rgbcolor='black'),
             line([BL,TL], rgbcolor='black')]

    proj = sum(circles) + sum(centres) + sum(sing_pts) + sum(corners) + sum(lines)
    proj.set_axes_range(-Xmax, Xmax, -Ymax, Ymax)
    return proj

def is_redundant(P, alphas):
    """Return True iff P is strictly covered by any of the hemispheres
    S_a for a in the list alphas.
    """
    return any(is_under(P,a)==1 for a in alphas)


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
    xalphas = sum([[translate_cusp(a,t) for t in
                    [a+b*w for a in [-1,0,1] for b in [-1,0,1]]] for a in alphas],
                  [])
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
        for k in range(max(i,j)+1, n):
            if i!=k and j!=k and {i,k} in ij_list and {j,k} in ij_list:
                ak = xalphas[k]
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
    return -half<x and x<= half and -half<f*y and f*y<=half

def cusp_in_rectangle(a):
    return in_rectangle(to_k(a))

def in_rectangle(a):
    f = 1 + a.parent().disc()%2
    return xy_in_rectangle(xy_coords(a), f)

def is_alpha_surrounded(a, alist, slist, debug=False):
    """Given alist, a candidate list of principal cusps, and slist, a
    complete list of singular points (excluding oo), tests whether the
    disc S_a is contained in the union of the translates of the S_b
    for b in alist, apart from any singular points on the boundary.

    Returns either (True, xlist) with xlist a list of any translates
    needed, or (False, None)
    """
    k = alist[0].number_field()
    w = k.gen()
    xlist = set()

    # convert a to a point with radius
    A = cusp_to_point(a)
    if debug:
        print("A = {}".format(A))
    # extract the relevant singular points, if any:
    Slist = [to_k(s) for s in slist]
    Slist = [s for s in Slist if is_under([s,0], a)==0]
    if debug:
        print("Slist = {}".format(Slist))
    # extend the candidate list by including offsets:
    offsets = [a+b*w for a in [-1,0,1] for b in [-1,0,1]]
    alistx = sum([[translate_cusp(b,t) for t in offsets] for b in alist], [])
    # extract the relevant alphas:
    Alist = [cusp_to_point(a) for a in alistx]
    Alist = [a for a in Alist if circles_intersect(A, a)]
    # sort these by slope:
    Alist.sort(key=lambda a: slope(a[0], A[0]))
    if debug:
        print("Alist (sorted) = {}".format(Alist))
        print("Relative slopes: {}".format([slope(a[0],A[0]) for a in Alist]))

    def test(A1, A2, detail=False):
        if not in_first_half(A1[0], A2[0], A[0]):
            if detail:
                print("test 1 fails (not close together)")
            return False, None, None
        if not circles_intersect(A1, A2):
            if detail:
                print("test 2 fails (circles do not intersect)")
            return False, None, None
        t = bi_inter(cusp(A1[0]), cusp(A2[0]))
        if t is None:
            if detail:
                print("test 3 fails (intersection does not cover boundary)")
            return False, None, None
        c, r2 = t
        if detail:
            print("c = {}, r2 = {}".format(c, r2))
        x,y = xy_coords(c)
        if -half<x and x<= half and -half<y and y<=half:
            x1 = A1[0]
            x2 = A2[0]
        else:
            x1 = x2 = None
        if r2==0:
            return (c in Slist), x1, x2
        else:
            return is_under(t, a), x1, x2

    for i, b in enumerate(Alist):
        flag, a1, a2 = test(Alist[i-1], b, debug)
        if flag:
            if a1 and cusp(a1) not in alist:
                xlist.add(a1)
            if a2 and cusp(a2) not in alist:
                xlist.add(a2)
        else:
            if debug:
                print("test fails for {} and {}".format(Alist[i-1], Alist[i]))
                print("slopes {} and {}".format(slope(Alist[i-1][0], A[0]),slope(Alist[i][0], A[0])))
                test(Alist[i-1], b, True)
            return False, None
    return True, xlist


def are_alphas_surrounded(alist, slist, debug=False):
    """Given alist, a candidate list of principal cusps, and slist, a
    complete list of singular points (excluding oo), tests whether every
    disc S_a is contained in the union of the translates of the S_b
    for b in alist, apart from any singular points on the boundary.

    Returns either (True, xlist) with xlist a list of any translates
    needed, or (False, None)
    """
    xlist = set()
    for a in alist:
        flag, xlist1 = is_alpha_surrounded(a, alist, slist, debug)
        if not flag:
            return False, None
        if debug:
            print("{} is ok".format(a))
        xlist = xlist.union(xlist1)
    xlist = [cusp(x) for x in xlist]
    return True, xlist

def elements_of_norm(k, n):
    return iter(k.elements_of_norm(n))

def elements_of_norm_upto(k, n, start=1):
    return chain(*(iter(k.elements_of_norm(n)) for n in range(start, n+1)))

def reduced_numerators(s):
    K = s.parent()
    one = K.ideal(1)
    for r in K.ideal(s).residues():
        if K.ideal(r,s)==one:
            yield r

def principal_cusps_iter(k, maxnorm_s):
    for s in elements_of_norm_upto(k, maxnorm_s):
        for r in reduced_numerators(s):
            yield cusp(r/s)

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
    return [cusp(a) for a in alphas]

def cong_mod(r1, r2, s):
    return ((r1-r2)/s).is_integral()

def find_edge_pairs(alphas, debug=False):
    from utils import add_two_alphas, add_four_alphas, nf

    k = nf(alphas[0])
    Ireps = smallest_ideal_class_representatives(k)
    A1 = [a for a in alphas if a.denominator()==1]
    A2 = [a for a in alphas if a.denominator()==2]
    A3 = [a for a in alphas if a.denominator()==3]
    A = [a for a in alphas if a.denominator() not in [1,2,3]]
    S = list(set(k(a.denominator()) for a in A))
    S.sort(key = lambda z: z.norm())

    new_alphas = []
    M_alphas = []
    pluspairs = []
    minuspairs = []
    fours = []

    for s in S:
        if not ispos(s):
            s = -s
            sgn = -1
        else:
            sgn = +1
        if debug:
            print("s = {}".format(s))
        As = [a for a in A if makepos(a.denominator())==s]
        for a in As:
            if debug:
                print("  a = {}".format(a))
            r = sgn*k(a.numerator())
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
    return A1, A2, A3, new_alphas, M_alphas, pluspairs, minuspairs, fours

