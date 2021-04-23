
def fundamental_domain(k):
    """
    For an imaginary quadratic field K, the function returns
    two lists, S and V. S is the list of hemispheres that form
    the floor of a fundamental domain. V is the list of vertices
    coming from the intersection of 3 or more hemispheres.
    """ 
    S = initial_list(k)
    V, S = clean_list(k, S)
    V, S = Swan(k, S, V)
    return S, V

#----- Plotting functions --------------------------

def plot_FunDomain_projection(k, S, V, Full=True):
    sk = singular_points_in_F(k)
    if Full == True:
        BS = []
        BV = []
        for s in S:
            BS.append(s)
            if s[0][0]>0:
                new_s = (k((-s[0][0],s[0][1])), s[1])
                BS.append(new_s)
        for v in V + [(s[0], s[1], 0) for s in sk]:
            BV.append(v)
            if v[0]>0:
                new_v = (-v[0], v[1], v[2])
                BV.append(new_v)
    else:
        BS = S
        BV = V + [(s[0], s[1], 0) for s in sk]
    L = []
    D = k.absolute_generator().norm()
    w = RR(D).sqrt()
    for s in BS:
        centre = (s[0][0], s[0][1]*w)
        rad = 1/RR(s[1]).sqrt()
        L.append(circle(centre, rad, aspect_ratio=1, \
                 fill=True, alpha = 0.2, rgbcolor = 'blue'))
    for v in BV:
        P = (v[0], v[1]*w)
        L.append(point(P, rgbcolor = 'red', pointsize=22))
    proj = sum(L)
    if ZZ(D).mod(4) == 3:
        b = 4
    else:
        b = 2
    if Full == True:
        proj.set_axes_range(-0.5, 0.5, 0, w/b)
    else:
        proj.set_axes_range(0, 0.5, 0, w/b)
    return proj

def plot_Bianchi_diagram(k, S):
    BS = []
    for s in S:
        BS.append(s)
        if s[0][0]>0:
            new_s = (k((-s[0][0],s[0][1])), s[1])
            BS.append(new_s)
    X, Y, Z = var('X, Y, Z')
    L = []
    D = k.absolute_generator().norm()
    w = RR(D).sqrt()
    b = 2
    if ZZ(D).mod(4)==3:
        b = 4
    Yrange = ceil(w/b)
    for s in BS:
        c = k(s[0])
        eq = (X - c[0])^2 + (Y - c[1]*w)^2 + Z^2 - 1/s[1]
        L.append(implicit_plot3d(eq, (Y, 0, Yrange ),  (X, -0.5, 0.5), \
                           (Z, 0, 1), plot_points=60, aspect_ratio=1)) 
    return sum(L)

#---------------------------------------------------------------------
#----- Finding an initial list of hemispheres ------------------------

#----- Auxiliar functions -----------------------------
@cached_function
def elements_of_norm(k,i):
    """
    Returns a list of elements of k of norm i, for k=Q(sqrt(-d)) imaginary
    quadratic field with d>3. The list contains both 'x' and '-x'.
    """
    L = []
    for l in k.elements_of_norm(i):
        L.append(l)
        L.append(-l)
    return L

def list_of_bounded_elems(k, minB, maxB):
    """
    Returns a list of elements of k with norm in the interval [minB, maxB]
    (as before k=Q(sqrt(-d)) imaginary quadratic field with d>3).
    """
    L = []
    for i in range(minB + 1,maxB + 1):
        L = L + elements_of_norm(k,i)
    return L


def in_circle(k, P, s, OnlyInside=False):
    """
    Returns True if the point P = (P0, P1) (coordinates in k basis) is 
    inside the circle s (where s = [alpha, i] with alpha in k the centre of
    the circle, i in ZZ the square of the inverse of the radius).
    
    Optional parameter 'OnlyInside' to specify if we want strict inclusion.
    """
    rad = 1/s[1] #this rad is the square of the radius, and is rational
    if OnlyInside:
        try:
            return (k(P)-s[0]).norm() < rad 
        except TypeError:
            D = k.absolute_generator().norm()
            return ((s[0][0]-P[0])^2 + (s[0][1]*sqrt(D) - P[1])^2)< rad
    else:
        try:
            return (k(P)-s[0]).norm() <= rad 
        except TypeError:
            D = k.absolute_generator().norm()
            return ((s[0][0]-P[0])^2 + (s[0][1]*sqrt(D) - P[1])^2)<= rad

def find_circles(k, P, S, OnlyInside=False):
    """
    Given 'P' a point and 'S' a list of circles,returns a list of the
    circles in S that contain the point 'P'.
    Optional argument 'OnlyInside' to specify strict inclusion.
    """
    S0 = []
    for s in S:
        if in_circle(k, P, s, OnlyInside):
            S0.append(s)
    return S0

def hem_covered(k, s, L):
    """
    Given 's' hemisphere, 'L' list of hemispheres, returns 'True' if 's' is 
    contained in one of the hemispheres of 'L'.
    """
    rad = 1/sqrt(s[1])
    for l in L:
        #we first check if the centre of 's' is contained in the hemisphere
        #'l' of L (in fact in the projection of 'l' on the plane).
        if in_circle(k, s[0], l, OnlyInside=True):
            rad_l = 1/sqrt(l[1])
            if (s[0] - l[0]).norm() <= (rad - rad_l)^2:
                return True
    return False

#----- Main functions to list hemispheres --------------------

def initial_list(k):
    """
    Returns a list of hemispheres that cover the fundamental rectangle on the 
    floor of the hyperbolic 3-space.
    """
    dk = k.discriminant()
    if dk.mod(4) == 1:
        D = dk.abs()
        b = 4 #fund region is [0, 1/2] x[0, sqrt(D)/b]
    else:
        D = (dk/4).abs()
        b = 2
    #we initialize the list of hemispheres with the only principal
    #hemispheres of radius 1
    if b==2:
        L = [(k(0), 1), (k(1), 1), (k((0,1)),1)]
    else:
        L = [(k(0), 1), (k(1), 1), (k((1/2,1/2)),1)]
    num = list_of_bounded_elems(k, 0, 1 + D)
    if k.class_number()==1: #for h=1 I use a different recursion
        L = check_rectangle(k, (0, 0), 1/2, 1/b, num, L)
    else:
        i = 2
        Sing = singular_points_in_F(k)
        check = False
        Pcov = []
        while check==False: #we loop over the size of the radius
            Li, num = list_for_radius_i(k, i, num, L)
            if Li:
                L = L + Li
                if not check:
                    check = True
                    for P in Sing:
                        if not (P in Pcov):
                            if find_circles(k, P, Li):
                                checkP = check_cover_for_sing_P(k, P, L)
                            else: #not adding new circles through P
                                checkP = False
                            if checkP:
                                Pcov.append(P)
                            check = check and checkP
            i = i + 1 
        R = list_rectangles(k, Sing, L)
        L = check_rectangles(k, i, num, R, L) #breadth-first recursion check
    return L

def list_for_radius_i(k, i, num, S, R = None):
    """
    Given a list of hemispheres 'S' and an integer 'i', return a list of 
    hemispheres of radius sqrt(1/i) which are not covered by any of the bigger
    hemispheres already in 'S' and that intersect the rectangle 'R'.

    By default 'R' is the fundamental rectangle.
    Updates the list of possible numerators (for the centres, given by cusps)
    """
    L = []
    dk = k.discriminant()
    if dk.mod(4) == 1:
        D = dk.abs()
        b = 4 #fund region is [0, 1/2] x[0, sqrt(D)/b]
    else:
        D = (dk/4).abs()
        b = 2
    if not R:
        R = ((0, 0), 1/2, 1/b)
    w = RR(D).sqrt()
    rad = RR(1/i).sqrt() #same radius since all denominators have same norm=i
    den = elements_of_norm(k,i)
    den = [l for l in den if l[1]>0 or (l[1]==0 and l[0]>=0)]
    n = k(num[-1]).norm()
    if n < i*(1 + D):
        num = num + list_of_bounded_elems(k, n, i*(1 + D))
    for y in den: #we add all hemispheres of radius rad...
        for x in num:
            alpha = x/y
            if (-1 <= alpha[0] <= 1) and (0 <= alpha[1] <=1): 
                if k.ideal(x, y).norm()==1:
                    if (-rad + R[0][0]) <= alpha[0] <= (R[0][0] + R[1] + rad):
                        if (-rad + R[0][1]*w) <= alpha[1]*w <= \
                           (rad + (R[0][1] + R[2])*w):
                            if not hem_covered(k, (alpha, i), S):
                                L.append((alpha, i))
                                S = S + L
    return L, num

#----- Functions concerning singular points ---------------------

# Set of singular points ----------------------------------------

def singular_points(k):
    """
    Returns a list of representatives for the singular points modulo 
    translations by elements of the ring of integers of the field.

    Algorithm: Swan.
    """
    L = []
    dk = k.discriminant()
    if dk.mod(4) == 1:
        D = dk.abs()
    else:
        D = (dk/4).abs()
    if ZZ(D).mod(4) == 3:
        s = 4
        step = 2
    else:
        s = 2
        step = 1
    while (3/4)*s^2<=D:
        for r in range(-s/2 + 1, s/2 + 1):
            if s^2<=(r^2 + D):
                if (step*s).divides(r^2 + D):
                    s0 = ZZ(s/step)
                    for p in range(s0):
                        if ZZ(p).gcd(s0)==1:
                            L.append((p*r/s, p/s))
        s = s + step
    return L

def singular_points_in_F(k):
    """
    Returns a list with the singular points which are on the floor of the 
    fundamental region.
    """
    L = singular_points(k)
    LinF = []
    dk = k.discriminant()
    if dk.mod(4)==1:
        h = 1/4
        d = 2
    else:
        h = 1/2
        d = 1
    for P in L: #we test if P or a translated P is inside Fk
        sign = 1
        t_rang = xrange(P[0] - 1/2, P[0] + 1)
        if P[0]<0:
            sign = -1
            t_rang = xrange(sign*(P[0]), sign*(P[0] - 1/2) + 1)
        for k in range(-d*(h - P[1]), d*P[1] + 1):
            if 0< (P[1] - k/d) <= h:
                for t in t_rang:
                    if k/d <= (P[0] - sign*t) <= (1/2 + k/d):
                        LinF.append((P[0] - sign*t - k/d, P[1] - k/d))
    return LinF

#----- Checking covering of singular points (all is 2-dimensional)----

def check_cover_for_sing_P(k, P, S):
    """
    Returns True if the point P is totally surrounded by circles (through P).
    """
    D = k.absolute_generator().norm()
    S0 = find_circles(k, P, S)
    P = k(P)
    # first we find the list of angles around which the point P is covered
    L = []
    for s in S0:
        if P[1]==s[0][1]:
            m = oo
        else:
            m = - (P[0]-s[0][0])/((P[1]-s[0][1])*sqrt(D))
        zeta = arctan(m)
        if m==0:
            if P[1]>s[0][1]:
                L.append((pi, 2*pi))
            else:
                L.append((0, pi))
        else:
            if m>0:
                if P[0]<s[0][0]:
                    L = L + [(pi + zeta, 2*pi + 1), (-1, zeta)]
                else:
                    L.append((zeta, zeta + pi))
            else:
                if P[0]<s[0][0]:
                    L = L + [(2*pi + zeta, 2*pi +1), (-1, pi + zeta)]
                else:
                    L.append((pi + zeta, 2*pi + zeta))

    # now we check if our list of angles covers all of [0, 2*pi]
    bmax = 0
    while bmax<=2*pi:
        cov = 0
        for j in L:
            if j[0]<bmax and j[1]>bmax:
                bmax = j[1]
                cov = 1
        if cov==0:
            return False
    return True

# Intersection of 'principal' circles surrounding a singular point:

def intersection_point(k, P, c1, c2):
    """
    P singular point, c1, c2 circles through P.
    Returns other point of intersection of c1 and c2 (different from P)
    """
    D = k.absolute_generator().norm()
    rad1 = 1/c1[1] #square of the radius in fact
    rad2 = 1/c2[1]
    P = k(P)
    a = c2[0][0] - c1[0][0]
    b = c2[0][1] - c1[0][1]
    c = 1/2*(rad1 - rad2 + c2[0][0]^2 - c1[0][0]^2 + \
        (c2[0][1]^2 - c1[0][1]^2)*D)
    if a.is_zero():
        y = c/(b*D) 
        B = -2*c1[0][0]
        x = -B - P[0]
        return (x, y)
    else:
        A = D*(b/a)^2 + 1
        B = 2*(c1[0][0]*b/a - b*c/a^2 - c1[0][1])
        y = -B/A - P[1]
        return (1/a*(c-b*y*D) , y)

def intersection_singular_circles(k, P, c1, c2):
    """
    P singular point, c1, c2 circles through P.
    Checks if c1 and c2 intersect and not only touch at P, and then
    returns the intersection point.
    """
    rad1 = 1/c1[1]
    rad2 = 1/c2[1]
    # we are only interested in circles that intersect in more than
    # one point (=sing point)
    if k(c1[0] - c2[0]).norm() < rad1 + rad2 + 2*sqrt(rad1*rad2):
        return intersection_point(k, P, c1, c2)
    else:
        return ()

def intersection_sing_point(k, P, S):
    """
    P singular point, S list of circles.
    Returns all points of intersection (except trivial P) between the
    circles of S that go through P.
    """
    s = find_circles(k, P, S)
    L = []
    Ldone = []
    for c1 in s:
        Ldone.append(c1)
        for c2 in s:
            if not (c2 in Ldone):
                l = intersection_singular_circles(k, P, c1, c2)
                if l and not l==P: #catch wrong intersec point
                    L.append(l)
    return L

#----- Rectangles around singular points ---------------------

def rect_around_sing_point(k, sk, S):
    """
    For sk singular point, returns a rectangle around sk which we can
    guarantee is covered by the circles in the set S.
    The rectangle is given in the form (Point, width, height).
    """
    dk = k.discriminant()
    if dk.mod(4) == 1:
        h = 1/4 #fund region is [0, 1/2] x[0, sqrt(D)*h]
    else:
        h = 1/2
    L = intersection_sing_point(k, sk, S)
    # now find intersection point inside fundamental region at minimum
    # distance to sk
    minP = 0
    P = ()
    for u in L:
        if 0<=u[0]<=1/2 and 0<=u[1]<=h:
            if minP==0:
                P = u
                minP = (k(sk) - k(P)).norm()
            else:
                d = (k(sk) - k(u)).norm()
                if d < minP:
                    minP = d
                    P = u
    # if all intersection points fall out of the fundamental region:
    if not P:
        maxP1 = max([v[1] for v in L])
        minP1 = min([v[1] for v in L])
        if maxP1 > h:
            return ((0, minP1), 1/2, h - minP1)
        else:
            if minP1 > 0:
                return ((0, minP1), 1/2, maxP1 - minP1)
            else:
                return ((0, 0), 1/2, maxP1)

    # in most cases at least one of the intersection points is inside the
    # fundamental region; then we find the the rectangle as follows:
    wP = 2*(sk[0] - P[0])
    hP = 2*(sk[1] - P[1])
    if wP > 0:
        r1 = P[0]
        # adjust the rectangle in case crosses out of fundamental region
        if (r1 + wP) > 1/2:
            wP = 1/2 - r1
    else:
        wP = -wP
        r1 = P[0] - wP
        if r1 < 0: # adjusting rectangle
            wP = r1 + wP
            r1 = 0
    if hP > 0:
        r2 = P[1]
        if (r2 + hP) > h:
            hP = h - r2
    else:
        hP = -hP
        r2 = P[1] - hP
        if r2 < 0:
            hP = r2 + hP
            r2 = 0
    return (r1, r2), wP, hP

#----- Checking rectangles -----------------------------------------
#----- Checking if a list of circles (hemispheres) covers all of F_K

def rectangle_cov(k, r, S):
    """
    Checks if rectangle r is covered by the circles in S.
    Returns 0 if any corner of r is not covered, 1 if all corners are covered
    (by different circles) and 2 if all corners covered by one circle from S
    (so the whole of r is covered).
    """
    c = find_circles(k, r[0], S)
    if not c:
        return 0 # = 0: need more hemispheres, v is not covered at all!
    vv = [(r[0][0] + r[1], r[0][1]), (r[0][0], r[0][1] + r[2]), \
          (r[0][0] + r[1], r[0][1] + r[2])]
    check = 2
    for u in vv:
        if check == 1: #two vertices at least not covered by same hemisphere
            if not find_circles(k, u, S):
                return 0 #that is, this particular vertex not covered at all!
        else:
            c_u = [s for s in c if in_circle(k, u, s)]
            if not c_u: #vertex u not in any circle s containing v =>c_u empty
                check = 1
                if not find_circles(k, u, S): #no circle contains v
                    return 0
            else:
                c = c_u
    return check #if check=2, whole rectangle covered by one circle

def check_rectangle(k, v, w, h, num, S):
    """
    Checks if rectangle of width 'w', height 'h' is covered by the circles of 
    the set S.
    Adds more hemispheres to S until the rectangle is covered.
    We only use this function for fields with class number 1.
    """
    check = rectangle_cov(k, (v, w, h), S)
    if check==2:
        return S
    if check==0:
        i = S[-1][1] #last radius in S
        while check==0:
            Li, num = list_for_radius_i(k, i + 1, num, S)
            if Li: #we have added something new
                S = S + Li
                check = rectangle_cov(k, (v, w, h), S)
            i = i + 1
        if check == 2:
            return S
    if check==1:
        R0 = rectangle_subdivision(k, (v, w, h))
        for r in R0:
            Sr = check_rectangle(k, r[0], r[1], r[2], num, S)
            S = Sr
        return S

def check_rectangles(k, i, num, R, S):
    """
    Checks if all rectangles in the list R are covered by the circles
    of the set S. Adds more hemispheres to S until all rectangles are
    covered. We use this function for fields with class number > 1,
    where previous (depth) recursive algorithm doesn't work in general.
    """
    if not R:
        return S
    R0 = []
    for r in R:
        check = rectangle_cov(k, r, S)
        if check == 1:
            R0 = R0 + rectangle_subdivision(k, r)
        if check == 0:
            norm_d = i
            while check==0:
                #print "checking r = ", r
                Li, num = list_for_radius_i(k, norm_d + 1, num, S, r)
                if Li: #we have added something new
                    #print "adding ", Li
                    S = S + Li
                    check = rectangle_cov(k, r, S)
                norm_d = norm_d + 1
            if check == 1:
                R0 = R0 + rectangle_subdivision(k, r)
    return check_rectangles(k, i, num, R0, S)

def rectangle_subdivision(k, r):
    """
    Returns list of rectangles (roughly squares) that gives a subdivision
    of rectangle 'r'
    """
    R = []
    sqD = RR(k.absolute_generator().norm()).sqrt()
    P, wP, hP = r
    if hP*sqD > wP:
        t = 2
        th = t*floor(hP*sqD/wP)
    else:
        th = 2
        t = 2*floor(wP/(hP*sqD))
    rw = wP/t
    rh = hP/th
    for j in range(th):
        for i in range(t):
            u1 = P[0] + i*rw
            u2 = P[1] + j*rh
            R.append(((u1, u2), rw, rh))
    return R

def list_rectangles(k, sk, S):
    """
    Returns a list of rectangles that contain all the areas of fundamental
    region that need checking (excludes surroundings of singular points).
    """
    dk = k.discriminant()
    if dk.mod(4) == 1:
        D = dk.abs()
        h = 1/4 #fund region is [0, 1/2] x[0, sqrt(D)*h]
    else:
        D = (dk/4).abs()
        h = 1/2
    R = []
    R = R + subdivide_rect(k, (0, 0), 1/2, h, sk, S)
    return R

def subdivide_rect(k, P, wP, hP, sk, S):
    """
    P, wP, hP a rectangle; sk set of singular points; S set of circles.
    Subdivides rectangle (P, wP, hP) into rectangles that don't contain 
    singular points.
    """
    R = []
    L = [(rect_around_sing_point(k, s, S)) for s in sk]
    num, sing = number_sing(k, P, wP, hP, L)
    if num == 0:
        R.append((P, wP, hP))
    if num == 1:
        R = R + one_subdiv(k, P, wP, hP, sing[0], S)
    if num > 1:
        R0 = rectangle_subdivision(k, (P, wP, hP))
        for r in R0:
            R = R + subdivide_rect(k, r[0], r[1], r[2], sk, S)
    return R

def number_sing(k, p, wp, hp, Lr):
    """
    p, wp, hp a rectangle; Lr list of neighborhoods of singular points.
    Returns number of singular points and list of rectangles around them
    that cut the given rectangle (p, wp, hp)
    """
    n = 0
    L = [] # list of the neighborhoods inside our main rectangle
    for r in Lr:
        if (p[0]<=r[0][0]<(p[0]+wp)) or (p[0]<(r[0][0] + r[1])<(p[0]+wp)):
            if (p[1]<r[0][1]<(p[1]+hp)) or (p[1]<(r[0][1] + r[2])<(p[1]+hp)):
                n = n + 1
                L.append(r)
    return n, L

def one_subdiv(k, P, wP, hP, rk, S):
    """
    Return subdivision of a rectangle (P, wP, hP) which includes only one
    singular rectangule rk.
    """
    R = []
    u, wu, hu = rk[0], rk[1], rk[2]
    w0 = 0
    if u[0] > P[0]:
        w0 = u[0] - P[0]
        R.append(((P[0], P[1]), w0, hP))
    if (u[0] + wu) < (P[0] + wP):
        w1 = (u[0] + wu) - (P[0] + w0)
        R.append(((P[0] + w0 + w1, P[1]), (P[0] + wP) - (u[0] + wu), hP))
    else:
        w1 = (P[0] + wP) - (P[0] + w0)
    if u[1] > P[1]:
        R.append(((P[0] + w0, P[1]), w1, u[1] - P[1]))
    if (u[1] + hu) < (P[1] + hP):
        R.append(((P[0] + w0, u[1] + hu), w1, (P[1] + hP) - (u[1] + hu)))
    return R

#--------------------------------------------------------
#----- Cleaning list of redundant hemispheres -----------

#----- Intersections of 3 hemispheres ----------------

def circles_intersection_check(k, c1, c2):
    """
    Checks if circles c1, c2 intersect; returns True or False.
    We assume circles are not included one inside the other.
    """
    rad1 = 1/c1[1]
    rad2 = 1/c2[1]
    c_distance = sqrt(k(c1[0]-c2[0]).norm())
    if RR(c_distance) < RR(sqrt(rad1)) + RR(sqrt(rad2)):
        return True
    else:
        return False

def hemispheres_intersection(k, s1, s2, s3):
    """
    Checks if hemispheres s1, s2, s3 truly intersect (assumes that they
    intersect pairwise), and if any of them happens to be covered by the
    union of the other two. Returns redundant hemisphere (if any) and
    point of intersection (if any).
    """    
    a1 = s2[0][0] - s1[0][0]
    b1 = s2[0][1] - s1[0][1]
    a2 = s3[0][0] - s1[0][0]
    b2 = s3[0][1] - s1[0][1]
    if (b1.is_zero() and b2.is_zero()) or (a1.is_zero() and a2.is_zero()) \
          or (b2*a1==b1*a2):
        s = check_redundant_hem(k, s1, s2, s3)
        return s, ()
    D = k.absolute_generator().norm()
    rad1 = 1/s1[1]
    rad2 = 1/s2[1]
    rad3 = 1/s3[1]
    c1 = 1/2*(rad1 - rad2 + s2[0][0]^2 - s1[0][0]^2 \
                          + (s2[0][1]^2 - s1[0][1]^2)*D)
    c2 = 1/2*(rad1 - rad3 + s3[0][0]^2 - s1[0][0]^2 \
                          + (s3[0][1]^2 - s1[0][1]^2)*D) 

    if a1.is_zero(): #then a2 nonzero
        yP = c1/(b1*D)
        xP = 1/a2*(c2 - b2*D*yP)
    else: #a1 not zero
        yP = (c2 - c1*a2/a1)/(D*(b2 - a2*b1/a1))
        xP = 1/a1*(c1 - b1*D*yP)
    if not all([in_circle(k, (xP, yP), s) for s in [s1, s2, s3]]):
        return (), ()
    tP = rad1 - (k((xP, yP)) - s1[0]).norm()
    return (), (xP, yP, tP)

def check_redundant_hem(k, s1, s2, s3):
    """
    Given hemispheres s1, s2, s3 with centres in same line, the function
    checks if any of the 3 hemispheres is covered by the other two (by
    checking the position of intersection lines).
    """
    D = k.absolute_generator().norm()
    L = [((s[0][0], s[0][1]),s[1]) for s in [s1, s2, s3]]
    L.sort()
    s1, s2, s3 = L[0], L[1], L[2]
    b12 = s2[0][1] - s1[0][1]
    b23 = s3[0][1] - s2[0][1]
    rad1 = 1/s1[1]
    rad2 = 1/s2[1]
    rad3 = 1/s3[1]
    c12 = 1/2*(rad1 - rad2 + s2[0][0]^2 - s1[0][0]^2 \
                           + (s2[0][1]^2 - s1[0][1]^2)*D)
    c23 = 1/2*(rad2 - rad3 + s3[0][0]^2 - s2[0][0]^2 \
                           + (s3[0][1]^2 - s2[0][1]^2)*D) 

    if b12.is_zero(): #horizontal line
        a12 = s2[0][0] - s1[0][0]
        a23 = s3[0][0] - s1[0][0]
        if c23/a23 <= c12/a12:
            return (k((s2[0][0], s2[0][1])), s2[1])
        else:
            return ()
    a12 = s2[0][0] - s1[0][0]
    if (-a12/b12) > 0: 
        if c23/b23 >= c12/b12:
            return (k((s2[0][0], s2[0][1])), s2[1])
        else:
            return ()
    else:
        if c23/b23 <= c12/b12:
            return (k((s2[0][0], s2[0][1])), s2[1])
        else:
            return()
    return ()

# the function that 'cleans' the list of unnecessary hemispheres

def clean_list(k, S):
    """
    Cleans list of hemispheres S of unnecessary ones.
    Returns the list of points of 3-intersection and the list of hemispheres.
    """
    dk = k.discriminant()
    if dk.mod(4) == 1:
        h = 1/4 #fund region is [0, 1/2] x[0, sqrt(D)*h]
    else:
        h = 1/2   
    Sing = singular_points_in_F(k)
    L = S[:]
    IntPoints = []
    for s1 in L:
        check = True
        checkSing = False
        if s1 in S:
            nP = 0
            Ldone = [s1]
            for s2 in L:
                if s2 in S:
                    Ldone.append(s2)
                    if (not(s2==s1)) and circles_intersection_check(k, s1, s2):
                        for s3 in [s for s in L if s in S and not s in Ldone]:
                            if circles_intersection_check(k, s1, s3) and \
                                       circles_intersection_check(k, s2, s3):
                                red, P = hemispheres_intersection(k, s1, s2, s3)
                                if red and red in S:
                                    S.remove(red)
                                elif P: #no redundant hemispheres
                                    if P in IntPoints:
                                        check = False
                                    elif P[2]==0:
                                        if (P[0], P[1]) in Sing:
                                            checkSing = True
                                    else:
                                        sP = find_hemispheres(k, P, S)
                                        if not sP: #nothing covers P
                                            check = False #we need s1
                                            IntPoints.append(P)
                                        else:
                                            nP = nP + 1
            if checkSing and nP==0:
                check = False
            if check and (s1 in S):
                S.remove(s1)
    return IntPoints, S

def find_hemispheres(k, P, S, OnSurface=False):
    """
    Finds hemispheres in S that cover the point P.
    """
    L = []
    sP = find_circles(k, (P[0], P[1]), S)
    for s in sP:
        t = 1/s[1] - (k((P[0], P[1])) - s[0]).norm()
        if OnSurface:
            if t==P[2]:
                L.append(s)
        else:
            if t>P[2]:
                L.append(s)
    return L

#--------------------------------------------------------
#----- Swan's algorithm ---------------------------------

def find_principal_hemisphere(k, P, maxNorm=None):
    """
    Finds principal hemisphere that covers P, and that does it with maximum
    height over P.
    """
    den = list_of_bounded_elems(k, 0, 1/P[2])
    den = [l for l in den if l[1]>0 or (l[1]==0 and l[0]>=0)]
    maxnormmu = k(den[-1]).norm()
    if maxNorm:
        if maxnormmu > maxNorm:
            den = [l for l in den if l.norm()<=maxNorm]
            maxnormmu = maxNorm

    z = k((P[0], P[1]))
    B = (1 + RR((z*den[-1]).norm()).sqrt())^2
    lam = [0] + list_of_bounded_elems(k, 0, B.ceil())
    tmax = 0
    sP = () # the hemisphere covering P

    for mu in den:
        normmu = k(mu).norm()
        if 1/normmu < tmax:
            return sP
        Bmu = (((RR(1 - normmu*P[2])).sqrt() \
               + RR((z*mu).norm()).sqrt())^2).ceil()
        for l in [l for l in lam if k(l).norm()<Bmu]:
            if k.ideal(l, mu).norm()==1:
                if (mu*z - l).norm() + normmu*P[2] < 1:
                    t = 1/normmu - (z - l/mu).norm()
                    if t > tmax:
                        tmax = t
                        sP = (k(l/mu), normmu)
    return sP

def Swan(k, S, V=None, maxNorm=None):
    """
    Swan's algorithm.
    Returns the list of hemispheres that give the floor of the fundamental
    region, and a list of all the intersection points (they are potential
    vertices).
    """
    LVchecked = []
    LVremove = []
    while True:
        check = True
        V = clear_inter_list(k, V)
        for v in V:
            if check and not v in LVchecked:
                if 0<=len(find_hemispheres(k, v, S, OnSurface=True))<3:
                    #the point shouldnt really be here
                    LVchecked.append(v) 
                    LVremove.append(v)
                else:
                    #print "find hemisph for v", v
                    sv = find_principal_hemisphere(k, v, maxNorm)
                    if sv:
                        v0 = v
                        if not sv in S: 
                            S.append(sv) 
                        check = False
                    else:
                        LVchecked.append(v)
        if check:
            break
        #print "checking point v0", v0
        V, S = adjust_list(k, v0, V, S)
    V = [v for v in V if not v in LVremove]
    return V, S

def clear_inter_list(k, V):
    """
    Clears the list of intersection points of symmetries (to avoid checking
    things twice).
    """
    dk = k.discriminant()
    if dk.mod(4) == 1:
        h = 1/4
    else:
        h = 1/2
    L = V[:]
    for v in V:
        if v[1]<=h:
            if v[0]<0:
                if [u for u in V if u[0]==-v[0] and u[1]==v[1]]:
                    L.remove(v)
            if v[0]>1/2:
                if [u for u in V if u[0]==1 - v[0] and u[1]==v[1]]:
                    L.remove(v)
        else:
            if 0<=v[0]<=1/2:
                if h == 1/2: # we have symmetry
                    if [u for u in V if u[0]==v[0] and u[1]==2*h - v[1]]:
                        L.remove(v)
                else: # we have glide reflection
                    if [u for u in V if u[0]==(1/2-v[0]) and \
                                        u[1]==2*h - v[1]]:
                        L.remove(v)
            else:
                L.remove(v)
    return L

def adjust_list(k, v0, V, S):
    """
    After adding a hemisphere that covers the point v0, we check again our
    list S for new redundancies resulting from that addition.
    """
    dk = k.discriminant()
    if dk.mod(4) == 1:
        h = 1/4 #fund region is [0, 1/2] x[0, sqrt(D)*h]
    else:
        h = 1/2
    Sv = [s for s in S if in_circle(k, k((v0[0], v0[1])), s)]
    LRed = []
    Ldone = []
    NewPoints = []
    for s1 in Sv:
        nP = 0 #number of intersection points for s1 that are covered
        checkSing = False
        check = True
        if s1 in S:
            Ldone = [s1]
            for s2 in S:
                Ldone.append(s2)
                if (not (s2==s1)) and circles_intersection_check(k, s1, s2):
                    for s3 in [s for s in S if s in S and not s in Ldone]:
                        if circles_intersection_check(k, s1, s3) and \
                                        circles_intersection_check(k, s2, s3):
                            red, P = hemispheres_intersection(k, s1, s2, s3)
                            if red:
                                LRed.append(red)
                            elif P: #no redundant hemispheres and there is int
                                if P in NewPoints:
                                    check = False
                                elif P[2]==0:
                                    if 0<=P[0]<=1/2 and 0<=P[1]<=h:
                                        checkSing = True
                                else:
                                    sP = find_hemispheres(k, P, S)
                                    if not sP: #no spheres covering P
                                        check = False #we need to keep s1
                                        if not P in V:
                                            NewPoints.append(P)
                                    else:
                                        nP = nP + 1
                                        #nP = number of covered points
                                        if P in V:
                                        #P not true vertex
                                            V.remove(P)
        if checkSing and nP==0:
        #all other intersection points except singular one are uncovered
        #this check is necessary for spheres that only intersect at sing point
            check = False
        if check:
            #S.remove(s1)
            LRed.append(s1)
    for s in LRed:
        if s in S:
            S.remove(s)
    V = V + NewPoints
    return V, S
