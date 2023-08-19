# -*- coding: utf-8 -*-


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
        if s[1] in ZZ:
            rad = 1/RR(s[1]).sqrt()
        else:
            rad = RR(s[1]).sqrt()
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
        if s[1] in ZZ:
            rad = 1/s[1]
        else:
            rad = s[1]
        eq = (X - c[0])^2 + (Y - c[1]*w)^2 + Z^2 - rad
        L.append(implicit_plot3d(eq, (Y, 0, Yrange ),  (X, -0.5, 0.5), \
                           (Z, 0, 1), plot_points=60, aspect_ratio=1)) 
    return sum(L)


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
    if s[1] in ZZ:
        rad = 1/s[1] #this rad is the square of the radius, and is a rational
    else:
        rad = s[1]
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
    if s[1] in ZZ:
        rad = 1/sqrt(s[1])
    else:
        rad = sqrt(s[1])
    for l in L:
        #we first check if the new "hemisphere" 's' is already in L
        if s[0]==l[0]:
            return True
        #we then check if the centre of 's' is contained in the hemisphere
        #'l' of L (in fact in the projection of 'l' on the plane).
        if in_circle(k, s[0], l, OnlyInside=True):
            if l[1] in ZZ:
                rad_l = 1/sqrt(l[1])
            else:
                rad_l = sqrt(l[1])
            if (s[0] - l[0]).norm() <= (rad - rad_l)^2:
                return True
    return False

# plotting functions:

def plot_circles(S, k):
    L = []
    D = k.absolute_generator().norm()
    w = RR(D).sqrt()
    for s in S:
        centre = (s[0][0], s[0][1]*w)
        if s[1] in ZZ:
            rad = 1/RR(s[1]).sqrt()
        else:
            rad = RR(s[1]).sqrt()
        L.append(circle(centre, rad, aspect_ratio=1, \
                 fill=True, alpha = 0.2, rgbcolor = 'blue'))
    return sum(L)

def plot_sing_points(S, k):
    L = []
    D = k.absolute_generator().norm()
    w = RR(D).sqrt()
    for s in S:
        P = (s[0], s[1]*w)
        L.append(point(P, rgbcolor = 'red', pointsize=22))
    return sum(L)

def plot_rect(v, w, h, k):
    D = k.absolute_generator().norm()
    vv = [v, (v[0] + w, v[1]), (v[0] + w, v[1] + h), (v[0], v[1] + h)] 
    P = [(v[0], v[1]*RR(D).sqrt()) for v in vv]
    P.append(P[0])
    return line(P, rgbcolor = hue(0.75), thickness=1.5)
    #return list_plot(P, rgbcolor = 'green', pointsize=22)

def plot_spheres(S, k):
    X, Y, Z = var('X, Y, Z')
    L = []
    dk = k.discriminant()
    if dk.mod(4) == 1:
        D = dk.abs()
    else:
        D = (dk/4).abs()
    sqrtD = RR(D).sqrt()
    for s in S:
        c = k(s[0])
        if s[1] in ZZ:
            rad2 = 1/s[1]
        else:
            rad2 = s[1]
        eq = (X - c[0])^2 + (Y - c[1]*sqrtD)^2 + Z^2 - 1/s[1]
        L.append(implicit_plot3d(eq, (X, -1, 1), (Y, 0, 2 ), (Z, 0, 2), aspect_ratio=1))#, opacity=0.5))
    return sum(L)

# Main hemispheres enumerating functions:

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
    #we initialize the list of hemispheres with the only principal hemispheres 
    #with radius 1
    if b==2:
        L = [(k(0), 1), (k(1), 1), (k((0,1)),1)]
    else:
        L = [(k(0), 1), (k(1), 1), (k((1/2,1/2)),1)]
    num = list_of_bounded_elems(k, 0, 1 + D) 
    i = 1 # last radius used (should be i=2 for check_rectangle)
    # now breadth-first recursion check
    L = check_rectangles(k, i, num, [((0, 0), 1/2, 1/b)], L)
    return L

def list_for_radius_i(k, i, num, S, R = None):
    """
    Given a list of hemispheres 'S' and an integer 'i', return a list of 
    hemispheres of radius sqrt(1/i) which are not covered by any of the bigger
    hemispheres already in 'S' and that intersect the rectangle 'R'.

    By default 'R' is the fundamental rectangle.
    I update the list of possible numerators (for the centres, given by cusps)
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
                else:
                    Idxy = k.ideal(x, y)
                    if (not Idxy.is_principal()) and (Idxy^2).is_principal():
                        radD = RR(Idxy.norm()/i).sqrt()
                        if (-radD + R[0][0]) <= alpha[0] <= (R[0][0] + R[1] + radD):
                            if (-radD + R[0][1]*w) <= alpha[1]*w <= \
                            (radD + (R[0][1] + R[2])*w):
                                if not hem_covered(k, (alpha, i), S):
                                    L.append((alpha, Idxy.norm()/i))
                                    S = S + L

    return L, num


# Functions concerning singular points:
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

# Checking covers ---------------------------------------
# Checking if a list of circles (hemispheres) covers the whole F_K

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
                if not find_circles(k, u, S): #no circle contains the vertex v
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
        if not i in ZZ:
            i = k(S[-1][0]).denominator()
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
    Checks if all rectangles in the list R are covered by the circles of the 
    set S. Adds more hemispheres to S until all rectangles are covered.
    We use this function for fields with class number > 1, where previous 
    (depth) recursive algorithm doesn't work in general.
    """
    if not R:
        return S
    R0 = []
    for r in R:
        check = rectangle_cov(k, r, S)
        if check == 1:
            R0 = R0 + rectangle_subdivision(k, r)
        if check == 0:
            #i = S[-1][1] #last radius in S
            #if not i in ZZ:
            #    i = k(S[-1][0]).denominator()
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
    Returns list of rectangles (roughly squares) that gives a subdivision of 
    rectangle 'r'
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


# Cleaning initial list

# 1: Intersections of 3 hemispheres:

def circles_intersection_check(k, c1, c2):
    """
    Checks if circles c1, c2 intersect; returns True or False.
    We assume circles are not included one inside the other.
    """
    if c1[1] in ZZ:
        rad1 = 1/c1[1]
    else:
        rad1 = c1[1]
    if c2[1] in ZZ:
        rad2 = 1/c2[1]
    else:
        rad2 = c2[1]
    c_distance = sqrt(k(c1[0]-c2[0]).norm())
    if RR(c_distance) < RR(sqrt(rad1)) + RR(sqrt(rad2)):
        return True
    else:
        return False

def hemispheres_intersection(k, s1, s2, s3):
    """
    Checks if hemispheres s1, s2, s3 truly intersect (assumes that they
    intersect pairwise), and if any of them happens to be covered by the union
    of the other two.
    Returns redundant hemisphere (if any) and point of intersection (if any).
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
    if s1[1] in ZZ:
        rad1 = 1/s1[1]
    else:
        rad1 = s1[1]
    if s2[1] in ZZ:
        rad2 = 1/s2[1]
    else:
        rad2 = s2[1]
    if s3[1] in ZZ:
        rad3 = 1/s3[1]
    else:
        rad3 = s3[1]
    c1 = 1/2*(rad1 - rad2 + s2[0][0]^2 - s1[0][0]^2 + (s2[0][1]^2 - s1[0][1]^2)*D)
    c2 = 1/2*(rad1 - rad3 + s3[0][0]^2 - s1[0][0]^2 + (s3[0][1]^2 - s1[0][1]^2)*D) 

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
    Given hemispheres s1, s2, s3 with centres in same line, the function checks
    if any of the 3 hemispheres is covered by the other two (by checking the 
    position of intersection lines).
    """
    D = k.absolute_generator().norm()
    L = [((s[0][0], s[0][1]),s[1]) for s in [s1, s2, s3]]
    L.sort()
    s1, s2, s3 = L[0], L[1], L[2]
    b12 = s2[0][1] - s1[0][1]
    b23 = s3[0][1] - s2[0][1]
    if s1[1] in ZZ:
        rad1 = 1/s1[1]
    else:
        rad1 = s1[1]
    if s2[1] in ZZ:
        rad2 = 1/s2[1]
    else:
        rad2 = s2[1]
    if s3[1] in ZZ:
        rad3 = 1/s3[1]
    else:
        rad3 = s3[1]
    c12 = 1/2*(rad1 - rad2 + s2[0][0]^2 - s1[0][0]^2 + (s2[0][1]^2 - s1[0][1]^2)*D)
    c23 = 1/2*(rad2 - rad3 + s3[0][0]^2 - s2[0][0]^2 + (s3[0][1]^2 - s2[0][1]^2)*D) 

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

# the function that cleans the list of unnecessary hemispheres

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
    L = S[:]
    IntPoints = []
    for s1 in L:
        check = True
        if s1 in S:
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
                                    elif P[2]>0:
                                        if 0<=P[0]<=1/2 and 0<=P[1]<=h:
                                            sP = find_hemispheres(k, P, S)
                                            if not sP: #nothing covers P
                                                check = False #we need s1
                                                IntPoints.append(P)
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
        if s[1] in ZZ:
            radsq = 1/s[1]
        else:
            radsq = s[1]
        t = radsq - (k((P[0], P[1])) - s[0]).norm()
        if OnSurface:
            if t==P[2]:
                L.append(s)
        else:
            if t>P[2]:
                L.append(s)
    return L

# Now Swan's algorithm:

def find_principal_hemisphere(k, P):
    """
    Finds a principal hemisphere that covers P with maximum height over P.
    """
    M = RR(k.discriminant().abs()/3).sqrt()
    den = list_of_bounded_elems(k, 0, M/P[2])
    den = [l for l in den if l[1]>0 or (l[1]==0 and l[0]>=0)]
    z = k((P[0], P[1]))
    maxnormmu = k(den[-1]).norm()
    B = ((RR(M - k(den[0]).norm()*P[2]/M)).sqrt() + RR((z*den[-1]).norm()).sqrt())^2
    lam = [0] + list_of_bounded_elems(k, 0, B.ceil())
    tmax = 0
    sP = ()
    for mu in den:
        normmu = k(mu).norm()
        if M/normmu < tmax:
            return sP
        Bmu = (((RR(M - normmu*P[2])).sqrt() + RR((z*mu).norm()).sqrt())^2).ceil()
        for l in [l for l in lam if k(l).norm()<Bmu]:
            if k.ideal(l, mu).norm()==1:
                if (mu*z - l).norm() + normmu*P[2] < 1:
                    t = 1/normmu - (z - l/mu).norm()
                    if t > tmax:
                        tmax = t
                        sP = (k(l/mu), normmu)
            else:
                Idlm = k.ideal(l, mu)
                if (not Idlm.is_principal()) and (Idlm^2).is_principal():
                    if (z -k(l/mu)).norm() + P[2] < Idlm.norm()/normmu:
                        t = Idlm.norm()/normmu - (z - l/mu).norm()
                        if t > tmax:
                            tmax = t
                            sP = (k(l/mu), Idlm.norm()/normmu)
    return sP

def Swan(k, S, V=None):
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
                    sv = find_principal_hemisphere(k, v)
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
                    if [u for u in V if u[0]==(1/2-v[0]) and u[1]==2*h - v[1]]:
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
                            elif P: #no redundant hemispheres
                                if P in NewPoints:
                                    check = False
                                elif P[2]>0:
                                    sP = find_hemispheres(k, P, S)
                                    if not sP: #no spheres covering P
                                        check = False #we need to keep s1
                                        if not P in V:
                                            NewPoints.append(P)
                                    elif P in V: #P not true vertex
                                        V.remove(P)
        if check:
            #S.remove(s1)
            LRed.append(s1)
    for s in LRed:
        if s in S:
            S.remove(s)
    V = V + NewPoints
    return V, S

