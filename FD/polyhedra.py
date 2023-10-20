##############################################
#
# Standard polyhedra (as graphs)
#
##############################################

from sage.all import Graph, parametric_plot3d, sin, cos, pi, CC

Tetrahedron = Graph([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
                    name="tetrahedron")

Cube = Graph([(0, 1), (0, 3), (0, 7), (1, 2), (1, 6), (2, 4), (2, 7),
              (3, 5), (3, 6), (4, 5), (4, 6), (5, 7)],
             name="cube")

SquarePyramid = Graph([(0, 1), (0, 2), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (3, 4)],
                    name = "square pyramid")

HexagonalPrism = Graph([(0, 1), (0, 5), (0, 9), (1, 2), (1, 10),
                        (2, 4), (2, 5), (3, 5), (3, 7), (3, 9),
                        (4, 8), (4, 10), (6, 7), (6, 9), (6, 11),
                        (7, 8), (8, 11), (10, 11)],
                       name = "hexagonal prism")

HexagonalCap = Graph([(1,2),(2,3),(3,4),(4,5),(5,6),(6,1),
              (1,7),(2,7),(3,8),(4,8),(5,9),(6,9),
              (7,8),(8,9),(9,7)],
             name="hexagonal cap")

TruncatedTetrahedron = Graph([(0, 6), (0, 8), (0, 10), (1, 5), (1, 8), (1, 11),
                  (2, 7), (2, 9), (2, 11), (3, 4), (3, 5), (3, 9),
                  (4, 5), (4, 6), (6, 10), (7, 9), (7, 10), (8, 11)],
                 name = "truncated tetrahedron")

TriangularPrism = Graph([(1,2),(2,3),(3,1),(4,5),(5,6),(6,4),(1,4),(2,5),(3,6)],
              name="triangular prism")

Cuboctahedron = Graph([(0, 2), (0, 6), (0, 7), (0, 11),
                       (1, 4), (1, 8), (1, 9), (1, 11),
                       (2, 6), (2, 8), (2, 10),
                       (3, 5), (3, 6), (3, 9), (3, 10),
                       (4, 5), (4, 7), (4, 9),
                       (5, 6), (5, 7), (7, 11), (8, 10), (8, 11),
                       (9, 10)],
                      name = "cuboctahedron")

Octahedron = Graph( [(0, 1), (0, 2), (0, 4), (0, 5), (1, 2), (1, 3),
                     (1, 4), (2, 3), (2, 5), (3, 4), (3, 5), (4, 5)],
                    name = "octahedron")

# Some compund polyhedra encountered:

# = cube with an inner diagonal face, also two triangular prisms glued along a square face
# eg d=42
SlicedCube = Graph([(0, 1), (0, 2), (0, 3), (1, 0), (1, 4), (1, 5), (6, 3), (6, 2), (6, 7),
                    (7, 5), (7, 4), (7, 6), (3, 6), (3, 0), (3, 2), (3, 4), (5, 7), (5, 1),
                    (5, 4), (5, 2), (2, 6), (2, 0), (2, 3), (2, 5), (4, 7), (4, 1), (4, 5), (4, 3)],
                   name = "sliced cube")

# Two hexagonal prisms glued along a hexagon
# eg d=57
DoubleHexagonalPrism = Graph([(0, 1), (0, 2), (0, 3), (4, 2), (4, 5), (4, 6), (7, 5), (7, 1), (7, 8),
                              (1, 0), (1, 7), (1, 9), (2, 4), (2, 0), (2, 10), (5, 7), (5, 4), (5, 11),
                              (3, 12), (3, 0), (3, 10), (3, 9), (8, 13), (8, 7), (8, 9), (8, 11),
                              (6, 14), (6, 4), (6, 11), (6, 10), (10, 15), (10, 2), (10, 3), (10, 6),
                              (9, 16), (9, 1), (9, 8), (9, 3), (11, 17), (11, 5), (11, 6), (11, 8),
                              (12, 15), (12, 3), (12, 16), (13, 16), (13, 8), (13, 17), (14, 17),
                              (14, 6), (14, 15), (15, 12), (15, 10), (15, 14), (16, 13), (16, 9),
                              (16, 12), (17, 14), (17, 11), (17, 13)],
                             name = "double hexagonal prism")

# Half a stella octangula, also a square prism with two tetrahedra glued on non-adjacent triangles
#  (includes a vertex of degree 6)
# eg d=29, 53
HalfStar = Graph( [(0, 1), (0, 2), (0, 3), (4, 5), (4, 6), (4, 3), (1, 0), (1, 2), (1, 5),
                   (1, 3), (5, 4), (5, 6), (5, 1), (5, 3), (2, 6), (2, 0), (2, 1), (2, 3),
                   (6, 2), (6, 4), (6, 5), (6, 3), (3, 4), (3, 0), (3, 2), (3, 5), (3, 1), (3, 6)],
                  name = "half star")

# hexagonal prism with 3 triangular prisms on alternate square faces
# eg d=33, 85
TentedHexPrism = Graph( [(0, 1), (0, 2), (0, 3), (0, 4), (5, 6), (5, 7), (5, 8), (5, 1), (9, 4), (9, 10),
                         (9, 11), (9, 6), (1, 0), (1, 12), (1, 13), (1, 5), (4, 9), (4, 14), (4, 15), (4, 0),
                         (6, 5), (6, 16), (6, 17), (6, 9), (12, 3), (12, 13), (12, 1), (12, 8), (14, 11),
                         (14, 15), (14, 4), (14, 3), (16, 8), (16, 17), (16, 6), (16, 11), (3, 12), (3, 0),
                         (3, 2), (3, 14), (8, 16), (8, 5), (8, 7), (8, 12), (11, 14), (11, 9), (11, 10),
                         (11, 16), (15, 2), (15, 4), (15, 14), (13, 7), (13, 1), (13, 12), (17, 10), (17, 6),
                         (17, 16), (2, 15), (2, 3), (2, 0), (7, 13), (7, 8), (7, 5), (10, 17), (10, 11), (10, 9)],
                        name = "tented hexagonal prism")

# dipyramid: two tetrahedra
# eg d=38,403
Dipyramid = Graph( [(0, 1), (0, 2), (0, 3), (3, 4), (3, 1), (3, 0), (3, 2), (1, 0), (1, 4), (1, 2),
                    (1, 3), (2, 0), (2, 1), (2, 4), (2, 3), (4, 3), (4, 1), (4, 2)],
                   name = "dipyramid")

# triangular prism + square pyramid
# eg d=41,131,179
PrismPyramid = Graph( [(0, 1), (0, 2), (0, 3), (0, 4), (2, 0), (2, 1), (2, 3), (2, 5), (3, 4), (3, 0),
                       (3, 2), (3, 6), (4, 6), (4, 0), (4, 3), (4, 1), (6, 4), (6, 3), (6, 5), (1, 5),
                       (1, 2), (1, 0), (1, 4), (5, 1), (5, 2), (5, 6)],
                      name = "triangular prism plus square pyramid")

# tetrahedron + square pyramid
# (includes a vertex of degree 5)
# eg d=143
TetraPyramid = Graph( [(0, 1), (0, 2), (0, 3), (0, 4), (4, 1), (4, 2), (4, 5), (4, 0), (2, 3), (2, 5),
                       (2, 1), (2, 0), (2, 4), (1, 0), (1, 4), (1, 2), (5, 3), (5, 2), (5, 4), (3, 5),
                       (3, 2), (3, 0)],
                      name = "tetrahedron plus square pyramid")

# Two cubes
# eg d=209
DoubleCube = Graph( [(1, 0), (0, 5), (6, 0), (2, 1), (3, 1), (4, 1), (2, 7), (2, 8), (3, 5), (7, 3),
                     (9, 3), (6, 4), (8, 4), (9, 4), (10, 5), (10, 6), (11, 7), (11, 8), (10, 9), (9, 11)],
                    name = "DoubleCube")

# Cube and two tents
# eg d=209
DoubleTentedCube = Graph([(1, 0), (3, 0), (5, 0), (2, 1), (3, 1), (4, 1), (6, 2), (7, 2), (2, 8),
                          (3, 8), (9, 3), (5, 4), (6, 4), (9, 4), (9, 5), (6, 10), (6, 11), (7, 8),
                          (7, 10), (11, 8), (9, 11), (10, 11)],
                         name = "DoubleTentedCube")


# 10 triangles + 2 squares
# (includes a vertex of degree 6)
# eg d=101
Unnamed_1 = Graph( [(0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (1, 5), (1, 6), (1, 4), (6, 5), (6, 7),
                    (6, 1), (6, 3), (5, 6), (5, 7), (5, 4), (5, 1), (3, 2), (3, 8), (3, 0), (3, 6),
                    (2, 3), (2, 4), (2, 8), (2, 0), (7, 4), (7, 5), (7, 6), (7, 8), (4, 5), (4, 7),
                    (4, 8), (4, 2), (4, 0), (4, 1), (8, 2), (8, 4), (8, 3), (8, 7)],
                   name = "unnamed #1")

# Two split cubes glued to a third cube in the middle
# eg d=149
Unnamed_2 = Graph([(1, 0), (4, 0), (5, 0), (2, 1), (3, 1), (3, 2), (2, 4), (2, 6), (5, 3), (3, 6),
                   (7, 3), (5, 4), (4, 14), (11, 5), (5, 14), (14, 6), (15, 6), (7, 8), (7, 11),
                   (7, 12), (7, 15), (9, 8), (8, 12), (9, 10), (9, 11), (11, 10), (12, 10), (13, 10),
                   (11, 13), (12, 15), (14, 13), (13, 15)],
                  name = "unnamed #2")

# 4 squares, 2 triangles
# eg d=471, 519, 615, 663, 759, 807, 903, 951

Unnamed_3 = Graph( [(1, 0), (2, 0), (0, 4), (2, 1), (3, 1), (2, 5), (2, 6), (3, 4), (3, 5), (4, 6), (5, 6)],
                  name = "unnamed #3")

Unknown = Graph(name="unknown")

# To add a new type of polyhedron, given one G of the new type:
#
# V = G.vertices()
# E = [tuple(V.index(ei) for ei in e[:2])  for e in G.edges()]
# assert G.is_isomorphic(Graph(E))
# E
#
# now add this code with the list of pairs E as first argument:
# NewPoly = Graph(..., name="new name")
#

all_polys = (Tetrahedron, Cube, Octahedron,
             TriangularPrism, SquarePyramid, HexagonalPrism,
             HexagonalCap, TruncatedTetrahedron, Cuboctahedron,
             SlicedCube, DoubleHexagonalPrism, HalfStar, TentedHexPrism,
             Dipyramid, PrismPyramid, TetraPyramid, DoubleCube, DoubleTentedCube,
             Unnamed_1, Unnamed_2, Unnamed_3, Unknown)

def poly_type(pol):
    for G in all_polys:
        if pol.is_isomorphic(G):
            return G.name()
    return "unknown"

def all_poly_types(pols):
    return [poly_type(pol) for pol in pols]

def poly_types(pols):
    n_unknown = sum(1 for pol in pols if poly_type(pol)=='unknown')
    return dict([(std_pol.name(), sum(1 for pol in pols if pol.is_isomorphic(std_pol)))
                 for std_pol in all_polys] + [('unknown', n_unknown)])

def vertical_halfline(x0,y0, h=1):
    return parametric_plot3d([x0,y0,lambda x:x], (0,h))

def vertical_semicircle(x1,y1,x2,y2):
    x0 = (x1+x2)/2
    y0 = (y1+y2)/2
    r = ((x1-x2)**2 + (y1-y2)**2).sqrt() /2
    return parametric_plot3d(
        [lambda t: x0+(x1-x0)*cos(t),
         lambda t: y0+(y1-y0)*cos(t),
         lambda t: r*sin(t)],
        (0,pi)
    )

def plot_edge(e, k, h=1):
    emb = next(e for e in k.embeddings(CC) if e(k.gen()).imag()>0)
    from utils import to_k
    if e[0].is_infinity():
        x0,y0 = emb(to_k(e[1],k))
        return vertical_halfline(x0,y0,h)
    if e[1].is_infinity():
        x0,y0 = emb(to_k(e[0],k))
        return vertical_halfline(x0,y0,h)
    (x1,y1), (x2,y2) = (emb(to_k(c,k)) for c in e)
    return vertical_semicircle(x1,y1,x2,y2)

def poly_display(P, k, h=0):
    from utils import cusp_from_string
    print(poly_type(P))
    V = [cusp_from_string(v, k) for v in P.vertices()]
    print(f"  vertices: {V}")
    E = [[cusp_from_string(v, k) for v in e[:2]] for e in P.edges()]
    print(f"  edges: {E}")
    F = [[[cusp_from_string(v, k) for v in e] for e in f] for f in P.faces()]
    print(f"  faces: {F}")
    if h:
        return sum([plot_edge(e, k, h) for e in E])
