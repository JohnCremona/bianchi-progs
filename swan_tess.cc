// FILE SWAN_TESS.CC: implementation of tessellation-finding from output of Swan's algorithm

#include "swan_utils.h"
#include "swan_alphas.h"
#include "swan_tess.h"

// Return all singular polyhedra
vector<POLYHEDRON>
singular_polyhedra(const CuspList& sigmas, const CuspList& alphas, int verbose)
{
  int i, j, k, n;
  RatQuad infty = RatQuad::infinity();
  vector<POLYHEDRON> polys;

  vector<pair<RatQuad,H3point>> sRlist; // list of (sigma, R) with R a fake corner at sigma

  // find and cache neighbours and fake corners for each finite sigma:
  map<RatQuad, CuspList, RatQuad_comparison> sigma_nbrs;
  for ( const auto& s : sigmas )
    {
      if (s.is_infinity())
        continue;
      auto alist = sorted_neighbours(s, alphas);
      sigma_nbrs[s] = alist;
      n = alist.size();
      for (i=0; i<n; i++)
        {
          RatQuad a = alist[i], b = alist[(i+1)%n];
          H3point R = bi_inter(a, b);
          sRlist.push_back({s,R});
        }
    }
  if (verbose)
    cout << " from " << sigmas.size()-1 << " finite sigmas we have " << sRlist.size() << " (sigma,R) pairs" << endl;

  // local function to test for being fundamental or oo:
  Quad x;
  auto is_fund = [alphas, &x](const RatQuad& a) {
    return a.is_infinity() || cusp_index_with_translation(a, alphas, x)>=0;
  };

  // Now process these {sigma,R} pairs; each gives a polyhedron with
  // one singular vertex (sigma) and we check off other pairs which
  // will give an SL2-congruent polyhedron
  vector<int> flags(sRlist.size(), 0);
  j = 0;
  for (const auto& sR : sRlist)
    {
      if (flags[j]) {j++; continue;}
      flags[j] = 1;
      vector<int> orbit = {j};
      RatQuad s = sR.first;
      H3point R = sR.second;
      if (verbose)
        cout<<" - processing #"<<j<<" (sigma,R) = ("<<s<<","<<R<<")"<<endl;
      CuspList alist = covering_hemispheres(R);
      POLYHEDRON P;
      // The vertices are oo, s, and the a in alist for which s is on S_a
      P.vertices.push_back(infty); // s will be added later, so now we
                                   // only have the principal vertices
      for (const auto& a : alist)
        {
          if ((a-s).norm()==radius_squared(a))
            P.vertices.push_back(a);
        }
      if (verbose)
        cout<<" - vertices are " << P.vertices <<endl;
      // the singular edges are {s,a} for all principal vertices a
      for (const auto& a : P.vertices)
        {
          P.edges.push_back({s,a});
          P.edges.push_back({a,s});
        }
      // the principal edges are {a1,a2} iff M_a1(a2) is fundamental;
      // we normalize the M_a so that M_a(s) is a reduced singular
      // point; then for each such a we can flag that the pair
      // {M_a(s), M_a(R)} is dealt with
      for (const auto& a : P.vertices)
        {
          mat22 M = Malpha(a, s, sigmas, i, k);
          for ( const auto& b : P.vertices)
            {
              if ((a!=b) && is_fund(M(b)))
                P.edges.push_back({a,b});
            }
          // check off flags k where {M(s), M(R)}=sRlist[k], given that M(s)=sigmas[i]
          pair<RatQuad,H3point> sR2 = {sigmas[i], M(R)};
          auto search = std::find(sRlist.begin(), sRlist.end(), sR2);
          assert (search!=sRlist.end());
          k = std::distance(sRlist.begin(), search);
          if ((k!=j) && (flags[k]==0))
            {
              flags[k] = 1;
              orbit.push_back(k);
            }
        }
      // add s to the vertices
      P.vertices.push_back(s);
      // fill in the faces
      fill_faces(P, verbose>1);
      // finished
      if (verbose)
        cout << "Constructed one singular " << poly_name(P)<<" from orbit "<<orbit<<endl;
      polys.push_back(P);
      j++;
    }
  if (verbose)
    {
      n = polys.size();
      if (n==1)
        cout<<" 1 singular polyhedron constructed"<<endl;
      else
        cout<<n<<" singular polyhedra constructed"<<endl;
      // for (const auto& P: polys)
      //   cout << poly_name(P) <<endl;
    }
  return polys;
}

// given vertices and edges, fill in faces:
void fill_faces(POLYHEDRON& P, int verbose)
{
  // create ends: map from v to list of w with vw an edge
  map<RatQuad, CuspList, RatQuad_comparison> ends;
  for (const auto& e : P.edges)
    ends[e.alpha()].push_back(e.beta());
  int nde=P.edges.size(); // number of directed edges; will decrease to 0
  if (verbose)
    {
      cout << "filling in faces of a polyhedron with "<<nde<<" directed edges" << endl;
      for (const auto& v : P.vertices)
        cout << v <<" --> " << ends[v] << endl;
    }

  RatQuad v, w, x, y;
  vector<CuspList> badstarts; // will hold any invalid v->w->x found
  CuspList face;

  while (nde) // continue until no directed edges are left
    {
      // find directed path v->w->x
      int found = 0;
      for ( const auto& vi: P.vertices)
        {
          if (found) break;
          for ( const auto& wi: ends[vi])
            {
              if (found) break;
              for ( const auto& xi: ends[wi])
                {
                  if (found) break;
                  if (xi!=vi)
                    {
                      v = vi; w = wi; x = xi;
                      face = {v,w,x};
                      found = std::find(badstarts.begin(), badstarts.end(), face) == badstarts.end();
                      if (found)
                        badstarts.push_back(face);
                    }
                }
            }
        }
      if (verbose)
        cout << " - trying "<<v<<"-->"<<w<<"-->"<<x<<endl;

      // Check that no vertices y are on the "wrong" side,
      // i.e. sign_im_cr(v, w, x, yi), and collect those with
      // sign_im_cr(v, w, x, yi)=0 as a new face:
      int good = 1;
      for ( const auto& yi : P.vertices )
        {
          if (!good)
            break;
          if ((yi==v)||(yi==w)||(yi==x))
            continue;
          int side = sign_im_cr(v, w, x, yi);
          if (side==0)
            face.push_back(yi);
          good = (side>=0);
        }
      if (!good)
        continue;

      // now we have a genuine oriented face

      // CuspList sface = face;
      // std::sort(sface.begin(), sface.end(), RatQuad_cmp);
      // if (std::find(sorted_faces.begin(), sorted_faces.end(), sface) != sorted_faces.end())
      //   continue; // duplicate
      // sorted_faces.push_back(sface);

      int nface = face.size();
      if (verbose)
        cout << " - found a new face with "<<nface<<" sides (unsorted):"<<face<<endl;

      // delete the relevant directed edges, starting with v->w->x
      if (verbose>1)
        cout<<" removing directed edge from "<<v<<" to "<<w<<endl;
      ends[v].erase(std::remove(ends[v].begin(), ends[v].end(), w), ends[v].end());
      nde--;
      if (verbose>1)
        cout<<" removing directed edge from "<<w<<" to "<<x<<endl;
      ends[w].erase(std::remove(ends[w].begin(), ends[w].end(), x), ends[w].end());
      nde--;
      // now face[2] is x; for j>=2 we reorder the face[k] for k>j if
      // necessary so that the edge from face[j] to face[j+1] is one
      // of the directed edges from face[j]:
      for (int j=2; j<nface; j++)
        {
          RatQuad x = face[j], y;
          if (verbose>1) cout<<"j="<<j<<", x="<<x<<endl;
          if (j==nface-1)
            {
              y = face[0]; // wrap around
              if (verbose>1) cout<<" k="<<0<<", y="<<y<<endl;
            }
          else
            {
              // find k>j such that x=face[j] --> y=face[k] is a directed edge
              for (int k=j+1; k<nface; k++)
                {
                  y = face[k];
                  if (verbose>1) cout<<" k="<<k<<", y="<<y<<endl;
                  if (std::find(ends[x].begin(), ends[x].end(), y) == ends[x].end())
                    {
                      if (verbose>1) cout<<" y IS NOT in ends[x] = "<<ends[x]<<endl;
                      continue;
                    }
                  if (verbose>1) cout<<" y IS in ends[x] = "<<ends[x]<<endl;
                  // now x-->y is a directed edge
                  if (k > j+1) // swap face[j+1] and face[k]
                    {
                      if (verbose>1)
                        cout<<"swapping vertices "<<k<<" ("<<y<<") and "<<j+1<<" ("<<face[j+1]<<")\n";
                      face[k]   = face[j+1];
                      face[j+1] = y;
                    }
                  break; // from the loop over k
                }
            }
          // delete y from ends[x]
          if (verbose>1)
            cout<<" removing directed edge from "<<x<<" to "<<y<<endl;
          ends[x].erase(std::remove(ends[x].begin(), ends[x].end(), y), ends[x].end());
          nde--;
          if (verbose>1)
            cout<<"Now the face is "<<face<<endl;
        }
      P.faces.push_back(face);
      if (verbose)
        cout << " we now have "<<P.faces.size()<<" faces; "<<nde << " directed edges remain" <<endl;
    }
}

// return a polyhedron from the j'th corner Plist[j]
POLYHEDRON
principal_polyhedron(int j, const CuspList& alphas, const H3pointList& Plist,
                     vector<int>& flags, int verbose)
{
  if (verbose)
    cout << " - using corner #"<<j<<" = "<<Plist[j]<<"...\n";
  H3point P = Plist[j];
  flags[j] = 1;
  vector<int> orbit = {j};

  POLYHEDRON poly;
  RatQuad infty = RatQuad::infinity();

  // Find all a with P on S_a; these & oo are the vertices of the polyhedron
  CuspList alist = covering_hemispheres(P);
  poly.vertices = alist;
  poly.vertices.insert(poly.vertices.begin(), infty);

  int nverts = poly.vertices.size();
  if (verbose)
    cout << " - polyhedron has "<< nverts << " vertices ("<<alist.size()
         <<" S_a go through P: "<<alist<<")"<<endl;

  // local function to test for being fundamental or oo:
  Quad x;
  auto is_fund = [alphas, &x](const RatQuad& a) {
    return a.is_infinity() || cusp_index_with_translation(a, alphas, x)>=0;
  };

  int i, k;
  for ( const auto& a : alist)
    {
      // First add edges {oo,a} and {a,oo} for a fundamental:
      if (is_fund(a))
        {
          poly.edges.push_back({infty,a});
          poly.edges.push_back({a,infty});
        }
      // Then add edges {a,b} for finite a,b, when M_a(b) is fundamental.
      // NB this is equivalent to M_b(a) fundamental, so {b,a} will also be added.
      mat22 M = Malpha(a, P, Plist, i, k);
      for ( const auto& b : alist)
        {
          if ((a!=b) && is_fund(M(b)))
            poly.edges.push_back({a,b});
        }
      // check off flags i, j where M(P)=Plist[i], u*M(P)=Plist[k]
      if ((i!=j) && (flags[i]==0))
        {
          flags[i] = 1;
          orbit.push_back(i);
        }
      // if ((k!=j) && (k!=i) && (flags[k]==0))
      //   {
      //     flags[k] = 1;
      //     orbit.push_back(k);
      //   }
    }
  int nedges = poly.edges.size()/2;
  int nfaces = 2+nedges-nverts; // Euler's formula!
  if (verbose)
    {
      cout << " - orbit " << orbit << " of size " << orbit.size() << endl;
      cout << " - polyhedron has (V,E,F)=("<<nverts<<","<<nedges<<","<<nfaces<<"):\n"; //<<poly << endl;
      cout << " - now filling in face data..."<<endl;
    }
  fill_faces(poly, verbose>1);
  if (verbose)
    {
      cout << "After filling in faces, polyhedron has "<<poly.faces.size()<<" faces:\n"<<poly.faces << endl;
      cout << "VEF(poly) = " << VEF(poly) <<endl;
      cout << "VEFx(poly) = " << VEFx(poly) <<endl;
      cout << " -- it is a "<<poly_name(poly)<<endl;
    }
  return poly;
}

// return a list of all principal polyhedra
vector<POLYHEDRON>
principal_polyhedra(const CuspList& alphas, int verbose)
{
  if (verbose)
    cout<<"Finding principal polyhedra from "<<alphas.size()<<" alphas..."<<flush;
  auto Plist = triple_intersections(alphas, verbose);
  int n = Plist.size();
  if (verbose)
    cout<<" number of corners is "<<n<<endl;
  vector<int>flags(n, 0);

  vector<POLYHEDRON> polys;
  for (int j=0; j<n; j++)
    {
      if (!flags[j])
        polys.push_back(principal_polyhedron(j, alphas, Plist, flags, verbose));
    }
  if (verbose)
    {
      cout<<polys.size()<<" principal polyhedra constructed"<<endl;
    }
  return polys;
}

// return a list of all polyhedra, principal and singular
vector<POLYHEDRON>
all_polyhedra(const CuspList& alphas, const CuspList& sigmas, int verbose)
{
  vector<POLYHEDRON> polys = principal_polyhedra(alphas, verbose);
  vector<POLYHEDRON> spolys = singular_polyhedra(alphas, sigmas, verbose);
  polys.insert(polys.end(), spolys.begin(), spolys.end());
  return polys;
}

// Given a polygon with either all vertices principal or just one
// singular; if there's a singular vertex, rotate to put it at the
// end; then apply M in SL(2,OK) taking the first two vertices to [a,oo]
// with a in alist (reduced fundamental alphas).  Set sing to 1 iff singular.

//#define DEBUG_NORMALISE
CuspList normalise_polygon( const CuspList& face, const CuspList& alphas, const CuspList& sigmas, int& sing)
{
#ifdef DEBUG_NORMALISE
  cout<<"normalising face "<<face<<endl;
#endif
  int n = face.size();
  CuspList newface = face;

  // local function to test for being finite singular:
  Quad x;
  auto is_sing = [sigmas, &x](const RatQuad& a) {
    return (!a.is_infinity()) && cusp_index_with_translation(a, sigmas, x)>=0;
  };

  // count the number of singular vertices:
  sing = std::count_if(face.begin(), face.end(), is_sing);

  // all vertices should be principal, expect for triangles where at most one can be singular:
  assert ((sing==0) || ((sing==1)&&(n==3)));
#ifdef DEBUG_NORMALISE
  cout<<" face "<<face<<" has "<<n<<" vertices of which "<<sing<<" are singular "<<endl;
#endif
  if (sing) // rotate until the (3rd and last) vertex is the singular one
    {
      if (is_sing(face[0])) // move 1st to end
        std::rotate(newface.begin(), newface.begin() + 1, newface.end());
      else
        if (is_sing(face[1])) // move 1st two to end
          std::rotate(newface.begin(), newface.begin() + 2, newface.end());
      assert (is_sing(newface[2]));
#ifdef DEBUG_NORMALISE
      cout<<"after moving singular vertex to the end, face is "<<newface<<endl;
#endif
    }

  int j, k; // not used here
  mat22 M = Malpha(newface[1], newface[0], alphas, j, k);
#ifdef DEBUG_NORMALISE
  cout<<"Applying matrix M = "<<M<<endl;
#endif
  std::transform(newface.begin(), newface.end(), newface.begin(),
                 [M](RatQuad a) {return M(a);});
  assert (newface[1]==RatQuad::infinity());
  assert (std::find(alphas.begin(), alphas.end(), newface[0]) != alphas.end());
#ifdef DEBUG_NORMALISE
  cout<<"after applying M, face is "<<newface<<endl;
#endif
  return newface;
}

// Given a polygon, move the first n (default 1) vertices to the end:
CuspList rotate_polygon( const CuspList& face, int n)
{
  CuspList newface = face;
  std::rotate(newface.begin(), newface.begin() + n, newface.end());
  return newface;
}

// Given a polygon, reverse the order of its vertices:
CuspList reverse_polygon( const CuspList& face)
{
  CuspList newface = face;
  std::reverse(newface.begin(), newface.end());
  return newface;
}

// Given a polygon, negate its vertices (i.e. apply a transformation in GL2 not SL2):
CuspList negate_polygon( const CuspList& face)
{
  CuspList newface = face; // assignment just to set the size
  std::transform(face.begin(), face.end(),
                 newface.begin(),
                 [](const RatQuad& v){return fundunit*v;});
  return newface;
}

// extract all oriented faces up to SL2-equivalence, returning lists of
// (0) principal triangles {a1, oo, a2} with a1 reduced fundamental
// (1) principal squares   {a1, oo, a2, a3} with a1 reduced fundamental
// (2) principal hexagons  {a1, oo, a2, a3, a4, a5, a6} with a1 reduced fundamental
// (3) singular triangles  {a, oo, s} with a reduced fundamental, s singular

// Omit the following:
// hexagon in a hexagonal cap (9,15,4,3,1)
// square in square pyramid (5,8,4,1,0), half-star (7,14,8,1,0), tetra+square pyr (6,11,6,1,0)
// as are redundant for 1-homology.
// Ditto the aaa-triangle in an aaas tetrahedron
vector<vector<CuspList>> get_faces( const vector<POLYHEDRON>& all_polys,
                                    const CuspList& alphas, const CuspList& sigmas,
                                    int verbose)
{
  // list of reduced polygons, up to GL2-equivalence
  vector<CuspList> aaa_triangles, aas_triangles, squares, hexagons;
  // extended list, including rotations and reflections
  vector<CuspList> aaa_triangle_copies, aas_triangle_copies, square_copies, hexagon_copies;
  int nTignored=0, nQignored=0, nHignored=0; // count faces ignored as redundant

  // local function to test for being finite singular:
  Quad x;
  auto is_sing = [sigmas, &x](const RatQuad& a) {
    return (!a.is_infinity()) && cusp_index_with_translation(a, sigmas, x)>=0;
  };

  for ( const auto& P : all_polys )
    {
      // count the number of singular vertices:
      int Psing = std::count_if(P.vertices.begin(), P.vertices.end(), is_sing);
      if (verbose)
        {
          cout<<"\nprocessing a "<<poly_name(P);
          if (Psing) cout<<" (an aaas tetrahedron)";
          cout<<endl;
        }

      auto vefx = VEFx(P);
      for ( const auto& F : P.faces )
      {
        int redundant = 0;
        int n = F.size();
        if (verbose)
          cout<<"\nprocessing face "<<F<<endl;
        if ((n==4) && (vefx[3]==1))
          {
            if (verbose)
              cout << " - ignoring the square face of a " << poly_name(P) << endl;
            nQignored++;
            redundant = 1;
          }
        if ((n==6) && (vefx[4]==1))
          {
            if (verbose)
              cout << " - ignoring the hexagon face of a " << poly_name(P) << endl;
            nHignored++;
            redundant = 1;
          }
        int Fsing;
        auto face = normalise_polygon(F, alphas, sigmas, Fsing);
        if (verbose)
          {
            cout<<" - after normalising, face "<< face;
            if (Fsing) cout<<" (singular)";
            cout<<endl;
          }
        if (Psing && !Fsing)
          {
            if (verbose)
              cout << " - ignoring the aaa face of an aaas " << poly_name(P) << endl;
            nTignored++;
            redundant = 1;
          }
        vector<CuspList>& face_copies = (n==3? (Fsing? aas_triangle_copies: aaa_triangle_copies) :
                                         (n==4? square_copies :
                                          hexagon_copies));
        vector<CuspList>& faces = (n==3? (Fsing? aas_triangles : aaa_triangles) :
                                   (n==4? squares :
                                    hexagons));
        if (std::find(face_copies.begin(), face_copies.end(), face) == face_copies.end())
          { // new face
            if (!redundant)
              {
                faces.push_back(face);
                if (verbose) cout << "Keeping face "<<face<<endl;
              }
            face_copies.push_back(face);
            auto rface = normalise_polygon(reverse_polygon(face), alphas, sigmas, Fsing);
            auto mface = normalise_polygon(negate_polygon(face), alphas, sigmas, Fsing);
            auto mrface = normalise_polygon(negate_polygon(rface), alphas, sigmas, Fsing);
            face_copies.push_back(rface);
            face_copies.push_back(mface);
            face_copies.push_back(mrface);
            if (!Fsing)  // apply (n-1) rotations to face and rface
              {
                for (int i=1; i<n; i++)
                  {
                    face = rotate_polygon(face);
                    face_copies.push_back(normalise_polygon(face, alphas, sigmas, Fsing));
                    rface = rotate_polygon(rface);
                    face_copies.push_back(normalise_polygon(rface, alphas, sigmas, Fsing));
                    mface = rotate_polygon(mface);
                    face_copies.push_back(normalise_polygon(mface, alphas, sigmas, Fsing));
                    mrface = rotate_polygon(mrface);
                    face_copies.push_back(normalise_polygon(mrface, alphas, sigmas, Fsing));
                  }
              }
          }
        else
          {
            if (verbose) cout <<"Ignoring face "<<face<<" as we already have a congruent one"<<endl;
          }
      }
    }
  // if (verbose)
    {
      cout<<"After processing "<<all_polys.size()<<" polyhedra, we have\n";
      cout<<aaa_triangles.size()<<" aaa-triangles ("<<nTignored<<" were ignored as redundant)\n";
      cout<<aas_triangles.size()<<" aas-triangles\n";
      cout<<squares.size()<<" squares ("<<nQignored<<" were ignored as redundant)\n";
      cout<<hexagons.size()<<" hexagons ("<<nHignored<<" were ignored as redundant)\n";
    }
  return {aaa_triangles, squares, hexagons, aas_triangles};
}

// Return j' so that M_a(a)=alphas[j'] where a=alphas[j]
int alpha_index_flip(int j, const CuspList& alphas)
{
  int jd;
  mat22 M = Malpha(alphas[j], alphas, jd);
  return jd;
}

string encode_int_list(char type, const vector<INT> data)
{
  ostringstream s;
  s << Quad::d << " " << type;
  for (const INT& a : data) s << " " << a;
  return s.str();
}

// NB POLYGON is typedef'd to pair<vector<int>, vector<Quad>>

// Return string for POLYGON representing an aaa-triangle, aas-triangle, quadrilateral or hexagon
string polygon_string(const POLYGON& P, int sing)
{
  vector<INT> data;
  for ( const int i : P.first)
    data.push_back(INT(i));
  for ( const Quad& u : P.second)
    {
      data.push_back(u.re());
      data.push_back(u.im());
    }
  int n = P.first.size(); // = 3, 4 or 6
  char type = (n==3? (sing? 'U' : 'T') : (n==4? 'Q' : 'H'));
  return encode_int_list(type, data);
}

// For any face, return the string which encodes it in the geodata files
string face_encode(const CuspList& face, const CuspList& alphas, const CuspList& sigmas)
{
  int sing;
  POLYGON P = make_polygon(face, alphas, sigmas, sing);
  return polygon_string(P, sing);
}

// Convert an actual polygon (aaa- or aas-triangle, quadrilateral or
// hexagon) as list of vertices to the POLYGON {{i,j,k},{u,...}}, setting
// sing to 1 for an aas-triangle, else to 0
POLYGON make_polygon(const CuspList& face, const CuspList& alphas, const CuspList& sigmas, int& sing)
{
  int n = face.size();
  sing = 0; // will be set to 1 if n==3 and it's an aas-triangle
  if (n==3) return make_triangle(face, alphas, sigmas, sing);
  if (n==4) return make_quadrilateral(face, alphas);
  if (n==6) return make_hexagon(face, alphas);
  assert ( 0 && "invalid face in make_polygon");
  return {{},{}};
}

// Convert an actual aaa- or aas-triangle as list of vertices [a,oo,b]
// or [a,oo,s] to a POLYGON {{i,j,k},{u}}, setting sing to 1 for an
// aas-triangle, else to 0
POLYGON make_triangle(const CuspList& T, const CuspList& alphas, const CuspList& sigmas, int& sing)
{
  assert (T[1]==RatQuad::infinity());
  Quad u, x;
  int i, j, k;
  i = cusp_index(T[0], alphas);
  assert (i>=0);
  mat22 M = Malpha(T[0]);
  j = cusp_index_with_translation(T[2], alphas, u);
  sing = (j<0);
  if (sing) // must be an aas-triangle
    {
      j = cusp_index_with_translation(T[2], sigmas, u);
      k = cusp_index_with_translation(M(T[2]), sigmas, x);
    }
  else
    {
      k = cusp_index_with_translation(M(T[2]), alphas, x);
    }
  assert (j>=0);
  assert (k>=0);
  return {{i,j,k}, {u}};
}

// Convert an actual quadrilateral as list of vertices [a,oo,b,c] to
// the POLYGON {{i,j,k},{x,y,z}}
POLYGON make_quadrilateral(const CuspList& Q, const CuspList& alphas)
{
  assert (Q[1]==RatQuad::infinity());
  int i,id,j,jd,k,kd,l;
  Quad x, y, z;
  i = cusp_index(Q[0], alphas);
  assert (i>=0);
  mat22 Mi = Malpha(Q[0], alphas, id); // id not used
  jd = cusp_index_with_translation(Q[2], alphas, z);
  assert (j>=0);
  j = alpha_index_flip(jd, alphas);
  l = cusp_index_with_translation(Mi(Q[3]), alphas, y);
  assert (l>=0);
  mat22 Mjd = Malpha(alphas[jd], alphas, id); // id not used
  kd = cusp_index_with_translation(Mjd(Q[3]-z), alphas, x);
  k = alpha_index_flip(kd, alphas);

  mat22 Mj = Malpha(alphas[j], alphas, id); // id not used
  mat22 Mk = Malpha(alphas[k], alphas, id); // id not used
  mat22 Ml = Malpha(alphas[l], alphas, id); // id not used
  RatQuad alpha1 = x + Mk.image_oo();
  RatQuad alpha2 = y + Ml.preimage_oo();
  mat22 M = Mi*mat22::Tmat(z)*Mj;
  assert (M(alpha1) == alpha2);
  return {{i,j,k,l}, {x,y,z}};
}

// Convert an actual hexagon as list of vertices [a_i, oo, a_j, b_2,
// gamma, b_1] to the POLYGON {{i,j,k,l,m,n},{u,x1,y1,x2,y2}}
POLYGON make_hexagon(const CuspList& H, const CuspList& alphas)
{
  assert (H[1]==RatQuad::infinity());
  int i,id,j,k,l,m,n;
  Quad u,x1,y1,x2,y2;
  i = cusp_index(H[0], alphas);
  assert (i>=0);
  mat22 Mi = Malpha(alphas[i], alphas, id); // id not used
  k = cusp_index_with_translation(Mi(H[5]), alphas, x1);
  assert (k>=0);
  mat22 Mk = Malpha(alphas[k], alphas, id); // id not used
  mat22 M = Mk*mat22::Tmat(-x1)*Mi;
  m = cusp_index_with_translation(M(H[4]), alphas, y1);
  j = cusp_index_with_translation(H[2], alphas, u);
  assert (j>=0);
  mat22 Mj = Malpha(alphas[j], alphas, id); // id not used
  M = Mj*mat22::Tmat(-u);
  l = cusp_index_with_translation(M(H[3]), alphas, x2);
  assert (l>=0);
  mat22 Ml = Malpha(alphas[l], alphas, id); // id not used
  M = Ml*mat22::Tmat(-x2)*Mj*mat22::Tmat(-u);
  n = cusp_index_with_translation(M(H[4]), alphas, y2);

  // check for consistency
  mat22
    M1 = Mi.inverse()*mat22::Tmat(x1)*Mk.inverse()*mat22::Tmat(y1),
    M2 = mat22::Tmat(u)*Mj.inverse()*mat22::Tmat(x2)*Ml.inverse()*mat22::Tmat(y2);
  RatQuad
    gamma1 = M1(alphas[m]),
    gamma2 = M2(alphas[n]);
  assert (gamma1==gamma2);

  return {{i,j,k,l,m,n}, {u,x1,y1,x2,y2}};
}

void output_faces( const vector<vector<CuspList>>& aaa_squ_hex_aas,
                   const CuspList& alphas, const CuspList& sigmas,
                   int to_file, int to_screen)
{
  if (Quad::is_Euclidean)
    {
      if (to_screen)
        cout << "No face output as field is Euclidean" << endl;
      return;
    }
  const auto& aaa_triangles = aaa_squ_hex_aas[0];
  const auto& squares       = aaa_squ_hex_aas[1];
  const auto& hexagons      = aaa_squ_hex_aas[2];
  const auto& aas_triangles = aaa_squ_hex_aas[3];
  vector<CuspList> faces;
  // We do not output the triangle {0,oo,1},
  // and for the fields 19, 43, 67, 163 also not
  // the triangles {oo,w/2,(w-1)/2} or {oo,w/2,(w+1)/2}
  CuspList tri0 = {RatQuad(0), RatQuad::infinity(), RatQuad(1)};
  CuspList tri1 = {RatQuad(Quad::w,TWO), RatQuad::infinity(), RatQuad(Quad::w-1,TWO)};
  CuspList tri2 = {RatQuad(Quad::w,TWO), RatQuad::infinity(), RatQuad(Quad::w+1,TWO)};
  long d = Quad::d;
  for (const auto& T : aaa_triangles)
    {
      CuspList T2 = reverse_polygon(T);
      if ((T==tri0) || (T2==tri0))
        {
          // if (to_screen) cout << "No T line for " << T << " as it is a universal one" << endl;
          continue;
        }
      if ((d==19)||(d==43)||(d==67)||(d==163))
        {
          if ((T==tri1) || (T==tri2) || (T2==tri1) || (T2==tri2))
            {
              // if (to_screen) cout << "No T line for " << T << " as it is a standard one for d=16,43,67,163" << endl;
              continue;
            }
        }
      faces.push_back(T);
    }
  // if (to_screen) cout<<"aaa-triangles to be output as T lines: "<<faces<<endl;
  faces.insert(faces.end(), aas_triangles.begin(), aas_triangles.end());
  faces.insert(faces.end(), squares.begin(), squares.end());
  faces.insert(faces.end(), hexagons.begin(), hexagons.end());
  if (to_screen) cout<<"all faces to be output: "<<faces<<endl;

  ofstream geodata;
  stringstream ss;
  if (to_file)
    {
      ss << "geodata_" << Quad::d << ".dat";
      geodata.open(ss.str().c_str(), ios_base::app);
    }
  int nlines=0;
  for (const auto& face: faces)
    {
      string s = face_encode(face, alphas, sigmas);
      nlines++;
      if (to_file)
        {
          geodata << s << endl;
        }
      if (to_screen)
        {
          cout << s << endl;
        }
    }
  if (to_file)
    geodata.close();
  if (to_screen)
    cout << nlines << " lines output" <<endl;
}