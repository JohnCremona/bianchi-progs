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
              if (verbose)
                cout << " - checking off (s,R) pair #"<<k<<endl;
            }
        }
      // add s to the vertices
      P.vertices.push_back(s);
      // fill in the faces
      fill_faces(P, verbose);
      // finished
      if (verbose)
        cout << "Constructed one singular polyhedron: "<<poly_name(P) << endl;
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
      for (const auto& P: polys)
        cout << poly_name(P) <<endl;
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
  vector<CuspList> starters; // store the {v,w,x} to be considered
  // find edges v->w->x
  RatQuad v, w, x, y;
  for ( const auto& vi: P.vertices)
    {
      for ( const auto& wi: ends[vi])
        {
          for ( const auto& xi: ends[wi])
            {
              if (xi!=vi)
                starters.push_back({vi,wi,xi});
            }
        }
    }
  if (verbose)
    {
      cout << " - found "<<starters.size()<<" v->w->x triples" << endl;
      // cout << starters <<endl;
    }
  vector<CuspList> sorted_faces;
  for (const auto& vwx : starters)
    {
      if (nde==0)
        break;

      v = vwx[0]; w = vwx[1]; x = vwx[2];
      // Make sure v->w and w->x are still in play
      if (std::find(ends[v].begin(), ends[v].end(), w) == ends[v].end())
        continue;
      if (std::find(ends[w].begin(), ends[w].end(), x) == ends[w].end())
        continue;
      if (verbose)
        cout << " - trying "<<v<<"-->"<<w<<"-->"<<x<<endl;
      CuspList face = vwx, unders, overs;
      for ( const auto& yi : P.vertices )
        {
          if ((yi==v)||(yi==w)||(yi==x))
            continue;
          int side = sign_im_cr(v, w, x, yi);
          if (side>0)
            overs.push_back(yi);
          else
            {
              if (side<0)
                unders.push_back(yi);
              else
                face.push_back(yi);
            }
        }
      if (!(unders.empty())) // || overs.empty()))
        continue;
      // now we have a genuine (external) face
      int nface = face.size();
      CuspList sface = face;
      std::sort(sface.begin(), sface.end(), RatQuad_cmp);
      if (std::find(sorted_faces.begin(), sorted_faces.end(), sface) != sorted_faces.end())
        continue; // duplicate

      sorted_faces.push_back(sface);
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
          if (verbose)
            cout << " - checking off corner #"<<i<<endl;
        }
      if ((k!=j) && (k!=i) && (flags[k]==0))
        {
          flags[k] = 1;
          if (verbose)
            cout << " - checking off corner #"<<k<<endl;
        }
    }
  int nedges = poly.edges.size()/2;
  int nfaces = 2+nedges-nverts; // Euler's formula!
  if (verbose)
    {
      cout << " - polyhedron has (V,E,F)=("<<nverts<<","<<nedges<<","<<nfaces<<"):\n"; //<<poly << endl;
      cout << " - now filling in face data..."<<endl;
    }
  fill_faces(poly, verbose);
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
    cout<<"Finding principal polyhedra..."<<flush;
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
      for (const auto& P: polys)
        cout << poly_name(P) <<endl;
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
// end; then apply M in SL(2,OK) taking the first two vertices to oo
// and a in alist (reduced fundamental alphas).  Set sing to 1 iff singular.
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
  mat22 M = Malpha(newface[0], newface[1], alphas, j, k);
#ifdef DEBUG_NORMALISE
  cout<<"Applying matrix M = "<<M<<endl;
#endif
  std::transform(newface.begin(), newface.end(), newface.begin(),
                 [M](RatQuad a) {return M(a);});
  assert (newface[0]==RatQuad::infinity());
  assert (std::find(alphas.begin(), alphas.end(), newface[1]) != alphas.end());
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

// extract all oriented faces up to SL2-equivalence, returning lists of
// (0) principal triangles {oo, a1, a2} with a1 reduced fundamental
// (1) principal squares   {oo, a1, a2, a3} with a1 reduced fundamental
// (2) principal hexagons  {oo, a1, a2, a3, a4, a5, a6} with a1 reduced fundamental
// (3) singular triangles  {oo, a, s} with a reduced fundamental, s singular
vector<vector<CuspList>> get_faces( const vector<POLYHEDRON>& all_polys,
                                    const CuspList& alphas, const CuspList& sigmas,
                                    int verbose)
{
  // list of reduced polygons, up to GL2-equivalence
  vector<CuspList> aaa_triangles, aas_triangles, squares, hexagons;
  // extended list, including rotations and reflections
  vector<CuspList> aaa_triangle_copies, aas_triangle_copies, square_copies, hexagon_copies;
  for ( const auto& P : all_polys )
    for ( const auto& F : P.faces )
      {
        if (verbose)
          cout<<"\nprocessing face "<<F<<endl;
        int sing;
        auto face = normalise_polygon(F, alphas, sigmas, sing);
        if (verbose)
          {
            cout<<" - after normalising, face "<< face;
            if (sing) cout<<" (singular)";
            cout<<endl;
          }
        int n = face.size();
        vector<CuspList>& face_copies = (n==3? (sing? aas_triangle_copies: aaa_triangle_copies) :
                                         (n==4? square_copies :
                                          hexagon_copies));
        vector<CuspList>& faces = (n==3? (sing? aas_triangles : aaa_triangles) :
                                   (n==4? squares :
                                    hexagons));
        if (std::find(face_copies.begin(), face_copies.end(), face) == face_copies.end())
          { // new face
            faces.push_back(face);
            face_copies.push_back(face);
            auto rface = normalise_polygon(reverse_polygon(face), alphas, sigmas, sing);
            face_copies.push_back(rface);
            if (!sing)  // apply (n-1) rotations to face and rface
              {
                for (int i=1; i<n; i++)
                  {
                    face = rotate_polygon(face);
                    face_copies.push_back(normalise_polygon(face, alphas, sigmas, sing));
                    rface = rotate_polygon(rface);
                    face_copies.push_back(normalise_polygon(rface, alphas, sigmas, sing));
                  }
              }
          }
      }
  if (verbose)
    {
      cout<<"After processing "<<all_polys.size()<<" polyhedra, we have\n";
      cout<<aaa_triangles.size()<<" aaa-triangles\n";
      cout<<aas_triangles.size()<<" aas-triangles\n";
      cout<<squares.size()<<" squares\n";
      cout<<hexagons.size()<<" hexagons\n";
    }
  return {aaa_triangles, squares, hexagons, aas_triangles};
}

/////////////////////////////////////////////////////////////////////////////////
#if(0) // obsolete functions giving singular polyhedra decomposed into tetrahedra

// return a tetrahedron (list of triangles) from a list of 4 vertices,
// arranged so that each is the first vertex of one triangle
POLYHEDRON tetrahedron(const CuspList& V)
{
  assert (V.size()==4);
  POLYHEDRON P;
  P.vertices = V;
  P.edges = {
    {V[0],V[1]}, {V[1],V[2]}, {V[2],V[0]},
    {V[1],V[3]}, {V[3],V[2]}, {V[2],V[1]},
    {V[2],V[3]}, {V[3],V[0]}, {V[0],V[2]},
    {V[3],V[1]}, {V[1],V[0]}, {V[0],V[3]}
  };
  P.faces = {
    {V[0], V[1], V[2]},
    {V[1], V[3], V[2]},
    {V[2], V[3], V[0]},
    {V[3], V[1], V[0]}
  };
  return P;
}

// return a list of tetrahedra (i.e. lists of 4 cusps (oo, sigmas[j],
// a1, a2) with a1,a2 fundamental); as a side-effect set flags[j]=1
// and flags[j']=1 where M_a_i(sigma[j])=sigma[j'] for i=1,2, for
// each.
vector<POLYHEDRON>
singular_tetrahedra(int j, const CuspList& sigmas, const CuspList& alphas, vector<int>& flags, int verbose)
{
  flags[j] = 1;
  RatQuad sigma = sigmas[j], infty = RatQuad::infinity();
  if (verbose)
    cout << " - finding singular tetrahedra for sigmas["<<j<<"]="<<sigma<<"..."<<flush;
  auto alist = sorted_neighbours(sigma, alphas);
  int n = alist.size();
  if (verbose) cout<<" constructing "<<n<<endl;
  vector<POLYHEDRON> tetras;
  tetras.reserve(n);

  int i;
  for (i=0; i<n; i++)
    {
      CuspList V = {infty, sigma, alist[(i==0? n-1 : i-1)], alist[i]};
      tetras.push_back(tetrahedron(V));
    }

  for ( const auto& a : alist)
    {
      Malpha(a, sigma, sigmas, i);
      if ((i!=j) && (flags[i]==0))
        {
          flags[i] = 1;
          if (verbose)
            cout << " - checking off sigmas["<<i<<"]="<<sigmas[i]<<endl;
        }
    }
  return tetras;
}

// return a list of all singular tetrahedra
vector<POLYHEDRON>
singular_tetrahedra(const CuspList& sigmas, const CuspList& alphas, int verbose)
{
  int n = sigmas.size();
  if (verbose)
    cout << "Finding singular tetrahedra for "<<n-1<<" finite singular points" << endl;
  vector<int> flags(n, 0); // will set to 1 as each sigma is accounted for
  flags[0] = 1; // since sigmas[0]=oo
  vector<POLYHEDRON> tetras;
  for (int j=0; j<n; j++)
    {
      if (flags[j]) continue;
      vector<POLYHEDRON> tetras1 = singular_tetrahedra(j, sigmas, alphas, flags, verbose);
      if (verbose)
        {
          cout << "sigma = "<<sigmas[j]<<":\n";
          for ( const auto& tetra: tetras1 )
            cout << " "<<tetra << endl;
        }
      tetras.insert(tetras.end(), tetras1.begin(), tetras1.end());
    }
  return tetras;
}

// Return the polyhedron with vertices oo, s=sigmas[j], a1, a2, ...,
// an with ai fundamental such that s is on S_ai.

// As a side-effect set flags[j]=1 and flags[j']=1 where
// M_a_i(sigma[j])=sigma[j'] for each i.

// The polyhedron is an n-dipyramid with poles at oo and sigmas[j]; it
// could be decomposed into n tetrahedra with vertices oo,s,ai,a{i+1}.

POLYHEDRON
singular_polyhedron(int j, const CuspList& sigmas, const CuspList& alphas, vector<int>& flags, int verbose)
{
  flags[j] = 1;
  RatQuad sigma = sigmas[j], infty = RatQuad::infinity();
  if (verbose)
    cout << " - finding singular polyhedron for sigmas["<<j<<"]="<<sigma<<"..."<<flush;
  auto alist = sorted_neighbours(sigma, alphas);
  int n = alist.size();

  POLYHEDRON P;

  // set vertices
  P.vertices = {infty, sigma};
  P.vertices.insert(P.vertices.end(), alist.begin(), alist.end());

  // set edges and faces
  int i, k;
  for (i=0; i<n; i++)
    {
      RatQuad a = alist[i], b = alist[(i+1)%n];
      P.edges.push_back({infty,a});
      P.edges.push_back({a,infty});
      P.edges.push_back({sigma,a});
      P.edges.push_back({a,sigma});
      P.edges.push_back({a,b});
      P.edges.push_back({b,a});
      P.faces.push_back({infty, a, b});
      P.faces.push_back({b, a, sigma});
    }

  // check off this sigma and others in its orbit
  for ( const auto& a : alist)
    {
      Malpha(a, sigma, sigmas, i, k);
      if ((i!=j) && (flags[i]==0))
        {
          flags[i] = 1;
          if (verbose)
            cout << " - checking off sigmas["<<i<<"]="<<sigmas[i]<<endl;
        }
      if ((k!=j) && (k!=i) && (flags[k]==0))
        {
          flags[k] = 1;
          if (verbose)
            cout << " - checking off sigmas["<<k<<"]="<<sigmas[k]<<endl;
        }
    }
  assert ((int)P.vertices.size()==n+2);
  assert ((int)P.edges.size()==6*n);
  assert ((int)P.faces.size()==2*n);
  return P;
}

#endif
