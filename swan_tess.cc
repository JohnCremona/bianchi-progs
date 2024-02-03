// FILE SWAN_TESS.CC: implementation of tessellation-finding from output of Swan's algorithm

#include "swan_utils.h"
#include "swan_alphas.h"
#include "swan_tess.h"

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
      if (!(unders.empty() || overs.empty()))
        continue;
      // now we have a genuine (external) face
      int nface = face.size();
      CuspList sface = face;
      std::sort(sface.begin(), sface.end(), RatQuad_cmp);
      if (std::find(sorted_faces.begin(), sorted_faces.end(), sface) != sorted_faces.end())
        continue; // duplicate

      sorted_faces.push_back(sface);
      P.faces.push_back(face);
      if (verbose)
        cout << " - found a new face with "<<nface<<" sides :"<<face<<endl;

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
          if (j==nface-1)
            {
              y = face[0]; // wrap around
              if (verbose>1)
                cout<<" removing directed edge from "<<x<<" to "<<y<<endl;
              ends[x].erase(std::remove(ends[x].begin(), ends[x].end(), y), ends[x].end());
              nde--;
              break;
            }
          for (int k=j+1; k<nface; k++)
            {
              y = face[k];
              if (std::find(ends[x].begin(), ends[x].end(), y) != ends[x].end())
                {
                  // delete y from ends[x]
                  if (verbose>1)
                    cout<<" removing directed edge from "<<x<<" to "<<y<<endl;
                  ends[x].erase(std::remove(ends[x].begin(), ends[x].end(), y), ends[x].end());
                  nde--;
                  // swap face[j+1] and face[k] if k>j+1
                  if (k!=j+1)
                    std::swap(face[k], face[j+1]);
                  break;
                }
            }
          if (verbose>1)
            cout<<"Now the face is "<<face<<endl;
        }
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
      if (flags[j]) continue;
      POLYHEDRON poly = principal_polyhedron(j, alphas, Plist, flags, verbose);
      polys.push_back(poly);
    }
  if (verbose)
    {
      cout<<polys.size()<<" principal polyhedra constructed"<<endl;
      for (const auto& P: polys)
        cout << poly_name(P) <<endl;
    }
  return polys;
}
