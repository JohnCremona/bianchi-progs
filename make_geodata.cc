// FILE MAKE_GEODATA.CC: creates a geodata file for a field if does
// not already exist, using Swan's algorithm + tessellation code

#include "geometry.h"
#include "swan.h"

#define VERBOSE 0 // verbose setting to use if not overridden locally

int main ()
{
  int verbose = VERBOSE;
  long D;
  cerr <<"Verbose? "; cin>>verbose; cerr<<endl;
  cerr << "Enter field (square-free d or absolute discriminant D): " << flush;  cin >> D;
  cerr << endl;

  long d = (D%4==0? D/4: D);
  long maxpnorm=100;
  Quad::field(d,maxpnorm);

  if (verbose)
    Quad::displayfield(cout);
  else
    cout << "Field Q(sqrt("<<-d<<"))\tdiscriminant = -"<<D<<endl;

  // See if precomputed geometry data already exists

  if (verbose)
    cout << "Reading data from geodata/geodata_"<<d<<".dat if possible, else creating from scratch..." <<endl;
  Quad::setup_geometry("geodata", verbose);
  if (verbose)
    cout << "done."<<endl;
  cout << Quad::SD.get_alphas().size()<<" alphas, "<<Quad::SD.get_sigmas().size()<<" sigmas"<<endl;
  cout << Quad::SD.T_faces.size() << " principal triangles" << endl;
  cout << Quad::SD.U_faces.size() << " singular triangles" << endl;
  cout << Quad::SD.Q_faces.size() << " squares" << endl;
  cout << Quad::SD.H_faces.size() << " hexagons" << endl;
}
