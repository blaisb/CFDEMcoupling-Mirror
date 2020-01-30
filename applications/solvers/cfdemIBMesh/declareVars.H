// Create dictionnary and read properties
IOdictionary immersedBoundaryProperties
(
    IOobject
    (
    "immersedBoundaryProperties",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
    )
);

word objectType(immersedBoundaryProperties.lookup("objectType"));
if (objectType != "mesh")
{
    FatalError<< "Currently only mesh is supported for body type!" << abort(FatalError); 
}

Info << "Object type : " << objectType << " -- detected" << endl;
Info << "Immersed Boundary method subDictionaries have been created" << endl;

// *********** Mesh ************************
// Get Names of STL Meshes

fileNameList surfaceName(immersedBoundaryProperties.lookup("surfaceName"));
int nSurf = surfaceName.size();

Info << surfaceName << endl;

//*********** Motion ***********************
// Get Motion properties
scalarList omega(immersedBoundaryProperties.lookup("omega"));
scalarList volumeMesh(immersedBoundaryProperties.lookup("volumeMesh"));
vectorList axis(immersedBoundaryProperties.lookup("axis"));
vectorList centerOfRotation(immersedBoundaryProperties.lookup("centerOfRotation"));
vectorList linearVel(immersedBoundaryProperties.lookup("linearVelocity"));

boolList invariant(immersedBoundaryProperties.lookup("invariant"));
boolList enableVolumeMeshCorrection(immersedBoundaryProperties.lookup("enableVolumeMeshCorrection"));

//************ IBM parameters *************
scalar alpha(readScalar(immersedBoundaryProperties.lookup("alpha")));
int logInterval=int(readScalar(immersedBoundaryProperties.lookup("logInterval")));
int logCount=0;

Info << omega << endl;
Info << volumeMesh << endl;
Info << centerOfRotation << endl;
Info << invariant << endl;
Info << linearVel << endl;
Info << axis << endl;

//******** Surface instantiation ***********
autoPtr<triSurface>* surf;
autoPtr<pointField>* points;
autoPtr<triSurfaceSearch> querySurf(NULL);
int NrParticleCells=mesh.nPoints();

surf = new autoPtr<triSurface>[nSurf];
points = new autoPtr<pointField>[nSurf];

// Halo and body list for force post-processing and volume correction
autoPtr<labelList>* bodyCells(new autoPtr<labelList>[nSurf]);
autoPtr<labelList>* haloCells(new autoPtr<labelList>[nSurf]);
List<label> nBodyCells(nSurf);
List<label> nHaloCells(nSurf);

// surfaces attribution 
for (int i=0;i<nSurf;i++)
{
    surf[i].reset(new triSurface(surfaceName[i]));
    points[i].reset(new pointField(surf[i]->points()));

    // Write information about the STL meshes
    surf[i]->writeStats(Info);

    // Allocate the lists
    haloCells[i].reset(new labelList(NrParticleCells));
    bodyCells[i].reset(new labelList(NrParticleCells));
}

//********** Local calculations for motion ******************

scalar dt   = runTime.deltaT().value();
int nVel=0;
vector lever(0,0,0);

// ROTATION
// define omega in rad.s^-1
vectorList omegaV(nSurf);
for (int i=0; i<nSurf; i++)
{
    omegaV[i] = omega[i] * axis[i] * 2.*PI/60.;
}

// define rotation tensor
tensorList T(nSurf);
for (int i =0 ; i<nSurf; i++)
{
    vector omegaT(axis[i]*omega[i]/60.*2.*PI*dt);
    Info << "omegaT : " << omegaT << endl;
    vector n1(1,0,0);
    vector n2(Foam::cos(-omegaT[2]),Foam::sin(omegaT[2]),0);
    T[i]=rotationTensor(n1, n2);
    Info << "T : " << T[i] << endl;
}

// TRANSLATION
vectorList translation(nSurf);
for (int i=0 ; i<nSurf; i++)
{
translation[i] = linearVel[i] * dt;
}

//********* Tools for volume rescaling of the geometry *****************

// No need to allocate less than the number of mesh point
// saves one free parameter also
bool isFound=false;

// Volume correction utilities
scalar totalVol;
scalar haloVol;

// Create a cell volume field. This is necessary since volume of cell is not
// a volScalarField
volScalarField cellVolume
(
 IOobject
 (
  "cellVolume",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
 ),
 mesh,
 dimensionedScalar("zero",dimVolume,0.0)
);

cellVolume.ref() = mesh.V(); 


