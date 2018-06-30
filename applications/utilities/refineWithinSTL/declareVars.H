// Create dictionnary and read properties
IOdictionary immersedBoundaryProperties
(
    IOobject
    (
    "refineWithinSTLProperties",
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

//******** Surface instantiation ***********
autoPtr<triSurface>* surf;
autoPtr<pointField>* points;
autoPtr<triSurfaceSearch> querySurf(NULL);

surf = new autoPtr<triSurface>[nSurf];
points = new autoPtr<pointField>[nSurf];

// surfaces attribution 
for (int i=0;i<nSurf;i++)
{
    surf[i].reset(new triSurface(surfaceName[i]));
    points[i].reset(new pointField(surf[i]->points()));

    // Write information about the STL meshes
    surf[i]->writeStats(Info);
}


