// this could be used to calculate the force on the IB

vector drag(0,0,0);
vector viscDrag(0,0,0);
vector pressureDrag(0,0,0);
Info << "calc drag force:"<< endl;

scalar rho=1;


volVectorField h1=rho*nu*fvc::laplacian(U);
h1.ref() *= U.mesh().V();
volVectorField h2=-fvc::grad(p);
h2.ref() *= U.mesh().V();

volVectorField h3=inside*f; // is it body * f or inside * f ---> it should actually be inside *f
h3.ref() *= U.mesh().V();

viscDrag = gSum(h1);
pressureDrag = gSum(h2);
drag = gSum(h3);

Info << "viscous drag =" << viscDrag << endl;
Info << "pressure drag =" << pressureDrag << endl;
Info << "constraint drag =" << drag << endl;

drag = vector(0.,0.,0.);
vector torque(0.,0.,0.);

logCount++;
if(logCount==logInterval)
{

    for (int i =0 ; i< nSurf ; i++)
    {
        for(int j=0 ; j<nBodyCells[i] ; j++)
        {
            drag += body[bodyCells[i]()[j]]*f[bodyCells[i]()[j]]*U.mesh().V()[bodyCells[i]()[j]];
            lever =  mesh.C()[bodyCells[i]()[j]]-centerOfRotation[i];
            torque += mesh.V()[bodyCells[i]()[j]] * (f[bodyCells[i]()[j]]^lever) ;//* body[bodyCells[i]()[j]] ;
        }
        reduce(drag, sumOp<vector>());
        reduce(torque, sumOp<vector>());
        Info << "---> Drag for " << surfaceName[i] << " is " << drag << endl;
        Info << "---> Torque is " << surfaceName[i] << " is " << torque << endl;

        dragLog[i] << runTime.timeName()  << " " << drag[0] << " "<< drag[1] << " " << drag[2] <<"\n" << std::flush;
        torqueLog[i] << runTime.timeName() << " " << torque[0] << " " << torque[1] << " " << torque[2] <<"\n" << std::flush;
    }

    logCount=0;
}



