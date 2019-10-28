// calc pos of particle and initialize the solid body motion field Uo

// find initial position after a restart
if(firstTime)
{
    pos = pos+timeAllreadyRun*vel;
}

// find new pos of TS
pos = pos+runTime.deltaT().value()*vel;

// find cells the particle is in and correct the body and Uo fields
if(mag(vel) > SMALL || firstTime) // re-locate cells only if object is moving
{
    firstTime = false;
    nrParticleCells=0;
    forAll(U,cellI)
    {
        dist = mag(U.mesh().C()[cellI]-pos);

        body[cellI]=0.;
        // Reset Up vector
        for (int i=0 ; i < 2 ; i++)
        {
            Up[i]=0.;
        }

        if(dist<rad)
        {
            particleCells[nrParticleCells] = cellI;
            nrParticleCells++;
            lever = U.mesh().C()[cellI]-pos;
            Up = vel + (omega^lever);
            body[cellI]=1.;
        }
        Uo[cellI]= Up;    
    }
}
Uo.correctBoundaryConditions();


// info about cells found
Info << "particle position "<< pos << endl;
Info << "particle is in "<< nrParticleCells << " , cells" << endl;
/*for(int i=0;i<nrParticleCells;i++)
{
    Info << "cell = " << particleCells[i] << endl;
}*/
