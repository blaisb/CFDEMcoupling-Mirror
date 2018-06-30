/* ----------------------------------------------------------------------
    Drag accumulator using various drag models
    Created by Bruno Blais
    URPEI, Polytechnique Montreal
    Bruno.Blais@polymtl.ca
    
    Based on the work of:
     
    LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
    Transfer Simulations

    LIGGGHTS is part of the CFDEMproject
    www.liggghts.com | www.cfdem.com
    
    Christoph Kloss, christoph.kloss@cfdem.com
    Copyright 2009-2012 JKU Linz
    Copyright 2012-     DCS Computing GmbH, Linz

    LIGGGHTS is based on LAMMPS
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    This software is distributed under the GNU General Public License.

    See the README file in the top-level directory.

    Functionnality :
    The drag accumulator implements numerous drag model that can then be
    calculated from the DEM side
        RongDrag Model
        RongDrag relaxed --> Rong Drag model but with a saturation on the
                                void fraction dependency
        singleDrag      ---> drag for a single particle. No void fraction
        oneWayDrag      ---> Only DEM feels the drag, no two-way coupling

------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "comm.h"
#include "math.h"
#include "math_const.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "fix_cfd_coupling_force_accumulator.h"
#include "fix_property_atom.h"

#define SMALL 1e-6

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingForceAccumulator::FixCfdCouplingForceAccumulator(LAMMPS *lmp, int narg, char **arg) :
    FixCfdCouplingForce(lmp,narg,arg),
    fix_Ksl_(0),
    fix_uf_(0),
    fix_dragAcc_(0)
{

    int iarg = 3;

    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;
        if(strcmp(arg[iarg],"RongDragRelaxed")==0)
        {
            if (narg<iarg+3)
            {
		error->fix_error(FLERR,this,"Insufficient number of argument for RongDragRelaxed fix\n arg1: rho fluid arg2: mu fluid  arg3: minimum void fraction\n");
            }
            dragModel_="RongDragRelaxed";
            iarg++;
            rhof_= atof(arg[iarg]);
            iarg++;
            muf_= atof(arg[iarg]);
            iarg++;
            voidMin_=atof(arg[iarg]);
            fprintf(screen,"Minimal void fraction is : %5.5e \n",voidMin_);
	}
        else if(strcmp(arg[iarg],"RongDrag")==0 )
	{
	    if (narg<iarg+2)
	    {
		error->fix_error(FLERR,this,"Insufficient number of argument for RongDrag fix\n arg1: rho fluid arg2: mu fluid \n");
	    }
	    dragModel_="RongDrag";
	
	    iarg++;
	    rhof_= atof(arg[iarg]);
	    iarg++;
	    muf_= atof(arg[iarg]);
	}
        else if(strcmp(arg[iarg],"singleDrag")==0)
        {
	    if (narg<iarg+2)
	    {
		error->fix_error(FLERR,this,"Insufficient number of argument for singleDrag fix\n arg1: rho fluid arg2: mu fluid \n");
	    }
	    dragModel_="singleDrag";

	    iarg++;
	    rhof_= atof(arg[iarg]);
	    iarg++;
	    muf_= atof(arg[iarg]);
        }
        else if(strcmp(arg[iarg],"oneWayDrag")==0)
        {
 	    if (narg<iarg+2)
	    {
		error->fix_error(FLERR,this,"Insufficient number of argument for oneWayDrag fix\n arg1: rho fluid arg2: mu fluid \n");
	    }
	    dragModel_="oneWayDrag";

	    iarg++;
	    rhof_= atof(arg[iarg]);
	    iarg++;
	    muf_= atof(arg[iarg]);
        }
	else
	{
	    error->fix_error
	    (FLERR,this,"The accumulator fix has not been implemented for drag formulations other than:\n -RongDrag\n -RongDragRelaxed\n-singleDrag\n-oneWayDrag\n");
	}
    }
    
    nevery = 1;
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingForceAccumulator::~FixCfdCouplingForceAccumulator()
{
}

/* ---------------------------------------------------------------------- */
int FixCfdCouplingForceAccumulator::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceAccumulator::post_create()
{
    // do mother class init w/o dragforce
    FixCfdCouplingForce::post_create();

    // register Ksl
    if(!fix_Ksl_)
    {
        const char* fixarg[9];
        fixarg[0]="Ksl";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="Ksl";
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_Ksl_ = modify->add_fix_property_atom(9,(char**)fixarg,style);
    }

    // register uf
    if(!fix_uf_)
    {
        const char* fixarg[11];
        fixarg[0]="uf";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="uf";
        fixarg[4]="vector"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fix_uf_ = modify->add_fix_property_atom(11,(char**)fixarg,style);
    }

    // register dragAcc
    if(!fix_dragAcc_)
    {
        const char* fixarg[11];
        fixarg[0]="dragAcc";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="dragAcc";
        fixarg[4]="vector"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fix_dragAcc_ = modify->add_fix_property_atom(11,(char**)fixarg,style);
    }

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceAccumulator::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_Ksl_) modify->delete_fix("Ksl");
    if(unfixflag && fix_uf_) modify->delete_fix("uf");
    if(unfixflag && fix_dragAcc_) modify->delete_fix("dragAcc");
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceAccumulator::init()
{
    FixCfdCouplingForce::init();

    int nlocal = atom->nlocal;
    double **dragAcc = fix_dragAcc_->array_atom;


    // Additional force accumulation vector to push to OF
    fix_coupling_->add_push_property("dragAcc","vector-atom");

    // values to come from OF
    fix_coupling_->add_pull_property("Ksl","scalar-atom");
    fix_coupling_->add_pull_property("uf","vector-atom");
   
    //Reinitialise dragAcc
    for (int i = 0 ; i < nlocal ; i++)
    {
	dragAcc[i][0]=0.;
	dragAcc[i][1]=0.;
	dragAcc[i][2]=0.;
    }

    deltaT_ = update->dt * force->ftm2v;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceAccumulator::post_force(int vflag)
{

    double **dragAcc = fix_dragAcc_->array_atom;
    int nlocal = atom->nlocal;
    // If first time step of coupling erase the drag accumulator 
    // so that the accumulation can start from the current iteration
    // and then transfer to CFD

    if (update->ntimestep==(update->firststep+1))
    {
        //std::cout <<"Clearing drag accumulator" << std::endl;
        for (int i =0 ; i<nlocal ; i++)
        {
            vectorZeroize3D(dragAcc[i]);
        }
    }

    //Rong drag implementation
    if (strcmp("RongDrag",dragModel_)==0)
    {
        RongDrag();
    }
    //Rong drag with void fraction saturation
    else if (strcmp("RongDragRelaxed",dragModel_)==0)
    {
        RongDragRelaxed();
    }
    //Rong drag implementation
    else if (strcmp("singleDrag",dragModel_)==0)
    {
        singleDrag();
    }
    //One way coupling implementation
    else if (strcmp("oneWayDrag",dragModel_)==0)
    {
        oneWayDrag();
    }

}
/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceAccumulator::end_of_step()
{
    // Empty object for now
}


void FixCfdCouplingForceAccumulator::RongDrag()
{
    //Implementation of the RongDrag() calculation
    //double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double *r = atom->radius;
    //double *rmass = atom->rmass;
    //double *mass = atom->mass;

    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double *Ksl = fix_Ksl_->vector_atom; // Ksl is the (1.-voidfraction) in the accumulator, this makes initialization easier
    double **uf = fix_uf_->array_atom;
    double **dragAcc = fix_dragAcc_->array_atom;
    double **dragforce = fix_dragforce_->array_atom;

    double Rep, Cd, Cd0, beta, voidfraction;
    double frc[3];

    double magUr; // norm of relative velocity

    for (int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {

            //if (rmass) m = rmass[i];
            //else	m=mass[type[i]];

            magUr = sqrt(
                    (uf[i][0] - v[i][0]) * (uf[i][0] - v[i][0]) +
                    (uf[i][1] - v[i][1]) * (uf[i][1] - v[i][1]) +
                    (uf[i][2] - v[i][2]) * (uf[i][2] - v[i][2]) 
                    );

            voidfraction= MAX(1.-Ksl[i],0.);

            //Calculate particulate Reynolds
            Rep = 2.*r[i] * voidfraction*magUr*rhof_/muf_ + SMALL;

            //Calculate Cd0
            Cd0 =(0.63+4.8/(sqrt(Rep))) * (0.63+4.8/(sqrt(Rep))) ;

            beta = 2.65 * (voidfraction + 1.) -
                (5.3-3.5*voidfraction)*voidfraction*voidfraction
                *exp(-0.5*pow((1.5-log10(Rep)),2));

            Cd = 0.5 * rhof_ * MathConst::MY_PI * r[i] * r[i] * Cd0 * magUr * pow(voidfraction,2.-beta); 
            for(int dirI=0;dirI<3;dirI++)
            {
                //Drag is calculated explicitly in the DEM
                frc[dirI] = Cd * (uf[i][dirI] - v[i][dirI]);
            }

          // add force to accumulation vector
          vectorAdd3D(dragAcc[i],frc,dragAcc[i]);
	
          // add force due to drag to force balance on particle
          vectorAdd3D(f[i],frc,f[i]);

          // add other forces (pressure, Archimedes, etc.)
          vectorAdd3D(f[i],dragforce[i],f[i]);
     
          // add up forces for post-proc
          vectorAdd3D(dragforce_total,frc,dragforce_total);
      }
    }
}


void FixCfdCouplingForceAccumulator::RongDragRelaxed()
{
    /**********************************************************************************
        Implementation of a drag calculation that can take into account a staturation
        of the void fraction for stability reasons
    ***********************************************************************************/
    double **v = atom->v;
    double **f = atom->f;
    double *r = atom->radius;

    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double *Ksl = fix_Ksl_->vector_atom; // Ksl is the (1.-voidfraction) in the accumulator, this makes initialization easier
    double **uf = fix_uf_->array_atom;
    double **dragAcc = fix_dragAcc_->array_atom;
    double **dragforce = fix_dragforce_->array_atom;

    double Rep, Cd, Cd0, beta, voidfraction;
    double frc[3];

    double magUr; // norm of relative velocity

    for (int i = 0; i < nlocal; i++)
    {
          if (mask[i] & groupbit)
        {

            magUr = sqrt(
                    (uf[i][0] - v[i][0]) * (uf[i][0] - v[i][0]) +
                    (uf[i][1] - v[i][1]) * (uf[i][1] - v[i][1]) +
                    (uf[i][2] - v[i][2]) * (uf[i][2] - v[i][2]) 
                    );

            voidfraction= MAX(1.-Ksl[i],voidMin_);

            //Calculate particulate Reynolds
            Rep = 2.*r[i] * voidfraction*magUr*rhof_/muf_ + SMALL;

            //Calculate Cd0
            Cd0 =(0.63+4.8/(sqrt(Rep))) * (0.63+4.8/(sqrt(Rep))) ;

            beta = 2.65 * (voidfraction + 1.) -
                (5.3-3.5*voidfraction)*voidfraction*voidfraction
                *exp(-0.5*pow((1.5-log10(Rep)),2));

            Cd = 0.5 * rhof_ * MathConst::MY_PI * r[i] * r[i] * Cd0 * magUr * pow(voidfraction,2.-beta); 
            for(int dirI=0;dirI<3;dirI++)
            {
                //Drag is calculated explicitly in the DEM
                frc[dirI] = Cd * (uf[i][dirI] - v[i][dirI]);
            }

            // add force to accumulation vector
            vectorAdd3D(dragAcc[i],frc,dragAcc[i]);

            // add force due to drag to force balance on particle
            vectorAdd3D(f[i],frc,f[i]);

            // add other forces (pressure, Archimedes, etc.)
            vectorAdd3D(f[i],dragforce[i],f[i]);

            // add up forces for post-proc
            vectorAdd3D(dragforce_total,frc,dragforce_total);
      }
  }
}



void FixCfdCouplingForceAccumulator::singleDrag()
{
    /**************************************************************************
        Implementation of a drag calculation that does not use void fraction 
        and does not account for hindering in the drag calculation
    **************************************************************************/
    double **v = atom->v;
    double **f = atom->f;
    double *r = atom->radius;

    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double **uf = fix_uf_->array_atom;
    double **dragAcc = fix_dragAcc_->array_atom;
    double **dragforce = fix_dragforce_->array_atom;

    double Rep, Cd, Cd0, sRep;
    double frc[3];

    double magUr; // norm of relative velocity

    for (int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {

            magUr = sqrt(
                    (uf[i][0] - v[i][0]) * (uf[i][0] - v[i][0]) +
                    (uf[i][1] - v[i][1]) * (uf[i][1] - v[i][1]) +
                    (uf[i][2] - v[i][2]) * (uf[i][2] - v[i][2]) 
                    );

            //Calculate particulate Reynolds
            Rep = 2.*r[i] *magUr*rhof_/muf_ + SMALL;
            sRep=sqrt(Rep);

            //Calculate Cd0
            Cd0 =(0.63+4.8/(sRep)) * (0.63+4.8/(sRep)) ;

            Cd = 0.5 * rhof_ * MathConst::MY_PI * r[i] * r[i] * Cd0 * magUr; 
            for(int dirI=0;dirI<3;dirI++)
            {
                //Drag is calculated explicitly in the DEM
                frc[dirI] = Cd * (uf[i][dirI] - v[i][dirI]);
            }

            // add force to accumulation vector
            vectorAdd3D(dragAcc[i],frc,dragAcc[i]);

            // add force due to drag 
            vectorAdd3D(f[i],frc,f[i]);

            // add other forces (pressure, Archimedes, etc.)
            vectorAdd3D(f[i],dragforce[i],f[i]);

            // add up forces for post-proc
            vectorAdd3D(dragforce_total,frc,dragforce_total);
      }
    }
}

void FixCfdCouplingForceAccumulator::oneWayDrag()
{
    //Implementztion of a drag calculation that does not use void fraction hindering
    //double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double *r = atom->radius;

    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double **uf = fix_uf_->array_atom;
    double **dragforce = fix_dragforce_->array_atom;

    double Rep, Cd, Cd0, sRep;
    double frc[3];

    double magUr; // norm of relative velocity

    for (int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            magUr = sqrt(
                    (uf[i][0] - v[i][0]) * (uf[i][0] - v[i][0]) +
                    (uf[i][1] - v[i][1]) * (uf[i][1] - v[i][1]) +
                    (uf[i][2] - v[i][2]) * (uf[i][2] - v[i][2]) 
                    );

            //Calculate particulate Reynolds
            Rep = 2.*r[i] *magUr*rhof_/muf_ + SMALL;
            sRep=sqrt(Rep);

            //Calculate Cd0
            Cd0 =(0.63+4.8/(sRep)) * (0.63+4.8/(sRep)) ;

            Cd = 0.5 * rhof_ * MathConst::MY_PI * r[i] * r[i] * Cd0 * magUr; 
            for(int dirI=0;dirI<3;dirI++)
            {
                //Drag is calculated explicitly in the DEM
                frc[dirI] = Cd * (uf[i][dirI] - v[i][dirI]);
            }

            // add force due to drag 
            vectorAdd3D(f[i],frc,f[i]);

            // add other forces (pressure, Archimedes, etc.)
            vectorAdd3D(f[i],dragforce[i],f[i]);

            // add up forces for post-proc
            vectorAdd3D(dragforce_total,frc,dragforce_total);
      }
    }
}

