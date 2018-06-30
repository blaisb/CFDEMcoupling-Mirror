/* ----------------------------------------------------------------------
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
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(couple/cfd/force/accumulator,FixCfdCouplingForceAccumulator)

#else

#ifndef LMP_FIX_CFD_COUPLING_FORCE_ACCUMULATOR_H
#define LMP_FIX_CFD_COUPLING_FORCE_ACCUMULATOR_H

#include <stdlib.h>
#include "string.h"
#include "fix_cfd_coupling_force.h"


namespace LAMMPS_NS {

class FixCfdCouplingForceAccumulator : public FixCfdCouplingForce  {
 public:
  FixCfdCouplingForceAccumulator(class LAMMPS *, int, char **);
  ~FixCfdCouplingForceAccumulator();
  void post_create();
  void pre_delete(bool unfixflag);

  int setmask();
  virtual void init();
  void post_force(int);
  void end_of_step();

  void RongDrag();
  void RongDragRelaxed();
  void singleDrag();
  void oneWayDrag();

 protected:
  double deltaT_;
  char const * dragModel_;
  double muf_;
  double rhof_;
  double voidMin_; //Minimal void fraction that is used in drag calculation
  class FixPropertyAtom* fix_Ksl_;
  class FixPropertyAtom* fix_uf_;
  class FixPropertyAtom* fix_dragAcc_;
};

}

#endif
#endif
