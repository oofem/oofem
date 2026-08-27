/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef latticedirichletcoupling_h
#define latticedirichletcoupling_h

#include "boundarycondition.h"
#include "intarray.h"

///@name Input fields for LatticeDirichletCoupling
//@{
#define _IFT_LatticeDirichletCoupling_Name "latticedirichletcoupling"
#define _IFT_LatticeDirichletCoupling_couplingelements "couplingelements"
//@}

namespace oofem {
class Dof;
class TimeStep;

/**
 * Dirichlet-type hydro-mechanical coupling boundary condition (Approach 2 of
 * Grassl, Fahy, Gallipoli, Wheeler 2015). Prescribes the pore-fluid pressure
 * P_f at a transport lattice node to the distance-weighted average of the
 * (compression-only) normal stress of coupled mechanical lattice elements,
 * read from the mechanical slave problem of a StaggeredProblem.
 */
class LatticeDirichletCoupling : public BoundaryCondition
{
protected:
    /// Mechanical (slave problem) elements whose normal stress drives the pressure.
    IntArray couplingElements;

    // Explicit (one-step-lagged) staggered coupling: the pressure prescribed at
    // step n uses the mechanical stress committed at step n-1. Dirichlet coupling
    // is mechanical to transport, so the mechanical slave (prob2) is solved before
    // the transport problem within a step (sm-first); its current stress is thus
    // already step n, and these cache the value so give() can return step n-1's.
    // Set on the (target) time, because the transport field applies Dirichlet
    // values through the give(.., double time) overload, not the TimeStep one.
    double cachedValue = 0.;
    double previousValue = 0.;
    double cachedTime = -1.;

public:
    LatticeDirichletCoupling(int n, Domain *d) : BoundaryCondition(n, d) { }

    void initializeFrom(const std::shared_ptr<InputRecord> &ir) override;

    double give(Dof *dof, ValueModeType mode, double time) override;

    const char *giveClassName() const override { return "LatticeDirichletCoupling"; }
    const char *giveInputRecordName() const override { return _IFT_LatticeDirichletCoupling_Name; }

protected:
    /// Distance-weighted average of the (compression-only) current normal stress of the coupled mechanical elements.
    double computeCouplingValue(Dof *dof);
};
} // namespace oofem
#endif // latticedirichletcoupling_h
