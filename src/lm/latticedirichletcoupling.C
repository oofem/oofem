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

#include "lm/latticedirichletcoupling.h"
#include "lm/Elements/latticestructuralelement.h"
#include "domain.h"
#include "engngm.h"
#include "staggeredproblem.h"
#include "dof.h"
#include "dofmanager.h"
#include "timestep.h"
#include "floatarray.h"
#include "intarray.h"
#include "classfactory.h"

#include <cmath>

namespace oofem {
REGISTER_BoundaryCondition(LatticeDirichletCoupling);

void
LatticeDirichletCoupling :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    BoundaryCondition :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, this->couplingElements, _IFT_LatticeDirichletCoupling_couplingelements);
}

double
LatticeDirichletCoupling :: give(Dof *dof, ValueModeType mode, double time)
{
    // One-step-lagged (explicit) staggered coupling: return the value computed
    // from the mechanical stress of the previous step. The mechanical slave
    // (prob2) is solved before the transport problem within a step (sm-first —
    // Dirichlet couples mechanical to transport, so the mechanical stress must
    // exist when the transport reads it), so its current stress is already this
    // step's; the status returns the previous step's. The first step returns 0.
    if ( time != this->cachedTime ) {
        this->previousValue = this->cachedValue;
        this->cachedValue = this->computeCouplingValue(dof);
        this->cachedTime = time;
    }
    return this->previousValue;
}

double
LatticeDirichletCoupling :: computeCouplingValue(Dof *dof)
{
    // Reach across to the mechanical slave problem of the driving StaggeredProblem.
    // coupledModels.at(1) is the Dirichlet coupling slot (the mechanical slave index).
    EngngModel *master = this->giveDomain()->giveEngngModel()->giveMasterEngngModel();
    StaggeredProblem *sp = dynamic_cast< StaggeredProblem * >(master);
    if ( !sp ) {
        OOFEM_ERROR("master engng model is not a StaggeredProblem (coupling requires staggered solve)");
    }
    IntArray coupledModels;
    sp->giveCoupledModels(coupledModels);
    if ( coupledModels.giveSize() < 1 || coupledModels.at(1) == 0 ) {
        OOFEM_ERROR("mechanical coupling slot (coupling field, slot 1) is not set");
    }
    Domain *smDomain = sp->giveSlaveProblem( coupledModels.at(1) )->giveDomain(1);

    const FloatArray &nodeCoords = dof->giveDofManager()->giveCoordinates();

    double nominator = 0., denominator = 0.;
    FloatArray gpCoords;
    for ( int i = 1; i <= this->couplingElements.giveSize(); i++ ) {
        LatticeStructuralElement *el = static_cast< LatticeStructuralElement * >
            ( smDomain->giveElement( this->couplingElements.at(i) ) );

        el->giveGpCoordinates(gpCoords);

        // Read the temp (current, not-yet-committed) normal stress: in the
        // staggered flow the mechanical material status is not committed between
        // the sub-problem solves, so the committed value is still zero here. The
        // one-step lag is supplied by the cache in give().
        double normalStress = el->giveTempNormalStress();

        // Compression only: tensile normal stress does not generate pore pressure here.
        if ( normalStress > 0. ) {
            normalStress = 0.;
        }

        double dx = gpCoords.at(1) - nodeCoords.at(1);
        double dy = gpCoords.at(2) - nodeCoords.at(2);
        double distance = std::sqrt(dx * dx + dy * dy);

        nominator   += distance * normalStress;
        denominator += distance;
    }

    return ( denominator != 0. ) ? nominator / denominator : 0.;
}
} // namespace oofem
