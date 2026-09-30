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

#include "lm/latticeneumanncoupling.h"
#include "domain.h"
#include "engngm.h"
#include "staggeredproblem.h"
#include "node.h"
#include "dof.h"
#include "floatarray.h"
#include "intarray.h"
#include "unknownnumberingscheme.h"
#include "classfactory.h"
#include "dofiditem.h"

#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace oofem {
REGISTER_BoundaryCondition(LatticeNeumannCoupling);

void
LatticeNeumannCoupling :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    ActiveBoundaryCondition :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, this->smNodes, _IFT_LatticeNeumannCoupling_smnodes);
    IR_GIVE_FIELD(ir, this->tmNodes, _IFT_LatticeNeumannCoupling_tmnodes);
    IR_GIVE_FIELD(ir, this->directionVector, _IFT_LatticeNeumannCoupling_direction);

    if ( this->smNodes.giveSize() != this->tmNodes.giveSize() ) {
        OOFEM_ERROR("smnodes (%d) and tmnodes (%d) must have the same size",
                    this->smNodes.giveSize(), this->tmNodes.giveSize());
    }
}

Node *
LatticeNeumannCoupling :: giveCoupledTransportNode(int tmNodeId)
{
    // Reach across to the transport slave problem of the driving StaggeredProblem
    // coupledModels.at(2) is the
    // Neumann coupling slot (the transport slave index).
    EngngModel *master = this->giveDomain()->giveEngngModel()->giveMasterEngngModel();
    StaggeredProblem *sp = dynamic_cast< StaggeredProblem * >(master);
    if ( !sp ) {
        OOFEM_ERROR("master engng model is not a StaggeredProblem (coupling requires staggered solve)");
    }
    IntArray coupledModels;
    sp->giveCoupledModels(coupledModels);
    if ( coupledModels.giveSize() < 2 || coupledModels.at(2) == 0 ) {
        OOFEM_ERROR("transport coupling slot (coupling field, slot 2) is not set");
    }
    return sp->giveSlaveProblem( coupledModels.at(2) )->giveDomain(1)->giveNode(tmNodeId);
}

void
LatticeNeumannCoupling :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                         CharType type, ValueModeType mode,
                                         const UnknownNumberingScheme &s, FloatArray *eNorms,
                                         void *lock)
{
    if ( type != ExternalForcesVector ) {
        return;
    }

    IntArray dofIdArray = {
        D_u, D_v, R_w
    };
    IntArray loc;
    FloatArray fext;

    for ( int pos = 1; pos <= this->smNodes.giveSize(); ++pos ) {
        Node *smNode = this->giveDomain()->giveNode( this->smNodes.at(pos) );
        Node *tmNode = this->giveCoupledTransportNode( this->tmNodes.at(pos) );

        // Fluid pressure at the coupled transport node
        double pf = tmNode->giveDofWithID(P_f)->giveUnknown(mode, tStep);

        // In-plane distance between the mechanical node and its transport counterpart.
        const FloatArray &cm = smNode->giveCoordinates();
        const FloatArray &ct = tmNode->giveCoordinates();
        double dx = ct.at(1) - cm.at(1);
        double dy = ct.at(2) - cm.at(2);
        double distance = std::sqrt(dx * dx + dy * dy);

        // f = P_f * distance * direction (pressure transferred to a nodal force).
        fext = this->directionVector;
        fext.times(pf * distance);

        smNode->giveLocationArray(dofIdArray, loc, s);
#ifdef _OPENMP
        if ( lock ) omp_set_lock(static_cast< omp_lock_t * >(lock));
#endif
        answer.assemble(fext, loc);
#ifdef _OPENMP
        if ( lock ) omp_unset_lock(static_cast< omp_lock_t * >(lock));
#endif
    }
}
} // namespace oofem
