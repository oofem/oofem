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


#include "interactionpfemparticle.h"
#include "timestep.h"
#include "classfactory.h"
#include "floatmatrix.h"
#include "fluidstructureproblem.h"
#include "../../sm/EngineeringModels/structengngmodel.h"
#include "domain.h"
#include "dof.h"
#include "mathfem.h"


namespace oofem {
REGISTER_DofManager(InteractionPFEMParticle);
/**
 * Constructor. Creates a particle with number n, belonging to aDomain.
 */
InteractionPFEMParticle :: InteractionPFEMParticle(int n, Domain *aDomain) : PFEMParticle(n, aDomain), coupledNode(0)
{ }
//from hanging node

/**
 * Gets from the source line from the data file all the data of the receiver.
 */
IRResultType
InteractionPFEMParticle :: initializeFrom(InputRecord *ir)
{
    IRResultType result;							// Required by IR_GIVE_FIELD macro

    result = PFEMParticle :: initializeFrom(ir);

    IR_GIVE_OPTIONAL_FIELD(ir, coupledNode, _IFT_InteractionPFEMParticle_CoupledNode);

    return (result != IRRT_OK) ? result : IRRT_OK;
}

/**
 * Checks internal data consistency in node.
 */
int
InteractionPFEMParticle :: checkConsistency()
{
    return PFEMParticle :: checkConsistency();
}

void
InteractionPFEMParticle :: updateYourself(TimeStep *tStep)
{
    PFEMParticle :: updateYourself(tStep);
}

void
InteractionPFEMParticle :: givePrescribedUnknownVector(FloatArray &answer, const IntArray &dofIDArry,
                                                       ValueModeType mode, TimeStep *stepN)
{
    answer.resize( dofIDArry.giveSize() );


    FloatArray velocities;
    FluidStructureProblem *fsiProblem = this->giveFluidStructureMasterProblem();
    if (fsiProblem) {
        StructuralEngngModel *structuralProblem = this->giveStructuralProblem();
        if ( structuralProblem ) {
            int j = 1;
            if (fsiProblem->giveIterationNumber() < 1) {  
                for (int dofid: dofIDArry ) {
                    answer.at(j++) = this->giveDofWithID( dofid )->giveBcValue(mode, stepN);
                }
            } else {
                DofManager *dman = structuralProblem->giveDomain(1)->giveDofManager(coupledNode);
                //dman->giveUnknownVectorOfType(velocities, VelocityVector, VM_Velocity, stepN);
                //dman->giveUnknownVector(velocities, dofIDArry, VM_Velocity, stepN);
                dman->giveCompleteUnknownVector(velocities, VM_Velocity, stepN);
                for ( int dofid: dofIDArry) {
                    answer.at(dofid) = velocities.at( dofid );
                }
            }
        }
    }

    // Transform to global c.s.
    FloatMatrix L2G;
    if (this->computeL2GTransformation(L2G, dofIDArry)) {
        answer.rotatedWith(L2G, 'n');
    }
}

void
InteractionPFEMParticle::giveCoupledVelocities(FloatArray &answer, TimeStep *stepN)
{
    StructuralEngngModel* structuralProblem = this->giveStructuralProblem();
    if ( structuralProblem ) {
        DofManager *dman = structuralProblem->giveDomain(1)->giveDofManager(coupledNode);
        dman->giveCompleteUnknownVector(answer, VM_Velocity, stepN);
    }
}

void
InteractionPFEMParticle :: printOutputAt(FILE *stream, TimeStep *stepN)
{
    PFEMParticle :: printOutputAt(stream, stepN);
}

#ifdef __OOFEG
void InteractionPFEMParticle :: drawScalar(oofegGraphicContext &gc)
{
    PFEMParticle :: drawScalar(gc);
}
#endif

StructuralEngngModel*
InteractionPFEMParticle :: giveStructuralProblem()
{
    StructuralEngngModel *structuralProblem = NULL;
    FluidStructureProblem *fsiProblem = this->giveFluidStructureMasterProblem();
    if (fsiProblem)
    {
        for ( int i = 1; i <= fsiProblem->giveNumberOfSlaveProblems(); i++ ) {
            structuralProblem = dynamic_cast<StructuralEngngModel*>(fsiProblem->giveSlaveProblem(i));
        }
    }
    return structuralProblem;
}

FluidStructureProblem*
InteractionPFEMParticle :: giveFluidStructureMasterProblem()
{
    FluidStructureProblem *fsiProblem = dynamic_cast<FluidStructureProblem*>(domain->giveEngngModel()->giveMasterEngngModel());

    return fsiProblem;
}
} // end namespace oofem
