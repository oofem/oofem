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

#include "pfemelement.h"
#include "domain.h"
#include "timestep.h"
#include "node.h"
#include "dof.h"
#include "load.h"
#include "boundaryload.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "verbose.h"
#include "material.h"
#include "elementside.h"
#include "mathfem.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
PFEMElement :: PFEMElement(int n, Domain *aDomain) :
    FMElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{ }


PFEMElement :: ~PFEMElement()
// Destructor.
{ }

IRResultType
PFEMElement :: initializeFrom(InputRecord *ir)
{
    //const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    //IRResultType result;                               // Required by IR_GIVE_FIELD macro

    FMElement :: initializeFrom(ir);

    this->computeGaussPoints();
    return IRRT_OK;
}


void
PFEMElement ::  giveCharacteristicMatrix(FloatMatrix &answer,
                                         CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver according to mtrx
//
{
    FloatMatrix reducedAnswer;
    int size = this->computeNumberOfDofs();
    answer.resize(size, size);
    answer.zero();

    if ( mtrx == LumpedMassMatrix ) {
        this->computeDiagonalMassMtrx(reducedAnswer, tStep);
        answer.assemble(reducedAnswer, this->giveVelocityDofMask(), this->giveVelocityDofMask());
    } else if ( mtrx == StiffnessMatrix ) { //Tangent stiffness matrix set as default
        this->computeStiffnessMatrix(reducedAnswer, TangentStiffness, tStep);
        answer.assemble(reducedAnswer, this->giveVelocityDofMask(), this->giveVelocityDofMask());
    } else if ( mtrx == TangentStiffnessMatrix ) {
        this->computeStiffnessMatrix(reducedAnswer, TangentStiffness, tStep);
        answer.assemble(reducedAnswer, this->giveVelocityDofMask(), this->giveVelocityDofMask());
    } else if ( mtrx == PressureGradientMatrix ) {
        this->computeGradientMatrix(reducedAnswer, tStep);
        answer.assemble(reducedAnswer, this->giveVelocityDofMask(), this->givePressureDofMask());
    } else if ( mtrx == DivergenceMatrix ) {
        this->computeDivergenceMatrix(reducedAnswer, tStep);
        answer.assemble(reducedAnswer, this->givePressureDofMask(), this->giveVelocityDofMask());
    } else if ( mtrx == PressureLaplacianMatrix ) {
        this->computePressureLaplacianMatrix(reducedAnswer, tStep);
        answer.assemble(reducedAnswer, this->givePressureDofMask(), this->givePressureDofMask());
    } else {
        OOFEM_ERROR("giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");
    }
}


void
PFEMElement ::  giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                         TimeStep *tStep)
//
// returns characteristics vector of receiver according to requested type
//
{
    FloatArray reducedVector;
    answer.resize(this->computeNumberOfDofs());
    answer.zero();

    if ( mtrx == ExternalForcesVector ) {
        this->computeForceVector(reducedVector, tStep);
        answer.assemble(reducedVector, this->giveVelocityDofMask()); 
    } else if ( mtrx == PressureGradientVector ) {
        FloatMatrix g;
        this->giveCharacteristicMatrix(g, PressureGradientMatrix, tStep);
        FloatArray p;
        this->computeVectorOfPressures(VM_Total, tStep, reducedVector);
        p.resize(this->computeNumberOfDofs());
        p.assemble(reducedVector, this->givePressureDofMask());
        answer.beProductOf(g, p);
    } else if ( mtrx == MassVelocityVector ) {
        FloatMatrix m;
        FloatArray u;
        this->giveCharacteristicMatrix(m, LumpedMassMatrix, tStep);
        this->computeVectorOfVelocities(VM_Total, tStep, reducedVector);
        u.resize(this->computeNumberOfDofs());
        u.assemble(reducedVector, this->giveVelocityDofMask());
        answer.beProductOf(m, u);
    } else if ( mtrx == LaplaceVelocityVector ) {
        FloatMatrix l;
        FloatArray u;
        this->giveCharacteristicMatrix(l, StiffnessMatrix, tStep);
        this->computeVectorOfVelocities(VM_Total, tStep, reducedVector);
        u.resize(this->computeNumberOfDofs());
        u.assemble(reducedVector, this->giveVelocityDofMask());
        answer.beProductOf(l, u);
    } else if ( mtrx == MassAuxVelocityVector ) {
        FloatMatrix m;
        FloatArray u;
        this->giveCharacteristicMatrix(m, LumpedMassMatrix, tStep);
        this->computeVectorOf(VM_Intermediate, tStep, u);
        answer.beProductOf(m, u);
    } else if ( mtrx == DivergenceAuxVelocityVector ) {
        FloatMatrix d;
        this->giveCharacteristicMatrix(d, DivergenceMatrix, tStep);
        FloatArray u_star;
        this->computeVectorOf(VM_Intermediate, tStep, u_star);
        answer.beProductOf(d, u_star);
    } else if ( mtrx == DivergenceVelocityVector ) {
        FloatMatrix d;
        this->giveCharacteristicMatrix(d, DivergenceMatrix, tStep);
        FloatArray u;
        this->computeVectorOfVelocities(VM_Total, tStep, u);
        answer.beProductOf(d, u);
    } else if ( mtrx == LumpedMassMatrix ) {
        this->computeDiagonalMassMtrx(reducedVector, tStep);
        answer.assemble(reducedVector, this->giveVelocityDofMask());
    } else {
        OOFEM_ERROR("giveCharacteristicVector: Unknown Type of characteristic mtrx.");
    }
}

int
PFEMElement :: checkConsistency()
// no check at the moment
{
    int result = 1;

    return result;
}

void
PFEMElement :: updateInternalState(TimeStep *stepN)
{
    FloatArray stress;
    // force updating strains & stresses
    for ( auto &iRule: integrationRulesArray ) {
        for ( GaussPoint *gp: *iRule ) {
            computeDeviatoricStress(stress, gp, stepN);
        }
    }
}

void
PFEMElement :: printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
#ifdef __PARALLEL_MODE
    fprintf( file, "element %d [%8d] :\n", this->giveNumber(), this->giveGlobalNumber() );
#else
    fprintf(file, "element %d :\n", number);
#endif
    for ( auto &iRule: integrationRulesArray ) {
        iRule->printOutputAt(file, tStep);
    }
}


#ifdef __OOFEG
int
PFEMElement :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                       int node, TimeStep *atTime)
{
    int indx = 1;
    Node *n = this->giveNode(node);

    if ( type == IST_Velocity ) {
        answer.resize( this->giveSpatialDimension() );
        int dofindx;
        if ( ( dofindx = n->findDofWithDofId(V_u) ) ) {
            answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(VM_Total, atTime);
        }

        if ( ( dofindx = n->findDofWithDofId(V_v) ) ) {
            answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(VM_Total, atTime);
        }

        if ( ( dofindx = n->findDofWithDofId(V_w) ) ) {
            answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(VM_Total, atTime);
        }

        return 1;
    } else if ( type == IST_Pressure ) {
        int dofindx;
        if ( ( dofindx = n->findDofWithDofId(P_f) ) ) {
            answer.resize(1);
            answer.at(1) = n->giveDof(dofindx)->giveUnknown(VM_Total, atTime);
            return 1;
        } else {
            return 0;
        }
    } else {
        return Element :: giveInternalStateAtNode(answer, type, mode, node, atTime);
    }
}

#endif

// might be needed for OOFEG-plot
// int
// PFEMElement :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type)
// {
//     if ( ( type == IST_Velocity ) ) {
//         IntArray mask;
//         int indx = 1;
//         answer.resize(3);
//         this->giveElementDofIDMask(EID_MomentumBalance, mask);
//         if ( mask.findFirstIndexOf(V_u) ) {
//             answer.at(1) = indx++;
//         }
//
//         if ( mask.findFirstIndexOf(V_v) ) {
//             answer.at(2) = indx++;
//         }
//
//         if ( mask.findFirstIndexOf(V_w) ) {
//             answer.at(3) = indx++;
//         }
//
//         return 1;
//     } else if ( type == IST_Pressure ) {
//         answer.resize(1);
//         answer.at(1) = 1;
//         return 1;
//     } else {
//         return Element :: giveIntVarCompFullIndx(answer, type);
//     }
// }
} // end namespace oofem
