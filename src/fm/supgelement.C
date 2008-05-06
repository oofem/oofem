/* $Header: /home/cvs/bp/oofem/tm/src/transportelement.C,v 1.3.4.1 2004/04/05 15:19:53 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "supgelement.h"
#include "domain.h"
#include "timestep.h"
#include "node.h"
#include "dof.h"
#include "load.h"
#include "boundaryload.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "intarray.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "debug.h"
#include "verbose.h"

#include "elementside.h"
#include "mathfem.h"
#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <stdio.h>
#endif

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#include "conTable.h"
#endif

SUPGElement :: SUPGElement(int n, Domain *aDomain) :
    FMElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{ }


SUPGElement :: ~SUPGElement()
// Destructor.
{ }

IRResultType
SUPGElement :: initializeFrom(InputRecord *ir)
{
    //const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    //IRResultType result;                               // Required by IR_GIVE_FIELD macro

    FMElement :: initializeFrom(ir);
    this->computeGaussPoints();
    return IRRT_OK;
}



void
SUPGElement ::  giveCharacteristicMatrix(FloatMatrix &answer,
                                         CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver accordind to mtrx
//
{
    if ( mtrx == AccelerationTerm_MB ) {
        this->computeAccelerationTerm_MB(answer, tStep);
    } else if ( mtrx == AdvectionDerivativeTerm_MB )  {
        this->computeAdvectionDerivativeTerm_MB(answer, tStep);
    } else if ( mtrx == DiffusionDerivativeTerm_MB )  {
        this->computeDiffusionDerivativeTerm_MB(answer, TangentStiffness, tStep);
    } else if ( mtrx == SecantDiffusionDerivativeTerm_MB )  {
        this->computeDiffusionDerivativeTerm_MB(answer, SecantStiffness, tStep);
    } else if ( mtrx == TangentDiffusionDerivativeTerm_MB )  {
        this->computeDiffusionDerivativeTerm_MB(answer, TangentStiffness, tStep);
    } else if ( mtrx == InitialDiffusionDerivativeTerm_MB )  {
        this->computeDiffusionDerivativeTerm_MB(answer, ElasticStiffness, tStep);
    } else if ( mtrx == PressureTerm_MB )  {
        this->computePressureTerm_MB(answer, tStep);
    } else if ( mtrx == LinearAdvectionTerm_MC )  {
        this->computeLinearAdvectionTerm_MC(answer, tStep);
    } else if ( mtrx == AdvectionDerivativeTerm_MC )  {
        this->computeAdvectionDerivativeTerm_MC(answer, tStep);
    } else if ( mtrx == DiffusionDerivativeTerm_MC )  {
        this->computeDiffusionDerivativeTerm_MC(answer, tStep);
    } else if ( mtrx == AccelerationTerm_MC )  {
        this->computeAccelerationTerm_MC(answer, tStep);
    } else if ( mtrx == PressureTerm_MC )  {
        this->computePressureTerm_MC(answer, tStep);
    } else if ( mtrx == LSICStabilizationTerm_MB )  {
        this->computeLSICStabilizationTerm_MB(answer, tStep);
    } else {
        _error("giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");
    }

    return;
}


void
SUPGElement ::  giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                         TimeStep *tStep)
//
// returns characteristics vector of receiver according to requested type
//
{
    if ( mtrx == AdvectionTerm_MB ) {
        this->computeAdvectionTerm_MB(answer, tStep);
    } else if ( mtrx == DiffusionTerm_MB )  {
        this->computeDiffusionTerm_MB(answer, tStep);
    } else if ( mtrx == AdvectionTerm_MC )  {
        this->computeAdvectionTerm_MC(answer, tStep);
    } else if ( mtrx == BCRhsTerm_MB )  {
        this->computeBCRhsTerm_MB(answer, tStep);
    } else if ( mtrx == BCRhsTerm_MC )  {
        this->computeBCRhsTerm_MC(answer, tStep);
    } else if ( mtrx == DiffusionTerm_MC )  {
        this->computeDiffusionTerm_MC(answer, tStep);
    } else {
        _error("giveCharacteristicVector: Unknown Type of characteristic mtrx.");
    }

    return;
}




double
SUPGElement :: giveCharacteristicValue(CharType mtrx, TimeStep *tStep)
{
    if ( mtrx == CriticalTimeStep ) {
        return this->computeCriticalTimeStep(tStep);
    } else {
        _error("giveCharacteristicValue: Unknown Type of characteristic mtrx.");
    }

    return 0.0;
}



int
SUPGElement :: checkConsistency()
//
// check internal consistency
// mainly tests, whether material and crossSection data
// are safe for conversion to "Structural" versions
//
{
    int result = 1;
    /*
     * if (!this->giveMaterial()->testMaterialExtension(Material_TransportCapability)) {
     * _warning("checkConsistency : material without support for transport problems");
     * result =0;
     * }
     */
    return result;
}

void
SUPGElement :: updateInternalState(TimeStep *stepN)
{
    int i, j;
    IntegrationRule *iRule;
    FloatArray stress;

    // force updating strains & stresses
    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            computeDeviatoricStress(stress, iRule->getIntegrationPoint(j), stepN);
        }
    }
}

void
SUPGElement :: printOutputAt(FILE *file, TimeStep *stepN)
// Performs end-of-step operations.
{
    int i;

#ifdef __PARALLEL_MODE
    fprintf( file, "element %d [%8d] :\n", this->giveNumber(), this->giveGlobalNumber() );
#else
    fprintf(file, "element %d :\n", number);
#endif

    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        integrationRulesArray [ i ]->printOutputAt(file, stepN);
    }
}


#ifdef __OOFEG
int
SUPGElement :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                       int node, TimeStep *atTime)
{
    int indx = 1;
    Node *n = this->giveNode(node);

    if ( type == IST_Velocity ) {
        answer.resize( this->giveSpatialDimension() );
        int dofindx;
        if ( ( dofindx = n->findDofWithDofId(V_u) ) ) {
            answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(EID_MomentumBalance, VM_Total, atTime);
        }

        if ( ( dofindx = n->findDofWithDofId(V_v) ) ) {
            answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(EID_MomentumBalance, VM_Total, atTime);
        }

        if ( ( dofindx = n->findDofWithDofId(V_w) ) ) {
            answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(EID_MomentumBalance, VM_Total, atTime);
        }

        return 1;
    } else if ( type == IST_Pressure ) {
        int dofindx;
        if ( ( dofindx = n->findDofWithDofId(P_f) ) ) {
            answer.resize(1);
            answer.at(1) = n->giveDof(dofindx)->giveUnknown(EID_ConservationEquation, VM_Total, atTime);
            return 1;
        } else {
            return 0;
        }
    } else {
        return Element :: giveInternalStateAtNode(answer, type, mode, node, atTime);
    }
}

#endif

int
SUPGElement :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type)
{
    if ( ( type == IST_Velocity ) ) {
        IntArray mask;
        int indx = 1;
        answer.resize(3);
        this->giveElementDofIDMask(EID_MomentumBalance, mask);
        if ( mask.findFirstIndexOf(V_u) ) {
            answer.at(1) = indx++;
        }

        if ( mask.findFirstIndexOf(V_v) ) {
            answer.at(2) = indx++;
        }

        if ( mask.findFirstIndexOf(V_w) ) {
            answer.at(3) = indx++;
        }

        return 1;
    } else if ( type == IST_Pressure ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else {
        return Element :: giveIntVarCompFullIndx(answer, type);
    }
}

/*
 * void
 * SUPGElement::computeVectorOfPrescribed (EquationID ut, ValueModeType type, TimeStep* stepN, FloatArray& answer)
 * {
 * double scale;
 *
 * Element::computeVectorOfPrescribed (ut, type, stepN, answer);
 *
 * if (domain->giveEngngModel()->giveEquationScalingFlag()) {
 *  if (ut == EID_MomentumBalance) {
 *    scale = domain->giveEngngModel()->giveVariableScale(VST_Velocity);
 *  } else if (ut == EID_ConservationEquation) {
 *    scale = domain->giveEngngModel()->giveVariableScale(VST_Pressure);
 *  } else scale = 1.0;
 *  answer.times (1.0/scale);
 * }
 * }
 */
