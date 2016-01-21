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

#include "cbselement.h"
#include "dof.h"
#include "node.h"
#include "integrationrule.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "dynamicinputrecord.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
CBSElement :: CBSElement(int n, Domain *aDomain) :
    FMElement(n, aDomain)
{ }


CBSElement :: ~CBSElement()
{ }

IRResultType
CBSElement :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    IR_GIVE_OPTIONAL_FIELD(ir, boundarySides, _IFT_CBSElement_bsides);
    if ( !boundarySides.isEmpty() ) {
        IR_GIVE_FIELD(ir, boundaryCodes, _IFT_CBSElement_bcodes);
    }

    return FMElement :: initializeFrom(ir);
}


void
CBSElement :: giveInputRecord(DynamicInputRecord &input)
{
    FMElement :: giveInputRecord(input);
    if ( boundarySides.giveSize() > 0 ) {
        input.setField(this->boundarySides, _IFT_CBSElement_bsides);
        input.setField(this->boundaryCodes, _IFT_CBSElement_bcodes);
    }
}


void
CBSElement :: giveCharacteristicMatrix(FloatMatrix &answer,
                                       CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver according to mtrx
//
{
    if ( mtrx == MassMatrix ) {
        this->computeConsistentMassMtrx(answer, tStep);
    } else {
        OOFEM_ERROR("Unknown Type of characteristic mtrx.");
    }
}


void
CBSElement :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                       TimeStep *tStep)
//
// returns characteristics vector of receiver according to requested type
//
{
    if ( mtrx == LumpedMassMatrix ) {
        this->computeDiagonalMassMtrx(answer, tStep);
    }
    //else if (mtrx == PrescribedDensityRhsVector)
    //  this->computePrescribedTermsII (answer, mode, tStep);
    else {
        OOFEM_ERROR("Unknown Type of characteristic mtrx.");
    }
}


void
CBSElement :: computePrescribedTermsI(FloatArray &answer, TimeStep *tStep)
{
    FloatMatrix mass;
    FloatArray usp;
    this->computeConsistentMassMtrx(mass, tStep);
    this->computeVectorOfVelocities(VM_Incremental, tStep, usp);
    answer.beProductOf(mass, usp);
    answer.negated();
}

#if 0
void
CBSElement :: computePrescribedTermsII(FloatArray &answer, ValueModeType mode, TimeStep *tStep)
{
    FloatMatrix lhs;
    FloatArray usp;
    this->computePressureLhs(lhs, tStep);
    this->computeVectorOfPressures(mode, tStep, usp);
    answer.beProductOf(lhs, usp);
    answer.negated();
}
#endif

int
CBSElement :: checkConsistency()
//
// check internal consistency
// mainly tests, whether material and crossSection data
// are safe for conversion to "Structural" versions
//
{
    int result = 1;
    return result;
}


void
CBSElement :: updateInternalState(TimeStep *tStep)
{
    FloatArray stress;

    // force updating strains & stresses
    for ( auto &iRule: integrationRulesArray ) {
        for ( GaussPoint *gp: *iRule ) {
            computeDeviatoricStress(stress, gp, tStep);
        }
    }
}


#ifdef __OOFEG
int
CBSElement :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                      int node, TimeStep *tStep)
{
    int indx = 1;
    Node *n = this->giveNode(node);

    if ( type == IST_Velocity ) {
        answer.resize( this->giveSpatialDimension() );
        std::vector< Dof* >::const_iterator dofindx;
        if ( ( dofindx = n->findDofWithDofId(V_u) ) != n->end() ) {
            answer.at(indx++) = (*dofindx)->giveUnknown(VM_Total, tStep);
        } else if ( ( dofindx = n->findDofWithDofId(V_v) ) != n->end() ) {
            answer.at(indx++) = (*dofindx)->giveUnknown(VM_Total, tStep);
        } else if ( ( dofindx = n->findDofWithDofId(V_w) ) != n->end() ) {
            answer.at(indx++) = (*dofindx)->giveUnknown(VM_Total, tStep);
        }

        return 1;
    } else if ( type == IST_Pressure ) {
        auto dofindx = n->findDofWithDofId(P_f);
        if ( dofindx != n->end() ) {
            answer.resize(1);
            answer.at(1) = (*dofindx)->giveUnknown(VM_Total, tStep);
            return 1;
        } else {
            return 0;
        }
    } else {
        return Element :: giveInternalStateAtNode(answer, type, mode, node, tStep);
    }
}

#endif
} // end namespace oofem
