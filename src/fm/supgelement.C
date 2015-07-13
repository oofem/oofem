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

#include "supgelement.h"
#include "domain.h"
#include "load.h"
#include "gausspoint.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "fluiddynamicmaterial.h"
#include "fluidcrosssection.h"
#include "dynamicinputrecord.h"
#include "engngm.h"
#include "node.h"
#include "dof.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
SUPGElement :: SUPGElement(int n, Domain *aDomain) :
    FMElement(n, aDomain), t_supg(0), t_pspg(0), t_lsic(0)
{ }


SUPGElement :: ~SUPGElement()
// Destructor.
{ }

IRResultType
SUPGElement :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro


    IR_GIVE_OPTIONAL_FIELD(ir, boundarySides, _IFT_SUPGElement_bsides);
    if ( !boundarySides.isEmpty() ) {
        IR_GIVE_FIELD(ir, boundaryCodes, _IFT_SUPGElement_bcodes);
    }

    return FMElement :: initializeFrom(ir);
}


void
SUPGElement :: giveInputRecord(DynamicInputRecord &input)
{
    FMElement :: giveInputRecord(input);
    if ( !boundarySides.isEmpty() ) {
        input.setField(this->boundarySides, _IFT_SUPGElement_bsides);
        input.setField(this->boundaryCodes, _IFT_SUPGElement_bcodes);
    }
}


void
SUPGElement :: giveCharacteristicMatrix(FloatMatrix &answer,
                                        CharType type, TimeStep *tStep)
//
// returns characteristics matrix of receiver according to mtrx
//
{

    if ( type == TangentStiffnessMatrix ) {
            // stokes flow only
        double dscale = this->giveDomain()->giveEngngModel()->giveVariableScale(VST_Density);
        double uscale = this->giveDomain()->giveEngngModel()->giveVariableScale(VST_Velocity);

        IntArray vloc, ploc;
        FloatMatrix h;
        int size = this->computeNumberOfDofs();
        this->giveLocalVelocityDofMap(vloc);
        this->giveLocalPressureDofMap(ploc);
        answer.resize(size, size);
        answer.zero();

        //this->computeAccelerationTerm_MB(h, tStep);
        //answer.assemble(h, vloc);
        this->computeDiffusionDerivativeTerm_MB(h, TangentStiffness, tStep);
        answer.assemble(h, vloc);
        this->computePressureTerm_MB(h, tStep);
        answer.assemble(h, vloc, ploc);
        //this->computeLSICStabilizationTerm_MB(h, tStep);
        //h.times( alpha * tStep->giveTimeIncrement() * lscale / ( dscale * uscale * uscale ) );
        //answer.assemble(h, vloc);
        this->computeBCLhsTerm_MB(h, tStep);
        if ( h.isNotEmpty() ) {
            answer.assemble(h, vloc);
        }

        this->computeBCLhsPressureTerm_MB(h, tStep);
        if ( h.isNotEmpty() ) {
            answer.assemble(h, vloc, ploc);
        }

        // conservation eq part
        this->computeLinearAdvectionTerm_MC(h, tStep);
        h.times( 1.0 / ( dscale * uscale ) );
        answer.assemble(h, ploc, vloc);
        this->computeBCLhsPressureTerm_MC(h, tStep);
        if ( h.isNotEmpty() ) {
            answer.assemble(h, ploc, vloc);
        }

        this->computeDiffusionDerivativeTerm_MC(h, tStep);
        answer.assemble(h, ploc, vloc);
        this->computePressureTerm_MC(h, tStep);
        answer.assemble(h, ploc);
    } else {
      OOFEM_ERROR("giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");
    }
   
}


void
SUPGElement :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                        TimeStep *tStep)
//
// returns characteristics vector of receiver according to requested type
//
{
    if ( mtrx == ExternalForcesVector ) {
        // stokes flow
        IntArray vloc, ploc;
        FloatArray h;
        int size = this->computeNumberOfDofs();
        this->giveLocalVelocityDofMap(vloc);
        this->giveLocalPressureDofMap(ploc);
        answer.resize(size);
        answer.zero();
        this->computeBCRhsTerm_MB(h, tStep);
        answer.assemble(h, vloc);
        this->computeBCRhsTerm_MC(h, tStep);
        answer.assemble(h, ploc);
    }

#if 1
    else if ( mtrx == InternalForcesVector ) {
        // stokes flow
        IntArray vloc, ploc;
        FloatArray h;
        int size = this->computeNumberOfDofs();
        this->giveLocalVelocityDofMap(vloc);
        this->giveLocalPressureDofMap(ploc);
        answer.resize(size);
        answer.zero();
        //this->computeAdvectionTerm_MB(h, tStep); answer.assemble(h, vloc);
        this->computeAdvectionTerm_MC(h, tStep);
        answer.assemble(h, ploc);
        this->computeDiffusionTerm_MB(h, tStep);
        answer.assemble(h, vloc);
        this->computeDiffusionTerm_MC(h, tStep);
        answer.assemble(h, ploc);

        FloatMatrix m1;
        FloatArray v;
        // add lsic stabilization term
        //this->giveCharacteristicMatrix(m1, LSICStabilizationTerm_MB, tStep);
        //m1.times( lscale / ( dscale * uscale * uscale ) );
        this->computeVectorOfVelocities(VM_Total, tStep, v);
        //h.beProductOf(m1, v);
        //answer.assemble(h, vloc);
        this->computeLinearAdvectionTerm_MC(m1, tStep);
        //m1.times( 1. / ( dscale * uscale ) );
        h.beProductOf(m1, v);
        answer.assemble(h, ploc);

        // add pressure term
        this->computePressureTerm_MB(m1, tStep);
        this->computeVectorOfPressures(VM_Total, tStep, v);
        h.beProductOf(m1, v);
        answer.assemble(h, vloc);

        // pressure term
        this->computePressureTerm_MC(m1, tStep);
        this->computeVectorOfPressures(VM_Total, tStep, v);
        h.beProductOf(m1, v);
        answer.assemble(h, ploc);
    }
#endif
    else {
        OOFEM_ERROR("Unknown Type of characteristic mtrx.");
    }
}





void
SUPGElement :: computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray eps;

    // compute deviatoric strain
    this->computeDeviatoricStrain(eps, gp, tStep);
    // call material to compute stress
    static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial()->computeDeviatoricStressVector(answer, gp, eps, tStep);
}



void
SUPGElement :: computeBCLhsTerm_MB(FloatMatrix &answer, TimeStep *tStep)
{
    bcType boundarytype;
    int nLoads = 0;
    //bcType loadtype;
    FloatMatrix helpMatrix;
    // loop over boundary load array

    answer.clear();

    nLoads = this->giveBoundaryLoadArray()->giveSize() / 2;
    if ( nLoads ) {
        for ( int i = 1; i <= nLoads; i++ ) {
            int n = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
            int side = boundaryLoadArray.at(i * 2);
            Load *load = domain->giveLoad(n);
            boundarytype = load->giveType();
            if ( boundarytype == SlipWithFriction ) {
                this->computeSlipWithFrictionBCTerm_MB(helpMatrix, load, side, tStep);
                answer.add(helpMatrix);
            } else if ( boundarytype == PenetrationWithResistance ) {
                this->computePenetrationWithResistanceBCTerm_MB(helpMatrix, load, side, tStep);
                answer.add(helpMatrix);
            } else {
                // OOFEM_ERROR("unsupported load type class");
            }
        }
    }

    nLoads = this->giveBodyLoadArray()->giveSize();

    if ( nLoads ) {
        bcGeomType ltype;
        for ( int i = 1; i <= nLoads; i++ ) {
            Load *load = domain->giveLoad( bodyLoadArray.at(i) );
            ltype = load->giveBCGeoType();
            if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ReinforceBVT ) ) {
                this->computeHomogenizedReinforceTerm_MB(helpMatrix, load, tStep);
                answer.add(helpMatrix);
            }
        }
    }
}

void
SUPGElement :: computeBCLhsPressureTerm_MB(FloatMatrix &answer, TimeStep *tStep)
{
    bcType boundarytype;
    int nLoads = 0;
    //bcType loadtype;
    FloatMatrix helpMatrix;
    // loop over boundary load array
    answer.clear();

    nLoads = this->giveBoundaryLoadArray()->giveSize() / 2;

    if ( nLoads ) {
        for ( int i = 1; i <= nLoads; i++ ) {
            int n = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
            int side = boundaryLoadArray.at(i * 2);
            Load *load = domain->giveLoad(n);
            boundarytype = load->giveType();
            if ( boundarytype == OutFlowBC ) {
                this->computeOutFlowBCTerm_MB(helpMatrix, side, tStep);
                answer.add(helpMatrix);
            } else {
                //_warning("computeForceLoadVector : unsupported load type class");
            }
        }
    }
}

void
SUPGElement :: computeBCLhsPressureTerm_MC(FloatMatrix &answer, TimeStep *tStep)
{
    int nLoads = 0;
    //bcType loadtype;
    FloatMatrix helpMatrix;

    nLoads = this->giveBodyLoadArray()->giveSize();
    answer.clear();
    if ( nLoads ) {
        bcGeomType ltype;
        for ( int i = 1; i <= nLoads; i++ ) {
            Load *load  = domain->giveLoad( bodyLoadArray.at(i) );
            ltype = load->giveBCGeoType();
            if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ReinforceBVT ) ) {
                this->computeHomogenizedReinforceTerm_MC(helpMatrix, load, tStep);
                answer.add(helpMatrix);
            }
        }
    }
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
    return result;
}

void
SUPGElement :: updateInternalState(TimeStep *tStep)
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
SUPGElement :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
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
