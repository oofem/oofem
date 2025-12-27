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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "prescribeddispslipmultiple.h"
#include "classfactory.h"
#include "node.h"
#include "masterdof.h"
#include "element.h"
#include "feinterpol.h"
#include "feinterpol2d.h"
#include "gausspoint.h"
#include "sparsemtrx.h"
#include "xfem/xfemelementinterface.h"
#include "xfem/integrationrules/discsegintegrationrule.h"
#include "timestep.h"
#include "function.h"
#include "sparselinsystemnm.h"
#include "unknownnumberingscheme.h"
#include "engngm.h"
#include "mathfem.h"
#include "dynamicinputrecord.h"


namespace oofem {
REGISTER_BoundaryCondition(PrescribedDispSlipMultiple);

PrescribedDispSlipMultiple :: PrescribedDispSlipMultiple(int n, Domain *d) :
    ActiveBoundaryCondition(n, d),
    PrescribedDispSlipHomogenization()
{
}

PrescribedDispSlipMultiple :: ~PrescribedDispSlipMultiple()
{
}


void PrescribedDispSlipMultiple :: initializeFrom(InputRecord &ir)
{
    ActiveBoundaryCondition :: initializeFrom(ir);
    PrescribedDispSlipHomogenization :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, bcs, _IFT_PrescribedDispSlipMultiple_BCs);
}


void PrescribedDispSlipMultiple :: giveInputRecord(DynamicInputRecord &input)
{
    ActiveBoundaryCondition :: giveInputRecord(input);
    PrescribedDispSlipHomogenization :: giveInputRecord(input);
    input.setField(this->bcs, _IFT_PrescribedDispSlipMultiple_BCs);
}


void PrescribedDispSlipMultiple::setSlipField( const FloatArray &t )
{
    this->slipField = t;
    for ( int i : this->bcs ) {
        auto bc = dynamic_cast<PrescribedDispSlipHomogenization*>(this->giveDomain()->giveBc(i));
        bc->setSlipField(t);
    }
}


void PrescribedDispSlipMultiple::setDispGradient( const FloatArray &t )
{
    this->slipGradient = FloatMatrix::fromArray(t);
    for ( int i : this->bcs ) {
        auto bc = dynamic_cast<PrescribedDispSlipHomogenization*>(this->giveDomain()->giveBc(i));
        bc->setDispGradient(t);
    }
}


void PrescribedDispSlipMultiple::setSlipGradient( const FloatArray &t )
{
    this->slipGradient = FloatMatrix::fromArray(t);
    for ( int i : this->bcs ) {
        auto bc = dynamic_cast<PrescribedDispSlipHomogenization*>(this->giveDomain()->giveBc(i));
        bc->setSlipGradient(t);
    }
}


void PrescribedDispSlipMultiple::setCenterCoordinate( FloatArray &x )
{
    this->mCenterCoord = x;
    for ( int i : this->bcs ) {
        auto bc = dynamic_cast<PrescribedDispSlipHomogenization*>(this->giveDomain()->giveBc(i));
        bc->setCenterCoordinate(x);
    }
}


DofManager *PrescribedDispSlipMultiple::giveInternalDofManager( int i )
{
    for ( int j: this->bcs ) {
        auto bc = dynamic_cast<ActiveBoundaryCondition*>(this->giveDomain()->giveBc(j));
        if ( bc != nullptr ) {
            return bc->giveInternalDofManager(i);
        }
    }
    return nullptr;
}


void PrescribedDispSlipMultiple :: scale(double s)
{
    this->dispField.times(s);
    this->dispGradient.times(s);
    this->slipField.times(s);
    this->slipGradient.times(s);
}


void PrescribedDispSlipMultiple :: computeStress(FloatArray &sigma, TimeStep *tStep)
{
    sigma.clear();
    for ( int i : this->bcs ) {
        FloatArray tmp;
        auto bc = dynamic_cast<PrescribedDispSlipHomogenization*>(this->giveDomain()->giveBc(i));
        bc->computeStress(tmp, tStep);
        sigma.add(tmp);
    }
}


void PrescribedDispSlipMultiple::computeTransferStress( FloatArray &bStress, TimeStep *tStep )
{
    bStress.clear();
    for ( int i : this->bcs ) {
        FloatArray tmp;
        auto bc = dynamic_cast<PrescribedDispSlipHomogenization*>(this->giveDomain()->giveBc(i));
        bc->computeTransferStress(tmp, tStep);
        bStress.add(tmp);
    }
}


void PrescribedDispSlipMultiple::computeReinfStress( FloatArray &rStress, TimeStep *tStep )
{
    rStress.clear();
    for ( int i : this->bcs ) {
        FloatArray tmp;
        auto bc = dynamic_cast<PrescribedDispSlipHomogenization*>(this->giveDomain()->giveBc(i));
        bc->computeReinfStress(tmp, tStep);
        rStress.add(tmp);
    }
}


void PrescribedDispSlipMultiple :: computeTangent(FloatMatrix &tangent, TimeStep *tStep)
{
    OOFEM_ERROR("Not implemented")
}
} /* namespace oofem */
