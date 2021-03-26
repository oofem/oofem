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
 *               Copyright (C) 1993 - 2021   Borek Patzak
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

#include "prescribedgradientmultiple.h"
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

#ifdef _OPENMP
#include <omp.h>
#endif

namespace oofem {
REGISTER_BoundaryCondition(PrescribedGradientMultiple);

PrescribedGradientMultiple :: PrescribedGradientMultiple(int n, Domain *d) :
    ActiveBoundaryCondition(n, d),
    PrescribedGradientHomogenization()
{
}

PrescribedGradientMultiple :: ~PrescribedGradientMultiple()
{
}


void PrescribedGradientMultiple :: initializeFrom(InputRecord &ir)
{
    ActiveBoundaryCondition :: initializeFrom(ir);
    PrescribedGradientHomogenization :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, bcs, _IFT_PrescribedGradientMultiple_BCs);
}


void PrescribedGradientMultiple :: giveInputRecord(DynamicInputRecord &input)
{
    ActiveBoundaryCondition :: giveInputRecord(input);
    PrescribedGradientHomogenization :: giveInputRecord(input);
    input.setField(this->bcs, _IFT_PrescribedGradientMultiple_BCs);
}


void PrescribedGradientMultiple :: setPrescribedGradient( const FloatMatrix &t )
{
    this->mGradient = t;
    for ( int i : this->bcs ) {
        auto bc = dynamic_cast<PrescribedGradientHomogenization*>(this->giveDomain()->giveBc(i));
        bc->setPrescribedGradient(t);
    }
}


void PrescribedGradientMultiple :: setPrescribedGradientVoigt( const FloatArray &t )
{
    PrescribedGradientHomogenization ::  setPrescribedGradientVoigt(t);
    for ( int i : this->bcs ) {
        auto bc = dynamic_cast<PrescribedGradientHomogenization*>(this->giveDomain()->giveBc(i));
        bc->setPrescribedGradientVoigt(t);
    }
}


void PrescribedGradientMultiple::setCenterCoordinate( FloatArray &x )
{
    this->mCenterCoord = x;
    for ( int i : this->bcs ) {
        auto bc = dynamic_cast<PrescribedGradientHomogenization*>(this->giveDomain()->giveBc(i));
        bc->setCenterCoordinate(x);
    }
}


DofManager *PrescribedGradientMultiple::giveInternalDofManager( int i )
{
    for ( int j: this->bcs ) {
        auto bc = dynamic_cast<ActiveBoundaryCondition*>(this->giveDomain()->giveBc(j));
        if ( bc != nullptr ) {
            return bc->giveInternalDofManager(i);
        }
    }
    return nullptr;
}


void PrescribedGradientMultiple :: scale(double s)
{
    this->mGradient.times(s);
}


void PrescribedGradientMultiple :: computeField(FloatArray &sigma, TimeStep *tStep)
{
    FloatArray tmp;
    sigma.clear();
    for ( int i : this->bcs ) {
        auto bc = dynamic_cast<PrescribedGradientHomogenization*>(this->giveDomain()->giveBc(i));
        bc->computeField(tmp, tStep);
        sigma.add(tmp);
    }
}


void PrescribedGradientMultiple :: computeTangent(FloatMatrix &tangent, TimeStep *tStep)
{
    OOFEM_ERROR("Not implemented")
}

} /* namespace oofem */
