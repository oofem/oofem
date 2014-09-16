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

#include "neumannmomentload.h"
#include "classfactory.h"
#include "function.h"
#include "inputrecord.h"
#include "domain.h"
#include "set.h"
#include "element.h"
#include "feinterpol.h"
#include "gausspoint.h"

namespace oofem {
REGISTER_BoundaryCondition(NeumannMomentLoad);

IRResultType
NeumannMomentLoad :: initializeFrom(InputRecord *ir)
{
    BoundaryLoad :: initializeFrom(ir);
    if ( componentArray.giveSize() != nDofs ) {
        OOFEM_ERROR("componentArray size mismatch");
    }

    IRResultType result;

    IR_GIVE_FIELD(ir, g, _IFT_NeumannMomentLoad_Gradient);
    p = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, p, _IFT_NeumannMomentLoad_Constant);
    IR_GIVE_FIELD(ir, cset, _IFT_NeumannMomentLoad_CenterSet);

    xbar.resize(0);

    return result;
}

void
NeumannMomentLoad :: computeXbar()
{

    if (!xbar.isEmpty()) return;

    xbar.resize(this->giveDomain()->giveNumberOfSpatialDimensions());
    xbar.zero();

    celements = this->giveDomain()->giveSet(cset)->giveElementList();

    for ( auto elementID : celements ) {

        Element *thisElement = this->giveDomain()->giveElement(elementID);
        FEInterpolation *i = thisElement->giveInterpolation();

        IntegrationRule *iRule = i->giveIntegrationRule(2);

        for ( GaussPoint * gp: * iRule ) {
            FloatArray coord;
            FloatArray *lcoords = gp->giveNaturalCoordinates();
            double detJ = i->giveTransformationJacobian(*lcoords, FEIElementGeometryWrapper(thisElement));

            i->local2global(coord, *lcoords, FEIElementGeometryWrapper(thisElement));
            coord.times(gp->giveWeight()*fabs(detJ));

            printf("(%f, %f, %f), detJ = %f\n", coord.at(1), coord.at(2), coord.at(3), detJ);

            xbar.add(coord);
        }
    }

    xbar.printYourself();

}

void
NeumannMomentLoad :: computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode)
{
    // we overload general implementation on the boundary load level due
    // to implementation efficiency

    computeXbar();

    double factor;

    if ( ( mode != VM_Total ) && ( mode != VM_Incremental ) ) {
        OOFEM_ERROR("mode not supported");
    }

    // ask time distribution

    /*
     * factor = this -> giveTimeFunction() -> at(tStep->giveTime()) ;
     * if ((mode==VM_Incremental) && (!tStep->isTheFirstStep()))
     * //factor -= this->giveTimeFunction()->at(tStep->givePreviousStep()->giveTime()) ;
     * factor -= this->giveTimeFunction()->at(tStep->giveTime()-tStep->giveTimeIncrement()) ;
     */
    factor = this->giveTimeFunction()->evaluate(tStep, mode);
    answer = componentArray;
    answer.times(factor);
}

}
