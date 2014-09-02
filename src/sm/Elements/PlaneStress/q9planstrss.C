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

#include "Elements/PlaneStress/q9planstrss.h"
#include "fei2dquadbiquad.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(Q9PlaneStress2d);

FEI2dQuadBiQuad Q9PlaneStress2d :: interpolation(1, 2);

Q9PlaneStress2d :: Q9PlaneStress2d(int n, Domain *aDomain) :
    PlaneStressElement(n, aDomain), ZZNodalRecoveryModelInterface(this), NodalAveragingRecoveryModelInterface()
{
    numberOfDofMans  = 9;
    numberOfGaussPoints = 4;
}


FEInterpolation *
Q9PlaneStress2d :: giveInterpolation() const { return & interpolation; }


Interface *
Q9PlaneStress2d :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    }

    return NULL;
}


void
Q9PlaneStress2d :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                              InternalStateType type, TimeStep *tStep)
{
    if ( numberOfGaussPoints != 4 ) {
        return;
    }

    GaussPoint *gp;

    if ( node < 5 ) {
        int i = 0;
        switch ( node ) {
        case 1: i = 4;
            break;
        case 2: i = 2;
            break;
        case 3: i = 1;
            break;
        case 4: i = 3;
            break;
        }

        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(i - 1);
        this->giveIPValue(answer, gp, type, tStep);
    } else {
        int i1 = 0, i2 = 0;
        switch ( node ) {
        case 5: i1 = 4;
            i2 = 2;
            break;
        case 6: i1 = 2;
            i2 = 1;
            break;
        case 7: i1 = 1;
            i2 = 3;
            break;
        case 8: i1 = 3;
            i2 = 4;
            break;
        }

        FloatArray contrib;
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(i1 - 1);
        this->giveIPValue(contrib, gp, type, tStep);
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(i2 - 1);
        this->giveIPValue(answer, gp, type, tStep);
        answer.add(contrib);
        answer.times(0.5);
    }
}


} // end namespace oofem
