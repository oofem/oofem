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

#include "../sm/Elements/Interfaces/intelline2.h"
#include "../sm/CrossSections/structuralinterfacecrosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "fei2dlinequad.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Element(IntElLine2);

FEI2dLineQuad IntElLine2 :: interp(2, 2);


IntElLine2 :: IntElLine2(int n, Domain *aDomain) : IntElLine1(n, aDomain)
{
    numberOfDofMans = 6;
}


void
IntElLine2 :: computeNmatrixAt(GaussPoint *ip, FloatMatrix &answer)
{
    // Returns the modified N-matrix which multiplied with u give the spatial jump.
    FloatArray N;
    interp.evalN( N, * ip->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(2, 12);
    answer.zero();
    answer.at(1, 1) = answer.at(2, 2) = -N.at(1);
    answer.at(1, 3) = answer.at(2, 4) = -N.at(2);
    answer.at(1, 5) = answer.at(2, 6) = -N.at(3);

    answer.at(1, 7) = answer.at(2, 8) = N.at(1);
    answer.at(1, 9) = answer.at(2, 10) = N.at(2);
    answer.at(1, 11) = answer.at(2, 12) = N.at(3);
}


void
IntElLine2 :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        //integrationRulesArray[0] = new LobattoIntegrationRule (1,domain, 1, 2); ///@todo - should be able to decide
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 2);
        integrationRulesArray [ 0 ]->SetUpPointsOnLine(4, _2dInterface); ///@todo - should be a parameter with num of ip
    }
}

FEInterpolation *
IntElLine2 :: giveInterpolation() const
{
    return & interp;
}


IRResultType
IntElLine2 :: initializeFrom(InputRecord *ir)
{
    return IntElLine1 :: initializeFrom(ir);   
}

} // end namespace oofem
