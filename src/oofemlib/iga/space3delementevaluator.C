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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "space3delementevaluator.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "matresponsemode.h"
#include "crosssection.h"
#include "structuralcrosssection.h"
#include "mathfem.h"
#include "iga.h"

namespace oofem {
/* 3D Space Elements */
void Space3dStructuralElementEvaluator :: computeNMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    int i, nDofMan;
    FloatArray N;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();

    interp->evalN(N, * gp->giveCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ));

    if ( ( nDofMan = interp->giveNumberOfKnotSpanBasisFunctions( * ( gp->giveIntegrationRule()->giveKnotSpan() ) ) ) == 0 ) { // HUHU
        nDofMan = gp->giveElement()->giveNumberOfDofManagers();
    }

    answer.resize(3, nDofMan * 3);
    answer.zero();

    for ( i = 1; i <= nDofMan; i++ ) {
        answer.at(1, i * 3 - 2) = N.at(i);
        answer.at(2, i * 3 - 1) = N.at(i);
        answer.at(3, i * 3 - 0) = N.at(i);
    }
}

void Space3dStructuralElementEvaluator :: computeBMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    int i, nDofMan;
    //IntArray dofmanSubElementMask;
    FloatMatrix d;

    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    // this uses FEIInterpolation::nodes2coords - quite inefficient in this case (large num of dofmans)
    interp->evaldNdx(d, * gp->giveCoordinates(),
                     FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ));

    if ( ( nDofMan = interp->giveNumberOfKnotSpanBasisFunctions( * ( gp->giveIntegrationRule()->giveKnotSpan() ) ) ) == 0 ) { // HUHU
        nDofMan = gp->giveElement()->giveNumberOfDofManagers();
    }

    answer.resize(6, nDofMan * 3);
    answer.zero();

    for ( i = 1; i <= nDofMan; i++ ) {
        answer.at(1, i * 3 - 2) = d.at(i, 1);
        answer.at(2, i * 3 - 1) = d.at(i, 2);
        answer.at(3, i * 3 - 0) = d.at(i, 3);

        answer.at(4, 3 * i - 1) = d.at(i, 3);
        answer.at(4, 3 * i - 0) = d.at(i, 2);

        answer.at(5, 3 * i - 2) = d.at(i, 3);
        answer.at(5, 3 * i - 0) = d.at(i, 1);

        answer.at(6, 3 * i - 2) = d.at(i, 2);
        answer.at(6, 3 * i - 1) = d.at(i, 1);
    }
}

double Space3dStructuralElementEvaluator :: computeVolumeAround(GaussPoint *gp)
{
    double determinant, weight, volume;
    determinant = fabs( this->giveElement()->giveInterpolation()
                       ->giveTransformationJacobian(* gp->giveCoordinates(),
                                                    FEIIGAElementGeometryWrapper( this->giveElement(),
                                                                                 gp->giveIntegrationRule()->giveKnotSpan() )) );
    weight      = gp->giveWeight();
    volume      = determinant * weight;

    return volume;
}
} // end namespace oofem
