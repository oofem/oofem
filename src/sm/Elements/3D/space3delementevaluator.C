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

#include "../sm/Elements/3D/space3delementevaluator.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "matresponsemode.h"
#include "crosssection.h"
#include "mathfem.h"
#include "iga/iga.h"

namespace oofem {
/* 3D Space Elements */
void Space3dStructuralElementEvaluator :: computeNMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatArray N;
    Element *element = this->giveElement();
    FEInterpolation *interp = element->giveInterpolation();

    interp->evalN( N, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( element, gp->giveIntegrationRule()->giveKnotSpan() ) );

    answer.beNMatrixOf(N, 3);
}

void Space3dStructuralElementEvaluator :: computeBMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix d;
    Element *element = this->giveElement();
    FEInterpolation *interp = element->giveInterpolation();
    // this uses FEInterpolation::nodes2coords - quite inefficient in this case (large num of dofmans)
    interp->evaldNdx( d, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( element, gp->giveIntegrationRule()->giveKnotSpan() ) );


    answer.resize(6, d.giveNumberOfRows() * 3);
    answer.zero();

    for ( int i = 1; i <= d.giveNumberOfRows(); i++ ) {
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
    double determinant = fabs( this->giveElement()->giveInterpolation()
                              ->giveTransformationJacobian( gp->giveNaturalCoordinates(),
                                                           FEIIGAElementGeometryWrapper( this->giveElement(),
                                                                                        gp->giveIntegrationRule()->giveKnotSpan() ) ) );
    return determinant *gp->giveWeight();
}


void Space3dStructuralElementEvaluator :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    static_cast< StructuralCrossSection * >( this->giveElement()->giveCrossSection() )->giveRealStress_3d(answer, gp, strain, tStep);
}


void Space3dStructuralElementEvaluator :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    static_cast< StructuralCrossSection * >( this->giveElement()->giveCrossSection() )->giveStiffnessMatrix_3d(answer, rMode, gp, tStep);
}
} // end namespace oofem
