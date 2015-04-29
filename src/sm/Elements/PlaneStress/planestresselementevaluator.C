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

#include "../sm/Elements/PlaneStress/planestresselementevaluator.h"
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
void PlaneStressStructuralElementEvaluator :: computeNMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatArray N;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    interp->evalN( N, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    answer.beNMatrixOf(N, 2);
}

void PlaneStressStructuralElementEvaluator :: computeBMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix d;

    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    // this uses FEInterpolation::nodes2coords - quite inefficient in this case (large num of dofmans)
    interp->evaldNdx( d, gp->giveNaturalCoordinates(),
                     FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    answer.resize(3, d.giveNumberOfRows() * 2);
    answer.zero();

    for ( int i = 1; i <= d.giveNumberOfRows(); i++ ) {
        answer.at(1, i * 2 - 1) = d.at(i, 1);
        answer.at(2, i * 2 - 0) = d.at(i, 2);

        answer.at(3, 2 * i - 1) = d.at(i, 2);
        answer.at(3, 2 * i - 0) = d.at(i, 1);
    }
}


double PlaneStressStructuralElementEvaluator :: computeVolumeAround(GaussPoint *gp)
{
    double determinant, weight, thickness, volume;
    determinant = fabs( this->giveElement()->giveInterpolation()
                       ->giveTransformationJacobian( gp->giveNaturalCoordinates(),
                                                    FEIIGAElementGeometryWrapper( this->giveElement(),
                                                                                 gp->giveIntegrationRule()->giveKnotSpan() ) ) );
    weight      = gp->giveWeight();
    thickness   = this->giveElement()->giveCrossSection()->give(CS_Thickness, gp);
    volume      = determinant * weight * thickness;

    return volume;
}


void PlaneStressStructuralElementEvaluator :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    static_cast< StructuralCrossSection * >( this->giveElement()->giveCrossSection() )->giveRealStress_PlaneStress(answer, gp, strain, tStep);
}

void PlaneStressStructuralElementEvaluator :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    static_cast< StructuralCrossSection * >( this->giveElement()->giveCrossSection() )->giveStiffnessMatrix_PlaneStress(answer, rMode, gp, tStep);
}
} // end namespace oofem
