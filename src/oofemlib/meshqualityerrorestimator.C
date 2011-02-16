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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#include "meshqualityerrorestimator.h"
#include "element.h"
#include "elementgeometrytype.h"
#include "mathfem.h"

#include "gaussintegrationrule.h"
#include "gausspnt.h"
#include "feinterpol.h"
#include "fei2dtrlin.h"
#include "fei2dtrquad.h"
#include "fei2dquadlin.h"
#include "fei2dquadquad.h"
#include "fei3dtrlin.h"
//#include "fei3dtrquad.h"
#include "fei3dhexalin.h"
#include "fei3dhexaquad.h"


namespace oofem {

// Set up the FEI.
FEI2dTrLin MeshQualityErrorEstimator :: fei2dtrlin(1, 2);
FEI2dTrQuad MeshQualityErrorEstimator :: fei2dtrquad(1, 2);
FEI2dQuadLin MeshQualityErrorEstimator :: fei2dquadlin(1, 2);
FEI2dQuadQuad MeshQualityErrorEstimator :: fei2dquadquad(1, 2);


FEI3dTrLin MeshQualityErrorEstimator :: fei3dtetlin;
//FEI3dTrQuad MeshQualityErrorEstimator :: fei3dtetquad;
FEI3dHexaLin MeshQualityErrorEstimator :: fei3dhexalin;
FEI3dHexaQuad MeshQualityErrorEstimator :: fei3dhexaquad;


double MeshQualityErrorEstimator :: giveElementError(EE_ErrorType type, Element *elem, TimeStep *tStep)
{
    // A good plan would probably be to let the element determines its own quality if they implement some interface.
    // otherwise use a sane default.
    double error;

    switch (elem->giveGeometryType()) {
        case EGT_triangle_1:
            {
                GaussIntegrationRule ir(1, elem, 1, 3);
                ir.setUpIntegrationPoints(_Triangle, 1, _2dFlow);
                error = computeJacobianError(MeshQualityErrorEstimator::fei2dtrlin, ir, elem);
            }
            // Alternative;
            // error = this->computeTriangleRadiusError(elem);
            break;
        case EGT_triangle_2:
            {
                GaussIntegrationRule ir(1, elem, 1, 3);
                ir.setUpIntegrationPoints(_Triangle, 4, _2dFlow);
                error = computeJacobianError(MeshQualityErrorEstimator::fei2dtrquad, ir, elem);
            }
            break;
        // EGT_line_2, EGT_quad, EGT_quad2, tetra1, tetra2, hexa, hexa_2, Composite, unknown?
        // EGT_point EGT_line_1 EGT_unknown
        default:
            error = 0.0;
            break;
    }
    return error;
}

double MeshQualityErrorEstimator :: computeTriangleRadiusError(Element *elem)
{
    // Outside/inside circle radius fraction based for quality measurement.
    // Zero for a perfect triangle,
    double a, b, c;
    FloatArray *c1, *c2, *c3;
    c1 = elem->giveNode(1)->giveCoordinates();
    c2 = elem->giveNode(2)->giveCoordinates();
    c3 = elem->giveNode(3)->giveCoordinates();
    a = c1->distance(*c2);
    b = c1->distance(*c3);
    c = c2->distance(*c3);
    return a*b*c/((b+c-a)*(a+c-b)*(a+b-c)) - 1.0;
    // Reciprocal error would be;
    // (b+c-a)*(a+c-b)*(a+b-c)/(a*b*c - (b+c-a)*(a+c-b)*(a+b-c));
    // Which is safe except for when all points coincide, i.e. a = b = c = 0
}

double MeshQualityErrorEstimator :: computeJacobianError(FEInterpolation &fei, IntegrationRule &ir, Element *elem)
{
    double min_rcond = 1.0, rcond;
    FloatMatrix jac;

    FEIElementGeometryWrapper egw(elem);
    for (int i = 1; i < ir.getNumberOfIntegrationPoints(); ++i) {
        fei.giveJacobianMatrixAt(jac, *ir.getIntegrationPoint(i)->giveCoordinates(), FEIElementGeometryWrapper(elem));
        rcond = jac.computeReciprocalCondition()*sgn(jac.giveDeterminant()); // Signed rcond. as inverted mappings are particularly bad.
        if (rcond < min_rcond) {
            min_rcond = rcond;
        }
    }
    return min_rcond < 0.001 ? 1000.0 : 1.0/min_rcond; // Cap it to avoid overflow.
}

double MeshQualityErrorEstimator :: giveValue(EE_ValueType type, TimeStep *tStep)
{
    double error = 0.0, temp;
    for (int i = 1; i <= this->domain->giveNumberOfElements(); ++i) {
        temp = this->giveElementError(unknownET, this->domain->giveElement(i), tStep);
        if (temp > error) {
            error = temp;
        }
    }
    return error;
}

int MeshQualityErrorEstimator :: estimateError(EE_ErrorMode mode, TimeStep *tStep)
{
    return true;
}

IRResultType MeshQualityErrorEstimator :: initializeFrom(InputRecord *ir)
{
    return ErrorEstimator :: initializeFrom(ir);
}
} // end namespace oofem
