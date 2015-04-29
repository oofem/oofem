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

#include "meshqualityerrorestimator.h"
#include "element.h"
#include "elementgeometrytype.h"
#include "mathfem.h"
#include "node.h"
#include "integrationrule.h"
#include "feinterpol.h"
#include "gausspoint.h"

namespace oofem {
double MeshQualityErrorEstimator :: giveElementError(EE_ErrorType type, Element *elem, TimeStep *tStep)
{
    // A good plan would probably be to let the element determines its own quality if they implement some interface.
    // otherwise use a sane default.
    double error;
    FEInterpolation *fei = elem->giveInterpolation();
    IntegrationRule *ir = elem->giveDefaultIntegrationRulePtr();
    if ( fei && ir ) {
        error = this->computeJacobianError(* fei, * ir, elem);
    } else {
        switch ( elem->giveGeometryType() ) {
        case EGT_triangle_1:
            error = this->computeTriangleRadiusError(elem);
            break;
        case EGT_triangle_2:
            error = this->computeTriangleRadiusError(elem);
            break;
        default:
            error = 0.0;
            break;
        }
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
    a = c1->distance(* c2);
    b = c1->distance(* c3);
    c = c2->distance(* c3);
    return a * b * c / ( ( b + c - a ) * ( a + c - b ) * ( a + b - c ) ) - 1.0;
    // Reciprocal error would be;
    // (b+c-a)*(a+c-b)*(a+b-c)/(a*b*c - (b+c-a)*(a+c-b)*(a+b-c));
    // Which is safe except for when all points coincide, i.e. a = b = c = 0
}

double MeshQualityErrorEstimator :: computeJacobianError(FEInterpolation &fei, IntegrationRule &ir, Element *elem)
{
    double min_rcond = 1.0, rcond;
    FloatMatrix jac;

    for ( GaussPoint *gp: ir ) {
        fei.giveJacobianMatrixAt( jac, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(elem) );
        rcond = jac.computeReciprocalCondition() * sgn( jac.giveDeterminant() ); // Signed rcond. as inverted mappings are particularly bad.
        if ( rcond < min_rcond ) {
            min_rcond = rcond;
        }
    }
    return min_rcond < 1e-6 ? 1e6 : 1.0 / min_rcond; // Cap it to avoid overflow.
}

double MeshQualityErrorEstimator :: giveValue(EE_ValueType type, TimeStep *tStep)
{
    double error = 0.0, temp;
    for ( auto &elem : this->domain->giveElements() ) {
        temp = this->giveElementError(unknownET, elem.get(), tStep);
        if ( temp > error ) {
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
