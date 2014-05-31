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

#include "qspace.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "structuralcrosssection.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(QSpace);

FEI3dHexaQuad QSpace :: interpolation;

QSpace :: QSpace(int n, Domain *aDomain) : NLStructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(this)
{
    numberOfDofMans = 20;
}


IRResultType
QSpace :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = 27;

    this->matRotation = ir->hasField(_IFT_QSpace_materialCoordinateSystem);

    return this->NLStructuralElement :: initializeFrom(ir);
}


void
QSpace :: giveMaterialOrientationAt(FloatArray &x, FloatArray &y, FloatArray &z, const FloatArray &lcoords)
{
    FloatMatrix jac;
    FloatArray help;
    this->interpolation.giveJacobianMatrixAt( jac, lcoords, FEIElementGeometryWrapper(this) );
    x.beColumnOf(jac, 1); // This is {dx/dxi, dy/dxi, dz/dxi}
    x.normalize();
    help.beColumnOf(jac, 2);
    z.beVectorProductOf(x, help); // Normal to the xi-eta plane.
    z.normalize();
    y.beVectorProductOf(z, x);
}


void
QSpace :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    answer = {D_u, D_v, D_w};
}


MaterialMode
QSpace :: giveMaterialMode()
{
    return _3dMat;
}


void
QSpace :: computeStressVector(FloatArray &answer, const FloatArray &e, GaussPoint *gp, TimeStep *tStep)
{
    ///@todo This whole function should be placed in the "Space3dStructuralElementEvaluator".
    if ( this->matRotation ) {
        ///@todo This won't work properly with "useUpdatedGpRecord" (!)
        FloatArray x, y, z;
        FloatArray rotStrain, s;

        this->giveMaterialOrientationAt( x, y, z, * gp->giveCoordinates() );
        // Transform from global c.s. to material c.s.
#if 0
        rotStrain = {
            e(0) * x(0) * x(0) + e(5) * x(0) * x(1) + e(1) * x(1) * x(1) + e(4) * x(0) * x(2) + e(3) * x(1) * x(2) + e(2) * x(2) * x(2),
            e(0) * y(0) * y(0) + e(5) * y(0) * y(1) + e(1) * y(1) * y(1) + e(4) * y(0) * y(2) + e(3) * y(1) * y(2) + e(2) * y(2) * y(2),
            e(0) * z(0) * z(0) + e(5) * z(0) * z(1) + e(1) * z(1) * z(1) + e(4) * z(0) * z(2) + e(3) * z(1) * z(2) + e(2) * z(2) * z(2),
            2 * e(0) * y(0) * z(0) + e(4) * y(2) * z(0) + 2 * e(1) * y(1) * z(1) + e(3) * y(2) * z(1) + e(5) * ( y(1) * z(0) + y(0) * z(1) ) + ( e(4) * y(0) + e(3) * y(1) + 2 * e(2) * y(2) ) * z(2),
            2 * e(0) * x(0) * z(0) + e(4) * x(2) * z(0) + 2 * e(1) * x(1) * z(1) + e(3) * x(2) * z(1) + e(5) * ( x(1) * z(0) + x(0) * z(1) ) + ( e(4) * x(0) + e(3) * x(1) + 2 * e(2) * x(2) ) * z(2),
            2 * e(0) * x(0) * y(0) + e(4) * x(2) * y(0) + 2 * e(1) * x(1) * y(1) + e(3) * x(2) * y(1) + e(5) * ( x(1) * y(0) + x(0) * y(1) ) + ( e(4) * x(0) + e(3) * x(1) + 2 * e(2) * x(2) ) * y(2)
        };
#else
        rotStrain.resize(6);
        rotStrain.at(1) = e(0) * x(0) * x(0) + e(5) * x(0) * x(1) + e(1) * x(1) * x(1) + e(4) * x(0) * x(2) + e(3) * x(1) * x(2) + e(2) * x(2) * x(2);
        rotStrain.at(2) = e(0) * y(0) * y(0) + e(5) * y(0) * y(1) + e(1) * y(1) * y(1) + e(4) * y(0) * y(2) + e(3) * y(1) * y(2) + e(2) * y(2) * y(2);
        rotStrain.at(3) = e(0) * z(0) * z(0) + e(5) * z(0) * z(1) + e(1) * z(1) * z(1) + e(4) * z(0) * z(2) + e(3) * z(1) * z(2) + e(2) * z(2) * z(2);
        rotStrain.at(4) = 2 * e(0) * y(0) * z(0) + e(4) * y(2) * z(0) + 2 * e(1) * y(1) * z(1) + e(3) * y(2) * z(1) + e(5) * ( y(1) * z(0) + y(0) * z(1) ) + ( e(4) * y(0) + e(3) * y(1) + 2 * e(2) * y(2) ) * z(2);
        rotStrain.at(5) = 2 * e(0) * x(0) * z(0) + e(4) * x(2) * z(0) + 2 * e(1) * x(1) * z(1) + e(3) * x(2) * z(1) + e(5) * ( x(1) * z(0) + x(0) * z(1) ) + ( e(4) * x(0) + e(3) * x(1) + 2 * e(2) * x(2) ) * z(2);
        rotStrain.at(6) = 2 * e(0) * x(0) * y(0) + e(4) * x(2) * y(0) + 2 * e(1) * x(1) * y(1) + e(3) * x(2) * y(1) + e(5) * ( x(1) * y(0) + x(0) * y(1) ) + ( e(4) * x(0) + e(3) * x(1) + 2 * e(2) * x(2) ) * y(2);
#endif
        this->giveStructuralCrossSection()->giveRealStress_3d(s, gp, rotStrain, tStep);
#if 0
        answer = {
            s(0) * x(0) * x(0) + 2 * s(5) * x(0) * y(0) + s(1) * y(0) * y(0) + 2 * ( s(4) * x(0) + s(3) * y(0) ) * z(0) + s(2) * z(0) * z(0),
            s(0) * x(1) * x(1) + 2 * s(5) * x(1) * y(1) + s(1) * y(1) * y(1) + 2 * ( s(4) * x(1) + s(3) * y(1) ) * z(1) + s(2) * z(1) * z(1),
            s(0) * x(2) * x(2) + 2 * s(5) * x(2) * y(2) + s(1) * y(2) * y(2) + 2 * ( s(4) * x(2) + s(3) * y(2) ) * z(2) + s(2) * z(2) * z(2),
            y(2) * ( s(5) * x(1) + s(1) * y(1) + s(3) * z(1) ) + x(2) * ( s(0) * x(1) + s(5) * y(1) + s(4) * z(1) ) + ( s(4) * x(1) + s(3) * y(1) + s(2) * z(1) ) * z(2),
            y(2) * ( s(5) * x(0) + s(1) * y(0) + s(3) * z(0) ) + x(2) * ( s(0) * x(0) + s(5) * y(0) + s(4) * z(0) ) + ( s(4) * x(0) + s(3) * y(0) + s(2) * z(0) ) * z(2),
            y(1) * ( s(5) * x(0) + s(1) * y(0) + s(3) * z(0) ) + x(1) * ( s(0) * x(0) + s(5) * y(0) + s(4) * z(0) ) + ( s(4) * x(0) + s(3) * y(0) + s(2) * z(0) ) * z(1)
        };
#else
        answer.resize(6);
        answer.at(1) = s(0) * x(0) * x(0) + 2 * s(5) * x(0) * y(0) + s(1) * y(0) * y(0) + 2 * ( s(4) * x(0) + s(3) * y(0) ) * z(0) + s(2) * z(0) * z(0);
        answer.at(2) = s(0) * x(1) * x(1) + 2 * s(5) * x(1) * y(1) + s(1) * y(1) * y(1) + 2 * ( s(4) * x(1) + s(3) * y(1) ) * z(1) + s(2) * z(1) * z(1);
        answer.at(3) = s(0) * x(2) * x(2) + 2 * s(5) * x(2) * y(2) + s(1) * y(2) * y(2) + 2 * ( s(4) * x(2) + s(3) * y(2) ) * z(2) + s(2) * z(2) * z(2);
        answer.at(4) = y(2) * ( s(5) * x(1) + s(1) * y(1) + s(3) * z(1) ) + x(2) * ( s(0) * x(1) + s(5) * y(1) + s(4) * z(1) ) + ( s(4) * x(1) + s(3) * y(1) + s(2) * z(1) ) * z(2);
        answer.at(5) = y(2) * ( s(5) * x(0) + s(1) * y(0) + s(3) * z(0) ) + x(2) * ( s(0) * x(0) + s(5) * y(0) + s(4) * z(0) ) + ( s(4) * x(0) + s(3) * y(0) + s(2) * z(0) ) * z(2);
        answer.at(6) = y(1) * ( s(5) * x(0) + s(1) * y(0) + s(3) * z(0) ) + x(1) * ( s(0) * x(0) + s(5) * y(0) + s(4) * z(0) ) + ( s(4) * x(0) + s(3) * y(0) + s(2) * z(0) ) * z(1);
#endif
    } else {
        this->giveStructuralCrossSection()->giveRealStress_3d(answer, gp, e, tStep);
    }
}


void
QSpace :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    ///@todo This whole function should be placed in the "Space3dStructuralElementEvaluator".
    this->giveStructuralCrossSection()->giveStiffnessMatrix_3d(answer, rMode, gp, tStep);
    if ( this->matRotation ) {
        FloatArray x, y, z;
        FloatMatrix Q;

        this->giveMaterialOrientationAt( x, y, z, * gp->giveCoordinates() );

#if 0
        Q = {
            { x(0) * x(0), x(1) * x(1), x(2) * x(2), x(1) * x(2), x(0) * x(2), x(0) * x(1) },
            { y(0) * y(0), y(1) * y(1), y(2) * y(2), y(1) * y(2), y(0) * y(2), y(0) * y(1) },
            { z(0) * z(0), z(1) * z(1), z(2) * z(2), z(1) * z(2), z(0) * z(2), z(0) * z(1) },
            { 2 * y(0) * z(0), 2 * y(1) * z(1), 2 * y(2) * z(2), y(2) * z(1) + y(1) * z(2), y(2) * z(0) + y(0) * z(2), y(1) * z(0) + y(0) * z(1) },
            { 2 * x(0) * z(0), 2 * x(1) * z(1), 2 * x(2) * z(2), x(2) * z(1) + x(1) * z(2), x(2) * z(0) + x(0) * z(2), x(1) * z(0) + x(0) * z(1) },
            { 2 * x(0) * y(0), 2 * x(1) * y(1), 2 * x(2) * y(2), x(2) * y(1) + x(1) * y(2), x(2) * y(0) + x(0) * y(2), x(1) * y(0) + x(0) * y(1) }
        };
#else
        Q.resize(6, 6);
        Q.at(1, 1) = x(0) * x(0);
        Q.at(2, 1) = x(1) * x(1);
        Q.at(3, 1) = x(2) * x(2);
        Q.at(4, 1) = x(1) * x(2);
        Q.at(5, 1) = x(0) * x(2);
        Q.at(6, 1) = x(0) * x(1);

        Q.at(1, 2) = y(0) * y(0);
        Q.at(2, 2) = y(1) * y(1);
        Q.at(3, 2) = y(2) * y(2);
        Q.at(4, 2) = y(1) * y(2);
        Q.at(5, 2) = y(0) * y(2);
        Q.at(6, 2) = y(0) * y(1);

        Q.at(1, 3) = z(0) * z(0);
        Q.at(2, 3) = z(1) * z(1);
        Q.at(3, 3) = z(2) * z(2);
        Q.at(4, 3) = z(1) * z(2);
        Q.at(5, 3) = z(0) * z(2);
        Q.at(6, 3) = z(0) * z(1);

        Q.at(1, 4) = 2 * y(0) * z(0);
        Q.at(2, 4) = 2 * y(1) * z(1);
        Q.at(3, 4) = 2 * y(2) * z(2);
        Q.at(4, 4) = y(2) * z(1) + y(1) * z(2);
        Q.at(5, 4) = y(2) * z(0) + y(0) * z(2);
        Q.at(6, 4) = y(1) * z(0) + y(0) * z(1);

        Q.at(1, 5) = 2 * x(0) * z(0);
        Q.at(2, 5) = 2 * x(1) * z(1);
        Q.at(3, 5) = 2 * x(2) * z(2);
        Q.at(4, 5) = x(2) * z(1) + x(1) * z(2);
        Q.at(5, 5) = x(2) * z(0) + x(0) * z(2);
        Q.at(6, 5) = x(1) * z(0) + x(0) * z(1);

        Q.at(1, 6) = 2 * x(0) * y(0);
        Q.at(2, 6) = 2 * x(1) * y(1);
        Q.at(3, 6) = 2 * x(2) * y(2);
        Q.at(4, 6) = x(2) * y(1) + x(1) * y(2);
        Q.at(5, 6) = x(2) * y(0) + x(0) * y(2);
        Q.at(6, 6) = x(1) * y(0) + x(0) * y(1);
#endif
        //printf("original = "); answer.printYourself();
        answer.rotatedWith(Q, 't');
        //printf("rotated = "); answer.printYourself();
        //OOFEM_ERROR("DEBUG QUIT");
    }
}


double
QSpace :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant = fabs( this->interpolation.giveTransformationJacobian( * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ) );
    double weight      = gp->giveWeight();

    return ( determinant * weight );
}


double
QSpace :: giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane)
{
    double factor = cbrt( ( double ) gp->giveIntegrationRule()->giveNumberOfIntegrationPoints() );
    return this->giveLenghtInDir(normalToCrackPlane) / factor;
}


void
QSpace :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 6);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


void
QSpace :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [6x60] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
// B matrix  -  6 rows : epsilon-X, epsilon-Y, epsilon-Z, gamma-YZ, gamma-ZX, gamma-XY  :
{
    FloatMatrix dnx;

    this->interpolation.evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(6, 60);
    answer.zero();

    for ( int i = 1; i <= 20; i++ ) {
        answer.at(1, 3 * i - 2) = dnx.at(i, 1);
        answer.at(2, 3 * i - 1) = dnx.at(i, 2);
        answer.at(3, 3 * i - 0) = dnx.at(i, 3);

        answer.at(4, 3 * i - 1) = dnx.at(i, 3);
        answer.at(4, 3 * i - 0) = dnx.at(i, 2);

        answer.at(5, 3 * i - 2) = dnx.at(i, 3);
        answer.at(5, 3 * i - 0) = dnx.at(i, 1);

        answer.at(6, 3 * i - 2) = dnx.at(i, 2);
        answer.at(6, 3 * i - 1) = dnx.at(i, 1);
    }
}


void
QSpace :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    FloatMatrix dnx;

    this->interpolation.evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(9, 60);
    answer.zero();

    for ( int i = 1; i <= dnx.giveNumberOfRows(); i++ ) {
        answer.at(1, 3 * i - 2) = dnx.at(i, 1);     // du/dx
        answer.at(2, 3 * i - 1) = dnx.at(i, 2);     // dv/dy
        answer.at(3, 3 * i - 0) = dnx.at(i, 3);     // dw/dz
        answer.at(4, 3 * i - 1) = dnx.at(i, 3);     // dv/dz
        answer.at(7, 3 * i - 0) = dnx.at(i, 2);     // dw/dy
        answer.at(5, 3 * i - 2) = dnx.at(i, 3);     // du/dz
        answer.at(8, 3 * i - 0) = dnx.at(i, 1);     // dw/dx
        answer.at(6, 3 * i - 2) = dnx.at(i, 2);     // du/dy
        answer.at(9, 3 * i - 1) = dnx.at(i, 1);     // dv/dx
    }
}


// ******************************
// ***  Surface load support  ***
// ******************************

IntegrationRule *
QSpace :: GetSurfaceIntegrationRule(int approxOrder)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, this, 1, 1);
    int npoints = iRule->getRequiredNumberOfIntegrationPoints(_Square, approxOrder);
    iRule->SetUpPointsOnSquare(npoints, _Unknown);
    return iRule;
}

void
QSpace :: computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *sgp)
{
    FloatArray n;
    interpolation.surfaceEvalN( n, iSurf, * sgp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(n, 3);
}

void
QSpace :: giveSurfaceDofMapping(IntArray &answer, int iSurf) const
{
    IntArray nodes;
    const int ndofsn = 3;

    interpolation.computeLocalSurfaceMapping(nodes, iSurf);

    answer.resize(24);

    for ( int i = 1; i <= 8; i++ ) {
        answer.at(i * ndofsn - 2) = nodes.at(i) * ndofsn - 2;
        answer.at(i * ndofsn - 1) = nodes.at(i) * ndofsn - 1;
        answer.at(i * ndofsn) = nodes.at(i) * ndofsn;
    }
}

double
QSpace :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    double determinant, weight, volume;
    determinant = fabs( interpolation.surfaceGiveTransformationJacobian( iSurf, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ) );

    weight = gp->giveWeight();
    volume = determinant * weight;

    return volume;
}

void
QSpace :: computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iSurf)
{
    interpolation.surfaceLocal2global( answer, iSurf, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
}

int
QSpace :: computeLoadLSToLRotationMatrix(FloatMatrix &answer, int iSurf, GaussPoint *gp)
{
    // returns transformation matrix from
    // surface local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)

    // definition of local c.s on surface:
    // local z axis - perpendicular to surface, pointing outwards from element
    // local x axis - is in global xy plane (perpendicular to global z axis)
    // local y axis - completes the righ hand side cs.

    /*
     * OOFEM_ERROR("surface local coordinate system not supported");
     * return 1;
     */
    FloatArray gc(3);
    FloatArray h1(3), h2(3), nn(3), n(3);
    IntArray snodes(4);

    answer.resize(3, 3);
    answer.zero();

    this->interpolation.computeSurfaceMapping(snodes, dofManArray, iSurf);
    for ( int i = 1; i <= 4; i++ ) {
        gc.add( * domain->giveNode( snodes.at(i) )->giveCoordinates() );
    }

    gc.times(1. / 4.);
    // determine "average normal"
    for ( int i = 1; i <= 4; i++ ) {
        int j = ( i ) % 4 + 1;
        h1.beDifferenceOf(* domain->giveNode( snodes.at(i) )->giveCoordinates(), gc);
        h2.beDifferenceOf(* domain->giveNode( snodes.at(j) )->giveCoordinates(), gc);
        n.beVectorProductOf(h1, h2);
        if ( n.computeSquaredNorm() > 1.e-6 ) {
            n.normalize();
        }

        nn.add(n);
    }

    nn.times(1. / 4.);
    if ( nn.computeSquaredNorm() < 1.e-6 ) {
        answer.zero();
    }

    nn.normalize();
    for ( int i = 1; i <= 3; i++ ) {
        answer.at(i, 3) = nn.at(i);
    }

    // determine lcs of surface
    // local x axis in xy plane
    double test = fabs(fabs( nn.at(3) ) - 1.0);
    if ( test < 1.e-5 ) {
        h1.at(1) = answer.at(1, 1) = 1.0;
        h1.at(2) = answer.at(2, 1) = 0.0;
    } else {
        h1.at(1) = answer.at(1, 1) = answer.at(2, 3);
        h1.at(2) = answer.at(2, 1) = -answer.at(1, 3);
    }

    h1.at(3) = answer.at(3, 1) = 0.0;
    // local y axis perpendicular to local x,z axes
    h2.beVectorProductOf(nn, h1);
    for ( int i = 1; i <= 3; i++ ) {
        answer.at(i, 2) = h2.at(i);
    }

    return 1;
}

Interface *
QSpace :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    }

    OOFEM_LOG_INFO("Interface on Qspace element not supported");
    return NULL;
}

void
QSpace :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(20);
    for ( int i = 1; i <= 20; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
QSpace :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i <= 20; i++ ) {
        if ( this->giveNode(i)->giveNumber() == pap ) {
            found = 1;
        }
    }

    if ( found ) {
        answer.at(1) = pap;
    } else {
        OOFEM_ERROR("unknown node number %d", pap);
    }
}

int
QSpace :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return this->integrationRulesArray [ 0 ]->giveNumberOfIntegrationPoints();
}


SPRPatchType
QSpace :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_3dBiQuadratic;
}


void
QSpace :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    answer.clear();
    OOFEM_WARNING("IP values will not be transferred to nodes. Use ZZNodalRecovery instead (parameter stype 1)");
}

void
QSpace :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep)
{
    answer.clear();
}
} // end namespace oofem
