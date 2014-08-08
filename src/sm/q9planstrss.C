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

#include "q9planstrss.h"
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
    NLStructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(this), NodalAveragingRecoveryModelInterface()
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
Q9PlaneStress2d :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [3x18] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    FloatMatrix dnx;

    this->interpolation.evaldNdx( dnx, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(3, 18);
    answer.zero();

    for ( int i = 1; i <= 9; i++ ) {
        answer.at(1, 2 * i - 1) = dnx.at(i, 1);
        answer.at(2, 2 * i - 0) = dnx.at(i, 2);

        answer.at(3, 2 * i - 1) = dnx.at(i, 2);
        answer.at(3, 2 * i - 0) = dnx.at(i, 1);
    }
}

void
Q9PlaneStress2d :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver,
// evaluated at gp.
{
    FloatArray n;

    this->interpolation.evalN( n, iLocCoord, FEIElementGeometryWrapper(this) );

    answer.beNMatrixOf(n, 2);
}

IRResultType
Q9PlaneStress2d :: initializeFrom(InputRecord *ir)
{
    return Element :: initializeFrom(ir);
}

void
Q9PlaneStress2d :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

double
Q9PlaneStress2d :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant, weight, thickness, volume;
    determinant = fabs( this->interpolation.giveTransformationJacobian( * gp->giveNaturalCoordinates(),
                                                                       FEIElementGeometryWrapper(this) ) );
    weight      = gp->giveWeight();
    thickness   = this->giveCrossSection()->give(CS_Thickness, gp);
    volume      = determinant * weight * thickness;

    return volume;
}


void
Q9PlaneStress2d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v};
}


double
Q9PlaneStress2d :: giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane)
{
    if ( normalToCrackPlane.at(3) < 0.999999 ) { //ensure that characteristic length is in the plane of element
        return this->giveLenghtInDir(normalToCrackPlane) / sqrt( ( double ) gp->giveIntegrationRule()->giveNumberOfIntegrationPoints() );
    } else { //otherwise compute out-of-plane characteristic length from element area
        return sqrt( this->computeVolumeAreaOrLength() / ( double ) gp->giveIntegrationRule()->giveNumberOfIntegrationPoints() );
    }
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


void
Q9PlaneStress2d :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
{
    /*
     *
     * computes interpolation matrix for element edge.
     * we assemble locally this matrix for only nonzero
     * shape functions.
     * (for example only two nonzero shape functions for 2 dofs are
     * necessary for linear plane stress tringle edge).
     * These nonzero shape functions are then mapped to
     * global element functions.
     *
     * Using mapping technique will allow to assemble shape functions
     * without regarding particular side
     */

    FloatArray n(3);
    this->interpolation.edgeEvalN( n, iedge, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(2, 6);
    answer.zero();

    answer.at(1, 1) = n.at(1);
    answer.at(1, 3) = n.at(2);
    answer.at(1, 5) = n.at(3);
    answer.at(2, 2) = n.at(1);
    answer.at(2, 4) = n.at(2);
    answer.at(2, 6) = n.at(3);
}


void
Q9PlaneStress2d :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    IntArray eNodes(3);
    this->interpolation.computeLocalEdgeMapping(eNodes,  iEdge);

    answer.resize(6);
    for ( int i = 1; i <= 3; i++ ) {
        answer.at(i * 2 - 1) = eNodes.at(i) * 2 - 1;
        answer.at(i * 2) = eNodes.at(i) * 2;
    }
}

double
Q9PlaneStress2d ::   computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double result = this->interpolation.edgeGiveTransformationJacobian( iEdge, * gp->giveNaturalCoordinates(),
                                                                       FEIElementGeometryWrapper(this) );
    return result *gp->giveWeight();
}

void
Q9PlaneStress2d :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    this->interpolation.edgeLocal2global( answer, iEdge, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
}


int
Q9PlaneStress2d :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    // returns transformation matrix from
    // edge local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)
    //
    FloatArray normal(2);
    answer.resize(2, 2);
    answer.zero();

    this->interpolation.edgeEvalNormal( normal, iEdge, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.at(1, 1) = normal.at(2);
    answer.at(1, 2) = normal.at(1);
    answer.at(2, 1) = -normal.at(1);
    answer.at(2, 2) = normal.at(2);

    return 1;
}
} // end namespace oofem
