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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#include "sm/Elements/quad2platesubsoil.h"
#include "fei2dquadquad.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(Quad2PlateSubSoil);

FEI2dQuadQuad Quad2PlateSubSoil :: interp_quad(1, 2);

Quad2PlateSubSoil :: Quad2PlateSubSoil(int n, Domain *aDomain) :
    Quad1PlateSubSoil(n, aDomain)
{
    numberOfGaussPoints = 4;
    numberOfDofMans = 8;
}


FEInterpolation *
Quad2PlateSubSoil :: giveInterpolation() const { return & interp_quad; }


FEInterpolation *
Quad2PlateSubSoil :: giveInterpolation(DofIDItem id) const
{
    return & interp_quad;
}


void
Quad2PlateSubSoil :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this, 1, 5);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}



void
Quad2PlateSubSoil :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [3x8] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    FloatArray n;
    FloatMatrix dn;

    this->interp_quad.evaldNdx( dn, gp->giveNaturalCoordinates(),  FEIElementGeometryWrapper(this) );
    this->interp_quad.evalN( n, gp->giveNaturalCoordinates(),  FEIElementGeometryWrapper(this) );

    answer.resize(3, 8);
    answer.zero();

    ///@todo Check sign here
    for ( int i = 0; i < 8; ++i ) {
        answer(0, i) = n(i); // eps_z
        answer(1, i) = dn(i, 0); // gamma_xz
        answer(2, i) = dn(i, 1); // gamma_yz
    }
}


void
Quad2PlateSubSoil :: initializeFrom(InputRecord &ir)
{
    this->numberOfGaussPoints = 4;
    StructuralElement :: initializeFrom(ir);
}


void
Quad2PlateSubSoil :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(8);
    for ( int i = 1; i < 8; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
Quad2PlateSubSoil :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i <= 8; i++ ) {
        if ( pap == this->giveNode(i)->giveNumber() ) {
            found = 1;
        }
    }

    if ( found ) {
        answer.at(1) = pap;
    } else {
        OOFEM_ERROR("node %d not found on element %d", pap, this->giveNumber());
    }
}

void
Quad2PlateSubSoil ::computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the [1x8] displacement interpolation matrix {N}
{
    FloatArray N(8);
    giveInterpolation()->evalN(N, iLocCoord, FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(N, 1);
}

void
Quad2PlateSubSoil :: computeSurfaceNMatrix(FloatMatrix &answer, int boundaryID, const FloatArray &lcoords)
{
    if (boundaryID == 1) {
        this->computeNmatrixAt(lcoords, answer);
    } else {
        OOFEM_ERROR("computeSurfaceNMatrix: Only one surface is supported with id=1");
    }
}



} // end namespace oofem
