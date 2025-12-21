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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "sm/Elements/tria2platesubsoil.h"
#include "fei2dtrquad.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(Tria2PlateSubSoil);

FEI2dTrQuad Tria2PlateSubSoil :: interp_quad(1, 2);

Tria2PlateSubSoil :: Tria2PlateSubSoil(int n, Domain *aDomain) :
    Tria1PlateSubSoil(n, aDomain)
{
    numberOfGaussPoints = 4;
    numberOfDofMans = 6;
}


FEInterpolation *
Tria2PlateSubSoil :: giveInterpolation(DofIDItem id) const
{
    return & interp_quad;
}


FEInterpolation *
Tria2PlateSubSoil :: giveInterpolation() const { return & interp_quad; }


void
Tria2PlateSubSoil :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this, 1, 5);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}



void
Tria2PlateSubSoil :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [3x6] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    FloatArray n;
    FloatMatrix dn;

    this->interp_quad.evaldNdx( dn, gp->giveNaturalCoordinates(),  FEIElementGeometryWrapper(this) );
    this->interp_quad.evalN( n, gp->giveNaturalCoordinates(),  FEIElementGeometryWrapper(this) );

    answer.resize(3, 6);
    answer.zero();

    ///@todo Check sign here
    for ( int i = 0; i < 6; ++i ) {
      answer(0, i) = n(i); // eps_z
      answer(1, i) = dn(i, 0); // gamma_xz
      answer(2, i) = dn(i, 1); // gamma_yz
    }
}


//TODO ZZNodalRecoveryModel can not determine some values.  This is caused by sum of row entries is zero for (N^T)N matrix for vertices, yielding zero entries in the lumped form.

void
Tria2PlateSubSoil :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(3);
    pap.at(1) = this->giveNode(1)->giveNumber();
    pap.at(2) = this->giveNode(2)->giveNumber();
    pap.at(3) = this->giveNode(3)->giveNumber();
}

void
Tria2PlateSubSoil :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    answer.resize(3);
    if ( pap == this->giveNode(1)->giveNumber() ) {
        answer.at(1) = pap;
        answer.at(2) = this->giveNode(4)->giveNumber();
        answer.at(3) = this->giveNode(6)->giveNumber();
    } else if ( pap == this->giveNode(2)->giveNumber() ) {
        answer.at(1) = pap;
        answer.at(2) = this->giveNode(5)->giveNumber();
        answer.at(3) = this->giveNode(4)->giveNumber();
    } else if ( pap == this->giveNode(3)->giveNumber() ) {
        answer.at(1) = pap;
        answer.at(2) = this->giveNode(6)->giveNumber();
        answer.at(3) = this->giveNode(5)->giveNumber();
    } else {
        OOFEM_ERROR("node unknown");
    }
}

SPRPatchType
Tria2PlateSubSoil :: SPRNodalRecoveryMI_givePatchType()
{ 
    return SPRPatchType_2dquadratic;
}

void
Tria2PlateSubSoil ::computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the [1x6] displacement interpolation matrix {N}
{
    FloatArray N(6);
    giveInterpolation()->evalN(N, iLocCoord, FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(N, 1);
}

void
Tria2PlateSubSoil :: computeSurfaceNMatrix(FloatMatrix &answer, int boundaryID, const FloatArray &lcoords)
{
    if (boundaryID == 1) {
        this->computeNmatrixAt(lcoords, answer);
    } else {
        OOFEM_ERROR("computeSurfaceNMatrix: Only one surface is supported with id=1");
    }
}

} // end namespace oofem
