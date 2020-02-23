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

#include "sm/Elements/quad1platesubsoil.h"
#include "sm/Materials/structuralms.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "fei2dquadlin.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "load.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(Quad1PlateSubSoil);

FEI2dQuadLin Quad1PlateSubSoil :: interp_lin(1, 2);

Quad1PlateSubSoil :: Quad1PlateSubSoil(int n, Domain *aDomain) :
    StructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(this),
    SPRNodalRecoveryModelInterface()
{
    numberOfGaussPoints = 4;
    numberOfDofMans = 4;
}


FEInterpolation *
Quad1PlateSubSoil :: giveInterpolation() const { return & interp_lin; }


FEInterpolation *
Quad1PlateSubSoil :: giveInterpolation(DofIDItem id) const
{
    return & interp_lin;
}


void
Quad1PlateSubSoil :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this, 1, 5);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


void
Quad1PlateSubSoil :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
{
    OOFEM_ERROR("Body load not supported, use surface load instead");
}


void
Quad1PlateSubSoil :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [3x4] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    FloatArray n;
    FloatMatrix dn;

    this->interp_lin.evaldNdx( dn, gp->giveNaturalCoordinates(),  FEIElementGeometryWrapper(this) );
    this->interp_lin.evalN( n, gp->giveNaturalCoordinates(),  FEIElementGeometryWrapper(this) );

    answer.resize(3, 4);
    answer.zero();

    ///@todo Check sign here
    for ( int i = 0; i < 4; ++i ) {
        answer(0, i) = n(i); // eps_z
        answer(1, i) = dn(i, 0); // gamma_xz
        answer(2, i) = dn(i, 1); // gamma_yz
    }
}


void
Quad1PlateSubSoil :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    answer = this->giveStructuralCrossSection()->giveGeneralizedStress_PlateSubSoil(strain, gp, tStep);
}


void
Quad1PlateSubSoil :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer = this->giveStructuralCrossSection()->give2dPlateSubSoilStiffMtrx(rMode, gp, tStep);
}


void
Quad1PlateSubSoil :: initializeFrom(InputRecord &ir)
{
    this->numberOfGaussPoints = 4;
    StructuralElement :: initializeFrom(ir);
}


void
Quad1PlateSubSoil :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_w};
}


void
Quad1PlateSubSoil :: computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp)
{
    FloatArray u, v;
    u.beDifferenceOf( this->giveNode(2)->giveCoordinates(), this->giveNode(1)->giveCoordinates() );
    v.beDifferenceOf( this->giveNode(3)->giveCoordinates(), this->giveNode(1)->giveCoordinates() );

    answer.beVectorProductOf(u, v);
    answer.normalize();
}


double
Quad1PlateSubSoil :: giveCharacteristicLength(const FloatArray &normalToCrackPlane)
//
// returns receiver's characteristic length for crack band models
// for a crack formed in the plane with normal normalToCrackPlane.
//
{
    return this->giveCharacteristicLengthForPlaneElements(normalToCrackPlane);
}


double
Quad1PlateSubSoil :: computeVolumeAround(GaussPoint *gp)
{
    double weight = gp->giveWeight();
    double detJ = fabs( this->interp_lin.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    return detJ * weight;
}


void
Quad1PlateSubSoil :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver.
{
    OOFEM_ERROR("Mass matrix not provided");
}


int
Quad1PlateSubSoil :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_ShellForceTensor ){
        FloatArray help;
        answer.resize(6);
        help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStressVector();
        answer.at(1) = 0.0; // nx
        answer.at(2) = 0.0; // ny
        answer.at(3) = help.at(1); // nz
        answer.at(4) = help.at(3); // vyz
        answer.at(5) = help.at(2); // vxz
        answer.at(6) = 0.0; // vxy
        return 1;
    }
    return StructuralElement :: giveIPValue(answer, gp, type, tStep);
}

Interface *
Quad1PlateSubSoil :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    }

    return NULL;
}

void
Quad1PlateSubSoil :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(4);
    for ( int i = 1; i < 5; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
Quad1PlateSubSoil :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i < 5; i++ ) {
        if ( pap == this->giveNode(i)->giveNumber() ) {
            found = 1;
        }
    }

    if ( found ) {
        answer.at(1) = pap;
    } else {
        OOFEM_ERROR("node unknown");
    }
}

void
Quad1PlateSubSoil ::computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the [1x4] displacement interpolation matrix {N}
{
    FloatArray N(4);
    giveInterpolation()->evalN(N, iLocCoord, FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(N, 1);
}


void
Quad1PlateSubSoil ::computeSurfaceNMatrix(FloatMatrix &answer, int boundaryID, const FloatArray &lcoords)
{
    if (boundaryID == 1) {
        this->computeNmatrixAt(lcoords, answer);
    } else {
        OOFEM_ERROR("computeSurfaceNMatrix: Only one surface is supported with id=1");
    }
}



// void
// Quad1PlateSubSoil :: giveSurfaceDofMapping(IntArray &answer, int iSurf) const
// {
//     answer.resize(4);
//     answer.zero();
//     if ( iSurf == 1 ) {
//         for (int i = 1; i <= 4; i++) {
//             answer.at(i) = i;
//         }
//     } else {
//         OOFEM_ERROR("wrong surface number");
//     }
// }


// double
// Quad1PlateSubSoil :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
// {
//     return this->computeVolumeAround(gp);
// }


// int
// Quad1PlateSubSoil :: computeLoadLSToLRotationMatrix(FloatMatrix &answer, int isurf, GaussPoint *gp)
// {
//     return 0;
// }


} // end namespace oofem
