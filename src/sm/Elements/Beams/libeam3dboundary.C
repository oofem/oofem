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
 *               Copyright (C) 1993 - 2019   Borek Patzak
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

#include "sm/Elements/Beams/libeam3dboundary.h"
#include "sm/Materials/structuralms.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(LIBeam3dBoundary);

LIBeam3dBoundary :: LIBeam3dBoundary(int n, Domain *aDomain) : LIBeam3d(n, aDomain)
    // Constructor.
{
    numberOfDofMans     = 3;
    referenceNode       = 0;
    length              = 0.;
}


void
LIBeam3dBoundary :: initializeFrom(InputRecord &ir)
{
    LIBeam3d :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, referenceNode, _IFT_LIBeam3dBoundary_refnode);
    if ( referenceNode == 0 ) {
        OOFEM_ERROR("wrong reference node specified");
    }

    IR_GIVE_FIELD(ir, location, _IFT_LIBeam3dBoundary_location);
}


void
LIBeam3dBoundary :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( inode == 3 ) {
        answer = { E_xx, E_xy, E_xz, E_yx, E_yy, E_yz, E_zx, E_zy, E_zz };
    } else {
        answer = { D_u, D_v, D_w, R_u, R_v, R_w };
    }
}


void
LIBeam3dBoundary :: giveSwitches(IntArray &answer, int location) {
    answer.resize(3);
    int counter = 1;
    for ( int x = -1; x <  2; x++ ) {
        for ( int y = -1; y <  2; y++ ) {
            for ( int z = -1; z <  2; z++ ) {
                if ( !( z == 0 && y == 0 && x == 0 ) ) {
                    if ( counter == location ) {
                        answer(0) = x;
                        answer(1) = y;
                        answer(2) = z;
                    }
                    counter++;
                }
            }
        }
    }
    return;
}


void
LIBeam3dBoundary :: recalculateCoordinates(int nodeNumber, FloatArray &coords)
{
    FloatArray unitCellSize;
    unitCellSize.resize(3);
    unitCellSize.at(1) = this->giveNode(3)->giveCoordinate(1);
    unitCellSize.at(2) = this->giveNode(3)->giveCoordinate(2);
    unitCellSize.at(3) = this->giveNode(3)->giveCoordinate(3);

    IntArray switches;
    this->giveSwitches(switches, this->location.at(nodeNumber) );

    coords.resize(3);
    coords.at(1) = this->giveNode(nodeNumber)->giveCoordinate(1) + switches.at(1) * unitCellSize.at(1);
    coords.at(2) = this->giveNode(nodeNumber)->giveCoordinate(2) + switches.at(2) * unitCellSize.at(2);
    coords.at(3) = this->giveNode(nodeNumber)->giveCoordinate(3) + switches.at(3) * unitCellSize.at(3);

    return;
}


double
LIBeam3dBoundary :: computeLength()
// Returns the length of the receiver.
{
    double dx, dy, dz;
    FloatArray nodeA, nodeB;

    if ( length == 0. ) {
        recalculateCoordinates(1, nodeA);
        recalculateCoordinates(2, nodeB);
        dx      = nodeB.at(1) - nodeA.at(1);
        dy      = nodeB.at(2) - nodeA.at(2);
        dz      = nodeB.at(3) - nodeA.at(3);
        length  = sqrt(dx * dx + dy * dy + dz * dz);
    }

    return length;
}


int
LIBeam3dBoundary :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    double ksi, n1, n2;
    FloatArray coordsNodeA, coordsNodeB;

    ksi = lcoords.at(1);
    n1 = ( 1. - ksi ) * 0.5;
    n2 = ( 1. + ksi ) * 0.5;
    recalculateCoordinates(1, coordsNodeA);
    recalculateCoordinates(2, coordsNodeB);

    answer.resize(3);
    answer.at(1) = n1 * coordsNodeA.at(1) + n2 * coordsNodeB.at(1);
    answer.at(2) = n1 * coordsNodeA.at(2) + n2 * coordsNodeB.at(2);
    answer.at(3) = n1 * coordsNodeA.at(3) + n2 * coordsNodeB.at(3);

    return 1;
}


int
LIBeam3dBoundary :: giveLocalCoordinateSystem(FloatMatrix &answer)
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
    FloatArray lx(3), ly(3), lz(3), help(3);
    double length = this->computeLength();
    Node *refNode;
    FloatArray coordsNodeA, coordsNodeB;

    answer.resize(3, 3);
    answer.zero();
    refNode = this->giveReferenceNode(referenceNode);
    recalculateCoordinates(1, coordsNodeA);
    recalculateCoordinates(2, coordsNodeB);

    for ( int i = 1; i <= 3; i++ ) {
        lx.at(i) = ( coordsNodeB.at(i) - coordsNodeA.at(i) ) / length;
        help.at(i) = ( refNode->giveCoordinate(i) - coordsNodeA.at(i) );
    }

    lz.beVectorProductOf(lx, help);
    lz.normalize();
    ly.beVectorProductOf(lz, lx);
    ly.normalize();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, i) = lx.at(i);
        answer.at(2, i) = ly.at(i);
        answer.at(3, i) = lz.at(i);
    }

    return 1;
}


bool
LIBeam3dBoundary :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    FloatMatrix lcs;

    answer.resize(this->computeNumberOfDofs(), this->computeNumberOfDofs() );
    answer.zero();

    this->giveLocalCoordinateSystem(lcs);
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(i, j) = lcs.at(i, j);
            answer.at(i + 3, j + 3) = lcs.at(i, j);
            answer.at(i + 6, j + 6) = lcs.at(i, j);
            answer.at(i + 9, j + 9) = lcs.at(i, j);
        }
    }
    //do not rotate extra DOFs
    for ( int i = 13; i <= this->computeNumberOfDofs(); i++ ) {
        for ( int j = 13; j <= this->computeNumberOfDofs(); j++ ) {
            if ( i == j ) {
                answer.at(i, j) = 1.;
            }
        }
    }
    return true;
}


void
LIBeam3dBoundary :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
// Returns the stiffness matrix of the receiver, expressed in the global
// axes.
{
    //get the raw stiffness matrix in local cs
    FloatMatrix Korig, GtoL, R, Rt, T, Tt, TtK;
    StructuralElement :: computeStiffnessMatrix(Korig, rMode, tStep);

    //rotate to global cs
    this->giveRotationMatrix(GtoL);
    R.beSubMatrixOf(GtoL, 1, 12, 1, 12);
    Korig.rotatedWith(R);

    //transform
    this->computeTransformationMatrix(T, tStep);
    Tt.beTranspositionOf(T);
    TtK.beProductOf(Tt, Korig);
    answer.beProductOf(TtK, T);

    //rotate back to local cs
    Rt.beTranspositionOf(GtoL);
    answer.rotatedWith(Rt);
}


void
LIBeam3dBoundary :: computeTransformationMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    FloatArray unitCellSize;
    unitCellSize.resize(3);
    unitCellSize.at(1) = this->giveNode(3)->giveCoordinate(1);
    unitCellSize.at(2) = this->giveNode(3)->giveCoordinate(2);
    unitCellSize.at(3) = this->giveNode(3)->giveCoordinate(3);

    IntArray switches1, switches2;
    this->giveSwitches(switches1, this->location.at(1) );
    this->giveSwitches(switches2, this->location.at(2) );

    FloatMatrix k1, k2;
    k1.resize(6, 9);
    k2.resize(6, 9);

    for ( int i = 1; i <= 3; i++ ) {
        k1.at(i, 3 * i - 2) = unitCellSize.at(1) * switches1.at(1);
        k1.at(i, 3 * i - 1) = unitCellSize.at(2) * switches1.at(2);
        k1.at(i, 3 * i) = unitCellSize.at(3) * switches1.at(3);
    }

    for ( int i = 1; i <= 3; i++ ) {
        k2.at(i, 3 * i - 2) = unitCellSize.at(1) * switches2.at(1);
        k2.at(i, 3 * i - 1) = unitCellSize.at(2) * switches2.at(2);
        k2.at(i, 3 * i) = unitCellSize.at(3) * switches2.at(3);
    }

    answer.resize(12, 12);
    answer.beUnitMatrix();
    answer.resizeWithData(12, 21);

    answer.assemble(k1, { 1, 2, 3, 4, 5, 6 }, { 13, 14, 15, 16, 17, 18, 19, 20, 21 });
    answer.assemble(k2, { 7, 8, 9, 10, 11, 12 }, { 13, 14, 15, 16, 17, 18, 19, 20, 21 });
}


void
LIBeam3dBoundary :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray u, uG, dispVecG, dispVecL;
    FloatMatrix B, T, GtoL, GtoLT, R;
    this->computeVectorOf(VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    //rotate to global
    this->giveRotationMatrix(GtoL);
    GtoLT.beTranspositionOf(GtoL);
    R.beSubMatrixOf(GtoL, 1, 12, 1, 12);
    uG.beProductOf(GtoLT, u);
    //transform
    this->computeTransformationMatrix(T, tStep);
    dispVecG.beProductOf(T, uG);
    //rotate to local
    dispVecL.beProductOf(R, dispVecG);

    this->computeBmatrixAt(gp, B);
    answer.beProductOf(B, dispVecL);
}


void
LIBeam3dBoundary :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatMatrix B, T, Tt, GtoL, R, Rt;
    FloatArray vStress, vStrain, fintsub, fintG, fintGT;

    //get the raw internal force vector in local cs
    answer.clear();
    fintsub.resize(12);

    for ( auto &gp : * this->giveDefaultIntegrationRulePtr() ) {
        this->computeBmatrixAt(gp, B);

        if ( useUpdatedGpRecord == 1 ) {
            StructuralMaterialStatus *matStat = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );
            vStress = matStat->giveStressVector();
        } else {
            if ( !this->isActivated(tStep) ) {
                vStrain.resize(StructuralMaterial :: giveSizeOfVoigtSymVector(gp->giveMaterialMode() ) );
                vStrain.zero();
            }
            this->computeStrainVector(vStrain, gp, tStep);
            this->computeStressVector(vStress, vStrain, gp, tStep);
        }

        // Compute nodal internal forces at nodes as f = B^T*Stress dV
        double dV  = this->computeVolumeAround(gp);

        if ( vStress.giveSize() == 6 ) {
            FloatArray stressTemp;
            StructuralMaterial :: giveReducedSymVectorForm(stressTemp, vStress, gp->giveMaterialMode() );
            fintsub.plusProduct(B, stressTemp, dV);
        } else {
            fintsub.plusProduct(B, vStress, dV);
        }
    }

    //rotate to global cs
    this->giveRotationMatrix(GtoL);
    R.beSubMatrixOf(GtoL, 1, 12, 1, 12);
    Rt.beTranspositionOf(R);
    fintG.beProductOf(Rt, fintsub);
    //transform
    this->computeTransformationMatrix(T, tStep);
    Tt.beTranspositionOf(T);
    fintGT.beProductOf(Tt, fintG);
    //rotate back to local cs
    answer.beProductOf(GtoL, fintGT);
}


int
LIBeam3dBoundary :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_DisplacementVector ) {
        FloatArray u;
        FloatMatrix N;
        this->computeVectorOf(VM_Total, tStep, u);
        u.resizeWithValues(12);
        this->computeNmatrixAt(gp->giveSubPatchCoordinates(), N); //no special treatment needed for this interpolation
        answer.beProductOf(N, u);
        return 1;
    }
    return LIBeam3d :: giveIPValue(answer, gp, type, tStep);
}
} // end namespace oofem
