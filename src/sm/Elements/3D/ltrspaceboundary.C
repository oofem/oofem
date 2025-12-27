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

#include "sm/Elements/3D/ltrspaceboundary.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "fei3dtetlin.h"
#include "classfactory.h"
#include "Materials/structuralms.h"
#include "parametermanager.h"
#include "paramkey.h"


namespace oofem {
REGISTER_Element(LTRSpaceBoundary);

ParamKey LTRSpaceBoundary::IPK_LTRSpaceBoundary_location("location");

FEI3dTetLin LTRSpaceBoundary :: interpolation;

LTRSpaceBoundary :: LTRSpaceBoundary(int n, Domain *aDomain) :
    Structural3DElement(n, aDomain), NodalAveragingRecoveryModelInterface(), SpatialLocalizerInterface(this)
{
    numberOfDofMans  = 5;
    numberOfGaussPoints = 1;
}

void
LTRSpaceBoundary :: initializeFrom(InputRecord &ir, int priority)
{
    Structural3DElement :: initializeFrom(ir, priority);
    ParameterManager &ppm = this->giveDomain()->elementPPM;
    PM_UPDATE_PARAMETER(location, ppm, ir, this->number, IPK_LTRSpaceBoundary_location, priority) ;
}

void
LTRSpaceBoundary :: postInitialize()
{
    Structural3DElement :: postInitialize();
    ParameterManager &ppm = domain->elementPPM;
    PM_ELEMENT_ERROR_IFNOTSET(ppm, this->number, IPK_LTRSpaceBoundary_location) ;
}

Interface *
LTRSpaceBoundary :: giveInterface(InterfaceType interface)
{
    if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >( this );
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >( this );
    }

    return NULL;
}


FEInterpolation *
LTRSpaceBoundary :: giveInterpolation() const
{
    return & interpolation;
}

int
LTRSpaceBoundary :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    FEInterpolation *fei = this->giveInterpolation();
    FloatArray n;
    FEIElementGeometryWrapper cellgeo(this);

    fei->evalN(n, lcoords, cellgeo); //this interpolation doesn't use cell geometry to compute shape functions, so it's ok to have it here

    answer.clear();
    for ( int i = 1; i <= 4; i++ ) {
        if ( location.at(i) != 0 ) { //recalculate vertex coordinates
            IntArray switches;
            FloatArray vertexCoords;
            recalculateCoordinates(i, vertexCoords);

            answer.add(n.at(i), vertexCoords);
        } else {
            answer.add(n.at(i), cellgeo.giveVertexCoordinates(i) );
        }
    }

    return true;
}

void
LTRSpaceBoundary :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( inode == 5 ) {
        answer = { E_xx, E_xy, E_xz, E_yx, E_yy, E_yz, E_zx, E_zy, E_zz };
    } else {
        answer = { D_u, D_v, D_w };
    }
}

void
LTRSpaceBoundary :: giveSwitches(IntArray &answer, int location) {
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

double
LTRSpaceBoundary :: computeVolumeAround(GaussPoint *gp)
{
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, detJ, weight, volume;
    FloatArray v1, v2, v3, v4;

    recalculateCoordinates(1, v1);
    recalculateCoordinates(2, v2);
    recalculateCoordinates(3, v3);
    recalculateCoordinates(4, v4);

    x1 = v1.at(1);
    x2 = v2.at(1);
    x3 = v3.at(1);
    x4 = v4.at(1);
    y1 = v1.at(2);
    y2 = v2.at(2);
    y3 = v3.at(2);
    y4 = v4.at(2);
    z1 = v1.at(3);
    z2 = v2.at(3);
    z3 = v3.at(3);
    z4 = v4.at(3);

    detJ = ( ( x4 - x1 ) * ( y2 - y1 ) * ( z3 - z1 ) - ( x4 - x1 ) * ( y3 - y1 ) * ( z2 - z1 ) +
             ( x3 - x1 ) * ( y4 - y1 ) * ( z2 - z1 ) - ( x2 - x1 ) * ( y4 - y1 ) * ( z3 - z1 ) +
             ( x2 - x1 ) * ( y3 - y1 ) * ( z4 - z1 ) - ( x3 - x1 ) * ( y2 - y1 ) * ( z4 - z1 ) );

    if ( detJ <= 0.0 ) {
        OOFEM_ERROR("negative volume");
    }

    weight = gp->giveWeight();
    volume = detJ * weight;

    return volume;
}

void
LTRSpaceBoundary :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, detJ;
    FloatArray v1, v2, v3, v4;
    FloatMatrix dNdx;

    recalculateCoordinates(1, v1);
    recalculateCoordinates(2, v2);
    recalculateCoordinates(3, v3);
    recalculateCoordinates(4, v4);

    x1 = v1.at(1);
    x2 = v2.at(1);
    x3 = v3.at(1);
    x4 = v4.at(1);
    y1 = v1.at(2);
    y2 = v2.at(2);
    y3 = v3.at(2);
    y4 = v4.at(2);
    z1 = v1.at(3);
    z2 = v2.at(3);
    z3 = v3.at(3);
    z4 = v4.at(3);

    detJ = ( ( x4 - x1 ) * ( y2 - y1 ) * ( z3 - z1 ) - ( x4 - x1 ) * ( y3 - y1 ) * ( z2 - z1 ) +
             ( x3 - x1 ) * ( y4 - y1 ) * ( z2 - z1 ) - ( x2 - x1 ) * ( y4 - y1 ) * ( z3 - z1 ) +
             ( x2 - x1 ) * ( y3 - y1 ) * ( z4 - z1 ) - ( x3 - x1 ) * ( y2 - y1 ) * ( z4 - z1 ) );

    if ( detJ <= 0.0 ) {
        OOFEM_ERROR("negative volume");
    }

    dNdx.resize(4, 3);
    dNdx.at(1, 1) = -( ( y3 - y2 ) * ( z4 - z2 ) - ( y4 - y2 ) * ( z3 - z2 ) );
    dNdx.at(2, 1) = ( y4 - y3 ) * ( z1 - z3 ) - ( y1 - y3 ) * ( z4 - z3 );
    dNdx.at(3, 1) = -( ( y1 - y4 ) * ( z2 - z4 ) - ( y2 - y4 ) * ( z1 - z4 ) );
    dNdx.at(4, 1) = ( y2 - y1 ) * ( z3 - z1 ) - ( y3 - y1 ) * ( z2 - z1 );

    dNdx.at(1, 2) = -( ( x4 - x2 ) * ( z3 - z2 ) - ( x3 - x2 ) * ( z4 - z2 ) );
    dNdx.at(2, 2) = ( x1 - x3 ) * ( z4 - z3 ) - ( x4 - x3 ) * ( z1 - z3 );
    dNdx.at(3, 2) = -( ( x2 - x4 ) * ( z1 - z4 ) - ( x1 - x4 ) * ( z2 - z4 ) );
    dNdx.at(4, 2) = ( x3 - x1 ) * ( z2 - z1 ) - ( x2 - x1 ) * ( z3 - z1 );

    dNdx.at(1, 3) = -( ( x3 - x2 ) * ( y4 - y2 ) - ( x4 - x2 ) * ( y3 - y2 ) );
    dNdx.at(2, 3) = ( x4 - x3 ) * ( y1 - y3 ) - ( x1 - x3 ) * ( y4 - y3 );
    dNdx.at(3, 3) = -( ( x1 - x4 ) * ( y2 - y4 ) - ( x2 - x4 ) * ( y1 - y4 ) );
    dNdx.at(4, 3) = ( x2 - x1 ) * ( y3 - y1 ) - ( x3 - x1 ) * ( y2 - y1 );
    dNdx.times(1. / detJ);

    answer.resize(6, dNdx.giveNumberOfRows() * 3);
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, 3 * i - 2) = dNdx.at(i, 1);
        answer.at(2, 3 * i - 1) = dNdx.at(i, 2);
        answer.at(3, 3 * i - 0) = dNdx.at(i, 3);

        answer.at(5, 3 * i - 2) = answer.at(4, 3 * i - 1) = dNdx.at(i, 3);
        answer.at(6, 3 * i - 2) = answer.at(4, 3 * i - 0) = dNdx.at(i, 2);
        answer.at(6, 3 * i - 1) = answer.at(5, 3 * i - 0) = dNdx.at(i, 1);
    }
}

void
LTRSpaceBoundary :: computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray u;
    this->computeVectorOf({ D_u, D_v, D_w }, VM_Total, tStep, u); // solution vector
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }
    u.resizeWithValues(12);

    // Displacement gradient H = du/dX
    FloatMatrix B;
    this->computeBHmatrixAt(gp, B);
    answer.beProductOf(B, u);

    // Deformation gradient F = H + I
    answer.at(1) += 1.0;
    answer.at(2) += 1.0;
    answer.at(3) += 1.0;
}

void
LTRSpaceBoundary :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the [ 9 x (nno * 3) ] displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// BH matrix  -  9 rows : du/dx, dv/dy, dw/dz, dv/dz, du/dz, du/dy, dw/dy, dw/dx, dv/dx
{
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, detJ;
    FloatArray v1, v2, v3, v4;
    FloatMatrix dNdx;

    recalculateCoordinates(1, v1);
    recalculateCoordinates(2, v2);
    recalculateCoordinates(3, v3);
    recalculateCoordinates(4, v4);

    x1 = v1.at(1);
    x2 = v2.at(1);
    x3 = v3.at(1);
    x4 = v4.at(1);
    y1 = v1.at(2);
    y2 = v2.at(2);
    y3 = v3.at(2);
    y4 = v4.at(2);
    z1 = v1.at(3);
    z2 = v2.at(3);
    z3 = v3.at(3);
    z4 = v4.at(3);

    detJ = ( ( x4 - x1 ) * ( y2 - y1 ) * ( z3 - z1 ) - ( x4 - x1 ) * ( y3 - y1 ) * ( z2 - z1 ) +
             ( x3 - x1 ) * ( y4 - y1 ) * ( z2 - z1 ) - ( x2 - x1 ) * ( y4 - y1 ) * ( z3 - z1 ) +
             ( x2 - x1 ) * ( y3 - y1 ) * ( z4 - z1 ) - ( x3 - x1 ) * ( y2 - y1 ) * ( z4 - z1 ) );

    if ( detJ <= 0.0 ) {
        OOFEM_ERROR("negative volume");
    }

    dNdx.resize(4, 3);
    dNdx.at(1, 1) = -( ( y3 - y2 ) * ( z4 - z2 ) - ( y4 - y2 ) * ( z3 - z2 ) );
    dNdx.at(2, 1) = ( y4 - y3 ) * ( z1 - z3 ) - ( y1 - y3 ) * ( z4 - z3 );
    dNdx.at(3, 1) = -( ( y1 - y4 ) * ( z2 - z4 ) - ( y2 - y4 ) * ( z1 - z4 ) );
    dNdx.at(4, 1) = ( y2 - y1 ) * ( z3 - z1 ) - ( y3 - y1 ) * ( z2 - z1 );

    dNdx.at(1, 2) = -( ( x4 - x2 ) * ( z3 - z2 ) - ( x3 - x2 ) * ( z4 - z2 ) );
    dNdx.at(2, 2) = ( x1 - x3 ) * ( z4 - z3 ) - ( x4 - x3 ) * ( z1 - z3 );
    dNdx.at(3, 2) = -( ( x2 - x4 ) * ( z1 - z4 ) - ( x1 - x4 ) * ( z2 - z4 ) );
    dNdx.at(4, 2) = ( x3 - x1 ) * ( z2 - z1 ) - ( x2 - x1 ) * ( z3 - z1 );

    dNdx.at(1, 3) = -( ( x3 - x2 ) * ( y4 - y2 ) - ( x4 - x2 ) * ( y3 - y2 ) );
    dNdx.at(2, 3) = ( x4 - x3 ) * ( y1 - y3 ) - ( x1 - x3 ) * ( y4 - y3 );
    dNdx.at(3, 3) = -( ( x1 - x4 ) * ( y2 - y4 ) - ( x2 - x4 ) * ( y1 - y4 ) );
    dNdx.at(4, 3) = ( x2 - x1 ) * ( y3 - y1 ) - ( x3 - x1 ) * ( y2 - y1 );
    dNdx.times(1. / detJ);

    answer.resize(9, dNdx.giveNumberOfRows() * 3);
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, 3 * i - 2) = dNdx.at(i, 1);     // du/dx
        answer.at(2, 3 * i - 1) = dNdx.at(i, 2);     // dv/dy
        answer.at(3, 3 * i - 0) = dNdx.at(i, 3);     // dw/dz
        answer.at(4, 3 * i - 1) = dNdx.at(i, 3);     // dv/dz
        answer.at(7, 3 * i - 0) = dNdx.at(i, 2);     // dw/dy
        answer.at(5, 3 * i - 2) = dNdx.at(i, 3);     // du/dz
        answer.at(8, 3 * i - 0) = dNdx.at(i, 1);     // dw/dx
        answer.at(6, 3 * i - 2) = dNdx.at(i, 2);     // du/dy
        answer.at(9, 3 * i - 1) = dNdx.at(i, 1);     // dv/dx
    }
}

void
LTRSpaceBoundary :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray u, dispVec;
    FloatMatrix B, T;
    this->computeVectorOf(VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    this->computeTransformationMatrix(T, tStep);
    dispVec.beProductOf(T, u);

    this->computeBmatrixAt(gp, B);
    answer.beProductOf(B, dispVec);
}

void
LTRSpaceBoundary :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatMatrix B, T, Tt;
    FloatArray vStress, vStrain, fintsub;

    // zero answer will resize accordingly when adding first contribution
    answer.clear();
    fintsub.resize(12);

    for ( auto &gp : * this->giveDefaultIntegrationRulePtr() ) {
        StructuralMaterialStatus *matStat = static_cast< StructuralMaterialStatus * >( this->giveCrossSection()->giveMaterial(gp)->giveStatus(gp) );

        // Engineering (small strain) stress
        if ( nlGeometry == 0 ) {
            this->computeBmatrixAt(gp, B);
            if ( useUpdatedGpRecord == 1 ) {
                vStress = matStat->giveStressVector();
            } else {
                if ( !this->isActivated(tStep) ) {
                    vStrain.resize(StructuralMaterial :: giveSizeOfVoigtSymVector(gp->giveMaterialMode() ) );
                    vStrain.zero();
                }
                this->computeStrainVector(vStrain, gp, tStep);
                this->computeStressVector(vStress, vStrain, gp, tStep);
            }
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

    this->computeTransformationMatrix(T, tStep);
    Tt.beTranspositionOf(T);

    answer.beProductOf(Tt, fintsub);
}

void
LTRSpaceBoundary :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    FloatMatrix Korig, T, Tt, TtK;
    NLStructuralElement :: computeStiffnessMatrix(Korig, rMode, tStep);

    this->computeTransformationMatrix(T, tStep);
    Tt.beTranspositionOf(T);

    TtK.beProductOf(Tt, Korig);
    answer.beProductOf(TtK, T);
}

void
LTRSpaceBoundary :: computeTransformationMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    FloatArray unitCellSize;
    unitCellSize.resize(3);
    unitCellSize.at(1) = this->giveNode(5)->giveCoordinate(1);
    unitCellSize.at(2) = this->giveNode(5)->giveCoordinate(2);
    unitCellSize.at(3) = this->giveNode(5)->giveCoordinate(3);

    IntArray switches1, switches2, switches3, switches4;
    this->giveSwitches(switches1, this->location.at(1) );
    this->giveSwitches(switches2, this->location.at(2) );
    this->giveSwitches(switches3, this->location.at(3) );
    this->giveSwitches(switches4, this->location.at(4) );

    FloatMatrix k1, k2, k3, k4;
    k1.resize(3, 9);
    k2.resize(3, 9);
    k3.resize(3, 9);
    k4.resize(3, 9);

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

    for ( int i = 1; i <= 3; i++ ) {
        k3.at(i, 3 * i - 2) = unitCellSize.at(1) * switches3.at(1);
        k3.at(i, 3 * i - 1) = unitCellSize.at(2) * switches3.at(2);
        k3.at(i, 3 * i) = unitCellSize.at(3) * switches3.at(3);
    }

    for ( int i = 1; i <= 3; i++ ) {
        k4.at(i, 3 * i - 2) = unitCellSize.at(1) * switches4.at(1);
        k4.at(i, 3 * i - 1) = unitCellSize.at(2) * switches4.at(2);
        k4.at(i, 3 * i) = unitCellSize.at(3) * switches4.at(3);
    }

    answer.resize(12, 12);
    answer.beUnitMatrix();
    answer.resizeWithData(12, 21);

    answer.assemble(k1, { 1, 2, 3 }, { 13, 14, 15, 16, 17, 18, 19, 20, 21 });
    answer.assemble(k2, { 4, 5, 6 }, { 13, 14, 15, 16, 17, 18, 19, 20, 21 });
    answer.assemble(k3, { 7, 8, 9 }, { 13, 14, 15, 16, 17, 18, 19, 20, 21 });
    answer.assemble(k4, { 10, 11, 12 }, { 13, 14, 15, 16, 17, 18, 19, 20, 21 });
}

void
LTRSpaceBoundary :: recalculateCoordinates(int nodeNumber, FloatArray &coords)
{
    FloatArray unitCellSize;
    unitCellSize.resize(3);
    unitCellSize.at(1) = this->giveNode(5)->giveCoordinate(1);
    unitCellSize.at(2) = this->giveNode(5)->giveCoordinate(2);
    unitCellSize.at(3) = this->giveNode(5)->giveCoordinate(3);

    IntArray switches;
    this->giveSwitches(switches, this->location.at(nodeNumber) );

    coords.resize(3);
    coords.at(1) = this->giveNode(nodeNumber)->giveCoordinate(1) + switches.at(1) * unitCellSize.at(1);
    coords.at(2) = this->giveNode(nodeNumber)->giveCoordinate(2) + switches.at(2) * unitCellSize.at(2);
    coords.at(3) = this->giveNode(nodeNumber)->giveCoordinate(3) + switches.at(3) * unitCellSize.at(3);

    return;
}

int
LTRSpaceBoundary :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
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
    return Element :: giveIPValue(answer, gp, type, tStep);
}

double
LTRSpaceBoundary :: giveLengthInDir(const FloatArray &normalToCrackPlane)
{
    double maxDis, minDis;
    int nnode = giveNumberOfNodes() - 1; //don't take control node

    FloatArray coords(3);
    recalculateCoordinates(1, coords);
    minDis = maxDis = normalToCrackPlane.dotProduct(coords, coords.giveSize() );

    for ( int i = 2; i <= nnode; i++ ) {
        FloatArray coords(3);
        recalculateCoordinates(i, coords);
        double dis = normalToCrackPlane.dotProduct(coords, coords.giveSize() );
        if ( dis > maxDis ) {
            maxDis = dis;
        } else if ( dis < minDis ) {
            minDis = dis;
        }
    }

    return maxDis - minDis;
}

void
LTRSpaceBoundary :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                               InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp;

    if ( numberOfGaussPoints != 1 ) {
        answer.clear(); // for more gp's need to be refined
        return;
    }

    gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    giveIPValue(answer, gp, type, tStep);
}
} // end namespace oofem
