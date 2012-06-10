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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "qspace.h"
#include "node.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"

namespace oofem {
FEI3dHexaQuad QSpace :: interpolation;

QSpace :: QSpace(int n, Domain *aDomain) : StructuralElement(n, aDomain)
    // Constructor.
{
    numberOfDofMans = 20;
    //nn = numberOfDofMans; // number of nodes
    //nnsurf = 8;           // number of nodes on surface
    //ndofsn = 3;           // number of DOFs on node
}


IRResultType
QSpace :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    this->StructuralElement :: initializeFrom(ir);
    numberOfGaussPoints = 27;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_QSpace_nip, "nip"); // Macro

    if ( ( numberOfGaussPoints != 8 ) && ( numberOfGaussPoints != 14 ) && ( numberOfGaussPoints != 27 ) && ( numberOfGaussPoints != 64 ) ) {
        numberOfGaussPoints = 27;
    }

    // set - up Gaussian integration points
    this->computeGaussPoints();

    return IRRT_OK;
}


void
QSpace :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    answer.setValues(3, D_u, D_v, D_w);
}


double
QSpace :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    double determinant = this->interpolation.giveTransformationJacobian(* aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));
    double weight      = aGaussPoint->giveWeight();

    return ( determinant * weight );
}


int
QSpace :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    this->interpolation.local2global(answer, lcoords, FEIElementGeometryWrapper(this));
    return 1;
}

double
QSpace :: giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane)
{
    double factor = __OOFEM_POW( ( double ) this->numberOfGaussPoints, 1. / 3. );
    return this->giveLenghtInDir(normalToCrackPlane) / factor;
}


void
QSpace :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    numberOfIntegrationRules = 1;
    integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];
    integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 6);
    integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Cube, numberOfGaussPoints, _3dMat);
}


void
QSpace :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
    int i;
    FloatArray n(20);

    answer.resize(3, 60);
    answer.zero();

    this->interpolation.evalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    for ( i = 1; i <= 20; i++ ) {
        answer.at(1, 3 * i - 2) = n.at(i);
        answer.at(2, 3 * i - 1) = n.at(i);
        answer.at(3, 3 * i - 0) = n.at(i);
    }
}


void
QSpace :: computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint, int i) {
    OOFEM_ERROR("NLBmatrix not implemented");
}


void
QSpace :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the [6x60] strain-displacement matrix {B} of the receiver, eva-
// luated at aGaussPoint.
// B matrix  -  6 rows : epsilon-X, epsilon-Y, epsilon-Z, gamma-YZ, gamma-ZX, gamma-XY  :
{
    int i;
    FloatMatrix dnx;

    this->interpolation.evaldNdx(dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(6, 60);
    answer.zero();

    for ( i = 1; i <= 20; i++ ) {
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
QSpace :: computeBFmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer) {
    OOFEM_ERROR("BFmatrix not implemented");
}


// ******************************
// ***  Surface load support  ***
// ******************************

IntegrationRule *
QSpace :: GetSurfaceIntegrationRule(int approxOrder)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, this, 1, 1);
    int npoints = iRule->getRequiredNumberOfIntegrationPoints(_Square, approxOrder);
    iRule->setUpIntegrationPoints(_Square, npoints, _Unknown);
    return iRule;
}

void
QSpace :: computeSurfaceNMatrixAt(FloatMatrix &answer, GaussPoint *sgp)
{
    FloatArray n(8);
    interpolation.surfaceEvalN(n, * sgp->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(3, 24);
    answer.zero();

    for ( int i = 0; i < 8; i++ ) { // loop over surfaces
        for ( int j = 0; j < 3; j++ ) { // loop over DOFs
            answer(j, j + i * 3) = n(i);
        }
    }
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
    determinant = fabs( interpolation.surfaceGiveTransformationJacobian(iSurf, * gp->giveCoordinates(), FEIElementGeometryWrapper(this)) );

    weight      = gp->giveWeight();
    volume      = determinant * weight;

    return volume;
}

void
QSpace :: computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iSurf)
{
    interpolation.surfaceLocal2global(answer, iSurf, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));
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
     * _error ("computeLoadLSToLRotationMatrix: surface local coordinate system not supported");
     * return 1;
     */
    int i, j;
    FloatArray gc(3);
    FloatArray h1(3), h2(3), nn(3), n(3);
    IntArray snodes(4);

    answer.resize(3, 3);

    this->interpolation.computeSurfaceMapping(snodes, dofManArray, iSurf);
    for ( i = 1; i <= 4; i++ ) {
        gc.add( * domain->giveNode( snodes.at(i) )->giveCoordinates() );
    }

    gc.times(1. / 4.);
    // determine "average normal"
    for ( i = 1; i <= 4; i++ ) {
        j = ( i ) % 4 + 1;
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
        return 1;
    }

    nn.normalize();
    for ( i = 1; i <= 3; i++ ) {
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
    for ( i = 1; i <= 3; i++ ) {
        answer.at(i, 2) = h2.at(i);
    }

    return 1;
}

Interface *
QSpace :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return ( ZZNodalRecoveryModelInterface * ) this;
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return ( SPRNodalRecoveryModelInterface * ) this;
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return ( NodalAveragingRecoveryModelInterface * ) this;
    }

    OOFEM_LOG_INFO("Interface on Qspace element not supported");
    return NULL;
}

int
QSpace :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) || ( type == IST_DamageTensor ) ) {
        return 6;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}

void
QSpace :: ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint, InternalStateType type)
{
    int i;
    FloatArray n;
    this->interpolation.evalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    if ( this->giveIPValueSize(type, aGaussPoint) ) {
        answer.resize(1, 20);
    } else {
        return;
    }

    for ( i = 1; i <= 20; i++ ) {
        answer.at(1, i)  = n.at(i);
    }
}

int
QSpace :: SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    return ZZNodalRecoveryMI_giveDofManRecordSize(type);
}

void
QSpace :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    int i;

    pap.resize(20);
    for ( i = 1; i <= 20; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
QSpace :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int i, found = 0;
    answer.resize(1);

    for ( i = 1; i <= 20; i++ ) {
        if ( this->giveNode(i)->giveNumber() == pap ) {
            found = 1;
        }
    }

    if ( found ) {
        answer.at(1) = pap;
    } else {
        _error2("SPRNodalRecoveryMI_giveDofMansDeterminedByPatch: unknown node number %d", pap);
    }
}

int
QSpace :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return numberOfGaussPoints;
}


void
QSpace :: SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp)
{
    if ( this->computeGlobalCoordinates( coords, * gp->giveCoordinates() ) == 0 ) {
        _error("SPRNodalRecoveryMI_computeIPGlobalCoordinates: computeGlobalCoordinates failed");
    }
}

SPRPatchType
QSpace :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_3dBiQuadratic;
}



void
QSpace :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    int size = NodalAveragingRecoveryMI_giveDofManRecordSize(type);
    answer.resize(size);
    answer.zero();
    _warning("Qspace element: IP values will not be transferred to nodes. Use ZZNodalRecovery instead (parameter stype 1)");
}

void
QSpace :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}
} // end namespace oofem
