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


#include "qtrspace.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "cltypes.h"
#include "structuralms.h"
#include "mathfem.h"
#include "structuralcrosssection.h"
#include "fei3dtetquad.h"
#include "classfactory.h"

#include <cstdio>

namespace oofem {

REGISTER_Element( QTRSpace );
 
FEI3dTetQuad QTRSpace :: interpolation;

QTRSpace :: QTRSpace(int n, Domain *aDomain) : NLStructuralElement(n, aDomain)
{
    numberOfDofMans = 10;
}


IRResultType
QTRSpace :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    this->NLStructuralElement :: initializeFrom(ir);
    numberOfGaussPoints = 4;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, _IFT_Element_nip);

    if ( ( numberOfGaussPoints != 1 ) && ( numberOfGaussPoints != 4 )&& ( numberOfGaussPoints != 5 ) ) {
        numberOfGaussPoints = 4;
    }

    return IRRT_OK;
}


FEInterpolation*
QTRSpace :: giveInterpolation() const
{
    return &interpolation;
}


void
QTRSpace :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
// returns DofId mask array for inode element node.
// DofId mask array determines the dof ordering requsted from node.
// DofId mask array contains the DofID constants (defined in cltypes.h)
// describing physical meaning of particular DOFs.
{
    answer.resize(3);

    answer.at(1) = D_u;
    answer.at(2) = D_v;
    answer.at(3) = D_w;
}


double
QTRSpace :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    double determinant = fabs ((this->interpolation.giveTransformationJacobian(* aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this))));
    double weight = aGaussPoint->giveWeight();
    return ( determinant * weight);
}


void
QTRSpace :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    numberOfIntegrationRules = 1;
    integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];
    integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 6);
    this->giveCrossSection()->setupIntegrationPoints( *integrationRulesArray[0], numberOfGaussPoints, this );
}


MaterialMode
QTRSpace :: giveMaterialMode()
{
    return _3dMat;
}


void
QTRSpace :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
    FloatArray n(10);

    answer.resize(3, 30);
    answer.zero();

    this->interpolation.evalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    for ( int i = 1; i <= 10; i++ ) {
        answer.at(1, 3 * i - 2) = n.at(i);
        answer.at(2, 3 * i - 1) = n.at(i);
        answer.at(3, 3 * i - 0) = n.at(i);
    }
}

#if 1

#endif

void
QTRSpace :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the [6x30] strain-displacement matrix {B} of the receiver, eva-
// luated at aGaussPoint.
// B matrix  -  6 rows : epsilon-X, epsilon-Y, epsilon-Z, gamma-YZ, gamma-ZX, gamma-XY  :
{
    FloatMatrix dnx;

    this->interpolation.evaldNdx(dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(6, 30);
    answer.zero();

    for ( int i = 1; i <= 10; i++ ) {
        answer.at(1, 3 * i - 2) = dnx.at(1, i);
        answer.at(2, 3 * i - 1) = dnx.at(2, i);
        answer.at(3, 3 * i - 0) = dnx.at(3, i);

        answer.at(4, 3 * i - 1) = dnx.at(3, i);
        answer.at(4, 3 * i - 0) = dnx.at(2, i);

        answer.at(5, 3 * i - 2) = dnx.at(3, i);
        answer.at(5, 3 * i - 0) = dnx.at(1, i);

        answer.at(6, 3 * i - 2) = dnx.at(2, i);
        answer.at(6, 3 * i - 1) = dnx.at(1, i);

    }
}


void
QTRSpace :: computeBHmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
{
    FloatMatrix dnx;
  
    this->interpolation.evaldNdx(dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(9, 30);
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

Interface *
QTRSpace :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >( this );
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >( this );
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >( this );
    }

    OOFEM_LOG_INFO("Interface on QTRSpace element not supported");
    return NULL;
}

int
QTRSpace :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) || ( type == IST_DamageTensor ) ) {
        return 6;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}

int
QTRSpace :: SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    return ZZNodalRecoveryMI_giveDofManRecordSize(type);
}

void
QTRSpace :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(10);
    for ( int i = 1; i <= 10; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
QTRSpace :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i <= 10; i++ ) {
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
QTRSpace :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return numberOfGaussPoints;
}


void
QTRSpace :: SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp)
{
    if ( this->computeGlobalCoordinates( coords, * gp->giveCoordinates() ) == 0 ) {
        _error("SPRNodalRecoveryMI_computeIPGlobalCoordinates: computeGlobalCoordinates failed");
    }
}

SPRPatchType
QTRSpace :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_3dBiQuadratic;
}


void
QTRSpace :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    int size = NodalAveragingRecoveryMI_giveDofManRecordSize(type);
    answer.resize(size);
    answer.zero();
    _warning("QTRSpace element: IP values will not be transferred to nodes. Use ZZNodalRecovery instead (parameter stype 1)");
}

void
QTRSpace :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}

} // end namespace oofem
