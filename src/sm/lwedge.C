/* $Header: /home/cvs/bp/oofem/sm/src/lwedge.C,v 1.3.4.1 2012/01/018 23:15:47 bp Exp $ */
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

#include "lwedge.h"
#include "node.h"
#include "material.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "cltypes.h"
#include "structuralms.h"
#include "mathfem.h"
#include "structuralcrosssection.h"

#include <cstdio>

namespace oofem {
FEI3dWedgeLin LWedge :: interpolation;

LWedge :: LWedge(int n, Domain *aDomain) : NLStructuralElement(n, aDomain)
    // Constructor.
{
    numberOfDofMans = 6;
}


IRResultType
LWedge :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    this->NLStructuralElement :: initializeFrom(ir);
    numberOfGaussPoints = 2;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_LWedge_nip, "nip"); // Macro

    if ( ( numberOfGaussPoints != 2 ) && ( numberOfGaussPoints != 9 ) ) {
        numberOfGaussPoints = 2;
    }

    // set - up Gaussian integration points
    this->computeGaussPoints();

    return IRRT_OK;
}


void
LWedge :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
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
LWedge :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    double determinant = this->interpolation.giveTransformationJacobian(* aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));
    double weight      = aGaussPoint->giveWeight();

    return ( determinant * weight );
}


void
LWedge :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    numberOfIntegrationRules = 1;
    integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];
    integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 6);
    integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Wedge, numberOfGaussPoints, this->giveMaterialMode());
}


void
LWedge :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
    FloatArray n(6);

    answer.resize(3, 18);
    answer.zero();

    this->interpolation.evalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    for ( int i = 1; i <= 6; i++ ) {
        answer.at(1, 3 * i - 2) = n.at(i);
        answer.at(2, 3 * i - 1) = n.at(i);
        answer.at(3, 3 * i - 0) = n.at(i);
    }
}


void
LWedge :: computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint, int i) 
// Returns the [60x60] nonlinear part of strain-displacement matrix {B} of the receiver,
// evaluated at aGaussPoint

{
    int j, k, l;
    FloatMatrix dnx;

    // compute the derivatives of shape functions
    this->interpolation.evaldNdx(dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(18, 18);
    answer.zero();

    // put the products of derivatives of shape functions into the "nonlinear B matrix",
    // depending on parameter i, which is the number of the strain component
    if ( i <= 3 ) {
        for ( k = 0; k < 6; k++ ) {
            for ( l = 0; l < 3; l++ ) {
                for ( j = 1; j <= 18; j += 3 ) {
                    answer.at(k * 3 + l + 1, l + j) = dnx.at(k + 1, i) * dnx.at( ( j - 1 ) / 3 + 1, i );
                }
            }
        }
    } else if ( i == 4 )        {
        for ( k = 0; k < 6; k++ ) {
            for ( l = 0; l < 3; l++ ) {
                for ( j = 1; j <= 18; j += 3 ) {
                    answer.at(k * 3 + l + 1, l + j) = dnx.at(k + 1, 2) * dnx.at( ( j - 1 ) / 3 + 1, 3 ) + dnx.at(k + 1, 3) * dnx.at( ( j - 1 ) / 3 + 1, 2 );
                }
            }
        }
    } else if ( i == 5 )        {
        for ( k = 0; k < 6; k++ ) {
            for ( l = 0; l < 3; l++ ) {
                for ( j = 1; j <= 18; j += 3 ) {
                    answer.at(k * 3 + l + 1, l + j) = dnx.at(k + 1, 1) * dnx.at( ( j - 1 ) / 3 + 1, 3 ) + dnx.at(k + 1, 3) * dnx.at( ( j - 1 ) / 3 + 1, 1 );
                }
            }
        }
    } else if ( i == 6 )        {
        for ( k = 0; k < 6; k++ ) {
            for ( l = 0; l < 3; l++ ) {
                for ( j = 1; j <= 18; j += 3 ) {
                    answer.at(k * 3 + l + 1, l + j) = dnx.at(k + 1, 1) * dnx.at( ( j - 1 ) / 3 + 1, 2 ) + dnx.at(k + 1, 2) * dnx.at( ( j - 1 ) / 3 + 1, 1 );
                }
            }
        }
    }
}
MaterialMode
LWedge :: giveMaterialMode()
{
  if(this->nlGeometry > 1)
    return _3dMat_F;
  else
    return _3dMat;
}

void
LWedge :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the [6x60] strain-displacement matrix {B} of the receiver, eva-
// luated at aGaussPoint.
// B matrix  -  6 rows : epsilon-X, epsilon-Y, epsilon-Z, gamma-YZ, gamma-ZX, gamma-XY  :
{
    int i;
    FloatMatrix dnx;

    this->interpolation.evaldNdx(dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(6, 18);
    answer.zero();

    for ( i = 1; i <= 6; i++ ) {
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
LWedge :: computeBFmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
{
  int i, j;
  FloatMatrix dnx;
  
  this->interpolation.evaldNdx(dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

  answer.resize(9, 18);
  answer.zero();

    for ( i = 1; i <= 3; i++ ) { // 3 spatial dimensions
        for ( j = 1; j <= 6; j++ ) { // 6 nodes
            answer.at(3 * i - 2, 3 * j - 2) =
                answer.at(3 * i - 1, 3 * j - 1) =
                    answer.at(3 * i, 3 * j) = dnx.at(j, i); // derivative of Nj wrt Xi
        }
    }
}
    

Interface *
LWedge :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >( this );
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >( this );
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >( this );
    }

    OOFEM_LOG_INFO("Interface on Lwedge element not supported");
    return NULL;
}

int
LWedge :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) || ( type == IST_DamageTensor ) ) {
        return 6;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}

void
LWedge :: ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint, InternalStateType type)
{
    int i;
    FloatArray n;
    this->interpolation.evalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    if ( this->giveIPValueSize(type, aGaussPoint) ) {
        answer.resize(1, 6);
    } else {
        return;
    }

    for ( i = 1; i <= 6; i++ ) {
        answer.at(1, i)  = n.at(i);
    }
}

int
LWedge :: SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    return ZZNodalRecoveryMI_giveDofManRecordSize(type);
}

void
LWedge :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(6);
    for ( int i = 1; i <= 6; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
LWedge :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i <= 6; i++ ) {
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
LWedge :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return numberOfGaussPoints;
}


void
LWedge :: SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp)
{
    if ( this->computeGlobalCoordinates( coords, * gp->giveCoordinates() ) == 0 ) {
        _error("SPRNodalRecoveryMI_computeIPGlobalCoordinates: computeGlobalCoordinates failed");
    }
}

SPRPatchType
LWedge :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_3dBiQuadratic;
}



void
LWedge :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    int size = NodalAveragingRecoveryMI_giveDofManRecordSize(type);
    answer.resize(size);
    answer.zero();
    _warning("Qspace element: IP values will not be transferred to nodes. Use ZZNodalRecovery instead (parameter stype 1)");
}

void
LWedge :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}
} // end namespace oofem
