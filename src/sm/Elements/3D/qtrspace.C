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


#include "Elements/3D/qtrspace.h"
#include "Materials/structuralms.h"
#include "CrossSections/structuralcrosssection.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "cltypes.h"
#include "mathfem.h"
#include "fei3dtetquad.h"
#include "classfactory.h"

#include <cstdio>

namespace oofem {
REGISTER_Element(QTRSpace);

FEI3dTetQuad QTRSpace :: interpolation;

QTRSpace :: QTRSpace(int n, Domain *aDomain) : Structural3DElement(n, aDomain), ZZNodalRecoveryModelInterface(this)
{
    numberOfDofMans = 10;
}


IRResultType
QTRSpace :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = 4;
    return Structural3DElement :: initializeFrom(ir);
}


FEInterpolation *
QTRSpace :: giveInterpolation() const
{
    return & interpolation;
}


Interface *
QTRSpace :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    }

    OOFEM_LOG_INFO("Interface on QTRSpace element not supported");
    return NULL;
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
        OOFEM_ERROR("unknown node number %d", pap);
    }
}

int
QTRSpace :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return numberOfGaussPoints;
}

SPRPatchType
QTRSpace :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_3dBiQuadratic;
}


void
QTRSpace :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    answer.clear();
    OOFEM_WARNING("IP values will not be transferred to nodes. Use ZZNodalRecovery instead (parameter stype 1)");
}

} // end namespace oofem
