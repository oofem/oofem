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

#include "Elements/3D/qwedge.h"
#include "Materials/structuralms.h"
#include "CrossSections/structuralcrosssection.h"
#include "fei3dwedgequad.h"
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
#include "classfactory.h"

#include <cstdio>

namespace oofem {
REGISTER_Element(QWedge);

FEI3dWedgeQuad QWedge :: interpolation;

QWedge :: QWedge(int n, Domain *aDomain) : Structural3DElement(n, aDomain), ZZNodalRecoveryModelInterface(this)
{
    numberOfDofMans = 15;
}


IRResultType
QWedge :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = 9;

    this->matRotation = ir->hasField(_IFT_QWedge_materialCoordinateSystem);

    return Structural3DElement :: initializeFrom(ir);
}

FEInterpolation *
QWedge :: giveInterpolation() const
{
    return & interpolation;
}


void
QWedge :: giveMaterialOrientationAt(FloatArray &x, FloatArray &y, FloatArray &z, const FloatArray &lcoords)
{
    ///@todo This is subject to change. I'm not sure which is the best way to define a local c.s.
    FloatMatrix jac;
    FloatArray help;
    this->interpolation.giveJacobianMatrixAt( jac, lcoords, FEIElementGeometryWrapper(this) );
    x.beColumnOf(jac, 1); // This is {dx/dxi, dy/dxi, dz/dxi}
    x.normalize();
    help.beColumnOf(jac, 2);
    z.beVectorProductOf(x, help); // Normal to the xi-eta plane.
    z.normalize();
    y.beVectorProductOf(z, x);
}


void
QWedge :: computeStressVector(FloatArray &answer, const FloatArray &e, GaussPoint *gp, TimeStep *tStep)
{
    ///@todo This whole function should be placed in the "Space3dStructuralElementEvaluator".
    if ( this->matRotation ) {
        ///@todo This won't work properly with "useUpdatedGpRecord" (!)
        FloatArray x, y, z;
        FloatArray rotStrain, s;

        this->giveMaterialOrientationAt( x, y, z, gp->giveNaturalCoordinates() );
        // Transform from global c.s. to material c.s.
        rotStrain = {
            e(0) * x(0) * x(0) + e(5) * x(0) * x(1) + e(1) * x(1) * x(1) + e(4) * x(0) * x(2) + e(3) * x(1) * x(2) + e(2) * x(2) * x(2),
            e(0) * y(0) * y(0) + e(5) * y(0) * y(1) + e(1) * y(1) * y(1) + e(4) * y(0) * y(2) + e(3) * y(1) * y(2) + e(2) * y(2) * y(2),
            e(0) * z(0) * z(0) + e(5) * z(0) * z(1) + e(1) * z(1) * z(1) + e(4) * z(0) * z(2) + e(3) * z(1) * z(2) + e(2) * z(2) * z(2),
            2 * e(0) * y(0) * z(0) + e(4) * y(2) * z(0) + 2 * e(1) * y(1) * z(1) + e(3) * y(2) * z(1) + e(5) * ( y(1) * z(0) + y(0) * z(1) ) + ( e(4) * y(0) + e(3) * y(1) + 2 * e(2) * y(2) ) * z(2),
            2 * e(0) * x(0) * z(0) + e(4) * x(2) * z(0) + 2 * e(1) * x(1) * z(1) + e(3) * x(2) * z(1) + e(5) * ( x(1) * z(0) + x(0) * z(1) ) + ( e(4) * x(0) + e(3) * x(1) + 2 * e(2) * x(2) ) * z(2),
            2 * e(0) * x(0) * y(0) + e(4) * x(2) * y(0) + 2 * e(1) * x(1) * y(1) + e(3) * x(2) * y(1) + e(5) * ( x(1) * y(0) + x(0) * y(1) ) + ( e(4) * x(0) + e(3) * x(1) + 2 * e(2) * x(2) ) * y(2)
        };
        this->giveStructuralCrossSection()->giveRealStress_3d(s, gp, rotStrain, tStep);
        answer = {
            s(0) * x(0) * x(0) + 2 * s(5) * x(0) * y(0) + s(1) * y(0) * y(0) + 2 * ( s(4) * x(0) + s(3) * y(0) ) * z(0) + s(2) * z(0) * z(0),
            s(0) * x(1) * x(1) + 2 * s(5) * x(1) * y(1) + s(1) * y(1) * y(1) + 2 * ( s(4) * x(1) + s(3) * y(1) ) * z(1) + s(2) * z(1) * z(1),
            s(0) * x(2) * x(2) + 2 * s(5) * x(2) * y(2) + s(1) * y(2) * y(2) + 2 * ( s(4) * x(2) + s(3) * y(2) ) * z(2) + s(2) * z(2) * z(2),
            y(2) * ( s(5) * x(1) + s(1) * y(1) + s(3) * z(1) ) + x(2) * ( s(0) * x(1) + s(5) * y(1) + s(4) * z(1) ) + ( s(4) * x(1) + s(3) * y(1) + s(2) * z(1) ) * z(2),
            y(2) * ( s(5) * x(0) + s(1) * y(0) + s(3) * z(0) ) + x(2) * ( s(0) * x(0) + s(5) * y(0) + s(4) * z(0) ) + ( s(4) * x(0) + s(3) * y(0) + s(2) * z(0) ) * z(2),
            y(1) * ( s(5) * x(0) + s(1) * y(0) + s(3) * z(0) ) + x(1) * ( s(0) * x(0) + s(5) * y(0) + s(4) * z(0) ) + ( s(4) * x(0) + s(3) * y(0) + s(2) * z(0) ) * z(1)
        };
    } else {
        this->giveStructuralCrossSection()->giveRealStress_3d(answer, gp, e, tStep);
    }
}


void
QWedge :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    ///@todo This whole function should be placed in the "Space3dStructuralElementEvaluator".
    this->giveStructuralCrossSection()->giveStiffnessMatrix_3d(answer, rMode, gp, tStep);
    if ( this->matRotation ) {
        FloatArray x, y, z;
        FloatMatrix Q;

        this->giveMaterialOrientationAt( x, y, z, gp->giveNaturalCoordinates() );

        Q = {
            { x(0) * x(0), x(1) * x(1), x(2) * x(2), x(1) * x(2), x(0) * x(2), x(0) * x(1) },
            { y(0) * y(0), y(1) * y(1), y(2) * y(2), y(1) * y(2), y(0) * y(2), y(0) * y(1) },
            { z(0) * z(0), z(1) * z(1), z(2) * z(2), z(1) * z(2), z(0) * z(2), z(0) * z(1) },
            { 2 * y(0) * z(0), 2 * y(1) * z(1), 2 * y(2) * z(2), y(2) * z(1) + y(1) * z(2), y(2) * z(0) + y(0) * z(2), y(1) * z(0) + y(0) * z(1) },
            { 2 * x(0) * z(0), 2 * x(1) * z(1), 2 * x(2) * z(2), x(2) * z(1) + x(1) * z(2), x(2) * z(0) + x(0) * z(2), x(1) * z(0) + x(0) * z(1) },
            { 2 * x(0) * y(0), 2 * x(1) * y(1), 2 * x(2) * y(2), x(2) * y(1) + x(1) * y(2), x(2) * y(0) + x(0) * y(2), x(1) * y(0) + x(0) * y(1) }
        };
        //printf("original = "); answer.printYourself();
        answer.rotatedWith(Q, 't');
        //printf("rotated = "); answer.printYourself();
        //OOFEM_ERROR("DEBUG QUIT");
    }
}


Interface *
QWedge :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    }

    OOFEM_LOG_INFO("Interface on QWedge element not supported");
    return NULL;
}

void
QWedge :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(20);
    for ( int i = 1; i <= 20; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
QWedge :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i <= 20; i++ ) {
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
QWedge :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return integrationRulesArray [ 0 ]->giveNumberOfIntegrationPoints();
}


SPRPatchType
QWedge :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_3dBiQuadratic;
}


void
QWedge :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    answer.clear();
    OOFEM_WARNING("IP values will not be transferred to nodes. Use ZZNodalRecovery instead (parameter stype 1)");
}

} // end namespace oofem
