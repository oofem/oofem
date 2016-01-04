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

#include "../sm/Elements/Shells/mitc4.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/Materials/structuralmaterial.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "fei2dquadlin.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "../sm/CrossSections/variablecrosssection.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "load.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(MITC4Shell);

FEI2dQuadLin MITC4Shell :: interp_lin(1, 2);

MITC4Shell :: MITC4Shell(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(this),
    SPRNodalRecoveryModelInterface(), SpatialLocalizerInterface(this)
{
    numberOfDofMans = 4;
    numberOfGaussPoints = 8;
}


FEInterpolation *
MITC4Shell :: giveInterpolation() const { return & interp_lin; }


FEInterpolation *
MITC4Shell :: giveInterpolation(DofIDItem id) const
{
    return & interp_lin;
}


Interface *
MITC4Shell :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
    }

    return NULL;
}


void
MITC4Shell :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(numberOfDofMans);
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}


void
MITC4Shell :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i <= numberOfDofMans; i++ ) {
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
MITC4Shell :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return this->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints();
}


SPRPatchType
MITC4Shell :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_3dBiLin;
}


void
MITC4Shell :: computeGaussPoints()
// Sets up the array containing the eight Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 6) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


void
MITC4Shell ::  giveDirectorVectors(FloatArray &V1, FloatArray &V2, FloatArray &V3, FloatArray &V4)
{
    FloatArray *c1, *c2, *c3, *c4;
    V1.resize(3);
    V2.resize(3);
    V3.resize(3);
    V4.resize(3);
    c1 = this->giveNode(1)->giveCoordinates();
    c2 = this->giveNode(2)->giveCoordinates();
    c3 = this->giveNode(3)->giveCoordinates();
    c4 = this->giveNode(4)->giveCoordinates();
    V1.at(1) = this->giveCrossSection()->give(CS_DirectorVectorX, * c1, this, false);
    V1.at(2) = this->giveCrossSection()->give(CS_DirectorVectorY, * c1, this, false);
    V1.at(3) = this->giveCrossSection()->give(CS_DirectorVectorZ, * c1, this, false);

    V2.at(1) = this->giveCrossSection()->give(CS_DirectorVectorX, * c2, this, false);
    V2.at(2) = this->giveCrossSection()->give(CS_DirectorVectorY, * c2, this, false);
    V2.at(3) = this->giveCrossSection()->give(CS_DirectorVectorZ, * c2, this, false);

    V3.at(1) = this->giveCrossSection()->give(CS_DirectorVectorX, * c3, this, false);
    V3.at(2) = this->giveCrossSection()->give(CS_DirectorVectorY, * c3, this, false);
    V3.at(3) = this->giveCrossSection()->give(CS_DirectorVectorZ, * c3, this, false);

    V4.at(1) = this->giveCrossSection()->give(CS_DirectorVectorX, * c4, this, false);
    V4.at(2) = this->giveCrossSection()->give(CS_DirectorVectorY, * c4, this, false);
    V4.at(3) = this->giveCrossSection()->give(CS_DirectorVectorZ, * c4, this, false);

    V1.normalize();
    V2.normalize();
    V3.normalize();
    V4.normalize();
}



void
MITC4Shell ::  giveLocalDirectorVectors(FloatArray &V1, FloatArray &V2, FloatArray &V3, FloatArray &V4)
{
    FloatArray V1g, V2g, V3g, V4g;
    this->giveDirectorVectors(V1g, V2g, V3g, V4g);
    this->computeGtoLRotationMatrix();

    V1.beProductOf(GtoLRotationMatrix, V1g);
    V2.beProductOf(GtoLRotationMatrix, V2g);
    V3.beProductOf(GtoLRotationMatrix, V3g);
    V4.beProductOf(GtoLRotationMatrix, V4g);
}

void
MITC4Shell :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the [6x24] displacement interpolation matrix {N} of the receiver,
// evaluated at gp.
// Zeroes in rows 4, 5, 6.
{
    FloatArray h(4);

    interp_lin.evalN( h, iLocCoord,  FEIElementGeometryWrapper(this) );

    double a1, a2, a3, a4;
    this->giveThickness(a1, a2, a3, a4);

    FloatArray V1, V2, V3, V4;
    this->giveLocalDirectorVectors(V1, V2, V3, V4);

    FloatArray e2 = {
        0, 1, 0
    };

    FloatArray V11(3), V12(3), V13(3), V14(3), V21(3), V22(3), V23(3), V24(3);
    V11.beVectorProductOf(e2, V1);
    V11.normalize();
    V12.beVectorProductOf(e2, V2);
    V12.normalize();
    V13.beVectorProductOf(e2, V3);
    V13.normalize();
    V14.beVectorProductOf(e2, V4);
    V14.normalize();

    V21.beVectorProductOf(V1, V11);
    V22.beVectorProductOf(V2, V12);
    V23.beVectorProductOf(V3, V13);
    V24.beVectorProductOf(V4, V14);

    answer.resize(6, 24);
    answer.zero();

    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = h.at(1);
    answer.at(1, 4) = -iLocCoord.at(3) / 2.0 *a1 *h.at(1) * V21.at(1);
    answer.at(1, 5) =  iLocCoord.at(3) / 2.0 *a1 *h.at(1) * V11.at(1);
    answer.at(2, 4) = -iLocCoord.at(3) / 2.0 *a1 *h.at(1) * V21.at(2);
    answer.at(2, 5) =  iLocCoord.at(3) / 2.0 *a1 *h.at(1) * V11.at(2);
    answer.at(3, 4) = -iLocCoord.at(3) / 2.0 *a1 *h.at(1) * V21.at(3);
    answer.at(3, 5) =  iLocCoord.at(3) / 2.0 *a1 *h.at(1) * V11.at(3);

    answer.at(1, 7)  = answer.at(2, 8) = answer.at(3, 9) = h.at(2);
    answer.at(1, 10)  = -iLocCoord.at(3) / 2.0 *a2 *h.at(2) * V22.at(1);
    answer.at(1, 11) =  iLocCoord.at(3) / 2.0 *a2 *h.at(2) * V12.at(1);
    answer.at(2, 10)  = -iLocCoord.at(3) / 2.0 *a2 *h.at(2) * V22.at(2);
    answer.at(2, 11) =  iLocCoord.at(3) / 2.0 *a2 *h.at(2) * V12.at(2);
    answer.at(3, 10)  = -iLocCoord.at(3) / 2.0 *a2 *h.at(2) * V22.at(3);
    answer.at(3, 11) =  iLocCoord.at(3) / 2.0 *a2 *h.at(2) * V12.at(3);

    answer.at(1, 13) = answer.at(2, 14) = answer.at(3, 15) = h.at(3);
    answer.at(1, 16) = -iLocCoord.at(3) / 2.0 *a3 *h.at(3) * V23.at(1);
    answer.at(1, 17) =  iLocCoord.at(3) / 2.0 *a3 *h.at(3) * V13.at(1);
    answer.at(2, 16) = -iLocCoord.at(3) / 2.0 *a3 *h.at(3) * V23.at(2);
    answer.at(2, 17) =  iLocCoord.at(3) / 2.0 *a3 *h.at(3) * V13.at(2);
    answer.at(3, 16) = -iLocCoord.at(3) / 2.0 *a3 *h.at(3) * V23.at(3);
    answer.at(3, 17) =  iLocCoord.at(3) / 2.0 *a3 *h.at(3) * V13.at(3);

    answer.at(1, 19) = answer.at(2, 20) = answer.at(3, 21) = h.at(4);
    answer.at(1, 23) = -iLocCoord.at(3) / 2.0 *a4 *h.at(4) * V24.at(1);
    answer.at(1, 24) =  iLocCoord.at(3) / 2.0 *a4 *h.at(4) * V14.at(1);
    answer.at(2, 23) = -iLocCoord.at(3) / 2.0 *a4 *h.at(4) * V24.at(2);
    answer.at(2, 24) =  iLocCoord.at(3) / 2.0 *a4 *h.at(4) * V14.at(2);
    answer.at(3, 23) = -iLocCoord.at(3) / 2.0 *a4 *h.at(4) * V24.at(3);
    answer.at(3, 24) =  iLocCoord.at(3) / 2.0 *a4 *h.at(4) * V14.at(3);
}


void
MITC4Shell :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    SimpleCrossSection *cs = dynamic_cast< SimpleCrossSection * >( this->giveCrossSection() );
    cs->give3dDegeneratedShellStiffMtrx(answer, rMode, gp, tStep);
}


void
MITC4Shell :: giveNodeCoordinates(double &x1, double &x2, double &x3, double &x4,
                                  double &y1, double &y2, double &y3, double &y4,
                                  double &z1, double &z2, double &z3, double &z4)
{
    FloatArray nc1(3), nc2(3), nc3(3), nc4(3);

    this->giveLocalCoordinates( nc1, * ( this->giveNode(1)->giveCoordinates() ) );
    this->giveLocalCoordinates( nc2, * ( this->giveNode(2)->giveCoordinates() ) );
    this->giveLocalCoordinates( nc3, * ( this->giveNode(3)->giveCoordinates() ) );
    this->giveLocalCoordinates( nc4, * ( this->giveNode(4)->giveCoordinates() ) );

    x1 = nc1.at(1);
    x2 = nc2.at(1);
    x3 = nc3.at(1);
    x4 = nc4.at(1);

    y1 = nc1.at(2);
    y2 = nc2.at(2);
    y3 = nc3.at(2);
    y4 = nc4.at(2);

    z1 = nc1.at(3);
    z2 = nc2.at(3);
    z3 = nc3.at(3);
    z4 = nc4.at(3);
}

void
MITC4Shell :: giveLocalCoordinates(FloatArray &answer, FloatArray &global)
// Returns global coordinates given in global vector
// transformed into local coordinate system of the
// receiver
{
    FloatArray offset;
    // test the parametr
    if ( global.giveSize() != 3 ) {
        OOFEM_ERROR("cannot transform coordinates - size mismatch");
        exit(1);
    }

    this->computeGtoLRotationMatrix();

    offset = global;
    offset.subtract( * this->giveNode(1)->giveCoordinates() );
    answer.beProductOf(GtoLRotationMatrix, offset);
}

IRResultType
MITC4Shell :: initializeFrom(InputRecord *ir)
{
    return this->NLStructuralElement :: initializeFrom(ir);
}


void
MITC4Shell :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v, D_w, R_u, R_v, R_w
    };
}

double
MITC4Shell :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double detJ, weight;
    FloatMatrix jacobianMatrix(3, 3);

    FloatArray lcoords(3);
    lcoords.at(1) = gp->giveNaturalCoordinate(1);
    lcoords.at(2) = gp->giveNaturalCoordinate(2);
    lcoords.at(3) = gp->giveNaturalCoordinate(3);

    weight = gp->giveWeight();

    this->giveJacobian(gp, jacobianMatrix);

    detJ = jacobianMatrix.giveDeterminant();
    return detJ * weight;
}


void
MITC4Shell :: giveJacobian(GaussPoint *gp, FloatMatrix &jacobianMatrix)
// Returns the jacobianMatrix
{
    FloatArray h(4);
    FloatMatrix dn(4, 2);
    FloatMatrix dndx(4, 2);

    jacobianMatrix.resize(3, 3);
    // get node coordinates
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
    this->giveNodeCoordinates(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4);

    FloatArray lcoords(3);
    lcoords.at(1) = gp->giveNaturalCoordinate(1);
    lcoords.at(2) = gp->giveNaturalCoordinate(2);
    lcoords.at(3) = gp->giveNaturalCoordinate(3);
    double r3 = lcoords.at(3);

    // get local director vector
    FloatArray V1, V2, V3, V4;
    this->giveLocalDirectorVectors(V1, V2, V3, V4);

    // get thickness
    double a1, a2, a3, a4;
    this->giveThickness(a1, a2, a3, a4);

    interp_lin.evalN( h, lcoords,  FEIElementGeometryWrapper(this) );
    interp_lin.giveDerivatives(dn, lcoords);

    FloatArray hk1(4);
    // derivatives of interpolation functions
    // dh(r1,r2)/dr1
    hk1.at(1) =  dn.at(1, 1);
    hk1.at(2) =  dn.at(2, 1);
    hk1.at(3) =  dn.at(3, 1);
    hk1.at(4) =  dn.at(4, 1);

    FloatArray hk2(4);
    // dh(r1,r2)/dr2
    hk2.at(1) =  dn.at(1, 2);
    hk2.at(2) =  dn.at(2, 2);
    hk2.at(3) =  dn.at(3, 2);
    hk2.at(4) =  dn.at(4, 2);

    double h1, h2, h3, h4;
    // interpolation functions - h(r1,r2)
    h1 = h.at(1);
    h2 = h.at(2);
    h3 = h.at(3);
    h4 = h.at(4);

    // Jacobian Matrix
    jacobianMatrix.at(1, 1) = hk1.at(1) * x1 + hk1.at(2) * x2 + hk1.at(3) * x3 + hk1.at(4) * x4 + r3 / 2. * ( a1 * hk1.at(1) * V1.at(1) + a2 * hk1.at(2) * V2.at(1) + a3 * hk1.at(3) * V3.at(1) + a4 * hk1.at(4) * V4.at(1) );
    jacobianMatrix.at(1, 2) = hk1.at(1) * y1 + hk1.at(2) * y2 + hk1.at(3) * y3 + hk1.at(4) * y4 + r3 / 2. * ( a1 * hk1.at(1) * V1.at(2) + a2 * hk1.at(2) * V2.at(2) + a3 * hk1.at(3) * V3.at(2) + a4 * hk1.at(4) * V4.at(2) );
    jacobianMatrix.at(1, 3) = r3 / 2. * ( a1 * hk1.at(1) * V1.at(3) + a2 * hk1.at(2) * V2.at(3) + a3 * hk1.at(3) * V3.at(3) + a4 * hk1.at(4) * V4.at(3) );
    jacobianMatrix.at(2, 1) = hk2.at(1) * x1 + hk2.at(2) * x2 + hk2.at(3) * x3 + hk2.at(4) * x4 + r3 / 2. * ( a1 * hk2.at(1) * V1.at(1) + a2 * hk2.at(2) * V2.at(1) + a3 * hk2.at(3) * V3.at(1) + a4 * hk2.at(4) * V4.at(1) );
    jacobianMatrix.at(2, 2) = hk2.at(1) * y1 + hk2.at(2) * y2 + hk2.at(3) * y3 + hk2.at(4) * y4 + r3 / 2. * ( a1 * hk2.at(1) * V1.at(2) + a2 * hk2.at(2) * V2.at(2) + a3 * hk2.at(3) * V3.at(2) + a4 * hk2.at(4) * V4.at(2) );
    jacobianMatrix.at(2, 3) = r3 / 2. * ( a1 * hk2.at(1) * V1.at(3) + a2 * hk2.at(2) * V2.at(3) + a3 * hk2.at(3) * V3.at(3) + a4 * hk2.at(4) * V4.at(3) );
    jacobianMatrix.at(3, 1) =  1. / 2. * ( a1 * h1 * V1.at(1) + a2 * h2 * V2.at(1) + a3 * h3 * V3.at(1) + a4 * h4 * V4.at(1) );
    jacobianMatrix.at(3, 2) =  1. / 2. * ( a1 * h1 * V1.at(2) + a2 * h2 * V2.at(2) + a3 * h3 * V3.at(2) + a4 * h4 * V4.at(2) );
    jacobianMatrix.at(3, 3) =  1. / 2. * ( a1 * h1 * V1.at(3) + a2 * h2 * V2.at(3) + a3 * h3 * V3.at(3) + a4 * h4 * V4.at(3) );
}


void
MITC4Shell :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [6x20] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    FloatArray h(4);
    FloatMatrix jacobianMatrix(3, 3);
    FloatMatrix dn(4, 2);
    FloatMatrix dndx(4, 2);

    answer.resize(6, 24);
    // get node coordinates
    double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
    this->giveNodeCoordinates(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4);

    // get gp coordinates
    double r1 = gp->giveNaturalCoordinate(1);
    double r2 = gp->giveNaturalCoordinate(2);
    double r3 = gp->giveNaturalCoordinate(3);

    FloatArray lcoords(3);
    lcoords.at(1) = r1;
    lcoords.at(2) = r2;
    lcoords.at(3) = r3;

    // get director vector
    FloatArray V1, V2, V3, V4;
    this->giveLocalDirectorVectors(V1, V2, V3, V4);

    // get thickness
    double a1, a2, a3, a4;
    this->giveThickness(a1, a2, a3, a4);

    interp_lin.evalN( h, lcoords,  FEIElementGeometryWrapper(this) );
    interp_lin.giveDerivatives(dn, lcoords);

    FloatArray hk1(4);
    // derivatives of interpolation functions
    // dh(r1,r2)/dr1
    hk1.at(1) =  dn.at(1, 1);
    hk1.at(2) =  dn.at(2, 1);
    hk1.at(3) =  dn.at(3, 1);
    hk1.at(4) =  dn.at(4, 1);

    FloatArray hk2(4);
    // dh(r1,r2)/dr2
    hk2.at(1) =  dn.at(1, 2);
    hk2.at(2) =  dn.at(2, 2);
    hk2.at(3) =  dn.at(3, 2);
    hk2.at(4) =  dn.at(4, 2);

    FloatArray hkx(4), hky(4);

    // Jacobian Matrix

    this->giveJacobian(gp, jacobianMatrix);

    FloatMatrix inv(3, 3);
    FloatMatrix inv2(2, 2);
    inv.beInverseOf(jacobianMatrix);

    inv2.beSubMatrixOf(inv, 1, 2, 1, 2);
    dndx.beProductTOf(dn, inv2);

    hkx.at(1)  = dndx.at(1, 1);
    hkx.at(2)  = dndx.at(2, 1);
    hkx.at(3)  = dndx.at(3, 1);
    hkx.at(4)  = dndx.at(4, 1);

    hky.at(1)  = dndx.at(1, 2);
    hky.at(2)  = dndx.at(2, 2);
    hky.at(3)  = dndx.at(3, 2);
    hky.at(4)  = dndx.at(4, 2);

    double sb = 2 * inv.at(1, 1) * inv.at(3, 3);
    double sa = 2 * inv.at(1, 2) * inv.at(3, 3);
    double cb = 2 * inv.at(2, 1) * inv.at(3, 3);
    double ca = 2 * inv.at(2, 2) * inv.at(3, 3);

    FloatArray e2 = {
        0, 1, 0
    };

    FloatArray V11(3), V12(3), V13(3), V14(3), V21(3), V22(3), V23(3), V24(3);
    V11.beVectorProductOf(e2, V1);
    V11.normalize();
    V12.beVectorProductOf(e2, V2);
    V12.normalize();
    V13.beVectorProductOf(e2, V3);
    V13.normalize();
    V14.beVectorProductOf(e2, V4);
    V14.normalize();

    V21.beVectorProductOf(V1, V11);
    V22.beVectorProductOf(V2, V12);
    V23.beVectorProductOf(V3, V13);
    V24.beVectorProductOf(V4, V14);


    answer.zero();
    answer.at(4, 1) = 1. / 32. * ( ( a1 * V1.at(1) + a2 * V2.at(1) ) * ( cb * ( 1. + r2 ) ) + ( a1 * V1.at(1) + a4 * V4.at(1) ) * ( ca * ( 1. + r1 ) ) );
    answer.at(4, 2) = 1. / 32. * ( ( a1 * V1.at(2) + a2 * V2.at(2) ) * ( cb * ( 1. + r2 ) ) + ( a1 * V1.at(2) + a4 * V4.at(2) ) * ( ca * ( 1. + r1 ) ) );
    answer.at(4, 3) = 1. / 32. * ( ( a1 * V1.at(3) + a2 * V2.at(3) ) * ( cb * ( 1. + r2 ) ) + ( a1 * V1.at(3) + a4 * V4.at(3) ) * ( ca * ( 1. + r1 ) ) );
    answer.at(4, 4) = -a1 / 32. * ( ( V21.dotProduct({ x1 - x2, y1 - y2, z1 - z2 }) * ( cb * ( 1. + r2 ) ) ) + ( V21.dotProduct({ x1 - x4, y1 - y4, z1 - z4 }) * ( ca * ( 1. + r1 ) ) ) );
    answer.at(4, 5) = a1 / 32. * ( ( V11.dotProduct({ x1 - x2, y1 - y2, z1 - z2 }) * ( cb * ( 1. + r2 ) ) ) + ( V11.dotProduct({ x1 - x4, y1 - y4, z1 - z4 }) * ( ca * ( 1. + r1 ) ) ) );

    answer.at(4, 7) = 1. / 32. * ( -( a1 * V1.at(1) + a2 * V2.at(1) ) * ( cb * ( 1. + r2 ) ) + ( a2 * V2.at(1) + a3 * V3.at(1) ) * ( ca * ( 1. - r1 ) ) );
    answer.at(4, 8) = 1. / 32. * ( -( a1 * V1.at(2) + a2 * V2.at(2) ) * ( cb * ( 1. + r2 ) ) + ( a2 * V2.at(2) + a3 * V3.at(2) ) * ( ca * ( 1. - r1 ) ) );
    answer.at(4, 9) = 1. / 32. * ( -( a1 * V1.at(3) + a2 * V2.at(3) ) * ( cb * ( 1. + r2 ) ) + ( a2 * V2.at(3) + a3 * V3.at(3) ) * ( ca * ( 1. - r1 ) ) );
    answer.at(4, 10) = -a1 / 32. * ( ( V21.dotProduct({ x1 - x2, y1 - y2, z1 - z2 }) * ( cb * ( 1. + r2 ) ) ) + ( V21.dotProduct({ x2 - x3, y2 - y3, z2 - z3 }) * ( ca * ( 1. - r1 ) ) ) );
    answer.at(4, 11) = a1 / 32. * ( ( V11.dotProduct({ x1 - x2, y1 - y2, z1 - z2 }) * ( cb * ( 1. + r2 ) ) ) + ( V11.dotProduct({ x2 - x3, y2 - y3, z2 - z3 }) * ( ca * ( 1. - r1 ) ) ) );

    answer.at(4, 13) = 1. / 32. * ( -( a3 * V3.at(1) + a4 * V4.at(1) ) * ( cb * ( 1. - r2 ) ) - ( a2 * V2.at(1) + a3 * V3.at(1) ) * ( ca * ( 1. - r1 ) ) );
    answer.at(4, 14) = 1. / 32. * ( -( a3 * V3.at(2) + a4 * V4.at(2) ) * ( cb * ( 1. - r2 ) ) - ( a2 * V2.at(2) + a3 * V3.at(2) ) * ( ca * ( 1. - r1 ) ) );
    answer.at(4, 15) = 1. / 32. * ( -( a3 * V3.at(3) + a4 * V4.at(3) ) * ( cb * ( 1. - r2 ) ) - ( a2 * V2.at(3) + a3 * V3.at(3) ) * ( ca * ( 1. - r1 ) ) );
    answer.at(4, 16) = -a1 / 32. * ( ( V21.dotProduct({ x4 - x3, y4 - y3, z4 - z3 }) * ( cb * ( 1. - r2 ) ) ) + ( V21.dotProduct({ x2 - x3, y2 - y3, z2 - z3 }) * ( ca * ( 1. - r1 ) ) ) );
    answer.at(4, 17) = a1 / 32. * ( ( V11.dotProduct({ x4 - x3, y4 - y3, z4 - z3 }) * ( cb * ( 1. - r2 ) ) ) + ( V11.dotProduct({ x2 - x3, y2 - y3, z2 - z3 }) * ( ca * ( 1. - r1 ) ) ) );

    answer.at(4, 19) = 1. / 32. * ( ( a3 * V3.at(1) + a4 * V4.at(1) ) * ( cb * ( 1. - r2 ) ) - ( a1 * V1.at(1) + a4 * V4.at(1) ) * ( ca * ( 1. + r1 ) ) );
    answer.at(4, 20) = 1. / 32. * ( ( a3 * V3.at(2) + a4 * V4.at(2) ) * ( cb * ( 1. - r2 ) ) - ( a1 * V1.at(2) + a4 * V4.at(2) ) * ( ca * ( 1. + r1 ) ) );
    answer.at(4, 21) = 1. / 32. * ( ( a3 * V3.at(3) + a4 * V4.at(3) ) * ( cb * ( 1. - r2 ) ) - ( a1 * V1.at(3) + a4 * V4.at(3) ) * ( ca * ( 1. + r1 ) ) );
    answer.at(4, 22) = -a1 / 32. * ( ( V21.dotProduct({ x4 - x3, y4 - y3, z4 - z3 }) * ( cb * ( 1. - r2 ) ) ) + ( V21.dotProduct({ x1 - x4, y1 - y4, z1 - z4 }) * ( ca * ( 1. + r1 ) ) ) );
    answer.at(4, 23) = a1 / 32. * ( ( V11.dotProduct({ x4 - x3, y4 - y3, z4 - z3 }) * ( cb * ( 1. - r2 ) ) ) + ( V11.dotProduct({ x1 - x4, y1 - y4, z1 - z4 }) * ( ca * ( 1. + r1 ) ) ) );

    answer.at(5, 1) = 1. / 32. * ( ( a1 * V1.at(1) + a2 * V2.at(1) ) * ( sb * ( 1. + r2 ) ) + ( a1 * V1.at(1) + a4 * V4.at(1) ) * ( sa * ( 1. + r1 ) ) );
    answer.at(5, 2) = 1. / 32. * ( ( a1 * V1.at(2) + a2 * V2.at(2) ) * ( sb * ( 1. + r2 ) ) + ( a1 * V1.at(2) + a4 * V4.at(2) ) * ( sa * ( 1. + r1 ) ) );
    answer.at(5, 3) = 1. / 32. * ( ( a1 * V1.at(3) + a2 * V2.at(3) ) * ( sb * ( 1. + r2 ) ) + ( a1 * V1.at(3) + a4 * V4.at(3) ) * ( sa * ( 1. + r1 ) ) );
    answer.at(5, 4) = -a1 / 32. * ( ( V21.dotProduct({ x1 - x2, y1 - y2, z1 - z2 }) * ( sb * ( 1. + r2 ) ) ) + ( V21.dotProduct({ x1 - x4, y1 - y4, z1 - z4 }) * ( sa * ( 1. + r1 ) ) ) );
    answer.at(5, 5) = a1 / 32. * ( ( V11.dotProduct({ x1 - x2, y1 - y2, z1 - z2 }) * ( sb * ( 1. + r2 ) ) ) + ( V11.dotProduct({ x1 - x4, y1 - y4, z1 - z4 }) * ( sa * ( 1. + r1 ) ) ) );

    answer.at(5, 7) = 1. / 32. * ( -( a1 * V1.at(1) + a2 * V2.at(1) ) * ( sb * ( 1. + r2 ) ) + ( a2 * V2.at(1) + a3 * V3.at(1) ) * ( sa * ( 1. - r1 ) ) );
    answer.at(5, 8) = 1. / 32. * ( -( a1 * V1.at(2) + a2 * V2.at(2) ) * ( sb * ( 1. + r2 ) ) + ( a2 * V2.at(2) + a3 * V3.at(2) ) * ( sa * ( 1. - r1 ) ) );
    answer.at(5, 9) = 1. / 32. * ( -( a1 * V1.at(3) + a2 * V2.at(3) ) * ( sb * ( 1. + r2 ) ) + ( a2 * V2.at(3) + a3 * V3.at(3) ) * ( sa * ( 1. - r1 ) ) );
    answer.at(5, 10) = -a1 / 32. * ( ( V21.dotProduct({ x1 - x2, y1 - y2, z1 - z2 }) * ( sb * ( 1. + r2 ) ) ) + ( V21.dotProduct({ x2 - x3, y2 - y3, z2 - z3 }) * ( sa * ( 1. - r1 ) ) ) );
    answer.at(5, 11) = a1 / 32. * ( ( V11.dotProduct({ x1 - x2, y1 - y2, z1 - z2 }) * ( sb * ( 1. + r2 ) ) ) + ( V11.dotProduct({ x2 - x3, y2 - y3, z2 - z3 }) * ( sa * ( 1. - r1 ) ) ) );

    answer.at(5, 13) = 1. / 32. * ( -( a3 * V3.at(1) + a4 * V4.at(1) ) * ( sb * ( 1. - r2 ) ) - ( a2 * V2.at(1) + a3 * V3.at(1) ) * ( sa * ( 1. - r1 ) ) );
    answer.at(5, 14) = 1. / 32. * ( -( a3 * V3.at(2) + a4 * V4.at(2) ) * ( sb * ( 1. - r2 ) ) - ( a2 * V2.at(2) + a3 * V3.at(2) ) * ( sa * ( 1. - r1 ) ) );
    answer.at(5, 15) = 1. / 32. * ( -( a3 * V3.at(3) + a4 * V4.at(3) ) * ( sb * ( 1. - r2 ) ) - ( a2 * V2.at(3) + a3 * V3.at(3) ) * ( sa * ( 1. - r1 ) ) );
    answer.at(5, 16) = -a1 / 32. * ( ( V21.dotProduct({ x4 - x3, y4 - y3, z4 - z3 }) * ( sb * ( 1. - r2 ) ) ) + ( V21.dotProduct({ x2 - x3, y2 - y3, z2 - z3 }) * ( sa * ( 1. - r1 ) ) ) );
    answer.at(5, 17) = a1 / 32. * ( ( V11.dotProduct({ x4 - x3, y4 - y3, z4 - z3 }) * ( sb * ( 1. - r2 ) ) ) + ( V11.dotProduct({ x2 - x3, y2 - y3, z2 - z3 }) * ( sa * ( 1. - r1 ) ) ) );

    answer.at(5, 19) = 1. / 32. * ( ( a3 * V3.at(1) + a4 * V4.at(1) ) * ( sb * ( 1. - r2 ) ) - ( a1 * V1.at(1) + a4 * V4.at(1) ) * ( sa * ( 1. + r1 ) ) );
    answer.at(5, 20) = 1. / 32. * ( ( a3 * V3.at(2) + a4 * V4.at(2) ) * ( sb * ( 1. - r2 ) ) - ( a1 * V1.at(2) + a4 * V4.at(2) ) * ( sa * ( 1. + r1 ) ) );
    answer.at(5, 21) = 1. / 32. * ( ( a3 * V3.at(3) + a4 * V4.at(3) ) * ( sb * ( 1. - r2 ) ) - ( a1 * V1.at(3) + a4 * V4.at(3) ) * ( sa * ( 1. + r1 ) ) );
    answer.at(5, 22) = -a1 / 32. * ( ( V21.dotProduct({ x4 - x3, y4 - y3, z4 - z3 }) * ( sb * ( 1. - r2 ) ) ) + ( V21.dotProduct({ x1 - x4, y1 - y4, z1 - z4 }) * ( sa * ( 1. + r1 ) ) ) );
    answer.at(5, 23) = a1 / 32. * ( ( V11.dotProduct({ x4 - x3, y4 - y3, z4 - z3 }) * ( sb * ( 1. - r2 ) ) ) + ( V11.dotProduct({ x1 - x4, y1 - y4, z1 - z4 }) * ( sa * ( 1. + r1 ) ) ) );



    answer.at(1, 1) = hkx.at(1);
    answer.at(1, 4) = -r3 / 2. *a1 *hkx.at(1) * V21.at(1);
    answer.at(1, 5) = r3 / 2. *a1 *hkx.at(1) * V11.at(1);
    answer.at(1, 7) = hkx.at(2);
    answer.at(1, 10) = -r3 / 2. *a2 *hkx.at(2) * V22.at(1);
    answer.at(1, 11) = r3 / 2. *a2 *hkx.at(2) * V12.at(1);
    answer.at(1, 13) =  hkx.at(3);
    answer.at(1, 16) = -r3 / 2. *a3 *hkx.at(3) * V23.at(1);
    answer.at(1, 17) = r3 / 2. *a3 *hkx.at(3) * V13.at(1);
    answer.at(1, 19) = hkx.at(4);
    answer.at(1, 22) = -r3 / 2. *a4 *hkx.at(4) * V24.at(1);
    answer.at(1, 23) = r3 / 2. *a4 *hkx.at(4) * V14.at(1);

    answer.at(2, 2) = hky.at(1);
    answer.at(2, 4) = -r3 / 2. *a1 *hky.at(1) * V21.at(2);
    answer.at(2, 5) = r3 / 2. *a1 *hky.at(1) * V11.at(2);

    answer.at(2, 8) = hky.at(2);
    answer.at(2, 10) = -r3 / 2. *a2 *hky.at(2) * V22.at(2);
    answer.at(2, 11) = r3 / 2. *a2 *hky.at(2) * V12.at(2);

    answer.at(2, 14) = hky.at(3);
    answer.at(2, 16) = -r3 / 2. *a3 *hky.at(3) * V23.at(2);
    answer.at(2, 17) = r3 / 2. *a3 *hky.at(3) * V13.at(2);

    answer.at(2, 20) = hky.at(4);
    answer.at(2, 22) = -r3 / 2. *a4 *hky.at(4) * V24.at(2);
    answer.at(2, 23) = r3 / 2. *a4 *hky.at(4) * V14.at(2);

    answer.at(6, 1) = hky.at(1);
    answer.at(6, 2) = hkx.at(1);
    answer.at(6, 4) = -r3 / 2. * a1 * ( hkx.at(1) * V21.at(2) + hky.at(1) * V21.at(1) );
    answer.at(6, 5) = r3 / 2. * a1 * ( hky.at(1) * V11.at(1) + hky.at(1) * V11.at(2) );

    answer.at(6, 7) = hky.at(2);
    answer.at(6, 8) = hkx.at(2);
    answer.at(6, 10) = -r3 / 2. * a2 * ( hkx.at(2) * V22.at(2) + hky.at(2) * V22.at(1) );
    answer.at(6, 11) = r3 / 2. * a2 * ( hky.at(2) * V12.at(1) + hky.at(2) * V12.at(2) );

    answer.at(6, 13) = hky.at(3);
    answer.at(6, 14) = hkx.at(3);
    answer.at(6, 16) = -r3 / 2. * a3 * ( hkx.at(3) * V23.at(2) + hky.at(3) * V23.at(1) );
    answer.at(6, 17) = r3 / 2. * a3 * ( hky.at(3) * V13.at(1) + hky.at(3) * V13.at(2) );

    answer.at(6, 19) = hky.at(4);
    answer.at(6, 20) = hkx.at(4);
    answer.at(6, 22) = -r3 / 2. * a4 * ( hkx.at(4) * V24.at(2) + hky.at(4) * V24.at(1) );
    answer.at(6, 23) = r3 / 2. * a4 * ( hky.at(4) * V14.at(1) + hky.at(4) * V14.at(2) );
}

void
MITC4Shell :: giveThickness(double &a1, double &a2, double &a3, double &a4)
{
    FloatArray *c1, *c2, *c3, *c4;
    c1 = this->giveNode(1)->giveCoordinates();
    c2 = this->giveNode(2)->giveCoordinates();
    c3 = this->giveNode(3)->giveCoordinates();
    c4 = this->giveNode(4)->giveCoordinates();

    a1 = this->giveCrossSection()->give(CS_Thickness, * c1, this, false);
    a2 = this->giveCrossSection()->give(CS_Thickness, * c2, this, false);
    a3 = this->giveCrossSection()->give(CS_Thickness, * c3, this, false);
    a4 = this->giveCrossSection()->give(CS_Thickness, * c4, this, false);
}


const FloatMatrix *
MITC4Shell :: computeGtoLRotationMatrix()
// Returns the rotation matrix of the receiver of the size [3,3]
// coords(local) = T * coords(global)
//
// local coordinate (described by vector triplet e1',e2',e3') is defined as follows:
//
// e1'    : [N4-N3]    Ni - means i - th node
// help   : [N2-N3]
// e3'    : e1' x help
// e2'    : e3' x e1'
{
    if ( !GtoLRotationMatrix.isNotEmpty() ) {
        FloatArray e1, e2, e3;


        this->computeLocalBaseVectors(e1, e2, e3);

        GtoLRotationMatrix.resize(3, 3);

        for ( int i = 1; i <= 3; i++ ) {
            GtoLRotationMatrix.at(1, i) = e1.at(i);
            GtoLRotationMatrix.at(2, i) = e2.at(i);
            GtoLRotationMatrix.at(3, i) = e3.at(i);
        }
    }

    return & GtoLRotationMatrix;
}

void
MITC4Shell :: computeLocalBaseVectors(FloatArray &e1, FloatArray &e2, FloatArray &e3)
{
    FloatArray help;

    FloatArray coordA, coordB;

    // compute A - (node1+node4)/2
    coordB.beDifferenceOf( * this->giveNode(1)->giveCoordinates(), * this->giveNode(4)->giveCoordinates() );
    coordB.times(0.5);
    coordB.add( * this->giveNode(4)->giveCoordinates() );

    // compute B - (node2+node3)/2
    coordA.beDifferenceOf( * this->giveNode(2)->giveCoordinates(), * this->giveNode(3)->giveCoordinates() );
    coordA.times(0.5);
    coordA.add( * this->giveNode(3)->giveCoordinates() );

    // compute e1' = [B-A]
    e1.beDifferenceOf(coordB, coordA);



    // compute A - (node2+node1)/2
    coordB.beDifferenceOf( * this->giveNode(1)->giveCoordinates(), * this->giveNode(2)->giveCoordinates() );
    coordB.times(0.5);
    coordB.add( * this->giveNode(2)->giveCoordinates() );

    // compute B - (node3+node4)/2
    coordA.beDifferenceOf( * this->giveNode(4)->giveCoordinates(), * this->giveNode(3)->giveCoordinates() );
    coordA.times(0.5);
    coordA.add( * this->giveNode(3)->giveCoordinates() );

    // compute e1' = [B-A]
    help.beDifferenceOf(coordB, coordA);

    // let us normalize e1'
    e1.normalize();

    // compute e3' : vector product of e1' x help
    e3.beVectorProductOf(e1, help);
    e3.normalize();

    // now from e3' x e1' compute e2'
    e2.beVectorProductOf(e3, e1);
}

void
MITC4Shell :: computeLToDirectorRotationMatrix(FloatMatrix &answer1, FloatMatrix &answer2, FloatMatrix &answer3, FloatMatrix &answer4)
// Returns the rotation matrix of the reciever of the size [2,3]
// {alpha_i,beta_i} = Ti * {rotL_xi, rotL_yi, rotL_zi}
//      alpha_i, beta_i - rotations about the components of director vector at node i
//      r1_i, r2_i, r3_i, - rotations about local coordinates e1', e2', e3'
//
// local coordinate (described by vector triplet e1',e2',e3') is defined as follows:
//
// e1'    : [N4-N3]    Ni - means i - th node
// help   : [N2-N3]
// e3'    : e1' x help
// e2'    : e3' x e1'
{
    FloatArray e1, e2, e3;

    this->computeLocalBaseVectors(e1, e2, e3);

    // get director vector
    FloatArray V1, V2, V3, V4;
    this->giveDirectorVectors(V1, V2, V3, V4);

    FloatArray V11(3), V12(3), V13(3), V14(3), V21(3), V22(3), V23(3), V24(3);
    V11.beVectorProductOf(e2, V1);
    V11.normalize();
    V12.beVectorProductOf(e2, V2);
    V12.normalize();
    V13.beVectorProductOf(e2, V3);
    V13.normalize();
    V14.beVectorProductOf(e2, V4);
    V14.normalize();

    V21.beVectorProductOf(V1, V11);
    V22.beVectorProductOf(V2, V12);
    V23.beVectorProductOf(V3, V13);
    V24.beVectorProductOf(V4, V14);

    answer1.resize(3, 3);
    answer2.resize(3, 3);
    answer3.resize(3, 3);
    answer4.resize(3, 3);

    answer1.at(1, 1) = V11.dotProduct(e1);
    answer1.at(1, 2) = V11.dotProduct(e2);
    answer1.at(1, 3) = V11.dotProduct(e3);
    answer1.at(2, 1) = V21.dotProduct(e1);
    answer1.at(2, 2) = V21.dotProduct(e2);
    answer1.at(2, 3) = V21.dotProduct(e3);
    answer1.at(3, 1) = V1.dotProduct(e1);
    answer1.at(3, 2) = V1.dotProduct(e2);
    answer1.at(3, 3) = V1.dotProduct(e3);

    answer2.at(1, 1) = V12.dotProduct(e1);
    answer2.at(1, 2) = V12.dotProduct(e2);
    answer2.at(1, 3) = V12.dotProduct(e3);
    answer2.at(2, 1) = V22.dotProduct(e1);
    answer2.at(2, 2) = V22.dotProduct(e2);
    answer2.at(2, 3) = V22.dotProduct(e3);
    answer2.at(3, 1) = V2.dotProduct(e1);
    answer2.at(3, 2) = V2.dotProduct(e2);
    answer2.at(3, 3) = V2.dotProduct(e3);

    answer3.at(1, 1) = V13.dotProduct(e1);
    answer3.at(1, 2) = V13.dotProduct(e2);
    answer3.at(1, 3) = V13.dotProduct(e3);
    answer3.at(2, 1) = V23.dotProduct(e1);
    answer3.at(2, 2) = V23.dotProduct(e2);
    answer3.at(2, 3) = V23.dotProduct(e3);
    answer3.at(3, 1) = V3.dotProduct(e1);
    answer3.at(3, 2) = V3.dotProduct(e2);
    answer3.at(3, 3) = V3.dotProduct(e3);

    answer4.at(1, 1) = V14.dotProduct(e1);
    answer4.at(1, 2) = V14.dotProduct(e2);
    answer4.at(1, 3) = V14.dotProduct(e3);
    answer4.at(2, 1) = V24.dotProduct(e1);
    answer4.at(2, 2) = V24.dotProduct(e2);
    answer4.at(2, 3) = V24.dotProduct(e3);
    answer4.at(3, 1) = V4.dotProduct(e1);
    answer4.at(3, 2) = V4.dotProduct(e2);
    answer4.at(3, 3) = V4.dotProduct(e3);
}

bool
MITC4Shell :: computeGtoLRotationMatrix(FloatMatrix &answer)
// Returns the rotation matrix of the receiver of the size [20,24]
// r(local) = T * r(global)
// for one node (r written transposed): {u,v,w,alpha,beta} = T * {u,v,w,r1,r2,r3}

{
    this->computeGtoLRotationMatrix();

    answer.resize(24, 24);
    answer.zero();

    FloatMatrix LtoDir1, LtoDir2, LtoDir3, LtoDir4;
    this->computeLToDirectorRotationMatrix(LtoDir1, LtoDir2, LtoDir3, LtoDir4);

    for ( int i = 0; i <= 3; i++ ) {
        answer.setSubMatrix(GtoLRotationMatrix, i * 6 + 1, i * 6 + 1);
    }
    ////
    FloatMatrix help;

    help.beProductOf(LtoDir1, GtoLRotationMatrix);
    answer.setSubMatrix(help, 4, 4);

    help.beProductOf(LtoDir2, GtoLRotationMatrix);
    answer.setSubMatrix(help, 10, 10);

    help.beProductOf(LtoDir3, GtoLRotationMatrix);
    answer.setSubMatrix(help, 16, 16);

    help.beProductOf(LtoDir4, GtoLRotationMatrix);
    answer.setSubMatrix(help, 22, 22);

    return 1;
}


void
MITC4Shell :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveRealStress_3dDegeneratedShell(answer, gp, strain, tStep);
}


void
MITC4Shell :: giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep)
// returns characteristic tensor of the receiver at given gp and tStep
{
    answer.resize(3, 3);
    answer.zero();
    this->computeGtoLRotationMatrix();
    StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( this->giveStructuralCrossSection()->giveMaterial(gp) );

    if ( ( type == GlobalForceTensor ) ) {
        FloatArray stress, localStress, localStrain;
        this->computeStrainVector(localStrain, gp, tStep);
        this->computeStressVector(localStress, localStrain, gp, tStep);
        mat->transformStressVectorTo(stress,  GtoLRotationMatrix, localStress, false);

        answer.at(1, 1) = stress.at(1);
        answer.at(2, 2) = stress.at(2);
        answer.at(3, 3) = stress.at(3);
        answer.at(1, 2) = stress.at(4);
        answer.at(2, 1) = stress.at(4);
        answer.at(2, 3) = stress.at(5);
        answer.at(3, 2) = stress.at(5);
        answer.at(1, 3) = stress.at(6);
        answer.at(3, 1) = stress.at(6);
    } else if ( ( type == GlobalStrainTensor ) ) {
        FloatArray strain, localStrain;
        this->computeStrainVector(localStrain, gp, tStep);
        mat->transformStrainVectorTo(strain,  GtoLRotationMatrix, localStrain, false);

        answer.at(1, 1) = strain.at(1);
        answer.at(2, 2) = strain.at(2);
        answer.at(3, 3) = strain.at(3);
        answer.at(2, 3) = strain.at(4) / 2.;
        answer.at(3, 2) = strain.at(4) / 2.;
        answer.at(1, 3) = strain.at(5) / 2.;
        answer.at(3, 1) = strain.at(5) / 2.;
        answer.at(1, 2) = strain.at(6) / 2.;
        answer.at(2, 1) = strain.at(6) / 2.;
    } else {
        OOFEM_ERROR("unsupported tensor mode");
        exit(1);
    }
}

void
MITC4Shell :: printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
    FloatArray v;

    fprintf(file, "element %d (%8d):\n", this->giveLabel(), number);

    for ( GaussPoint *gp : *integrationRulesArray [ 0 ] ) {
        fprintf( file, "  GP 1.%d :", gp->giveNumber() );
        this->giveIPValue(v, gp, IST_ShellStrainTensor, tStep);
        fprintf(file, "  strains    ");
        for ( auto &val : v ) {
            fprintf(file, " %.4e", val);
        }

        this->giveIPValue(v, gp, IST_ShellForceTensor, tStep);
        fprintf(file, "\n              stresses   ");
        for ( auto &val : v ) {
            fprintf(file, " %.4e", val);
        }

        fprintf(file, "\n");
    }
}



int
MITC4Shell :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FloatMatrix globTensor;
    CharTensor cht;

    answer.resize(6);

    if (  type == IST_ShellStrainTensor ) {
        cht = GlobalStrainTensor;

        this->giveCharacteristicTensor(globTensor, cht, gp, tStep);

        answer.at(1) = globTensor.at(1, 1); //xx
        answer.at(2) = globTensor.at(2, 2); //yy
        answer.at(3) = globTensor.at(3, 3); //zz
        answer.at(4) = 2 * globTensor.at(2, 3); //yz
        answer.at(5) = 2 * globTensor.at(1, 3); //xz
        answer.at(6) = 2 * globTensor.at(2, 3); //yz

        return 1;
    } else if ( type == IST_ShellForceTensor ) {
        cht = GlobalForceTensor;

        this->giveCharacteristicTensor(globTensor, cht, gp, tStep);

        answer.at(1) = globTensor.at(1, 1); //xx
        answer.at(2) = globTensor.at(2, 2); //yy
        answer.at(3) = globTensor.at(3, 3); //zz
        answer.at(4) = globTensor.at(2, 3); //yz
        answer.at(5) = globTensor.at(1, 3); //xz
        answer.at(6) = globTensor.at(1, 2); //xy

        return 1;
    } else {
        return NLStructuralElement :: giveIPValue(answer, gp, type, tStep);
    }
}


bool
MITC4Shell :: computeLocalCoordinates(FloatArray &answer, const FloatArray &coords)
//converts global coordinates to local planar area coordinates,
//does not return a coordinate in the thickness direction, but
//does check that the point is in the element thickness
{
    // rotate the input point Coordinate System into the element CS
    FloatArray inputCoords_ElCS;
    std :: vector< FloatArray >lc(3);
    FloatArray llc;
    this->giveLocalCoordinates( inputCoords_ElCS, const_cast< FloatArray & >(coords) );
    for ( int _i = 0; _i < 4; _i++ ) {
        this->giveLocalCoordinates( lc [ _i ], * this->giveNode(_i + 1)->giveCoordinates() );
    }
    bool inplane = interp_lin.global2local( llc, inputCoords_ElCS, FEIVertexListGeometryWrapper(lc) ) > 0;
    answer.resize(2);
    answer.at(1) = inputCoords_ElCS.at(1);
    answer.at(2) = inputCoords_ElCS.at(2);
    GaussPoint _gp(NULL, 1, answer, 2.0, _2dPlate);
    // now check if the third local coordinate is within the thickness of element
    bool outofplane = ( fabs( inputCoords_ElCS.at(3) ) <= this->giveCrossSection()->give(CS_Thickness, & _gp) / 2. );

    return inplane && outofplane;
}



int
MITC4Shell :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    double l1 = lcoords.at(1);
    double l2 = lcoords.at(2);
    double l3 = lcoords.at(3);

    answer.resize(3);
    for ( int _i = 1; _i <= 3; _i++ ) {
        answer.at(_i) = l1 * this->giveNode(1)->giveCoordinate(_i) + l2 *this->giveNode(2)->giveCoordinate(_i) + l3 *this->giveNode(3)->giveCoordinate(_i);
    }
    return true;
}


int
MITC4Shell :: computeLoadGToLRotationMtrx(FloatMatrix &answer)
// Returns the rotation matrix of the receiver of the size [5,6]
// f(local) = T * f(global)
{
    this->computeGtoLRotationMatrix();

    answer.resize(6, 6);
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, i) = answer.at(4, i + 3) = GtoLRotationMatrix.at(1, i);
        answer.at(2, i) = answer.at(5, i + 3) = GtoLRotationMatrix.at(2, i);
        answer.at(3, i) = answer.at(6, i + 3) = GtoLRotationMatrix.at(3, i);
    }

    return 1;
}


void
MITC4Shell :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                         InternalStateType type, TimeStep *tStep)
{
    double x1 = 0.0, x2 = 0.0, x3 = 0.0, y = 0.0;
    FloatMatrix A(4, 4);
    FloatMatrix b, r;
    FloatArray val;
    double u, v, w;

    int size = 0;

    for ( GaussPoint *gp : *integrationRulesArray [ 0 ] ) {
        giveIPValue(val, gp, type, tStep);
        if ( size == 0 ) {
            size = val.giveSize();
            b.resize(4, size);
            r.resize(4, size);
            A.zero();
            r.zero();
        }

        const FloatArray &coord = gp->giveNaturalCoordinates();
        u = coord.at(1);
        v = coord.at(2);
        w = coord.at(3);

        A.at(1, 1) += 1;
        A.at(1, 2) += u;
        A.at(1, 3) += v;
        A.at(1, 4) += w;
        A.at(2, 1) += u;
        A.at(2, 2) += u * u;
        A.at(2, 3) += u * v;
        A.at(2, 4) += u * w;
        A.at(3, 1) += v;
        A.at(3, 2) += v * u;
        A.at(3, 3) += v * v;
        A.at(3, 4) += v * w;
        A.at(4, 1) += w;
        A.at(4, 2) += w * u;
        A.at(4, 3) += w * v;
        A.at(4, 4) += w * w;

        for ( int j = 1; j <= size; j++ ) {
            y = val.at(j);
            r.at(1, j) += y;
            r.at(2, j) += y * u;
            r.at(3, j) += y * v;
            r.at(4, j) += y * w;
        }
    }

    A.solveForRhs(r, b);

    switch ( node ) {
    case 1:
        x1 =  1.0;
        x2 =  1.0;
        x3 =  1.0;
        break;
    case 2:
        x1 = -1.0;
        x2 =  1.0;
        x3 =  1.0;
        break;
    case 3:
        x1 = -1.0;
        x2 = -1.0;
        x3 =  1.0;
        break;
    case 4:
        x1 =  1.0;
        x2 = -1.0;
        x3 =  1.0;
        break;
    default:
        OOFEM_ERROR("unsupported node");
    }

    answer.resize(size);
    for ( int j = 1; j <= size; j++ ) {
        answer.at(j) = b.at(1, j) + x1 *b.at(2, j) * x2 * b.at(3, j) * x3 * b.at(4, j);
    }
}


void
MITC4Shell :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
{
    IntArray edgeNodes;
    FloatArray n;

    this->interp_lin.edgeEvalN( n, iedge, gp->giveNaturalCoordinates(), FEIVoidCellGeometry() );
    this->interp_lin.computeLocalEdgeMapping(edgeNodes, iedge);

    answer.beNMatrixOf(n, 6);
}


void
MITC4Shell :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    if ( iEdge == 1 ) { // edge between nodes 1 2
        answer = {
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
        };
    } else if ( iEdge == 2 ) { // edge between nodes 2 3
        answer = {
            7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18
        };
    } else if ( iEdge == 3 ) { // edge between nodes 3 4
        answer = {
            13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24
        };
    } else if ( iEdge == 4 ) { // edge between nodes 4 1
        answer = {
            19, 20, 21, 22, 23, 24, 1, 2, 3, 4, 5, 6
        };
    } else {
        OOFEM_ERROR("wrong edge number");
    }
}


double
MITC4Shell :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    std :: vector< FloatArray >lc = {
        FloatArray(3), FloatArray(3), FloatArray(3), FloatArray(3)
    };
    this->giveNodeCoordinates( lc [ 0 ].at(1), lc [ 1 ].at(1), lc [ 2 ].at(1), lc [ 3 ].at(1),
                               lc [ 0 ].at(2), lc [ 1 ].at(2), lc [ 2 ].at(2), lc [ 3 ].at(2),
                               lc [ 0 ].at(3), lc [ 1 ].at(3), lc [ 2 ].at(3), lc [ 3 ].at(3) );


    double detJ = this->interp_lin.edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(), FEIVertexListGeometryWrapper(lc) );
    return detJ * gp->giveWeight();
}

/*
 * void
 * MITC4Shell :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
 * {
 * std :: vector< FloatArray > lc = {FloatArray(3), FloatArray(3), FloatArray(3), FloatArray(3)};
 * this->giveNodeCoordinates(lc[0].at(1), lc[1].at(1), lc[2].at(1), lc[3].at(1),
 *                          lc[0].at(2), lc[1].at(2), lc[2].at(2), lc[3].at(2),
 *                          lc[0].at(3), lc[1].at(3), lc[2].at(3), lc[3].at(3));
 *
 * FloatArray local;
 * this->interp_lin.edgeLocal2global( local, iEdge, gp->giveNaturalCoordinates(), FEIVertexListGeometryWrapper(lc)  );
 * local.resize(3);
 * local.at(3) = 0.;
 *
 * answer = local; // MITC4 - todo
 * // answer.beProductOf(this->lcsMatrix, local);
 * }
 */

int
MITC4Shell :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    FloatArray e1, e2, e3, xl, yl;
    this->computeLocalBaseVectors(e1, e2, e3);

    IntArray edgeNodes;
    this->interp_lin.computeLocalEdgeMapping(edgeNodes, iEdge);

    xl.beDifferenceOf( * this->giveNode( edgeNodes.at(2) )->giveCoordinates(), * this->giveNode( edgeNodes.at(1) )->giveCoordinates() );

    xl.normalize();
    yl.beVectorProductOf(e3, xl);

    answer.resize(6, 6);
    answer.zero();

    answer.at(1, 1) = answer.at(4, 4) = e1.dotProduct(xl);
    answer.at(1, 2) = answer.at(4, 5) = e1.dotProduct(yl);
    answer.at(1, 3) = answer.at(4, 6) = e1.dotProduct(e3);
    answer.at(2, 1) = answer.at(5, 4) = e2.dotProduct(xl);
    answer.at(2, 2) = answer.at(5, 5) = e2.dotProduct(yl);
    answer.at(2, 3) = answer.at(5, 6) = e2.dotProduct(e3);
    answer.at(3, 1) = answer.at(6, 4) = e3.dotProduct(xl);
    answer.at(3, 2) = answer.at(6, 5) = e3.dotProduct(yl);
    answer.at(3, 3) = answer.at(6, 6) = e3.dotProduct(e3);


    return 1;
}



void
MITC4Shell :: computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *sgp)
{
    const FloatArray &coords2 = sgp->giveNaturalCoordinates();
    FloatArray coords(3);
    coords.at(1) = coords2.at(1);
    coords.at(2) = coords2.at(2);
    coords.at(3) = 0.0;
    this->computeNmatrixAt(coords, answer);
}

void
MITC4Shell :: giveSurfaceDofMapping(IntArray &answer, int iSurf) const
{
    int i;
    answer.resize(24);
    answer.zero();
    if ( iSurf == 1 ) {
        for ( i = 1; i <= 24; i++ ) {
            answer.at(i) = i;
        }
    } else {
        OOFEM_ERROR("wrong surface number");
    }
}

IntegrationRule *
MITC4Shell :: GetSurfaceIntegrationRule(int approxOrder)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, this, 1, 1);
    int npoints = iRule->getRequiredNumberOfIntegrationPoints(_Square, approxOrder);
    iRule->SetUpPointsOnSquare(npoints, _Unknown);
    return iRule;
}



double
MITC4Shell :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    double detJ, weight;
    FloatMatrix jacobianMatrix(2, 2);
    FloatMatrix dn;
    FloatArray lcoords(2);
    lcoords.at(1) = gp->giveNaturalCoordinate(1);
    lcoords.at(2) = gp->giveNaturalCoordinate(2);
    FloatArray x(4), y(4), z(4);
    this->giveNodeCoordinates( x.at(1), x.at(2), x.at(3), x.at(4), y.at(1), y.at(2), y.at(3), y.at(4), z.at(1), z.at(2), z.at(3), z.at(4) );


    weight = gp->giveWeight();

    interp_lin.giveDerivatives(dn, lcoords);

    for ( int i = 1; i <= dn.giveNumberOfRows(); i++ ) {
        double xl = x.at(i);
        double yl = y.at(i);

        jacobianMatrix.at(1, 1) += dn.at(i, 1) * xl;
        jacobianMatrix.at(1, 2) += dn.at(i, 1) * yl;
        jacobianMatrix.at(2, 1) += dn.at(i, 2) * xl;
        jacobianMatrix.at(2, 2) += dn.at(i, 2) * yl;
    }

    detJ = jacobianMatrix.giveDeterminant();
    return detJ * weight;
}
} // end namespace oofem
