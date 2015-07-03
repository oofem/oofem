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

#include "../sm/Elements/Shells/quad1mindlinshell3d.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "../sm/Loads/constantpressureload.h"
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
#include "fei2dquadlin.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(Quad1MindlinShell3D);

FEI2dQuadLin Quad1MindlinShell3D :: interp(1, 2);
IntArray Quad1MindlinShell3D :: shellOrdering = { 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23};
IntArray Quad1MindlinShell3D :: drillOrdering = { 6, 12, 18, 24};

Quad1MindlinShell3D :: Quad1MindlinShell3D(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(this),
    SPRNodalRecoveryModelInterface(),
    lnodes(4)
{
    numberOfGaussPoints = 4;
    this->numberOfDofMans = 4;
    this->reducedIntegrationFlag = false;
}


Quad1MindlinShell3D :: ~Quad1MindlinShell3D()
{
}


FEInterpolation *
Quad1MindlinShell3D :: giveInterpolation() const
{
    return & interp;
}


FEInterpolation *
Quad1MindlinShell3D :: giveInterpolation(DofIDItem id) const
{
    return & interp;
}


void
Quad1MindlinShell3D :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 5) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
    ///@todo Deal with updated geometries and such.
    this->computeLCS();
}


void
Quad1MindlinShell3D :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
{
    // Only gravity load
    double dV, density;
    FloatArray forceX, forceY, forceZ, glob_gravity, gravity, n;

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        OOFEM_ERROR("unknown load type");
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(glob_gravity, tStep, mode);
    // Transform the load into the local c.s.
    gravity.beProductOf(this->lcsMatrix, glob_gravity); ///@todo Check potential transpose here.

    if ( gravity.giveSize() ) {
        for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {

            this->interp.evalN( n, gp->giveNaturalCoordinates(), FEIVoidCellGeometry() );
            dV = this->computeVolumeAround(gp) * this->giveCrossSection()->give(CS_Thickness, gp);
            density = this->giveStructuralCrossSection()->give('d', gp);

            forceX.add(density * gravity.at(1) * dV, n);
            forceY.add(density * gravity.at(2) * dV, n);
            forceZ.add(density * gravity.at(3) * dV, n);
        }

        answer.resize(24);
        answer.zero();

        answer.at(1)  = forceX.at(1);
        answer.at(2)  = forceY.at(1);
        answer.at(3)  = forceZ.at(1);

        answer.at(7)  = forceX.at(2);
        answer.at(8)  = forceY.at(2);
        answer.at(9)  = forceZ.at(2);

        answer.at(13) = forceX.at(3);
        answer.at(14) = forceY.at(3);
        answer.at(15) = forceZ.at(3);

        answer.at(19) = forceX.at(4);
        answer.at(20) = forceY.at(4);
        answer.at(21) = forceZ.at(4);
    } else {
        answer.clear();
    }
}


void
Quad1MindlinShell3D :: computeSurfaceLoadVectorAt(FloatArray &answer, Load *load,
                                                  int iSurf, TimeStep *tStep, ValueModeType mode)
{
    BoundaryLoad *surfLoad = static_cast< BoundaryLoad * >(load);
    if ( dynamic_cast< ConstantPressureLoad * >(surfLoad) ) { // Just checking the type of b.c.
        // EXPERIMENTAL CODE:
        FloatArray n, gcoords, pressure;

        answer.resize(24);
        answer.zero();

        for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
            double dV = this->computeVolumeAround(gp);
            this->interp.evalN( n, gp->giveNaturalCoordinates(), FEIVoidCellGeometry() );
            this->interp.local2global( gcoords, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
            surfLoad->computeValueAt(pressure, tStep, gcoords, mode);

            answer.at(3) += n.at(1) * pressure.at(1) * dV;
            answer.at(9) += n.at(2) * pressure.at(1) * dV;
            answer.at(15) += n.at(3) * pressure.at(1) * dV;
            answer.at(21) += n.at(4) * pressure.at(1) * dV;
        }
        // Second surface is the outside;
        if ( iSurf == 2 ) {
            answer.negated();
        }
    } else {
        OOFEM_ERROR("only supports constant pressure boundary load.");
    }
}


void
Quad1MindlinShell3D :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
    FloatArray n, ns;
    FloatMatrix dn, dns;
    const FloatArray &localCoords = gp->giveNaturalCoordinates();

    this->interp.evaldNdx( dn, localCoords, FEIVertexListGeometryWrapper(lnodes) );
    this->interp.evalN( n, localCoords,  FEIVoidCellGeometry() );

    answer.resize(8, 4 * 5);
    answer.zero();

    // enforce one-point reduced integration if requested
    if ( this->reducedIntegrationFlag ) {
        FloatArray lc(2);
        lc.zero(); // set to element center coordinates

        this->interp.evaldNdx( dns, lc, FEIVertexListGeometryWrapper(lnodes) );
        this->interp.evalN( ns, lc,  FEIVoidCellGeometry() );
    } else {
        dns = dn;
        ns = n;
    }


    // Note: This is just 5 dofs (sixth column is all zero, torsional stiffness handled separately.)
    for ( int i = 0; i < 4; ++i ) {
        ///@todo Check the rows for both parts here, to be consistent with _3dShell material definition
        // Part related to the membrane (columns represent coefficients for D_u, D_v)
        answer(0, 0 + i * 5) = dn(i, 0);//eps_x = du/dx
        answer(1, 1 + i * 5) = dn(i, 1);//eps_y = dv/dy
        answer(2, 0 + i * 5) = dn(i, 1);//gamma_xy = du/dy+dv/dx
        answer(2, 1 + i * 5) = dn(i, 0);

        // Part related to the plate (columns represent the dofs D_w, R_u, R_v)
        ///@todo Check sign here
        answer(3 + 0, 2 + 2 + i * 5) = dn(i, 0);// kappa_x = d(fi_y)/dx
        answer(3 + 1, 2 + 1 + i * 5) =-dn(i, 1);// kappa_y = -d(fi_x)/dy
        answer(3 + 2, 2 + 2 + i * 5) = dn(i, 1);// kappa_xy=d(fi_y)/dy-d(fi_x)/dx
        answer(3 + 2, 2 + 1 + i * 5) =-dn(i, 0);

        // shear strains
        answer(3 + 3, 2 + 0 + i * 5) = dns(i, 0);// gamma_xz = fi_y+dw/dx
        answer(3 + 3, 2 + 2 + i * 5) = ns(i);
        answer(3 + 4, 2 + 0 + i * 5) = dns(i, 1);// gamma_yz = -fi_x+dw/dy
        answer(3 + 4, 2 + 1 + i * 5) = -ns(i);
    }


#if 0
    // Experimental MITC4 support.
    // Based on "Short communication A four-node plate bending element based on mindling/reissner plate theory and a mixed interpolation"
    // KJ Bathe, E Dvorkin

    double x1, x2, x3, x4;
    double y1, y2, y3, y4;
    double Ax, Bx, Cx, Ay, By, Cy;

    double r = localCoords[0];
    double s = localCoords[1];

    x1 = lnodes[0][0];
    x2 = lnodes[1][0];
    x3 = lnodes[2][0];
    x4 = lnodes[3][0];

    y1 = lnodes[0][1];
    y2 = lnodes[1][1];
    y3 = lnodes[2][1];
    y4 = lnodes[3][1];

    Ax = x1 - x2 - x3 + x4;
    Bx = x1 - x2 + x3 - x4;
    Cx = x1 + x2 - x3 - x4;

    Ay = y1 - y2 - y3 + y4;
    By = y1 - y2 + y3 - y4;
    Cy = y1 + y2 - y3 - y4;

    FloatMatrix jac;
    this->interp.giveJacobianMatrixAt(jac, localCoords, FEIVertexListGeometryWrapper(lnodes) );
    double detJ = jac.giveDeterminant();

    double rz = sqrt( sqr(Cx + r*Bx) + sqr(Cy + r*By)) / ( 16 * detJ );
    double sz = sqrt( sqr(Ax + s*Bx) + sqr(Ay + s*By)) / ( 16 * detJ );

    // TODO: Not sure about this part (the reference is not explicit about these angles. / Mikael
    // Not sure about the transpose either.
    OOFEM_WARNING("The MITC4 implementation isn't verified yet. Highly experimental");
    FloatArray dxdr = {jac(0,0), jac(0,1)};
    dxdr.normalize();
    FloatArray dxds = {jac(1,0), jac(1,1)};
    dxds.normalize();

    double c_b = dxdr(0); //cos(beta);
    double s_b = dxdr(1); //sin(beta);
    double c_a = dxds(0); //cos(alpha);
    double s_a = dxds(1); //sin(alpha);

    // gamma_xz = "fi_y+dw/dx" in standard formulation
    answer(6, 2 + 5*0) = rz * s_b * ( (1+s)) - sz * s_a * ( (1+r));
    answer(6, 2 + 5*1) = rz * s_b * (-(1+s)) - sz * s_a * ( (1-r));
    answer(6, 2 + 5*2) = rz * s_b * (-(1-s)) - sz * s_a * (-(1-r));
    answer(6, 2 + 5*3) = rz * s_b * ( (1-s)) - sz * s_a * (-(1+r));

    answer(6, 3 + 5*0) = rz * s_b * (y2-y1) * 0.5 * (1+s) - sz * s_a * (y4-y1) * 0.5 * (1+r); // tx1
    answer(6, 4 + 5*0) = rz * s_b * (x1-x2) * 0.5 * (1+s) - sz * s_a * (x1-x4) * 0.5 * (1+r); // ty1

    answer(6, 3 + 5*1) = rz * s_b * (y2-y1) * 0.5 * (1+s) - sz * s_a * (y3-x2) * 0.5 * (1+r); // tx2
    answer(6, 4 + 5*1) = rz * s_b * (x1-x2) * 0.5 * (1+s) - sz * s_a * (x2-x3) * 0.5 * (1+r); // ty2

    answer(6, 3 + 5*2) = rz * s_b * (y3-y4) * 0.5 * (1-s) - sz * s_a * (y3-y2) * 0.5 * (1-r); // tx3
    answer(6, 4 + 5*2) = rz * s_b * (x4-x3) * 0.5 * (1-s) - sz * s_a * (x2-x3) * 0.5 * (1-r); // ty3

    answer(6, 3 + 5*3) = rz * s_b * (y3-y4) * 0.5 * (1-s) - sz * s_a * (y4-y1) * 0.5 * (1-r); // tx4
    answer(6, 4 + 5*3) = rz * s_b * (x4-x3) * 0.5 * (1-s) - sz * s_a * (x1-x4) * 0.5 * (1-r); // ty4

    // gamma_yz = -fi_x+dw/dy in standard formulation
    answer(7, 2 + 5*0) = - rz * c_b * ( (1+s)) + sz * c_a * ( (1+r));
    answer(7, 2 + 5*1) = - rz * c_b * (-(1+s)) + sz * c_a * ( (1-r));
    answer(7, 2 + 5*2) = - rz * c_b * (-(1-s)) + sz * c_a * (-(1-r));
    answer(7, 2 + 5*3) = - rz * c_b * ( (1-s)) + sz * c_a * (-(1+r));

    answer(7, 3 + 5*0) = - rz * c_b * (y2-y1) * 0.5 * (1+s) + sz * c_a * (y4-y1) * 0.5 * (1+r); // tx1
    answer(7, 4 + 5*0) = - rz * c_b * (x1-x2) * 0.5 * (1+s) + sz * c_a * (x1-x4) * 0.5 * (1+r); // ty1

    answer(7, 3 + 5*1) = - rz * c_b * (y2-y1) * 0.5 * (1+s) + sz * c_a * (y3-x2) * 0.5 * (1+r); // tx2
    answer(7, 4 + 5*1) = - rz * c_b * (x1-x2) * 0.5 * (1+s) + sz * c_a * (x2-x3) * 0.5 * (1+r); // ty2

    answer(7, 3 + 5*2) = - rz * c_b * (y3-y4) * 0.5 * (1-s) + sz * c_a * (y3-y2) * 0.5 * (1-r); // tx3
    answer(7, 4 + 5*2) = - rz * c_b * (x4-x3) * 0.5 * (1-s) + sz * c_a * (x2-x3) * 0.5 * (1-r); // ty3

    answer(7, 3 + 5*3) = - rz * c_b * (y3-y4) * 0.5 * (1-s) + sz * c_a * (y4-y1) * 0.5 * (1-r); // tx4
    answer(7, 4 + 5*3) = - rz * c_b * (x4-x3) * 0.5 * (1-s) + sz * c_a * (x1-x4) * 0.5 * (1-r); // ty4
#endif
}


void
Quad1MindlinShell3D :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveGeneralizedStress_Shell(answer, gp, strain, tStep);
}


void
Quad1MindlinShell3D :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->give3dShellStiffMtrx(answer, rMode, gp, tStep);
}


void
Quad1MindlinShell3D :: computeVectorOfUnknowns(ValueModeType mode, TimeStep *tStep, FloatArray &shell, FloatArray &drill)
{
    FloatArray tmp;
    this->computeVectorOf(mode, tStep, tmp);
    shell.beSubArrayOf(tmp, this->shellOrdering);
    drill.beSubArrayOf(tmp, this->drillOrdering);
}


void
Quad1MindlinShell3D :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray shellUnknowns, tmp;
    FloatMatrix b;
    /* Here we do compute only the "traditional" part of shell strain vector, the quasi-strain related to rotations is not computed */
    this->computeVectorOfUnknowns(VM_Total, tStep, shellUnknowns, tmp);

    this->computeBmatrixAt(gp, b);
    answer.beProductOf(b, shellUnknowns);
}


void
Quad1MindlinShell3D :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // We need to overload this for practical reasons (this 3d shell has all 9 dofs, but the shell part only cares for the first 8)
    // This elements adds an additional stiffness for the so called drilling dofs, meaning we need to work with all 9 components.
    FloatMatrix b;
    FloatArray n, strain, stress;
    FloatArray shellUnknowns, drillUnknowns;
    bool drillCoeffFlag = false;

    // Split this for practical reasons into normal shell dofs and drilling dofs
    this->computeVectorOfUnknowns(VM_Total, tStep, shellUnknowns, drillUnknowns);

    FloatArray shellForces, drillMoment;
    StructuralCrossSection *cs = this->giveStructuralCrossSection();

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        this->computeBmatrixAt(gp, b);
        double dV = this->computeVolumeAround(gp);
        double drillCoeff = cs->give(CS_DrillingStiffness, gp);

        if ( useUpdatedGpRecord ) {
            stress = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStressVector();
        } else {
            strain.beProductOf(b, shellUnknowns);
            cs->giveGeneralizedStress_Shell(stress, gp, strain, tStep);
        }
        shellForces.plusProduct(b, stress, dV);

        // Drilling stiffness is here for improved numerical properties
        if ( drillCoeff > 0. ) {
            this->interp.evalN( n, gp->giveNaturalCoordinates(), FEIVoidCellGeometry() );
            for ( int j = 0; j < 4; j++ ) {
                n(j) -= 0.25;
            }
            double dtheta = n.dotProduct(drillUnknowns);
            drillMoment.add(drillCoeff * dV * dtheta, n); ///@todo Decide on how to alpha should be defined.
            drillCoeffFlag = true;
        }
    }

    answer.resize(24);
    answer.zero();
    answer.assemble(shellForces, this->shellOrdering);

    if ( drillCoeffFlag ) {
        answer.assemble(drillMoment, this->drillOrdering);
    }
}


void
Quad1MindlinShell3D :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    // We need to overload this for practical reasons (this 3d shell has all 9 dofs, but the shell part only cares for the first 8)
    // This elements adds an additional stiffness for the so called drilling dofs, meaning we need to work with all 9 components.
    FloatMatrix d, b, db;
    FloatArray n;
    bool drillCoeffFlag = false;

    FloatMatrix shellStiffness, drillStiffness;

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        this->computeBmatrixAt(gp, b);
        double dV = this->computeVolumeAround(gp);
        double drillCoeff = this->giveStructuralCrossSection()->give(CS_DrillingStiffness, gp);

        this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);

        db.beProductOf(d, b);
        shellStiffness.plusProductSymmUpper(b, db, dV);

        // Drilling stiffness is here for improved numerical properties
        if ( drillCoeff > 0. ) {
            this->interp.evalN( n, gp->giveNaturalCoordinates(), FEIVoidCellGeometry() );
            for ( int j = 0; j < 4; j++ ) {
                n(j) -= 0.25;
            }
            drillStiffness.plusDyadSymmUpper(n, drillCoeff * dV);
            drillCoeffFlag = true;
        }
    }
    shellStiffness.symmetrized();

    answer.resize(24, 24);
    answer.zero();
    answer.assemble(shellStiffness, this->shellOrdering);

    if ( drillCoeffFlag ) {
        drillStiffness.symmetrized();
        answer.assemble(drillStiffness, this->drillOrdering);
    }
}


IRResultType
Quad1MindlinShell3D :: initializeFrom(InputRecord *ir)
{
    this->reducedIntegrationFlag = ir->hasField(_IFT_Quad1MindlinShell3D_ReducedIntegration);
    return NLStructuralElement :: initializeFrom(ir);
}


void
Quad1MindlinShell3D :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v, D_w, R_u, R_v, R_w};
}


void
Quad1MindlinShell3D :: computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp)
{
    FloatArray u, v;
    u.beDifferenceOf( * this->giveNode(2)->giveCoordinates(), * this->giveNode(1)->giveCoordinates() );
    v.beDifferenceOf( * this->giveNode(3)->giveCoordinates(), * this->giveNode(1)->giveCoordinates() );

    answer.beVectorProductOf(u, v);
    answer.normalize();
}


double
Quad1MindlinShell3D :: giveCharacteristicLength(const FloatArray &normalToCrackPlane)
{
    return this->giveLengthInDir(normalToCrackPlane);
}


double
Quad1MindlinShell3D :: computeVolumeAround(GaussPoint *gp)
{
    double detJ, weight;

    weight = gp->giveWeight();
    detJ = fabs( this->interp.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIVertexListGeometryWrapper(lnodes) ) );
    return detJ * weight;
}


void
Quad1MindlinShell3D :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver.
{
    double mass = 0.;

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        mass += this->computeVolumeAround(gp) * this->giveStructuralCrossSection()->give('d', gp);
    }

    answer.resize(24, 24);
    answer.zero();
    for ( int i = 0; i < 4; i++ ) {
        answer(i * 6 + 0, i * 6 + 0) = mass * 0.25;
        answer(i * 6 + 1, i * 6 + 1) = mass * 0.25;
        answer(i * 6 + 2, i * 6 + 2) = mass * 0.25;
    }
}


int
Quad1MindlinShell3D :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FloatArray help;
    answer.resize(6);
    if ( type == IST_ShellForceTensor || type == IST_ShellMomentumTensor ) {
        if ( type == IST_ShellForceTensor ) {
            help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStressVector();
        } else {
            help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStrainVector();
        }
        answer.at(1) = help.at(1); // nx
        answer.at(2) = help.at(2); // ny
        answer.at(3) = 0.0; // nz
        answer.at(4) = help.at(8); // vyz
        answer.at(5) = help.at(7); // vxy
        answer.at(6) = help.at(3); // vxy
        return 1;
    } else if ( type == IST_ShellMomentumTensor || type == IST_ShellCurvatureTensor ) {
        if ( type == IST_ShellMomentumTensor ) {
            help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStressVector();
        } else {
            help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStrainVector();
        }
        answer.at(1) = help.at(4); // mx
        answer.at(2) = help.at(5); // my
        answer.at(3) = 0.0;      // mz
        answer.at(4) = 0.0;      // mzy
        answer.at(5) = 0.0;      // mzx
        answer.at(6) = help.at(6); // mxy
        return 1;
    } else {
        return NLStructuralElement :: giveIPValue(answer, gp, type, tStep);
    }
}


void
Quad1MindlinShell3D :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
{
    IntArray edgeNodes;
    FloatArray n;

    this->interp.edgeEvalN( n, iedge, gp->giveNaturalCoordinates(), FEIVoidCellGeometry() );
    this->interp.computeLocalEdgeMapping(edgeNodes, iedge);

    answer.beNMatrixOf(n, 6);
}


void
Quad1MindlinShell3D :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    if ( iEdge == 1 ) { // edge between nodes 1 2
        answer = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    } else if ( iEdge == 2 ) { // edge between nodes 2 3
        answer = { 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
    } else if ( iEdge == 3 ) { // edge between nodes 3 4
        answer = {13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
    } else if ( iEdge == 4 ) { // edge between nodes 4 1
        answer = {19, 20, 21, 22, 23, 24, 1, 2, 3, 4, 5, 6};
    } else {
        OOFEM_ERROR("wrong edge number");
    }
}


double
Quad1MindlinShell3D :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double detJ = this->interp.edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(), FEIVertexListGeometryWrapper(lnodes) );
    return detJ *gp->giveWeight();
}


void
Quad1MindlinShell3D :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    FloatArray local;
    this->interp.edgeLocal2global( local, iEdge, gp->giveNaturalCoordinates(), FEIVertexListGeometryWrapper(lnodes)  );
    local.resize(3);
    local.at(3) = 0.;
    answer.beProductOf(this->lcsMatrix, local);
}


int
Quad1MindlinShell3D :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    double dx, dy, length;
    IntArray edgeNodes;
    Node *nodeA, *nodeB;

    answer.resize(3, 3);
    answer.zero();

    this->interp.computeLocalEdgeMapping(edgeNodes, iEdge);

    nodeA = this->giveNode( edgeNodes.at(1) );
    nodeB = this->giveNode( edgeNodes.at(2) );

    dx = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
    dy = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
    length = sqrt(dx * dx + dy * dy);

    /// @todo I haven't even looked at this code yet / Mikael
    answer.at(1, 1) = 1.0;
    answer.at(2, 2) = dx / length;
    answer.at(2, 3) = -dy / length;
    answer.at(3, 2) = -answer.at(2, 3);
    answer.at(3, 3) = answer.at(2, 2);

    return 1;
}


void
Quad1MindlinShell3D :: computeLCS()
{
    lcsMatrix.resize(3, 3); // Note! G -> L transformation matrix
    FloatArray e1, e2, e3, help;

    // compute e1' = [N2-N1]  and  help = [N4-N1]
    e1.beDifferenceOf( * this->giveNode(2)->giveCoordinates(), * this->giveNode(1)->giveCoordinates() );
    help.beDifferenceOf( * this->giveNode(4)->giveCoordinates(), * this->giveNode(1)->giveCoordinates() );
    e1.normalize();
    e3.beVectorProductOf(e1, help);
    e3.normalize();
    e2.beVectorProductOf(e3, e1);
    for ( int i = 1; i <= 3; i++ ) {
        this->lcsMatrix.at(1, i) = e1.at(i);
        this->lcsMatrix.at(2, i) = e2.at(i);
        this->lcsMatrix.at(3, i) = e3.at(i);
    }

    for ( int i = 1; i <= 4; i++ ) {
        this->lnodes [ i - 1 ].beProductOf( this->lcsMatrix, * this->giveNode(i)->giveCoordinates() );
    }
}


bool
Quad1MindlinShell3D :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    answer.resize(24, 24);
    answer.zero();

    for ( int i = 0; i < 4; i++ ) { // Loops over nodes
        // In each node, transform global c.s. {D_u, D_v, D_w, R_u, R_v, R_w} into local c.s.
        answer.setSubMatrix(this->lcsMatrix, 1 + i * 6, 1 + i * 6);     // Displacements
        answer.setSubMatrix(this->lcsMatrix, 1 + i * 6 + 3, 1 + i * 6 + 3); // Rotations
    }

    return true;
}

Interface *
Quad1MindlinShell3D :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    }

    return NULL;
}



void
Quad1MindlinShell3D :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(4);
    for ( int i = 1; i < 5; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
Quad1MindlinShell3D :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
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
} // end namespace oofem
