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
 *               Copyright (C) 1993 - 2026   Borek Patzak
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

#include "latticelink3dnl.h"
#include "../sm/Materials/LatticeMaterials/latticematstatus.h"
#include "node.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(LatticeLink3dNL);

// Build the 6x6 block-diagonal rotation from a 3x3 direction-cosine matrix.
static void
expandTo6(const FloatMatrix &Q, FloatMatrix &Q6)
{
    Q6.resize(6, 6);
    Q6.zero();
    for ( int i = 1; i <= 3; ++i ) {
        for ( int j = 1; j <= 3; ++j ) {
            Q6.at(i, j) = Q.at(i, j);
            Q6.at(i + 3, j + 3) = Q.at(i, j);
        }
    }
}


void
LatticeLink3dNL::computeCurrentFrame(FloatMatrix &Qcur, const FloatArray &spinB)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    FloatMatrix R;
    LatticeStructuralElement::computeGlobalRotationMatrix(R, spinB);

    // Rotate each reference local axis (rows of localCoordinateSystem, global coords) by spin_B.
    Qcur.resize(3, 3);
    for ( int i = 1; i <= 3; ++i ) {
        FloatArray e0(3), ei;
        for ( int j = 1; j <= 3; ++j ) {
            e0.at(j) = this->localCoordinateSystem.at(i, j);
        }
        ei.beProductOf(R, e0);
        for ( int j = 1; j <= 3; ++j ) {
            Qcur.at(i, j) = ei.at(j);
        }
    }
}


void
LatticeLink3dNL::computeRotatedArm(FloatArray &cA, const FloatArray &spinA)
{
    FloatArray coordA = this->giveNode(1)->giveCoordinates();
    FloatArray coordB = this->giveNode(2)->giveCoordinates();
    FloatArray rigidGlobal = coordB;
    rigidGlobal.subtract(coordA);

    FloatMatrix RA;
    LatticeStructuralElement::computeGlobalRotationMatrix(RA, spinA);
    cA.beProductOf(RA, rigidGlobal);
}


void
LatticeLink3dNL::computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    FloatArray u;
    this->computeVectorOf(VM_Total, tStep, u);

    FloatArray coordA = this->giveNode(1)->giveCoordinates();
    FloatArray coordB = this->giveNode(2)->giveCoordinates();
    FloatArray rigidGlobal = coordB;
    rigidGlobal.subtract(coordA);

    // Node A rotation: rotate the rigid arm by spin_A.
    FloatArray spinA(3);
    spinA.at(1) = u.at(4); spinA.at(2) = u.at(5); spinA.at(3) = u.at(6);
    FloatMatrix RA;
    LatticeStructuralElement::computeGlobalRotationMatrix(RA, spinA);
    FloatArray cA;
    cA.beProductOf(RA, rigidGlobal);

    // Jump (global). Node B enters only via its translation (its arm is zero).
    answer.resize(6);
    answer.at(1) = u.at(7) - u.at(1) - cA.at(1) + rigidGlobal.at(1);
    answer.at(2) = u.at(8) - u.at(2) - cA.at(2) + rigidGlobal.at(2);
    answer.at(3) = u.at(9) - u.at(3) - cA.at(3) + rigidGlobal.at(3);

    // Rotational strain: exact relative rotation log(R_B R_A^T).
    FloatArray spinB(3);
    spinB.at(1) = u.at(10); spinB.at(2) = u.at(11); spinB.at(3) = u.at(12);
    FloatMatrix RB, RAt, Rrel;
    LatticeStructuralElement::computeGlobalRotationMatrix(RB, spinB);
    RAt.beTranspositionOf(RA);
    Rrel.beProductOf(RB, RAt);
    FloatArray theta;
    LatticeStructuralElement::logRotationVector(Rrel, theta);
    answer.at(4) = theta.at(1);
    answer.at(5) = theta.at(2);
    answer.at(6) = theta.at(3);

    //Node B rotation: project onto the current bond frame.
    FloatMatrix Qcur, Q6;
    this->computeCurrentFrame(Qcur, spinB);
    expandTo6(Qcur, Q6);

    // global -> current local
    answer.rotatedWith(Q6, 'n');
}


void
LatticeLink3dNL::computeNLBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    FloatArray u;
    this->computeVectorOf(VM_Total, tStep, u);
    FloatArray spinA(3), cA;
    spinA.at(1) = u.at(4); spinA.at(2) = u.at(5); spinA.at(3) = u.at(6);
    this->computeRotatedArm(cA, spinA);

    // Global geometric B matrix with the rotated node-A arm; node B has no arm.
    answer.resize(6, 12);
    answer.zero();
    answer.at(1, 1) = -1.;
    answer.at(1, 5) = -cA.at(3);
    answer.at(1, 6) =  cA.at(2);
    answer.at(1, 7) = 1.;
    answer.at(2, 2) = -1.;
    answer.at(2, 4) =  cA.at(3);
    answer.at(2, 6) = -cA.at(1);
    answer.at(2, 8) = 1.;
    answer.at(3, 3) = -1.;
    answer.at(3, 4) = -cA.at(2);
    answer.at(3, 5) =  cA.at(1);
    answer.at(3, 9) = 1.;
    answer.at(4, 4) = -1.;
    answer.at(4, 10) = 1.;
    answer.at(5, 5) = -1.;
    answer.at(5, 11) = 1.;
    answer.at(6, 6) = -1.;
    answer.at(6, 12) = 1.;
}


void
LatticeLink3dNL::giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatArray u, strain, stress, s;
    this->computeVectorOf(VM_Total, tStep, u);
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    FloatArray spinA(3), spinB(3), cA;
    spinA.at(1) = u.at(4);
    spinA.at(2) = u.at(5);
    spinA.at(3) = u.at(6);
    spinB.at(1) = u.at(10);
    spinB.at(2) = u.at(11);
    spinB.at(3) = u.at(12);
    this->computeRotatedArm(cA, spinA);

    FloatMatrix Qcur, Q6;
    this->computeCurrentFrame(Qcur, spinB);
    expandTo6(Qcur, Q6);

    answer.resize(12);
    answer.zero();

    for ( GaussPoint *gp: * this->giveDefaultIntegrationRulePtr() ) {
        LatticeMaterialStatus *matStat = static_cast< LatticeMaterialStatus * >( gp->giveMaterialStatus() );
        if ( useUpdatedGpRecord == 1 ) {
            stress = matStat->giveLatticeStress();
        } else {
            this->computeStrainVector(strain, gp, tStep);
            this->computeStressVector(stress, strain, gp, tStep);
        }
        if ( stress.giveSize() == 0 ) {
            break;
        }

        // Resultants: rotate to global, scale by bond area (as in the linear link).
        s = stress;
        s.rotatedWith(Q6, 't');                // local -> global
        double area = this->computeVolumeAround(gp) / this->giveLength();
        s.times(area);

        answer.at(1)  += -s.at(1);
        answer.at(2)  += -s.at(2);
        answer.at(3)  += -s.at(3);
        answer.at(4)  +=  s.at(2) * cA.at(3) - s.at(3) * cA.at(2) - s.at(4);
        answer.at(5)  += -s.at(1) * cA.at(3) + s.at(3) * cA.at(1) - s.at(5);
        answer.at(6)  +=  s.at(1) * cA.at(2) - s.at(2) * cA.at(1) - s.at(6);
        answer.at(7)  +=  s.at(1);
        answer.at(8)  +=  s.at(2);
        answer.at(9)  +=  s.at(3);
        answer.at(10) +=  s.at(4);              // node B: no arm (c_B = 0)
        answer.at(11) +=  s.at(5);
        answer.at(12) +=  s.at(6);
    }

    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}



void
LatticeLink3dNL::computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    FloatArray u;
    this->computeVectorOf(VM_Total, tStep, u);

    FloatArray spinB(3);
    spinB.at(1) = u.at(10);
    spinB.at(2) = u.at(11);
    spinB.at(3) = u.at(12);

    FloatMatrix Qcur, Q6;
    this->computeCurrentFrame(Qcur, spinB);
    expandTo6(Qcur, Q6);

    answer.resize(12, 12);
    answer.zero();

    for ( GaussPoint *gp: * this->giveDefaultIntegrationRulePtr() ) {
        FloatMatrix b, bt, d, rT, dR, rTDR, db, contrib;
        this->computeNLBmatrixAt(gp, b, tStep);
        this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);   // material tangent (local)

        // K = B^T ( Q^T D Q ) B * area  (geometry from deformed config; material as returned).
        rT.beTranspositionOf(Q6);
        dR.beProductOf(d, Q6);
        rTDR.beProductOf(rT, dR);
        db.beProductOf(rTDR, b);
        bt.beTranspositionOf(b);
        contrib.beProductOf(bt, db);

        double area = this->computeVolumeAround(gp) / this->giveLength();
        contrib.times(area);
        answer.add(contrib);
    }
}
} // end namespace oofem
