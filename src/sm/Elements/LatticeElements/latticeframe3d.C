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

#include "domain.h"
#include "latticeframe3d.h"
#include "../sm/Materials/LatticeMaterials/latticematstatus.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatmatrixf.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "mathfem.h"
#include "latticestructuralelement.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"
#include "sm/CrossSections/latticecrosssection.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "../sm/Materials/structuralmaterial.h"
#endif

namespace oofem {
REGISTER_Element(LatticeFrame3d);

LatticeFrame3d::LatticeFrame3d(int n, Domain *aDomain) : LatticeStructuralElement(n, aDomain)
{
    numberOfDofMans = 2;
}

LatticeFrame3d::~LatticeFrame3d()
{}


void
LatticeFrame3d::computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the strain matrix of the receiver.
{
    //Assemble Bmatrix (used to compute strains and rotations)
    answer.resize(6, 12);
    answer.zero();

    this->length = computeLength();

    //Normal displacement jump in x-direction
    //First node
    answer.at(1, 1) = -1.;
    answer.at(1, 2) = 0.;
    answer.at(1, 3) = 0.;
    answer.at(1, 4) = 0.;
    answer.at(1, 5) = 0.;
    answer.at(1, 6) = 0.;
    //Second node
    answer.at(1, 7) = 1.;
    answer.at(1, 8) = 0.;
    answer.at(1, 9) = 0.;
    answer.at(1, 10) = 0.;
    answer.at(1, 11) = 0.;
    answer.at(1, 12) = 0.;

    //Shear displacement jump in y-plane
    //first node
    answer.at(2, 1) = 0.;
    answer.at(2, 2) = -1.;
    answer.at(2, 3) =  0.;
    answer.at(2, 4) = 0.;
    answer.at(2, 5) = 0;
    answer.at(2, 6) = -this->length*(1-this->s)/2.;
    //Second node
    answer.at(2, 7) = 0.;
    answer.at(2, 8) = 1.;
    answer.at(2, 9) =  0.;
    answer.at(2, 10) = 0.;
    answer.at(2, 11) = 0;
    answer.at(2, 12) = -this->length*(1+this->s)/2.;

    //Shear displacement jump in z-plane
    //first node
    answer.at(3, 1) = 0.;
    answer.at(3, 2) = 0.;
    answer.at(3, 3) = -1.;
    answer.at(3, 4) = 0.;
    answer.at(3, 5) = this->length*(1-this->s)/2.;
    answer.at(3, 6) = 0.;
    //Second node
    answer.at(3, 7) = 0.;
    answer.at(3, 8) = 0.;
    answer.at(3, 9) =  1.;
    answer.at(3, 10) = 0.;
    answer.at(3, 11) = this->length*(1+this->s)/2.;
    answer.at(3, 12) = 0.;

    //Rotation around x-axis
    //First node
    answer.at(4, 1) = 0.;
    answer.at(4, 2) = 0;
    answer.at(4, 3) = 0.;
    answer.at(4, 4) = -1.;
    answer.at(4, 5) = 0.;
    answer.at(4, 6) = 0.;
    //Second node
    answer.at(4, 7) = 0.;
    answer.at(4, 8) = 0.;
    answer.at(4, 9) = 0.;
    answer.at(4, 10) = 1.;
    answer.at(4, 11) = 0.;
    answer.at(4, 12) = 0.;

    //Rotation around y-axis
    //First node
    answer.at(5, 1) = 0.;
    answer.at(5, 2) = 0.;
    answer.at(5, 3) = 0.;
    answer.at(5, 4) = 0.;
    answer.at(5, 5) = -1.;
    answer.at(5, 6) = 0.;
    //Second node
    answer.at(5, 7) = 0.;
    answer.at(5, 8) = 0.;
    answer.at(5, 9) =  0.;
    answer.at(5, 10) = 0.;
    answer.at(5, 11) = 1.;
    answer.at(5, 12) = 0.;

    //Rotation around z-axis
    //First node
    answer.at(6, 1) = 0.;
    answer.at(6, 2) = 0.;
    answer.at(6, 3) = 0.;
    answer.at(6, 4) = 0.;
    answer.at(6, 5) = 0.;
    answer.at(6, 6) = -1.;
    //Second node
    answer.at(6, 7) = 0.;
    answer.at(6, 8) = 0.;
    answer.at(6, 9) =  0.;
    answer.at(6, 10) = 0.;
    answer.at(6, 11) = 0.;
    answer.at(6, 12) = 1.;

    return;
}

void
LatticeFrame3d::giveGPCoordinates(FloatArray &coords)
{
    coords.resize(3);
    coords = this->globalCentroid;
    return;
}

double
LatticeFrame3d::giveLength()
{
    return this->length;
}



void
LatticeFrame3d::computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer =  static_cast< LatticeCrossSection * >( this->giveCrossSection() )->give3dFrameStiffnessMatrix(rMode, gp, tStep);
}

void
LatticeFrame3d::computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    answer = static_cast< LatticeCrossSection * >( this->giveCrossSection() )->giveFrameForces3d(strain, gp, tStep);
}

int
LatticeFrame3d::computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    double ksi, n1, n2;

    ksi = lcoords.at(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    answer.resize(3);
    answer.at(1) = n1 * this->giveNode(1)->giveCoordinate(1) + n2 * this->giveNode(2)->giveCoordinate(1);
    answer.at(2) = n1 * this->giveNode(1)->giveCoordinate(2) + n2 * this->giveNode(2)->giveCoordinate(2);
    answer.at(3) = n1 * this->giveNode(1)->giveCoordinate(3) + n2 * this->giveNode(2)->giveCoordinate(3);

    return 1;
}


void
LatticeFrame3d::computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                       TimeStep *tStep)
// Computes numerically the stiffness matrix of the receiver.
{
    FloatMatrix d, bi, bj, bjt, dbj, dij;

    this->length = computeLength();

    answer.resize(12, 12);
    answer.zero();
    this->computeBmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), bj);
    this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);

    dbj.beProductOf(d, bj);
    dbj.times(1. / length);
    bjt.beTranspositionOf(bj);
    answer.beProductOf(bjt, dbj);

    return;
}

void LatticeFrame3d::computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    integrationRulesArray.resize(1);
    integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
    integrationRulesArray [ 0 ]->SetUpPointsOnLine(1, _3dLattice);
}



double LatticeFrame3d::giveArea() {
  FloatArray lc(1);
  return this->giveCrossSection()->give(CS_Area, lc, this);
}

double LatticeFrame3d::giveIy() {
  FloatArray lc(1);
  return this->giveCrossSection()->give(CS_InertiaMomentY, lc, this);
}

double LatticeFrame3d::giveIz() {
  FloatArray lc(1);
  return this->giveCrossSection()->give(CS_InertiaMomentZ, lc, this);
}

double LatticeFrame3d::giveIk() {
  FloatArray lc(1);
  return this->giveCrossSection()->give(CS_TorsionConstantX, lc, this);
}

double LatticeFrame3d::giveShearAreaY() {
  FloatArray lc(1);
  return this->giveCrossSection()->give(CS_ShearAreaY, lc, this);
}

double LatticeFrame3d::giveShearAreaZ() {
  FloatArray lc(1);
  return this->giveCrossSection()->give(CS_ShearAreaZ, lc, this);
}


void
LatticeFrame3d::giveInternalForcesVector(FloatArray &answer,
                                         TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatMatrix b, bt;
    FloatArray u, stress, strain;

    this->length = computeLength();
    
    this->computeVectorOf(VM_Total, tStep, u);

    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.clear();

    this->computeBmatrixAt(integrationRulesArray [ 0 ]->getIntegrationPoint(0), b);
    bt.beTranspositionOf(b);

    if ( useUpdatedGpRecord == 1 ) {
        LatticeMaterialStatus *lmatStat = dynamic_cast< LatticeMaterialStatus * >( integrationRulesArray [ 0 ]->getIntegrationPoint(0)->giveMaterialStatus() );
        stress = lmatStat->giveLatticeStress();
    } else {
        if ( !this->isActivated(tStep) ) {
            strain.zero();
        }
        strain.beProductOf(b, u);
	strain.times(1./this->length);
        this->computeStressVector(stress, strain, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);
    }

    answer.beProductOf(bt, stress);

    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}



bool
LatticeFrame3d::computeGtoLRotationMatrix(FloatMatrix &answer)
{
    FloatMatrix lcs;
    answer.resize(12, 12);
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

    return 1;
}


int
LatticeFrame3d::giveLocalCoordinateSystem(FloatMatrix &answer)
{
    FloatArray lx, ly, lz, help(3);
    Node *nodeA, *nodeB;
    nodeA = this->giveNode(1);
    nodeB = this->giveNode(2);

    lx.beDifferenceOf(nodeB->giveCoordinates(), nodeA->giveCoordinates() );
    lx.normalize();

    if ( this->referenceNode ) {
        Node *refNode = this->giveDomain()->giveNode(this->referenceNode);
        help.beDifferenceOf(refNode->giveCoordinates(), nodeA->giveCoordinates() );

        lz.beVectorProductOf(lx, help);
        lz.normalize();
    } else if ( this->zaxis.giveSize() > 0 ) {
        lz = this->zaxis;
        lz.add(lz.dotProduct(lx), lx);
        lz.normalize();
    } else {
        FloatMatrix rot(3, 3);
        double theta = referenceAngle * M_PI / 180.0;

        rot.at(1, 1) = cos(theta) + pow(lx.at(1), 2) * ( 1 - cos(theta) );
        rot.at(1, 2) = lx.at(1) * lx.at(2) * ( 1 - cos(theta) ) - lx.at(3) * sin(theta);
        rot.at(1, 3) = lx.at(1) * lx.at(3) * ( 1 - cos(theta) ) + lx.at(2) * sin(theta);

        rot.at(2, 1) = lx.at(2) * lx.at(1) * ( 1 - cos(theta) ) + lx.at(3) * sin(theta);
        rot.at(2, 2) = cos(theta) + pow(lx.at(2), 2) * ( 1 - cos(theta) );
        rot.at(2, 3) = lx.at(2) * lx.at(3) * ( 1 - cos(theta) ) - lx.at(1) * sin(theta);

        rot.at(3, 1) = lx.at(3) * lx.at(1) * ( 1 - cos(theta) ) - lx.at(2) * sin(theta);
        rot.at(3, 2) = lx.at(3) * lx.at(2) * ( 1 - cos(theta) ) + lx.at(1) * sin(theta);
        rot.at(3, 3) = cos(theta) + pow(lx.at(3), 2) * ( 1 - cos(theta) );

        help.at(3) = 1.0;         // up-vector
        // here is ly is used as a temp var
        if ( fabs(lx.dotProduct(help) ) > 0.999 ) {  // Check if it is vertical
            ly = {
                0., 1., 0.
            };
        } else {
            ly.beVectorProductOf(lx, help);
        }
        lz.beProductOf(rot, ly);
        lz.normalize();
    }

    ly.beVectorProductOf(lz, lx);
    ly.normalize();

    answer.resize(3, 3);
    answer.zero();
    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, i) = lx.at(i);
        answer.at(2, i) = ly.at(i);
        answer.at(3, i) = lz.at(i);
    }

    return 1;
}


void
LatticeFrame3d::giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v, D_w, R_u, R_v, R_w
    };
}

void
LatticeFrame3d::initializeFrom(InputRecord &ir)
{
    LatticeStructuralElement::initializeFrom(ir);

    referenceNode = 0;
    referenceAngle = 0;
    this->zaxis.clear();
    if ( ir.hasField(_IFT_LatticeFrame3d_zaxis) ) {
        IR_GIVE_FIELD(ir, this->zaxis, _IFT_LatticeFrame3d_zaxis);
    } else if ( ir.hasField(_IFT_LatticeFrame3d_refnode) ) {
        IR_GIVE_FIELD(ir, referenceNode, _IFT_LatticeFrame3d_refnode);
        if ( referenceNode == 0 ) {
            OOFEM_WARNING("wrong reference node specified. Using default orientation.");
        }
    } else if ( ir.hasField(_IFT_LatticeFrame3d_refangle) ) {
        IR_GIVE_FIELD(ir, referenceAngle, _IFT_LatticeFrame3d_refangle);
    } else {
        throw ValueInputException(ir, _IFT_LatticeFrame3d_zaxis, "axis, reference node, or angle not set");
    }

    this->s = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, s, _IFT_LatticeFrame3d_s);
}


double
LatticeFrame3d::computeLength()
{
    double dx, dy, dz;
    Node *nodeA, *nodeB;

    if ( length == 0. ) {
        nodeA   = this->giveNode(1);
        nodeB   = this->giveNode(2);
        dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
        dy      = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
        dz      = nodeB->giveCoordinate(3) - nodeA->giveCoordinate(3);
        length  = sqrt(dx * dx + dy * dy + dz * dz);
    }

    return length;
}



void
LatticeFrame3d::computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver. This expression is
// valid in both local and global axes.
{
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    double density = static_cast< LatticeCrossSection * >( this->giveCrossSection() )->give('d', gp);
    double halfMass = density * computeVolumeAround(gp) / 2.;
    answer.resize(12, 12);
    answer.zero();
    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = halfMass;
    answer.at(7, 7) = answer.at(8, 8) = answer.at(9, 9) = halfMass;
}
} // end namespace oofem
