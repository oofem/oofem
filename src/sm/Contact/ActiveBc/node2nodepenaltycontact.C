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


#include "sm/Contact/ActiveBc/node2nodepenaltycontact.h"
#include "set.h"
#include "domain.h"
#include "node.h"
#include "floatmatrix.h"
#include "unknownnumberingscheme.h"
#include "sparsemtrx.h"
#include "classfactory.h"

namespace oofem {
REGISTER_BoundaryCondition(Node2NodePenaltyContact);


IRResultType
Node2NodePenaltyContact :: initializeFrom(InputRecord *ir)
{
    IRResultType result;
    IR_GIVE_FIELD(ir, this->penalty, _IFT_Node2NodePenaltyContact_penalty);
    this->useTangent = ir->hasField(_IFT_Node2NodePenaltyContact_useTangent);


    IR_GIVE_FIELD(ir, this->masterSet, _IFT_Node2NodePenaltyContact_masterSet);
    IR_GIVE_FIELD(ir, this->slaveSet, _IFT_Node2NodePenaltyContact_slaveSet);


    return ActiveBoundaryCondition :: initializeFrom(ir);
}


void
Node2NodePenaltyContact :: assemble(SparseMtrx &answer, TimeStep *tStep,
                                    CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale)
{
    if ( !this->useTangent || type != TangentStiffnessMatrix ) {
        return;
    }


    FloatMatrix K;
    IntArray loc, c_loc;

    IntArray dofIdArray = {
        D_u, D_v
    };

    if ( masterSet.giveSize() != slaveSet.giveSize() ) {
        OOFEM_ERROR("Number of master nodes does not match number of slave nodes")
    }

    for ( int pos = 1; pos <= masterSet.giveSize(); ++pos ) {
        Node *masterNode = this->giveDomain()->giveNode(masterSet.at(pos) );
        Node *slaveNode = this->giveDomain()->giveNode(slaveSet.at(pos) );

        masterNode->giveLocationArray(dofIdArray, loc, r_s);
        slaveNode->giveLocationArray(dofIdArray, c_loc, c_s);
        loc.followedBy(c_loc);

        this->computeTangentFromContact(K, masterNode, slaveNode, tStep);
        answer.assemble(loc, K);
    }
}

void
Node2NodePenaltyContact :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                          CharType type, ValueModeType mode,
                                          const UnknownNumberingScheme &s, FloatArray *eNorms)
{
    if ( type != ExternalForcesVector ) {
        return;
    }

    //IntArray dofIdArray = {D_u, D_v, D_w};
    IntArray dofIdArray = {
        D_u, D_v
    };

    IntArray loc, c_loc;
    FloatArray fext;

    if ( masterSet.giveSize() != slaveSet.giveSize() ) {
        OOFEM_ERROR("Number of master nodes does not match number of slave nodes")
    }

    for ( int pos = 1; pos <= masterSet.giveSize(); ++pos ) {
        Node *masterNode = this->giveDomain()->giveNode(masterSet.at(pos) );
        Node *slaveNode = this->giveDomain()->giveNode(slaveSet.at(pos) );

        masterNode->giveLocationArray(dofIdArray, loc, s);
        slaveNode->giveLocationArray(dofIdArray, c_loc, s);
        this->computeExternalForcesFromContact(fext, masterNode, slaveNode, tStep);
        loc.followedBy(c_loc);
        answer.assemble(fext, loc);
    }
}






void
Node2NodePenaltyContact :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    IntArray r_loc, c_loc;
    rows.resize(masterSet.giveSize() );
    cols.resize(masterSet.giveSize() );
    IntArray dofIdArray = {
        D_u, D_v
    };

    for ( int pos = 1; pos <= masterSet.giveSize(); ++pos ) {
        Node *masterNode = this->giveDomain()->giveNode(masterSet.at(pos) );
        Node *slaveNode = this->giveDomain()->giveNode(slaveSet.at(pos) );

        masterNode->giveLocationArray(dofIdArray, r_loc, r_s);
        slaveNode->giveLocationArray(dofIdArray, c_loc, c_s);

        // column block
        rows [ pos - 1 ] = r_loc;
        cols [ pos - 1 ] = c_loc;
    }
}









void
Node2NodePenaltyContact :: computeTangentFromContact(FloatMatrix &answer, Node *masterNode, Node *slaveNode, TimeStep *tStep)
{
    double gap;
    FloatArray Nv;
    this->computeGap(gap, masterNode, slaveNode, tStep);
    this->computeNormalMatrixAt(Nv, masterNode, slaveNode, tStep);
    answer.beDyadicProductOf(Nv, Nv);
    answer.times(this->penalty);
    if ( gap > 0.0 ) {
        answer.zero();
    }
}

void
Node2NodePenaltyContact :: computeGap(double &answer, Node *masterNode, Node *slaveNode, TimeStep *tStep)
{
    FloatArray xs, xm, uS, uM;
    xs = * slaveNode->giveCoordinates();
    xm = * masterNode->giveCoordinates();
    FloatArray normal = xs - xm;
    double norm = normal.computeNorm();
    if ( norm < 1.0e-8 ) {
        OOFEM_ERROR("Couldn't compute normal between master node (num %d) and slave node (num %d), nodes are too close to each other.",
                    masterNode->giveGlobalNumber(), slaveNode->giveGlobalNumber() );
    } else {
        normal.times(1.0 / norm);
    }

    slaveNode->giveUnknownVector(uS, { D_u, D_v }, VM_Total, tStep, true);
    masterNode->giveUnknownVector(uM, { D_u, D_v }, VM_Total, tStep, true);
    xs.add(uS);
    xm.add(uM);
    FloatArray dx = xs - xm;
    answer = dx.dotProduct(normal);
}


void
Node2NodePenaltyContact :: computeNormalMatrixAt(FloatArray &answer, Node *masterNode, Node *slaveNode, TimeStep *TimeStep)
{
    FloatArray xs, xm;
    xs = * slaveNode->giveCoordinates();
    xm = * masterNode->giveCoordinates();
    FloatArray normal = xs - xm;
    double norm = normal.computeNorm();
    if ( norm < 1.0e-8 ) {
        OOFEM_ERROR("Couldn't compute normal between master node (num %d) and slave node (num %d), nodes are too close to each other.",
                    masterNode->giveGlobalNumber(), slaveNode->giveGlobalNumber() );
    } else {
        normal.times(1.0 / norm);
    }
    // The normal is not updated for node2node which is for small deformations only
    // C = {n -n}
    answer = {
        normal.at(1), normal.at(2),
        -normal.at(1), -normal.at(2)
    };
}


void
Node2NodePenaltyContact :: computeExternalForcesFromContact(FloatArray &answer, Node *masterNode, Node *slaveNode, TimeStep *tStep)
{
    double gap;
    this->computeGap(gap, masterNode, slaveNode, tStep);
    this->computeNormalMatrixAt(answer, masterNode, slaveNode, tStep);
    if ( gap < 0.0 ) {
        answer.times(penalty * gap);
    } else {
        answer.times(0);
    }
}
} // namespace oofem
