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


#include "sm/Contact/ActiveBc/node2nodecontact.h"
#include "set.h"
#include "domain.h"
#include "node.h"
#include "floatmatrix.h"
#include "unknownnumberingscheme.h"
#include "sparsemtrx.h"

namespace oofem {


  IRResultType Node2NodePenaltyContact :: initializeFrom(InputRecord *ir)
  {
    IRResultType result;
    IR_GIVE_FIELD(ir, this->penalty, _IFT_Node2NodePenaltyContact_penalty);
    this->useTangent = ir->hasField( _IFT_Node2NodePenaltyContact_useTangent);


    IR_GIVE_FIELD(ir, this->masterSet, _IFT_Node2NodePenaltyContact_masterSet);
    IR_GIVE_FIELD(ir, this->slaveSet, _IFT_Node2NodePenaltyContact_slaveSet);

    return ActiveBoundaryCondition :: initializeFrom(ir);
  }





  void Node2NodePenaltyContact :: assemble(SparseMtrx &answer, TimeStep *tStep,
					   CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
  {
    if ( !this->useTangent || type != TangentStiffnessMatrix ) {
      return;
    }


    FloatMatrix Kaa, Kab;
    IntArray r_loc, c_loc;

    Set *masterSet = this->giveDomain()->giveSet(this->masterSet);
    Set *slaveSet = this->giveDomain()->giveSet(this->slaveSet);

    const IntArray &masterNodeList = masterSet->giveNodeList();
    const IntArray &slaveNodeList = slaveSet->giveNodeList();

    IntArray dofIdArray = {D_u, D_v, D_w};

    if(masterNodeList.giveSize() != slaveNodeList.giveSize()) {
      OOFEM_ERROR("Number of master nodes does not match number of slave nodes")
	}

    for ( int pos = 1; pos <= masterNodeList.giveSize(); ++pos ) {
      Node *masterNode = this->giveDomain()->giveNode( masterNodeList.at(pos) );
      Node *slaveNode = this->giveDomain()->giveNode( slaveNodeList.at(pos) );

      masterNode->giveLocationArray(dofIdArray, r_loc, r_s);
      slaveNode->giveLocationArray(dofIdArray, c_loc, c_s);
      this->computeTangentFromContact(Kaa, masterNode, slaveNode, tStep);
      Kab = Kaa;
      Kab.times(-1);
      answer.assemble(r_loc, r_loc, Kaa);
      answer.assemble(c_loc, c_loc, Kaa);
      answer.assemble(r_loc, c_loc, Kab);
      answer.assemble(c_loc, r_loc, Kab);
    }
  }

  void Node2NodePenaltyContact :: assembleVector(FloatArray &answer, TimeStep *tStep,
						 CharType type, ValueModeType mode,
						 const UnknownNumberingScheme &s, FloatArray *eNorms)
  {
    if ( type != ExternalForcesVector ) {
      return;
    }

    IntArray dofIdArray = {D_u, D_v, D_w};
    IntArray r_loc, c_loc;
    FloatArray fm, fs;
    IntArray loc, masterdofids, bNodes;
    
    Set *masterSet = this->giveDomain()->giveSet(this->set);
    Set *slaveSet = this->giveDomain()->giveSet(this->slaveSet);

    const IntArray &masterNodeList = masterSet->giveNodeList();
    const IntArray &slaveNodeList = slaveSet->giveNodeList();

    if(masterNodeList.giveSize() != slaveNodeList.giveSize()) {
      OOFEM_ERROR("Number of master nodes does not match number of slave nodes")
	}

    for ( int pos = 1; pos <= masterNodeList.giveSize(); ++pos ) {
      Node *masterNode = this->giveDomain()->giveNode( masterNodeList.at(pos) );
      Node *slaveNode = this->giveDomain()->giveNode( slaveNodeList.at(pos) );

      masterNode->giveLocationArray(dofIdArray, r_loc, s);
      slaveNode->giveLocationArray(dofIdArray, c_loc, s);
      this->computeInternalForcesFromContact(fm, masterNode, slaveNode, tStep);
      fs = fm;
      fs.times(-1);
      answer.assemble(fm, r_loc);
      answer.assemble(fs, c_loc);
    }
    
  }
  
  void
  Node2NodePenaltyContact :: computeTangentFromContact(FloatMatrix &answer, Node* masterNode, Node *slaveNode, TimeStep *tStep)
  {
    double gap;
    FloatArray Nv;
    this->computeGap(gap, masterNode, slaveNode, tStep);
    this->computeNormalMatrixAt(Nv, masterNode, slaveNode, tStep);
    answer.beTProductOf(Nv,Nv);
    // this is the interface stiffness and should be obtained from that model
    answer.times( this->penalty);
    answer.negated();
    if( gap > 0.0 ) {
      answer.zero();
    }    
  } 

  void
  Node2NodePenaltyContact :: computeGap(double &answer,  Node* masterNode, Node *slaveNode, TimeStep *tStep)
  {
    FloatArray xs, xm, uS, uM;
    xs = *slaveNode->giveCoordinates();
    xm = *masterNode->giveCoordinates();
    FloatArray normal = xs-xm;
    double norm = normal.computeNorm();
    if ( norm < 1.0e-8 ) {
      OOFEM_ERROR("Couldn't compute normal between master node (num %d) and slave node (num %d), nodes are too close to each other.", 
		  masterNode->giveGlobalNumber(), slaveNode->giveGlobalNumber() );
    } else {
      normal.times(1.0/norm);
    }
    
    slaveNode->giveUnknownVector(uS, {D_u, D_v, D_w}, VM_Total, tStep, true);
    masterNode->giveUnknownVector(uM, {D_u, D_v, D_w}, VM_Total, tStep, true);
    xs.add(uS);
    xm.add(uM);
    FloatArray dx = xs-xm;
    answer = dx.dotProduct(normal);

  }


  void
  Node2NodePenaltyContact :: computeNormalMatrixAt( FloatArray &answer, Node* masterNode, Node *slaveNode, TimeStep *TimeStep)
  {
    FloatArray xs, xm;
    xs = *slaveNode->giveCoordinates();
    xm = *masterNode->giveCoordinates();
    FloatArray normal = xs-xm;
    double norm = normal.computeNorm();
    if ( norm < 1.0e-8 ) {
      OOFEM_ERROR("Couldn't compute normal between master node (num %d) and slave node (num %d), nodes are too close to each other.", 
		  masterNode->giveGlobalNumber(), slaveNode->giveGlobalNumber() );
    } else {
      normal.times(1.0/norm);
    }
    // The normal is not updated for node2node which is for small deformations only
    // C = {n -n}
    answer = {  normal.at(1),  normal.at(2),  normal.at(3),
		-normal.at(1), -normal.at(2), -normal.at(3) };
  }


  void
  Node2NodePenaltyContact :: computeInternalForcesFromContact(FloatArray &answer,  Node *masterNode, Node *slaveNode, TimeStep *tStep)
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

  
}
