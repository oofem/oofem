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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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


#include "contactbc.h"
#include "set.h"
#include "domain.h"
#include "node.h"
#include "floatmatrix.h"
#include "sparsemtrx.h"
#include "contactpair.h"
#include "contactsearch.h"
#include "unknownnumberingscheme.h"
#include "timestep.h"
#include "vtkxmlexportmodule.h"
#include "vtkbaseexportmodule.h"

namespace oofem {
 

void
ContactBoundaryCondition :: initForNewIteration(TimeStep *tStep, int iter)
{
  if(iter % updateEachNthIter == 0) {
    this->giveContactSearchAlgorithm()->updateContactPairs(tStep);
  }

}
  
 


void
ContactBoundaryCondition :: assemble(SparseMtrx &answer, TimeStep *tStep, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale, void *lock)
{
    if ( type != TangentStiffnessMatrix ) {
        return;
    }

    FloatMatrix K;
    IntArray loc, node_loc;

    //iterate over all pairs of nodes and segments
    const auto& contactPairs = getContactPairs();
    for(auto const &cp : contactPairs) {
      if(cp->inContact()) {
	this->computeTangentFromContact(K, cp.get(), tStep);
	this->giveLocationArray(loc, r_s, cp.get());
	if(K.giveNumberOfRows() && K.giveNumberOfColumns()) {
	  answer.assemble(loc, K);
	}
      }
    }
}


void
ContactBoundaryCondition :: assembleVector(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode, const UnknownNumberingScheme &s, FloatArray *eNorms, void *lock)
{
    if ( type != InternalForcesVector ) {
        return;
    }


    IntArray loc;
    FloatArray fint;
    //iterate over all pairs of nodes and segments
    const auto& contactPairs = getContactPairs();
    for(auto const &cp : contactPairs) {
      if(cp->inContact()) {
	this->computeInternalForcesFromContact(fint, cp.get(), tStep);
	this->giveLocationArray(loc, s, cp.get());
	if(fint.giveSize()) {
	  answer.assemble(fint, loc);
	}
      }
    }
}




void
ContactBoundaryCondition :: assembleExtrapolatedForces(FloatArray &answer, TimeStep *tStep)
{
  /*    if ( type != TangentStiffnessMatrix ) {
        return;
    }
  */
    FloatArray fext;
    FloatMatrix K;
    IntArray loc;

    //iterate over all pairs of nodes and segments
    const auto& contactPairs = getContactPairs();
    for(auto const &cp : contactPairs) {
      if(cp->inContact()) {
	this->computeTangentFromContact(K, cp.get(), tStep);
	if(K.giveNumberOfRows() && K.giveNumberOfColumns()) {
	  //cp->computeVectorOf(VM_Incremental, tStep, delta_u);
	  FloatArray delta_u;
	  cp->computeVectorOf(VM_Total, tStep, delta_u);
	  FloatArray tmp;
	  if ( tStep->isTheFirstStep() ) {
	    tmp = delta_u;
	    tmp.zero();
	  } else {
	    cp->computeVectorOf(VM_Total, tStep->givePreviousStep(), tmp);
	  }
	  delta_u.subtract(tmp);
	  fext.beProductOf(K, delta_u);
	  EModelDefaultEquationNumbering dn;
	  this->giveLocationArray(loc, dn, cp.get());
	  answer.assemble(fext,loc);
	}
      }
    }
}


  
void
ContactBoundaryCondition :: giveLocationArray(IntArray &loc, const UnknownNumberingScheme &ns, const ContactPair *cp) const
{
  cp->giveLocationArray(dofs, loc, ns);  
}


void ContactBoundaryCondition :: postInitialize()
{
  this->setupContactSearchAlgorithm();
  this->giveContactSearchAlgorithm()->createContactPairs();

}


void ContactBoundaryCondition :: updateYourself(TimeStep *tStep)
{

  const auto& contactPairs = getContactPairs();
  for(auto &cp : contactPairs) {
    cp->updateYourself(tStep);
  }

}



void
ContactBoundaryCondition :: giveExportData(std::vector< ExportRegion > &vtkPieces, FloatArray shift, TimeStep *tStep )
{
  
    vtkPieces.resize(1);
 
    const auto& contactPairs = getContactPairs();
    int numCells = contactPairs.size();
    const int numCellNodes  = 2; // linear line
    int nNodes = numCells * numCellNodes;
    //
    vtkPieces.at(0).setNumberOfCells(numCells);
    vtkPieces.at(0).setNumberOfNodes(nNodes);
    //
    int val    = 1;
    int offset = 0;
    IntArray nodes(numCellNodes);
    int nodeNum = 1;
    int iElement = 1;
    FloatArray nodeCoords(3);
    IntArray connectivity(2);
    for(auto const &cp : contactPairs) {
      if(cp->inContact()) {
	nodeCoords = cp->giveMasterContactPoint()->giveGlobalCoordinates();
	if(shift.giveSize()){
	  nodeCoords.add(shift);
	}
	vtkPieces.at(0).setNodeCoords(nodeNum, nodeCoords);
	connectivity.at(1) = val++;
	nodeNum++;
	nodeCoords = cp->giveSlaveContactPoint()->giveGlobalCoordinates();
	if(shift.giveSize()){
	  nodeCoords.add(-1.*shift);
	}
	vtkPieces.at(0).setNodeCoords(nodeNum, nodeCoords);
	connectivity.at(2) = val++;
	nodeNum++;
	
	vtkPieces.at(0).setConnectivity(iElement, connectivity);
	offset += 2;
	vtkPieces.at(0).setOffset(iElement, offset);
	vtkPieces.at(0).setCellType(iElement, 3);
	iElement++;
      } else {
	numCells--;
	nNodes -= 2;
	vtkPieces.at(0).setNumberOfCells(numCells);
	vtkPieces.at(0).setNumberOfNodes(nNodes);   
      }
    } 
}



  
  
      
  
} // namespace oofem















