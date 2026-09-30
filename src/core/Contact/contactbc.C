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
#include "datastream.h"
#include "contextioerr.h"
#include "contactpoint.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cstdlib>
#include <cstdio>

namespace oofem {

namespace {
// Diagnostic-only: dump active master-facet ownership per Newton iteration
// when OOFEM_CONTACT_TRACE_FACETS is set in the environment.  Used to confirm
// or refute contact-facet chattering (the closest master facet flipping every
// iteration) as the cause of a non-decaying Newton residual cycle; see
// doc/contact-improvement-handoff.md, 2026-07-28 co47 crossed-tubes entry.
void traceActiveContactFacets(const ContactBoundaryCondition &bc, TimeStep *tStep, int iter)
{
  static const bool enabled = std::getenv("OOFEM_CONTACT_TRACE_FACETS") != nullptr;
  if (!enabled || tStep == nullptr) {
    return;
  }

  const auto &pairs = bc.getContactPairs();
  int activeCount = 0;
  for (const auto &cp : pairs) {
    if (cp->hasActiveContact()) {
      activeCount++;
    }
  }
  std::fprintf(stderr, "[contact-trace] bc=%d step=%d iter=%d active=%d\n",
               bc.giveNumber(), tStep->giveNumber(), iter, activeCount);

  for (std::size_t i = 0; i < pairs.size(); ++i) {
    ContactPair *cp = pairs[i].get();
    if (!cp->hasActiveContact()) {
      continue;
    }
    auto *slavePoint = dynamic_cast<FEContactPoint *>(cp->giveSlaveContactPoint());
    auto *masterPoint = dynamic_cast<FEContactPoint *>(cp->giveMasterContactPoint());
    const int slaveElemId = slavePoint ? slavePoint->giveContactElementId() : -1;
    const int masterElemId = masterPoint ? masterPoint->giveContactElementId() : -1;
    const FloatArray *mlc = masterPoint ? &masterPoint->giveLocalCoordinates() : nullptr;
    const double xi1 = (mlc && mlc->giveSize() > 0) ? mlc->at(1) : 0.0;
    const double xi2 = (mlc && mlc->giveSize() > 1) ? mlc->at(2) : 0.0;
    std::fprintf(stderr,
        "[contact-trace]   pair=%zu slaveElem=%d masterElem=%d feature=%d/%d mlc=(% .6f,% .6f) gap=% .6e\n",
        i, slaveElemId, masterElemId,
        static_cast<int>(cp->giveCurrentMasterFeatureType()), cp->giveCurrentMasterFeatureIndex(),
        xi1, xi2, cp->giveNormalGap());
  }
}
} // unnamed namespace


void
ContactBoundaryCondition :: initForNewIteration(TimeStep *tStep, int iter)
{
  if(iter % updateEachNthIter == 0) {
    this->giveContactSearchAlgorithm()->updateContactPairs(tStep);
  }

  traceActiveContactFacets(*this, tStep, iter);
}
  
 


void
ContactBoundaryCondition :: assemble(SparseMtrx &answer, TimeStep *tStep, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale, void *lock)
{
    if ( type != TangentStiffnessMatrix ) {
        return;
    }

    FloatMatrix K;
    IntArray rowLoc, colLoc;

    //iterate over all pairs of nodes and segments
    const auto& contactPairs = getContactPairs();
    for(auto const &cp : contactPairs) {
      if(cp->hasMasterContact()) {
        K.clear();
	this->computeTangentFromContact(K, cp.get(), tStep);
	this->giveRowLocationArray(rowLoc, r_s, cp.get());
	this->giveColumnLocationArray(colLoc, c_s, cp.get());
	if(K.giveNumberOfRows() && K.giveNumberOfColumns()) {
          // A smooth/frictionless tangent depends only on the current master
          // and slave facets.  Frictional convective transport across a facet
          // boundary additionally depends on exclusive nodes of the committed
          // master facet and therefore returns a wider, rectangular tangent.
          // Select the column map from the tangent support actually produced.
          if (K.giveNumberOfColumns() != colLoc.giveSize()) {
            IntArray currentFacetColLoc;
            cp->giveRowLocationArray(dofs, currentFacetColLoc, c_s);
            if (K.giveNumberOfColumns() == currentFacetColLoc.giveSize()) {
              colLoc = currentFacetColLoc;
            } else {
              OOFEM_ERROR(
                "ContactBoundaryCondition: tangent has %d columns, but current "
                "and transition location arrays have %d and %d entries",
                K.giveNumberOfColumns(), currentFacetColLoc.giveSize(),
                colLoc.giveSize());
            }
          }
	          K.times(scale);
#ifdef _OPENMP
          if (lock) omp_set_lock(static_cast<omp_lock_t *>(lock));
#endif
	  answer.assemble(rowLoc, colLoc, K);
#ifdef _OPENMP
          if (lock) omp_unset_lock(static_cast<omp_lock_t *>(lock));
#endif
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


    IntArray loc, dofIds;
    FloatArray fint;
    //iterate over all pairs of nodes and segments
    const auto& contactPairs = getContactPairs();
    for(auto const &cp : contactPairs) {
      if(cp->hasMasterContact()) {
        fint.clear();
	this->computeInternalForcesFromContact(fint, cp.get(), tStep);
	this->giveRowLocationArray(loc, s, cp.get());
	if(fint.giveSize()) {
#ifdef _OPENMP
          if (lock) omp_set_lock(static_cast<omp_lock_t *>(lock));
#endif
	  answer.assemble(fint, loc);
          if (eNorms) {
            dofIds.resize(fint.giveSize());
            for (int i = 1; i <= fint.giveSize(); ++i) {
              dofIds.at(i) = dofs.at((i - 1) % dofs.giveSize() + 1);
            }
            eNorms->assembleSquared(fint, dofIds);
          }
#ifdef _OPENMP
          if (lock) omp_unset_lock(static_cast<omp_lock_t *>(lock));
#endif
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
      if(cp->hasMasterContact()) {
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
		  this->giveRowLocationArray(loc, dn, cp.get());
	  answer.assemble(fext,loc);
	}
      }
    }
}


  
void
ContactBoundaryCondition :: giveRowLocationArray(
    IntArray &loc, const UnknownNumberingScheme &ns, const ContactPair *cp) const
{
  cp->giveRowLocationArray(dofs, loc, ns);
}

void
ContactBoundaryCondition :: giveColumnLocationArray(
    IntArray &loc, const UnknownNumberingScheme &ns, const ContactPair *cp) const
{
  cp->giveColumnLocationArray(dofs, loc, ns);
}


void ContactBoundaryCondition :: postInitialize()
{
  this->setupContactSearchAlgorithm();
  this->giveContactSearchAlgorithm()->createContactPairs();
  // Establish projection history for bodies supported by initially touching
  // frictional interfaces, so their first loaded tangent is not mechanically
  // singular in the tangential directions.
  this->giveContactSearchAlgorithm()->updateContactPairs(nullptr);
  const auto &contactPairs = this->getContactPairs();
  for (auto &cp : contactPairs) {
    cp->updateYourself(nullptr);
  }
  this->giveContactSearchAlgorithm()->commitSearchState(nullptr);
}

void
ContactBoundaryCondition :: saveContext(DataStream &stream, ContextMode mode)
{
  ActiveBoundaryCondition :: saveContext(stream, mode);
  if ((mode & CM_Definition) && !stream.write(updateEachNthIter)) {
    THROW_CIOERR(CIO_IOERR);
  }
}

void
ContactBoundaryCondition :: restoreContext(DataStream &stream, ContextMode mode)
{
  ActiveBoundaryCondition :: restoreContext(stream, mode);
  if ((mode & CM_Definition) && !stream.read(updateEachNthIter)) {
    THROW_CIOERR(CIO_IOERR);
  }
}


void ContactBoundaryCondition :: updateYourself(TimeStep *tStep)
{

  const auto& contactPairs = getContactPairs();
  for(auto &cp : contactPairs) {
    cp->updateYourself(tStep);
  }
  this->giveContactSearchAlgorithm()->commitSearchState(tStep);

}



void
ContactBoundaryCondition :: giveExportData(std::vector< ExportRegion > &vtkPieces, FloatArray shift, TimeStep *tStep )
{
  
    vtkPieces.resize(1);
 
    const auto& contactPairs = getContactPairs();
    const int numCellNodes  = 2; // linear line
    int numCells = 0;
    for (auto const &cp : contactPairs) {
      if (cp->hasMasterContact()) {
        numCells++;
      }
    }
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
    FloatArray nodeCoords;
    Coordinates updatedNodeCoords;
    IntArray connectivity(2);
    for(auto const &cp : contactPairs) {
      if(cp->hasMasterContact()) {
	cp->giveMasterContactPoint()->giveUpdatedCoordinates(updatedNodeCoords, tStep);
	nodeCoords = updatedNodeCoords;
	nodeCoords.resizeWithValues(3);
	if(shift.giveSize()){
	  nodeCoords.add(shift);
	}
	vtkPieces.at(0).setNodeCoords(nodeNum, nodeCoords);
	connectivity.at(1) = val++;
	nodeNum++;
	cp->giveSlaveContactPoint()->giveUpdatedCoordinates(updatedNodeCoords, tStep);
	nodeCoords = updatedNodeCoords;
	nodeCoords.resizeWithValues(3);
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
      }
    } 
}



  
  
      
  
} // namespace oofem




