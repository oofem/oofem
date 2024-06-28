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

#include "domain.h"
#include "intarray.h"
#include "element.h"
#include "set.h"
#include "load.h"
#include "bodyload.h"
#include "boundaryload.h"
#include "nodalload.h"
#include "activebc.h"
#include "timestep.h"
#ifdef __SM_MODULE
#include "Loads/structtemperatureload.h"
#include "Loads/structeigenstrainload.h"
#endif
namespace oofem {

BCTracker :: BCTracker(Domain* d) {
  this->domain = d;
}


void
BCTracker::initialize() {
  this->elemList.clear();
  this->elemList.resize(domain->giveNumberOfElements());

  int nbc = domain->giveNumberOfBoundaryConditions();
  for ( int ibc = 1; ibc <= nbc; ++ibc ) {
    GeneralBoundaryCondition *bc = domain->giveBc(ibc);
    ActiveBoundaryCondition *abc;
    Load *load;

    if ( ( abc = dynamic_cast< ActiveBoundaryCondition * >(bc) ) ) {
      continue;
    } else if ( bc->giveSetNumber() && ( load = dynamic_cast< Load * >(bc) )) {
      BodyLoad *bodyLoad;
      BoundaryLoad *bLoad;
      Set *set = domain->giveSet( bc->giveSetNumber() );

      if ( ( bodyLoad = dynamic_cast< BodyLoad * >(bc) ) ) { // Body load:
        const IntArray &elements = set->giveElementList();
        for ( int ielem = 1; ielem <= elements.giveSize(); ++ielem ) {
          Entry entry (ibc, 0);
          this->elemList[elements.at(ielem)-1].push_back(entry);
        }
      } else if ( ( bLoad = dynamic_cast< BoundaryLoad * >(bc) ) ) { // Boundary load:
        const IntArray &boundaries = set->giveBoundaryList();
        for ( int ibnd = 1; ibnd <= boundaries.giveSize() / 2; ++ibnd ) {
          int eid = boundaries.at(ibnd * 2 - 1) ;
          int bid = boundaries.at(ibnd * 2);
          Entry entry (ibc, bid);
          this->elemList[eid-1].push_back(entry);
        }

        ///@todo Should we have a seperate entry for edge loads? Just sticking to the general "boundaryload" for now.
        const IntArray &edgeBoundaries = set->giveEdgeList();
        for ( int ibnd = 1; ibnd <= edgeBoundaries.giveSize() / 2; ++ibnd ) {
          int eid = edgeBoundaries.at(ibnd * 2 - 1) ;
          int bid = edgeBoundaries.at(ibnd * 2);
          Entry entry(ibc, bid);
          this->elemList[eid-1].push_back(entry);
        }
      }
#ifdef __SM_MODULE
      else if ( ( load = dynamic_cast< StructuralTemperatureLoad * >(bc) ) || ( load = dynamic_cast< StructuralEigenstrainLoad * >(bc) )) { // Body load:
        const IntArray &elements = set->giveElementList();
        for ( int ielem = 1; ielem <= elements.giveSize(); ++ielem ) {
          Entry entry (ibc, 0);
          this->elemList[elements.at(ielem)-1].push_back(entry);
        }
      }
#endif
    }
  }// end loop over BCs
}
  
const BCTracker::entryListType&
BCTracker::getElementRecords (int elem) {
  if (this->elemList.empty()) {
    this->initialize();
  }
  return this->elemList.at(elem-1);
}

} // end namespace oofem
