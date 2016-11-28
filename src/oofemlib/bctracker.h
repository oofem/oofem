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

#ifndef bctracker_h
#define bctracker_h

#include <list>
#include <vector>

namespace oofem {
class Domain;

/**
 * This class keeps track of applied boundary conditions on individual entities.
 * Some of the BCs can be applied via sets, in this case its BC and associated set
 * that keeps track which element or node is subjected to it.
 * In some applications, however, one needs the list of BCs applied on the component.
 * Cration and management of this relation is the task of this class.
 * Note: only element list is managed at present
 **/

class BCTracker {
 public:
  /// Helper class storing a sigle record for component (element, node, etc)
  struct Entry {
    /// bc number
    int bcNumber;
    /// boundari ID if required
    int boundaryId;

    Entry (int bc, int bid) {bcNumber=bc; boundaryId = bid;}
  };
  typedef std::list<Entry> entryListType;
  
 private:
  /// list keeping element entries
  std::vector<entryListType> elemList;
  /// Domain link
  Domain *domain;
  
 public:
  BCTracker (Domain *d);

  void initialize ();
  const entryListType& getElementRecords (int elem);
};
  
} // end namespace oofem
#endif


  
