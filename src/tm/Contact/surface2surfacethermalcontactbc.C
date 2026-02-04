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
 *               Copyright (C) 1993 - 2023   Borek Patzak
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


#include "surface2surfacethermalcontactbc.h"
#include "set.h"
#include "domain.h"
#include "node.h"
#include "floatmatrix.h"
#include "classfactory.h"

namespace oofem {
REGISTER_BoundaryCondition(Surface2SurfaceThermalContactBoundaryCondition);

void
Surface2SurfaceThermalContactBoundaryCondition :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    ThermalContactBoundaryCondition :: initializeFrom(ir);
    // set of master nodes
    IR_GIVE_FIELD(ir, this->slaveSurfaceNumber, _IFT_Surface2SurfaceThermalContactBoundaryCondition_slaveSurfaceNum);
    // set of slave nodes
    IR_GIVE_FIELD(ir, this->masterSurfaceNumber, _IFT_Surface2SurfaceThermalContactBoundaryCondition_masterSurfaceNum);
    //
    this->slaveContactSurface = dynamic_cast<ThermalFEContactSurface*>( this->giveDomain()->giveContactSurface(slaveSurfaceNumber ));
    this->masterContactSurface = dynamic_cast<ThermalFEContactSurface*>( this->giveDomain()->giveContactSurface(masterSurfaceNumber ));
    
    

}

void
Surface2SurfaceThermalContactBoundaryCondition :: setupContactSearchAlgorithm()
{
  this->contactSearchAlgorithm = std::make_unique<ContactSearchAlgorithm_Surface2FESurface_3d>(this->slaveContactSurface, this->masterContactSurface, domain);
}



void
Surface2SurfaceThermalContactBoundaryCondition :: giveLocationArrays(std::vector< IntArray > &rows, std::vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    //returns all possible combinations of dof that can theoretically be triggered by contact
    //of any segment with any node. Room for optimization aplenty...
   IntArray m_loc, s_loc;

    rows.resize(0);
    cols.resize(0);

    const auto& contactPairs = getContactPairs();
    for(auto const &cp : contactPairs) {
      if(cp->inContact()) {
	cp->giveSlaveContactPoint()->giveLocationArray( s_loc, dofs, c_s);
	cp->giveMasterContactPoint()->giveLocationArray( m_loc, dofs, r_s);
	// insert location arrays into the answer arrays
	rows.push_back(s_loc);
	cols.push_back(m_loc);
	rows.push_back(m_loc);
	cols.push_back(s_loc);
      }
    }
    /*
    IntArray m_loc, s_loc;

    int ncombinations = nodeSet.giveSize();
    rows.resize(ncombinations * 2);
    cols.resize(ncombinations * 2);

    int pos = 0;
    for(auto const &cp : contactPairs) {
      cp->giveMasterContactPoint()->giveLocationArray( m_loc, dofs, r_s);
      cp->giveSlaveContactPoint()->giveLocationArray( s_loc, dofs, c_s);
      if(m_loc.giveSize() && s_loc.giveSize()) {
	// insert location arrays into the answer arrays
	rows [ pos ]     = m_loc;
	cols [ pos ]     = s_loc;
	rows [ pos + 1 ] = s_loc;
	cols [ pos + 1 ] = m_loc;
	pos += 2;
      }
      }*/
}


 
void
Surface2SurfaceThermalContactBoundaryCondition :: initForNewIteration(TimeStep *tStep, int iter)
{
  if(!initializedPairs && iter % updateEachNthIter == 0) {
    ContactBoundaryCondition :: initForNewIteration(tStep, iter);
    initializedPairs = true;
  }

}
  
} // namespace oofem
