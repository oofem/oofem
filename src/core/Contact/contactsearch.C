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

#include "contactsearch.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "floatmatrix.h"
#include "node.h"
#include "domain.h"
#include "fecontactsurface.h"
#include "contactpoint.h"
#include "contactpair.h"
#include <limits>


namespace oofem {


  ContactSearchAlgorithm_Surface2FESurface :: ContactSearchAlgorithm_Surface2FESurface(FEContactSurface *scs, FEContactSurface *mcs, Domain *d, int sd) : ContactSearchAlgorithm(d)
{
  this->slaveContactSurface = scs;
  this->masterContactSurface = mcs;
  this->surface_dimension  = sd;
}


void
ContactSearchAlgorithm_Surface2FESurface :: createContactPairs()
{
  contactPairs.clear();
  for ( int slave_ce = 1; slave_ce <= this->slaveContactSurface->giveNumberOfContactElements(); slave_ce++ ) {
    for ( GaussPoint *gp : * this->slaveContactSurface->giveContactElement_InSet(slave_ce)->giveDefaultIntegrationRulePtr() ) {
      auto cp_slave = std::make_unique<FEContactPoint_Slave>(this->slaveContactSurface, this->slaveContactSurface->giveContactElement_InSet(slave_ce)->giveNumber(),  surface_dimension, gp);
      auto contactPair = std::make_unique<ContactPair>(std::move(cp_slave));
      contactPairs.emplace_back(std::move(contactPair));
    }
  }
}

  ContactSearchAlgorithm_Surface2FESurface_3d ::   ContactSearchAlgorithm_Surface2FESurface_3d(FEContactSurface *scs, FEContactSurface *mcs, Domain *d) :   ContactSearchAlgorithm_Surface2FESurface(scs, mcs, d, 2)
{
}

  
void
ContactSearchAlgorithm_Surface2FESurface_3d :: updateContactPairs(TimeStep *tStep)
{

  FloatArrayF<3> normalVector, tangentVector1, tangentVector2;
  for (auto &cp : contactPairs) {
      auto slavePoint = dynamic_cast<FEContactPoint_Slave*> (cp->giveSlaveContactPoint());
      //iterate over all contact elements
      int closestContactElementId = -1;
      FloatArray contactPointLocalCoordinates;
      double gap = 0;
      for ( int i = 1; i <= this->masterContactSurface->giveNumberOfContactElements(); i++ ) {
	auto contactElement = this->masterContactSurface->giveContactElement_InSet(i);
	// check if we are not testing the point with its own element
	if(slavePoint->giveContactElementId() == contactElement->giveNumber()) {
	  break;
	}
	///first check the contact element bounding box
	/// if the element is inside, do more detailed search
	  /// otherwise continue;
	  /* auto boundingBox = contactElement.giveBoundingBox();
	     if(!boundingBox.contains(node.giveCoordinates())) {
	       continue;
	     } else {*/         
	auto [inElement,localCoords, newGap,normal, t1, t2] = this->masterContactSurface->findContactPointInElement_3d(slavePoint, contactElement, tStep);
	if(inElement) {
	  //FloatArray globalCoords, nodeCoords;
	  //
	  //contactElement->giveInterpolation()->local2global( globalCoords, localCoords, FEIElementDeformedGeometryWrapper(contactElement, tStep) );
	  //contact when gap is negative !!!
	  if ( closestContactElementId == -1. || newGap < gap ) {
	    //larger gap/penetration found, update it
	    gap = newGap;
	    normalVector = normal;
	    tangentVector1 = -t1;
	    tangentVector2 = -t2;
	    contactPointLocalCoordinates = localCoords;
	    closestContactElementId = contactElement->giveNumber();
	  }
	}
      }
      
      //update contact pairs for each node
      if(closestContactElementId) {
	auto master_point  =  std::make_unique<FEContactPoint_Master>(this->masterContactSurface, closestContactElementId, 2, contactPointLocalCoordinates);
	cp->setMasterContactPoint(std::move(master_point));
	cp->setNormalGap(gap);
	cp->setNormalVector(normalVector);
	cp->setTangentVector1(tangentVector1);
	cp->setTangentVector2(tangentVector2);  
 
      }
  }
  
}

  

  ContactSearchAlgorithm_Surface2FESurface_2d ::   ContactSearchAlgorithm_Surface2FESurface_2d(FEContactSurface *scs, FEContactSurface *mcs, Domain *d) :   ContactSearchAlgorithm_Surface2FESurface(scs, mcs, d, 1)
{
}
  

void
ContactSearchAlgorithm_Surface2FESurface_2d :: updateContactPairs(TimeStep *tStep)
{

  FloatArrayF<2> normalVector, tangentVector1;
  for (auto &cp : contactPairs) {
      auto slavePoint = dynamic_cast<FEContactPoint_Slave*> (cp->giveSlaveContactPoint());
       //iterate over all contact elements
      int closestContactElementId = -1;
      FloatArray contactPointLocalCoordinates, contactPointGlobalCoordinates;
      double gap = std::numeric_limits<double>::infinity();
      for ( int i = 1; i <= this->masterContactSurface->giveNumberOfContactElements(); i++ ) {
	auto masterContactElement = this->masterContactSurface->giveContactElement_InSet(i);
	// check if we are not testing the point with its own element
	if(slavePoint->giveContactElementId() == masterContactElement->giveNumber()) {
	  continue;
	}

	///first check the contact element bounding box
	/// if the element is inside, do more detailed search
	  /// otherwise continue;
	  /* auto boundingBox = contactElement.giveBoundingBox();
	     if(!boundingBox.contains(node.giveCoordinates())) {
	       continue;
	     } else {*/         
	auto [inElement,localCoords, newGap,normal, t1] = this->masterContactSurface->findContactPointInElement_2d(slavePoint, masterContactElement, tStep);
	if(inElement) {
	  FloatArray globalCoords, nodeCoords;
	  //
	  masterContactElement->giveInterpolation()->local2global( globalCoords, localCoords, FEIElementDeformedGeometryWrapper(masterContactElement, tStep) );
	  //contact when gap is negative!!!
	  if ( closestContactElementId == -1. || newGap < gap ) {
	    //larger gap/penetration found, update it
	    gap = newGap;
	    normalVector = normal;
	    tangentVector1 = t1;
	    contactPointLocalCoordinates = localCoords;
	    contactPointGlobalCoordinates = globalCoords;
	    closestContactElementId = masterContactElement->giveNumber();
	  }
	}
      }
      
      //update contact pairs for each node
      if(closestContactElementId) {
	auto masterContactPoint =  std::make_unique<FEContactPoint_Master>(this->masterContactSurface, closestContactElementId, 1, contactPointLocalCoordinates);
	cp->setMasterContactPoint(std::move(masterContactPoint));
	cp->setNormalGap(gap);
	cp->setNormalVector(normalVector);
	cp->setTangentVector1(tangentVector1);
      }
  }
  
}


  
} // end namespace oofem
