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

#include "Contact/contactelement.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"

#include "mathfem.h"
#include "feinterpol.h"
#include "node.h"



/*
#include "entityrenumberingscheme.h"
#include "contextioerr.h"
#include "feinterpol1d.h"
#include "feinterpol2d.h"
#include "feinterpol3d.h"
#include "dofmanager.h"
#include "node.h"
#include "gausspoint.h"
#include "unknownnumberingscheme.h"
#include "cltypes.h"
#include <cstdio>
*/
namespace oofem {
  ContactElement :: ContactElement(int n, Domain *aDomain) :  Element(n, aDomain)
{

}


ContactElement :: ~ContactElement()
{
}

AABB
ContactElement :: computeAABB()
{
  AABB aabb;
  int nnode = this->giveNumberOfNodes();
  for (int i=1; i<=nnode; i++) {
      Node *node = this->giveNode(i);
      const auto& coords = node->giveCoordinates();
      double x = coords.at(1);
      double y = coords.at(2);
      double z = ( coords.size() > 2 ) ? coords.at(3) : 0.0;
      aabb.merge(x,y,z);
  }
  return aabb;
}

AABB
ContactElement :: computeUpdatedAABB(TimeStep *tStep)
{
  if (tStep == nullptr) {
      return this->computeAABB();
  }
  AABB aabb;
  int nnode = this->giveNumberOfNodes();
  for (int i = 1; i <= nnode; i++) {
      Node *node = this->giveNode(i);
      double x = node->giveUpdatedCoordinate(1, tStep);
      double y = node->giveUpdatedCoordinate(2, tStep);
      double z = (node->giveCoordinates().giveSize() > 2) ? node->giveUpdatedCoordinate(3, tStep) : 0.0;
      aabb.merge(x, y, z);
  }
  return aabb;
}

void
ContactElement :: setContactOutputState(const GaussPoint *gp, double normalGap,
                                        double pressure, int status)
{
  if (gp == nullptr) {
    return;
  }
  contactOutputStates[gp] = { normalGap, pressure, static_cast<double>(status) };
}

int
ContactElement :: giveIPValue(FloatArray &answer, GaussPoint *gp,
                              InternalStateType type, TimeStep *tStep)
{
  if (type != IST_ContactNormalGap
      && type != IST_ContactPressure
      && type != IST_ContactStatus) {
    return Element :: giveIPValue(answer, gp, type, tStep);
  }

  answer.resize(1);
  const auto state = contactOutputStates.find(gp);
  if (state == contactOutputStates.end()) {
    answer.at(1) = 0.0;
  } else if (type == IST_ContactNormalGap) {
    answer.at(1) = state->second.normalGap;
  } else if (type == IST_ContactPressure) {
    answer.at(1) = state->second.pressure;
  } else {
    answer.at(1) = state->second.status;
  }
  return 1;
}

  

} // end namespace oofem
