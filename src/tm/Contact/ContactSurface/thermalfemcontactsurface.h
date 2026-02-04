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



#pragma once
#include "contactsurface.h"
#include "set.h"

#define _IFT_ThermalFEContactSurface_elementSurfaceSetNumber "elementsurfaceset"


namespace oofem {
class ThermalContactSurface : public ContactSurface
{
protected:
  IntArray elementSurfaces;
  
public:
    ThermalContactSurface(int n, Domain *aDomain) : ContactSurface(n, aDomain) {; }
    ~ThermalContactSurface() {};

    void initializeFrom(InputRecord &ir) override;
    void updateYourself(TimeStep *tStep) override;
    // post initialization
    void postInitialize() override;
    virtual IntArray *giveClosestElementForNode(Node *n) {return &closestContactElement.at(n);}
  
protected:
    std::vector<std::uniq_ptr<ContactElement>> contactElements;

  /*
    enum class UpdateMode { UM_Never, UM_EachStep, UM_EveryIteration };
    UpdateMode updateMode;
  */
    //gives the closest edge to a given node in the form of an IntArray(elempos,edgepos)
    //only computes it again if it wasn't determined for this node in this solution step yet
    //(i.e. the array of known closest edges is reset after step convergence)
  //void give_closestElement_and_contactPointLocal(IntArray &closestEdge, FloatArray &contactPoint, Node *node, TimeStep *tStep);

    //computes the closest point on a given element to a given node. Answer is in parametric coordinates of element boundary
    //specific per type of elements used - needs to be overriden by child class
  //  virtual bool computeContactPoint(FloatArray &ksi, Node *node, StructuralElement *element, int elemedge, TimeStep *tStep) = 0;


  
};
} //end namespace oofem
