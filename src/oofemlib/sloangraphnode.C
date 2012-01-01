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
 *               Copyright (C) 1993 - 2011   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Modified and optimized by: Borek Patzak */
/* Author: Milan Jirasek */

#include "sloangraphnode.h"
#include "sloangraph.h"
#include "mathfem.h"

namespace oofem {
SloanGraphNode :: SloanGraphNode(SloanGraph *graph, int numOld) : neighborList()
{
    this->graph = graph;
    NumberOld  = numOld;
    NumberNew  = 0;
    nodeStatus = Inactive;
    Degree    = 0;
    Distance  = -1;
    Priority  = -1;
}

SloanGraphNode :: ~SloanGraphNode()
{
}

void SloanGraphNode :: addNeighbor(int neighbor)
{
    // check if neighbor already in list
    dynaList< int > :: iterator pos;
    for ( pos = this->neighborList.begin(); pos != this->neighborList.end(); ++pos ) {
        if ( * pos == neighbor ) {
            return;
        }
    }

    Degree++;
    this->neighborList.pushFront(neighbor);
}

int SloanGraphNode :: computeProfileHeight()
{
    int numberMin = NumberNew;
    dynaList< int > :: iterator pos;
    for ( pos = this->neighborList.begin(); pos != this->neighborList.end(); ++pos ) {
        numberMin = min( numberMin, this->graph->giveNode(* pos)->giveNewNumber() );
    }

    return ( NumberNew - numberMin + 1 );
}
} // end namespace oofem
