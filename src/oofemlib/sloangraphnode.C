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
{ }

void SloanGraphNode :: addNeighbor(int newNeighbor)
{
    // check if neighbor already in list
    for ( int neighbor: this->neighborList ) {
        if ( neighbor == newNeighbor ) {
            return;
        }
    }

    Degree++;
    this->neighborList.push_front(newNeighbor);
}

int SloanGraphNode :: computeProfileHeight()
{
    int numberMin = NumberNew;
    for ( int neighbor: this->neighborList ) {
        numberMin = min( numberMin, this->graph->giveNode(neighbor).giveNewNumber() );
    }

    return ( NumberNew - numberMin + 1 );
}
} // end namespace oofem
