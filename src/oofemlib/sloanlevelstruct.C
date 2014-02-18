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

#include "sloanlevelstruct.h"
#include "sloangraphnode.h"
#include "sloangraph.h"
#include "intarray.h"
#include "domain.h"

namespace oofem {
#define LEVEL_LIST_GROW_CHUNK 50

SloanLevelStructure :: ~SloanLevelStructure()
{
    destroyLevels();
}

void
SloanLevelStructure :: destroyLevels()
{
    Structure.clear();
}

int
SloanLevelStructure :: formYourself(int limitWidth)
{
    if ( Structure.isNotEmpty() ) {
        return 1;
    }

    int nnodes = Graph->giveDomain()->giveNumberOfDofManagers();
    IntArray nodalStatuses(nnodes);
    IntArray workLevel;
    IntArray *Level;

    Level = new IntArray;
    Level->followedBy(Root);
    // mark root
    nodalStatuses.at(Root) = 1;

    IntArray *PrevLevel;
    std :: list< int > *Neighbor;
    std :: list< int > :: iterator pos2;
    SloanGraphNode *Node;
    int i, PrevLevelWidth, CurrLevelWidth;

    while ( !Level->isEmpty() ) { /* loop over levels */
        Structure.put(Structure.giveSize() + 1, Level);
        PrevLevel = Level;
        /* start new level */
        PrevLevelWidth = PrevLevel->giveSize();
        /* loop over nodes on prev. level */
        workLevel.clear();
        CurrLevelWidth = 0;
        for ( i = 1; i <= PrevLevelWidth; i++ ) {
            Node = Graph->giveNode( PrevLevel->at(i) );
            Neighbor = Node->giveNeighborList();
            for ( pos2 = Neighbor->begin(); pos2 != Neighbor->end(); ++pos2 ) {
                if ( nodalStatuses.at(* pos2) == 0 ) {
                    workLevel.followedBy(* pos2, LEVEL_LIST_GROW_CHUNK);
                    nodalStatuses.at(* pos2) = 1;
                    if ( ( limitWidth > 0 ) && ( ++CurrLevelWidth > limitWidth ) ) {
                        this->destroyLevels();
                        return 0; // zero mean aborted assembly
                    }
                }
            }
        }

        Level = new IntArray(workLevel);
    }

    delete Level;
    return 1;
}

void
SloanLevelStructure :: computeDepth()
{
    this->formYourself();
    Depth = Structure.giveSize();
}

void
SloanLevelStructure :: computeWidth()
{
    Width = 0;
    int i, LevelWidth;
    for ( i = 1; i <= giveDepth(); i++ ) {
        LevelWidth = giveLevel(i)->giveSize();
        if ( Width < LevelWidth ) {
            Width = LevelWidth;
        }
    }
}

IntArray *
SloanLevelStructure :: giveLevel(int num)
{
    if ( Structure.isEmpty() ) {
        this->formYourself();
    }

    if ( num < 1 || num > giveDepth() ) {
        OOFEM_WARNING2("LevelStructureClass::give_level - out of bounds (%d)", num);
        return NULL;
    }

    return Structure.at(num);
}
} // end namespace oofem
