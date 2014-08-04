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
    if ( ! Structure.empty() ) {
        return 1;
    }

    int nnodes = Graph->giveDomain()->giveNumberOfDofManagers();
    IntArray nodalStatuses(nnodes);
    IntArray Level = {Root};

    // mark root
    nodalStatuses.at(Root) = 1;

    while ( !Level.isEmpty() ) { /* loop over levels */
        Structure.push_back(Level);
        /* start new level */
        /* loop over nodes on prev. level */
        Level.resize(0);
        int CurrLevelWidth = 0;
        for ( int inode: Structure.back() ) {
            for ( int n: Graph->giveNode( inode ).giveNeighborList() ) {
                if ( nodalStatuses.at(n) == 0 ) {
                    Level.followedBy(n, LEVEL_LIST_GROW_CHUNK);
                    nodalStatuses.at(n) = 1;
                    if ( ( limitWidth > 0 ) && ( ++CurrLevelWidth > limitWidth ) ) {
                        this->destroyLevels();
                        return 0; // zero mean aborted assembly
                    }
                }
            }
        }
    }

    return 1;
}

void
SloanLevelStructure :: computeDepth()
{
    this->formYourself();
    Depth = Structure.size();
}

void
SloanLevelStructure :: computeWidth()
{
    Width = 0;
    for ( int i = 1; i <= giveDepth(); i++ ) {
        int LevelWidth = giveLevel(i).giveSize();
        if ( Width < LevelWidth ) {
            Width = LevelWidth;
        }
    }
}

IntArray &
SloanLevelStructure :: giveLevel(int num)
{
    if ( Structure.empty() ) {
        this->formYourself();
    }

    if ( num < 1 || num > giveDepth() ) {
        OOFEM_ERROR("out of bounds (%d)", num);
    }

    return Structure [ num - 1 ];
}
} // end namespace oofem
