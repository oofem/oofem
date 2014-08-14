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

#include "../sm/FETISolver/fetiboundarydofman.h"
#include "error.h"

using namespace std;


namespace oofem {
FETIBoundaryDofManager :: FETIBoundaryDofManager() : partitions(), codeNumbers()
{
    globalNumber = 0;
    numberOfPartitions = 0;
    ndofs = 0;
}

FETIBoundaryDofManager :: FETIBoundaryDofManager(int num, int part, int ndof) : partitions(part), codeNumbers(part * ndof)
{
    globalNumber = num;
    numberOfPartitions = part;
    ndofs = ndof;
}



FETIBoundaryDofManager :: FETIBoundaryDofManager(const FETIBoundaryDofManager &src)
{
    globalNumber = src.globalNumber;
    numberOfPartitions = src.numberOfPartitions;
    ndofs = src.ndofs;
    referencePartition = src.referencePartition;
    partitions = src.partitions;
    codeNumbers = src.codeNumbers;
}


void
FETIBoundaryDofManager :: addPartition(int partitionNumber)
{
    IntArray partitionToAdd(1);

    partitionToAdd.at(1) = partitionNumber;
    partitions.followedBy(partitionToAdd);
    if ( numberOfPartitions != 0 ) {
        if ( partitionNumber < referencePartition ) {
            referencePartition = partitionNumber;
        }
    } else {
        referencePartition = partitionNumber;
    }

    numberOfPartitions++;
}

int
FETIBoundaryDofManager :: setCodeNumbers(int &equationCounter)
{
    int i, size;

    this->codeNumbers.resize( ( size = ( numberOfPartitions - 1 ) * ndofs ) );
    for ( i = 1; i <= size; i++ ) {
        codeNumbers.at(i) = ++equationCounter;
    }

    return equationCounter;
}


int
FETIBoundaryDofManager :: giveCodeNumber(int partition_num, int dof_num)
{
    int indx = 0, i;
    for ( i = 1; i <= numberOfPartitions; i++ ) {
        if ( partitions.at(i) != referencePartition ) {
            indx++;
        }

        if ( partitions.at(i) == partition_num ) {
            break;
        }
    }

    if ( ( indx == 0 ) || ( partition_num == referencePartition ) ) {
        return 0;
    }

    if ( ( dof_num < 1 ) || ( dof_num > ndofs ) ) {
        OOFEM_ERROR("bad dof_num requested");
    }

    return codeNumbers.at( ( indx - 1 ) * ndofs + dof_num );
}

int
FETIBoundaryDofManager :: giveCompleteLocationArray(int rank, IntArray &locationArray)
{
    int indx = 0, i;
    for ( i = 1; i <= numberOfPartitions; i++ ) {
        if ( partitions.at(i) != referencePartition ) {
            indx++;
        }

        if ( partitions.at(i) == rank ) {
            break;
        }
    }

    if ( ( indx == 0 ) || ( rank == referencePartition ) ) {
        return 0;
    }

    locationArray.resize(ndofs);
    for ( i = 1; i <= ndofs; i++ ) {
        locationArray.at(i) = codeNumbers.at( ( indx - 1 ) * ndofs + i );
    }

    return 1;
}
} // end namespace oofem
