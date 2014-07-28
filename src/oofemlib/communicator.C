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

#include "communicator.h"
#include "intarray.h"

#include <cstdarg>

#ifdef __USE_MPI
 #include <mpi.h>
#endif

namespace oofem {
CommunicatorBuff :: CommunicatorBuff(int s, CommBuffType t)
{
    this->size = s;

    if ( size ) {
        processCommBuffs = new ProcessCommunicatorBuff * [ size ];
        for ( int i = 0; i < size; i++ ) {
            processCommBuffs [ i ] = new ProcessCommunicatorBuff(t);
        }
    } else {
        processCommBuffs = NULL;
    }
}

CommunicatorBuff :: ~CommunicatorBuff()
{
    for ( int i = 0; i < size; i++ ) {
        if ( processCommBuffs [ i ] ) {
            delete(processCommBuffs [ i ]);
        }
    }

    if ( processCommBuffs ) {
        delete[] processCommBuffs;
    }
}

Communicator :: Communicator(EngngModel *emodel, CommunicatorBuff *b, int rank, int size, CommunicatorMode m)
{
    this->engngModel = emodel;
    this->rank = rank;
    this->size = size;
    this->mode = m;

    if ( size ) {
        processComms = new ProcessCommunicator * [ size ];
        for ( int i = 0; i < size; i++ ) {
            processComms [ i ] =
                new ProcessCommunicator(b->giveProcessCommunicatorBuff ( i ), i, mode);
        }
    } else {
        processComms = NULL;
    }
}

Communicator :: ~Communicator()
{
    int i = size;

    if ( size ) {
        while ( i-- ) {
            delete(processComms [ i ]);
        }

        delete[]  processComms;
    }
}

int
Communicator :: initExchange(int tag)
{
    int result = 1;
    for  ( int i = 0; i < size; i++ ) {
        result &= this->giveProcessCommunicator(i)->initExchange(tag);
    }

    return result;
}

int
Communicator :: finishExchange()
{
    int result = 1;
    for  ( int i = 0; i < size; i++ ) {
        result &= this->giveProcessCommunicator(i)->finishExchange();
    }

    return result;
}



int
Communicator :: initSend(int tag)
{
    int result = 1;
    for  ( int i = 0; i < size; i++ ) {
        result &= this->giveProcessCommunicator(i)->initSend(tag);
    }

    return result;
}

int
Communicator :: initReceive(int tag)
{
    int result = 1;
    for  ( int i = 0; i < size; i++ ) {
        result &= this->giveProcessCommunicator(i)->initReceive(tag);
    }

    return result;
}

void
Communicator :: clearBuffers()
{
    for  ( int i = 0; i < size; i++ ) {
        this->giveProcessCommunicator(i)->clearBuffers();
    }
}

std :: string
Communicator :: errorInfo(const char *func) const
{
    return std::string("Communicator::") + func + ", Rank: " + std::to_string(rank);
}
} // end namespace oofem
