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

#ifdef __PARALLEL_MODE

#include "processcomm.h"
#include "intarray.h"

#include "combuff.h"
#include "dyncombuff.h"

#ifdef __USE_MPI
 #include <mpi.h>
#endif

namespace oofem {
ProcessCommunicatorBuff :: ProcessCommunicatorBuff(CommBuffType t)
{
    if ( t == CBT_static ) {
        send_buff = new StaticCommunicationBuffer(MPI_COMM_WORLD);
        recv_buff = new StaticCommunicationBuffer(MPI_COMM_WORLD);
    } else {
        send_buff = new DynamicCommunicationBuffer(MPI_COMM_WORLD);
        recv_buff = new DynamicCommunicationBuffer(MPI_COMM_WORLD);
    }
}


ProcessCommunicator :: ProcessCommunicator(ProcessCommunicatorBuff *b, int rank, CommunicatorMode m) :
    toSend(), toReceive()
{
    this->rank = rank;
    this->pcBuffer =  b;
    this->mode = m;
}


int
ProcessCommunicator :: initSend(int tag)
{
    int result = 1;
    if ( !toSend.isEmpty() || ( this->mode == CommMode_Dynamic ) ) {
        //  fprintf (stderr, "\nNlDEIDynamicDomainComunicator :: initExchange: sending to %d",rank);
        result = giveProcessCommunicatorBuff()->initSend(this->rank, tag);
    } else {
        giveProcessCommunicatorBuff()->initSendBuff();
    }

    return result;
}


int
ProcessCommunicator :: initReceive(int tag)
{
    int result = 1;
    if ( !toReceive.isEmpty() || ( this->mode == CommMode_Dynamic ) ) {
        //  fprintf (stderr, "\nNlDEIDynamicDomainComunicator :: initExchange: recv from %d",rank);
        result &= giveProcessCommunicatorBuff()->initReceive(this->rank, tag);
    } else {
        giveProcessCommunicatorBuff()->initRecvBuff();
    }

    return result;
}


int
ProcessCommunicator :: initExchange(int tag)
{
    int result = 1;
    result &= initSend(tag);
    result &= initReceive(tag);

    return result;
}

int
ProcessCommunicator :: finishExchange()
{
    return waitCompletion();
}

void
ProcessCommunicator :: clearBuffers()
{
    giveProcessCommunicatorBuff()->init();
}

int
ProcessCommunicator :: sendCompleted()
{
    if ( !toSend.isEmpty() || ( this->mode == CommMode_Dynamic ) ) {
        return giveProcessCommunicatorBuff()->sendCompleted();
    } else {
        return 1;
    }
}

int
ProcessCommunicator :: receiveCompleted()
{
    if ( !toReceive.isEmpty() || ( this->mode == CommMode_Dynamic ) ) {
        return giveProcessCommunicatorBuff()->receiveCompleted();
    } else {
        return 1;
    }
}

int
ProcessCommunicator :: testCompletion()
{
    return ( sendCompleted() && receiveCompleted() );
}

int
ProcessCommunicator :: waitCompletion()
{
    while ( !testCompletion() ) {
        ;
    }

    return 1;
}
} // end namespace oofem
#endif
