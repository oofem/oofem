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

/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

#include "timestep.h"
#include "engngm.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {

TimeStep :: TimeStep(int n, EngngModel *e, int mn, double tt, double dt, StateCounterType counter, TimeDiscretizationType td) :
    eModel(e), targetTime(tt), intrinsicTime(tt), deltaT(dt), solutionStateCounter(counter), 
    number(n), version(0), mstepNumber(mn), timeDiscretization(td)
{
    // Target time and intrinsic time is the same in the constructor.
}

TimeStep :: TimeStep(EngngModel *e)
{
    eModel = e;
    deltaT = 0.0;
    targetTime = 0.0;
    intrinsicTime = 0.0;
    solutionStateCounter = 0;
    number = -1;
    version = 0;
    mstepNumber = 0;
}

TimeStep :: TimeStep(const TimeStep &src)
{
    eModel = src.eModel;
    targetTime = src.targetTime;
    intrinsicTime = src.intrinsicTime;
    deltaT = src.deltaT;
    solutionStateCounter = src.solutionStateCounter;
    number = src.number;
    version = src.version;
    mstepNumber = src.mstepNumber;
}

TimeStep &
TimeStep :: operator = ( const TimeStep & src )
{
    eModel = src.eModel;
    targetTime = src.targetTime;
    intrinsicTime = src.intrinsicTime;
    deltaT = src.deltaT;
    solutionStateCounter = src.solutionStateCounter;
    number = src.number;
    version = src.version;
    mstepNumber = src.mstepNumber;

    return * this;
}




TimeStep *TimeStep :: givePreviousStep()
// Not accepted in-line.
{
    if ( isTheCurrentTimeStep() ) {
        return eModel->givePreviousStep();
    } else {
        OOFEM_ERROR("TimeStep::givePreviousStep Could not return previous step of noncurrent step");
    }

    return NULL; // to make compiler happy
}


bool TimeStep :: isNotTheLastStep()
// Returns True if the time history contains steps after the receiver,
// else returns False.
{
    return  ( number != eModel->giveNumberOfSteps() );
}

bool TimeStep :: isTheFirstStep()
{
    // Returns True if the receiver is the first time step,
    // according to first step number
    // else returns False.

    return  ( number == eModel->giveNumberOfFirstStep() );
}


bool TimeStep :: isIcApply()
{
    // Returns True if the receiver is the  time step,
    // when Initial conditions apply
    // else returns False.

    return  ( number == eModel->giveNumberOfTimeStepWhenIcApply() );
}



bool TimeStep :: isTheCurrentTimeStep()
// Not accepted in-line.
{
    return this == eModel->giveCurrentStep();
}


contextIOResultType
TimeStep :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    // write step number
    if ( !stream->write(& number, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // write meta step number
    if ( !stream->write(& mstepNumber, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // write target time
    if ( !stream->write(& this->targetTime, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // write intrinsic time
    if ( !stream->write(& this->intrinsicTime, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // write deltaT
    if ( !stream->write(& this->deltaT, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // write solutionStateCounter
    if ( !stream->write(& this->solutionStateCounter, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // write timeDiscretization
    int tDiscretization = (int) timeDiscretization;
    if ( !stream->write(&tDiscretization, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // return result back
    return CIO_OK;
}

contextIOResultType
TimeStep :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    // read step number
    if ( !stream->read(& number, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // read meta step number
    if ( !stream->read(& mstepNumber, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // read target time
    if ( !stream->read(& this->targetTime, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // read intrinsic time
    if ( !stream->read(& this->intrinsicTime, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // read deltaT
    if ( !stream->read(& this->deltaT, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // read solutionStateCounter
    if ( !stream->read(& this->solutionStateCounter, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // read timeDiscretization
    int tDiscretization = 0;
    if ( !stream->read(& tDiscretization, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    timeDiscretization = ( TimeDiscretizationType ) tDiscretization;

    // return result back
    return CIO_OK;
}
} // end namespace oofem
