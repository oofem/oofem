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

#include "engngm.h"
#include "classfactory.h"
#include "util.h"
#include "inputrecord.h"
#include "datareader.h"

#include <cstring>

namespace oofem {

EngngModel *InstanciateProblem(DataReader *dr, problemMode mode, int contextFlag, EngngModel *_master, bool parallelFlag)
{
    const char *__proc = "InstanciateProblem"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                       // Required by IR_GIVE_FIELD macro
    EngngModel *problem;
    std::string problemName, dataOutputFileName, desc;

    dataOutputFileName = dr->giveOutputFileName();
    desc = dr->giveDescription();

    /* here we need copy of input record. The pointer returned by dr->giveInputRecord can (and will)
     * be updated as reading e-model components (nodes, etc). But we need this record being available
     * through the whole e-model instanciation
     */
    InputRecord *emodelir = dr->giveInputRecord(DataReader :: IR_emodelRec, 1)->GiveCopy();
    result = emodelir->giveRecordKeywordField(problemName); ///@todo Make this function robust, it can't be allowed to fail (the record keyword is not a normal field-id)
    if ( result != IRRT_OK ) {
        IR_IOERR("", __proc, "", emodelir, result);
    }

    problem = classFactory.createEngngModel(problemName.c_str(), 1, _master);
    if (!problem) {
        OOFEM_ERROR2("EngngModel::InstanciateProblem - Failed to construct engineering model if type \"%s\".\n", problemName.c_str());
    }
    problem->setProblemMode(mode);
    problem->setParallelMode(parallelFlag);

    if ( contextFlag ) {
        problem->setContextOutputMode(COM_Always);
    }

    problem->instanciateYourself(dr, emodelir, dataOutputFileName.c_str(), desc.c_str());

    delete emodelir;

    return problem;
}

} // end namespace oofem
