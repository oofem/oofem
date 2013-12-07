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

#include "rvematerial.h"
#include "oofemtxtdatareader.h"
#include "util.h"
#include "inputrecord.h"

#include <cstdio>
//#include <unistd.h>

namespace oofem {
IRResultType RVEMaterial :: initializeFrom(InputRecord *ir)
{
    printf("rvematerial initializing...\n");

    const char *__proc = "initializeFrom";
    IRResultType result;

    // Read filename and instanciate RVE problem

    IR_GIVE_FIELD(ir, this->rveFilename, _IFT_RVEMaterial_fileName);

    OOFEM_LOG_INFO( "************************** Instanciating microproblem from file %s\n", rveFilename.c_str() );
    OOFEMTXTDataReader drMicro( rveFilename.c_str() );

    this->rve = InstanciateProblem(& drMicro, _processor, 0);
    drMicro.finish();

    this->rve->setProblemScale(microScale);
    this->rve->setProblemScale(microScale);
    this->rve->checkProblemConsistency();
    this->rve->initMetaStepAttributes( this->rve->giveMetaStep(1) );
    this->rve->giveNextStep();
    this->rve->init();

    OOFEM_LOG_INFO("************************** Microproblem at %p instanciated \n", rve);

    SupressRVEoutput = 0;

    IR_GIVE_OPTIONAL_FIELD(ir, SupressRVEoutput, _IFT_RVEMaterial_supressoutput);

    return IRRT_OK;
}

void
RVEMaterial :: suppressStdout()
{
    //    if (SupressRVEoutput) {
    //        fgetpos(stdout, &stdoutPos);
    //        stdoutFID=dup(fileno(stdout));
    //        freopen(this->rveLogFilename.c_str(), "a", stdout);
    //    }
}

void RVEMaterial :: enableStdout()
{
    //    if (SupressRVEoutput) {
    //        fflush(stdout);
    //        dup2(stdoutFID, fileno(stdout));
    //        close (stdoutFID);
    //        clearerr(stdout);
    //        fsetpos(stdout, &stdoutPos);        /* for C9X */
    //    }
}
}
