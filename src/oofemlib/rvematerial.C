/*
 * RVEMaterial.C
 *
 *  Created on: Mar 26, 2010
 *      Author: carl
 */

#include "rvematerial.h"
#include "oofemtxtdatareader.h"
#include "util.h"
#include "inputrecord.h"

#include <cstdio>
#include <unistd.h>

namespace oofem {

IRResultType RVEMaterial :: initializeFrom(InputRecord *ir)
{
	printf("rvematerial initializing...\n");

    const char *__proc = "initializeFrom";  // Required by IR_GIVE_FIELD macro
    IRResultType result;                    // Required by IR_GIVE_FIELD macro

    // Read filename and instanciate RVE problem

    IR_GIVE_FIELD(ir, this->rveFilename, IFT_MicroMaterialFileName, "file");

    OOFEM_LOG_INFO("************************** Instanciating microproblem from file %s\n", rveFilename.c_str());
    OOFEMTXTDataReader drMicro(rveFilename.c_str());

    this->rve = InstanciateProblem(& drMicro, _processor, 0);
    drMicro.finish();

    this->rve->setProblemScale(microScale);
    this->rve->setProblemScale(microScale);
    this->rve->checkProblemConsistency();
    this->rve->initMetaStepAttributes( this->rve->giveMetaStep( 1 ) );
    this->rve->giveNextStep();
    this->rve->init();

    OOFEM_LOG_INFO("************************** Microproblem at %p instanciated \n", rve);

    /* Get BCType which is the type of boundary condition imposed on the rve
     *
     */
    //IR_GIVE_FIELD(ir, BCType, IFT_RVEMaterial_bctype, "bctype");

    SupressRVEoutput=0;

    IR_GIVE_OPTIONAL_FIELD(ir, SupressRVEoutput, IFT_RVEMaterial_supressoutput, "supressoutput"); // Macro

    return IRRT_OK;
}

void
RVEMaterial :: suppressStdout()
{
    if (SupressRVEoutput) {
        fgetpos(stdout, &stdoutPos);
        stdoutFID=dup(fileno(stdout));
        freopen(this->rveLogFilename.c_str(), "a", stdout);
    }
}

void RVEMaterial :: enableStdout()
{
    if (SupressRVEoutput) {
        fflush(stdout);
        dup2(stdoutFID, fileno(stdout));
        close (stdoutFID);
        clearerr(stdout);
        fsetpos(stdout, &stdoutPos);        /* for C9X */
    }
}

}
