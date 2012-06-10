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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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
#include "usrdefsub.h"
#include "util.h"
#include "inputrecord.h"
#include "datareader.h"

#ifndef __MAKEDEPEND
 #include <cstring>
 #include <iostream>
#endif

namespace oofem {
char *giveLineFromInput(FILE *inputStream, char *line, int len)
//
// reads one line from inputStream - for private use only.
//
{
    char *ptr;

    giveRawLineFromInput(inputStream, line, len);
    // convert line to lowercase
    for ( ptr = line; ( * ptr = tolower(* ptr) ); ptr++ ) {
        ;
    }

    return line;
}

char *giveRawLineFromInput(FILE *inputStream, char *line, int len)
//
// reads one line from inputStream - for private use only.
//
{
    char *_res;
    do {
        _res = fgets(line, len, inputStream);
        if ( _res == NULL ) {
            OOFEM_ERROR("giveRawLineFromInput : End of file encountered");
        }
    } while ( * line == '#' ); // skip comments

    return line;
}


void giveInputDataFileName(std::string &dataInputFileName)

{
    // Returns the name of the file containing the data of the problem.
    printf("Please enter the name of the input data file : \n");
    std::getline(std::cin, dataInputFileName);
}


EngngModel *InstanciateProblem(DataReader *dr, problemMode mode, int contextFlag, EngngModel *_master, bool parallelFlag)
{
    const char *__keyword, *__proc = "InstanciateProblem"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                                 // Required by IR_GIVE_FIELD macro
    EngngModel *problem;
    std::string problemName, dataOutputFileName, desc;

    InputRecord *ir = dr->giveInputRecord(DataReader :: IR_outFileRec, 1);
    __keyword = NULL;
    result = ir->giveField(dataOutputFileName, IFT_EngngModel_outfile, __keyword);
    if ( result != IRRT_OK ) {
        IR_IOERR("", __proc, IFT_EngngModel_outfile, "Output file record", ir, result);
    }

    ir->finish();

    ir = dr->giveInputRecord(DataReader :: IR_jobRec, 1);
    __keyword = NULL;
    result = ir->giveField(desc, IFT_EngngModel_probdescription, __keyword);

    /* here we need copy of input record. The pointer returned by dr->giveInputRecord can (and will)
     * be updated as reading e-model components (nodes, etc). But we need this record being available
     * through the whole e-model instanciation
     */
    InputRecord *emodelir = dr->giveInputRecord(DataReader :: IR_emodelRec, 1)->GiveCopy();
    result = emodelir->giveRecordKeywordField(problemName);
    //result = IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);
    if ( result != IRRT_OK ) {
        IR_IOERR("", __proc, IFT_EngngModel_probname, "", emodelir, result);
    }

    problem = CreateUsrDefEngngModelOfType(problemName.c_str(), 1, _master);
    if (!problem) {
        OOFEM_ERROR("EngngModel::InstanciateProblem - Failed to construct engineering model.\n");
    }
    problem->setProblemMode(mode);
    problem->setParallelMode(parallelFlag);

    if ( contextFlag ) {
        problem->setContextOutputMode(COM_Always);
    }

    problem->instanciateYourself(dr, emodelir, dataOutputFileName.c_str(), desc.c_str());
    //emodelir.finish();
    delete(emodelir);

    return problem;
}

#define oofem_tmpstring_len 1024
static char oofem_tmpstring [ oofem_tmpstring_len + 1 ];

char *oofem_tmpstr(const char *src)
{
    strncpy(oofem_tmpstring, src, oofem_tmpstring_len);
    oofem_tmpstring [ oofem_tmpstring_len ] = '\0';
    return oofem_tmpstring;
}
} // end namespace oofem
