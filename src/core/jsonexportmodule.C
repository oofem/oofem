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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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


#include "jsonexportmodule.h"
#include "engngm.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "timestep.h"
#include "classfactory.h"

namespace oofem {

REGISTER_ExportModule(JsonExportModule)

#define _TR

JsonExportModule :: JsonExportModule(int n, EngngModel *e) : OutputExportModule(n, e)
{
    _TR;
}

void
JsonExportModule :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    _TR;
    this->initializeSilent=true;
    OutputExportModule::initializeFrom(ir);
    ctx=std::make_unique<JsonContext>(giveOutputStream(),ordered_json{});
    ctx->print({{"what","GLOBAL"},{"startTime",ctime(& emodel->giveStartTime())},{"description",emodel->giveDescription()}});
}


void
JsonExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
    _TR;
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) return;
    _TR;
    emodel->printOutputAt_json(*ctx, tStep, nodeSets, elementSets);
    ctx->print({{"what","GLOBAL"},{"step",tStep->giveNumber()},{"duration",emodel->giveSolutionStepTime()}});
}

void
JsonExportModule :: terminate()
{
    _TR;
    int rhrs, rmin, rsec, uhrs, umin, usec;
    time_t endTime = time(NULL);
    emodel->giveAnalysisTime(rhrs, rmin, rsec, uhrs, umin, usec);
    ctx->print({{"what","GLOBAL"},{"endTime",ctime(& endTime)},{"realTime",{{"hrs",rhrs},{"min",rmin},{"sec",rsec}}},{"userTime",{{"hrs",uhrs},{"min",umin},{"sec",usec}}}});
}

} // end namespace oofem
