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

#include "problemsequence.h"
#include "inputrecord.h"
#include "classfactory.h"
#include "util.h"
#include "inputrecord.h"
#include "datastream.h"
#include "oofemtxtdatareader.h"

#include <memory>

namespace oofem {

REGISTER_EngngModel(ProblemSequence)

ProblemSequence :: ProblemSequence(int i, EngngModel *master): EngngModel(i, master) {}

ProblemSequence :: ~ProblemSequence() {}


void ProblemSequence :: solveYourself()
{
    for ( auto &emodel : emodelList ) {
        emodel->solveYourself();
        ///@todo Still lacking the all important code to connect the subsequent analysis!
        // Options:
        // 1. Use initial conditions (for both primary and internal fields!)
        // 2. Copy domains over to the subsequent analysis
        // 2 might be easier, but 1 has advantages outside of this type of analysis.
    }
}


int ProblemSequence :: instanciateYourself(DataReader &dr, InputRecord &ir, const char *outFileName, const char *desc)
{
    int result = EngngModel :: instanciateYourself(dr, ir, dataOutputFileName.c_str(), desc);
    ir.finish();

    for ( auto &s : inputStreamNames ) {
        OOFEMTXTDataReader dr( inputStreamNames [ i - 1 ] );
        std :: unique_ptr< EngngModel >prob( InstanciateProblem(dr, this->pMode, this->contextOutputMode, NULL) );
        emodelList.emplace_back(std :: move(prob));
        if ( prob ) {
            return 0;
        }
    }

    return result;
}


void ProblemSequence :: initializeFrom(InputRecord &ir)
{
    EngngModel :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, inputStreamNames, _IFT_ProblemSequence_engineeringModels);
}


int ProblemSequence :: checkProblemConsistency()
{
    int ret = 1;
    for (auto &emodel : emodelList) {
        ret *= emodel->checkConsistency();
    }
    return ret;
}


void ProblemSequence :: saveContext(DataStream &stream, ContextMode mode)
{
    EngngModel :: saveContext(stream, mode);

    stream.write(activeModel);

    for ( auto &emodel : emodelList ) {
        emodel->saveContext(stream, mode);
    }
}


void ProblemSequence :: restoreContext(DataStream &stream, ContextMode mode)
{
    EngngModel :: restoreContext(stream, mode);

    stream.read(activeModel);

    for ( auto &emodel : emodelList ) {
        emodel->restoreContext(stream, mode);
    }
}


#ifdef __OOFEG
void ProblemSequence :: drawYourself(oofegGraphicContext &gc) {}
void ProblemSequence :: drawElements(oofegGraphicContext &gc) {}
void ProblemSequence :: drawNodes(oofegGraphicContext &gc) {}

#endif

}
