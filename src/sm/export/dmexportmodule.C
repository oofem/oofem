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

#include "dmexportmodule.h"
#include "timestep.h"
#include "engngm.h"
#include "dofmanager.h"
#include "dof.h"
#include "classfactory.h"
#include "domain.h"

#include "../sm/EngineeringModels/structengngmodel.h" // JB

namespace oofem {
REGISTER_ExportModule(DofManExportModule)

DofManExportModule :: DofManExportModule(int n, EngngModel *e) : ExportModule(n, e)
{
    this->dofManList.clear();
}


DofManExportModule :: ~DofManExportModule()
{ }


IRResultType
DofManExportModule :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    // Read in dofMan's to export - defaults to all
    IR_GIVE_OPTIONAL_FIELD(ir, this->dofManList, _IFT_DofManExportModule_dmlist);

    return ExportModule :: initializeFrom(ir);
}


void
DofManExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }

    FILE *stream = this->giveOutputStream(tStep);
    fprintf(stream, "%% DofMan DataFile Version 1.1\n"); // Version 1.1 ??
    fprintf( stream, "%% Output for time %f\n", tStep->giveTargetTime() );

    DofManager *dm;
    double x, y, z, displacement;
    int domainIndex = 1;
    Domain *d  = emodel->giveDomain(domainIndex);
    FloatArray reactions;


    IntArray dofManMap, dofidMap, eqnMap;
    StructuralEngngModel *strEngMod = dynamic_cast< StructuralEngngModel * >(emodel);
    if ( strEngMod ) {
        strEngMod->buildReactionTable(dofManMap, dofidMap, eqnMap, tStep, domainIndex);
        strEngMod->computeReaction(reactions, tStep, 1);
    }



    int nTotaldm = d->giveNumberOfDofManagers();
    int ndmInList = this->dofManList.giveSize();

    int ndm = ( ndmInList > 0 ) ? ndmInList : nTotaldm;

    fprintf(stream, "%% Primary fields \n");
    for ( int idm = 1; idm <= ndm; idm++ ) {
        if ( ndm < nTotaldm ) {
            dm = d->giveDofManager( this->dofManList.at(idm) );
        } else {
            dm = d->giveDofManager(idm);
        }
        x = dm->giveCoordinate(1);
        y = dm->giveCoordinate(2);
        z = dm->giveCoordinate(3);
        fprintf(stream, "%d %g %g %g ", dm->giveNumber(), x, y, z);
        for ( Dof *dof: *dm ) {
            displacement = dof->giveUnknown(VM_Total, tStep);
            fprintf(stream, " %g", displacement);
        }
        fprintf(stream, "\n");
    }

#if 0
    fprintf(stream, "%% Reaction forces \n");
    for ( int idm = 1; idm <= ndm; idm++ ) {
        int dManNum = this->dofManList.at(idm);
        dm = d->giveDofManager(dManNum);
        fprintf(stream, "%d", dManNum);

        for ( int i = 1; i <= numRestrDofs; i++ ) {
            if ( dofManMap.at(i) == dManNum ) { // dofman has reaction
                //double reaction =  reactions.at( eqnMap.at(i) );
                double reaction =  reactions.at(i);
                fprintf(stream, " %g", reaction);
            }
        }
        fprintf(stream, "\n");
    }

    for ( int i = 1; i <= numRestrDofs; i++ ) {
        double reaction =  reactions.at(i);
        fprintf(stream, " %g", reaction);
    }
#endif
    fclose(stream);
}

FILE *
DofManExportModule :: giveOutputStream(TimeStep *tStep)
{
    FILE *answer;

    std :: string fileName = this->giveOutputBaseFileName(tStep) + ".dm";
    if ( ( answer = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR("failed to open file %s", fileName.c_str() );
    }

    return answer;
}
} // end namespace
