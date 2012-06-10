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

#include "dmexportmodule.h"
#include "timestep.h"
#include "engngm.h"
#include "dofmanager.h"
#include "dof.h"

namespace oofem {
DofManExportModule :: DofManExportModule(int n, EngngModel *e) : ExportModule(n, e)
{}


DofManExportModule :: ~DofManExportModule()
{}


IRResultType
DofManExportModule :: initializeFrom(InputRecord *ir)
{
    ExportModule :: initializeFrom(ir);
    return IRRT_OK;
}


void
DofManExportModule :: doOutput(TimeStep *tStep)
{
    if ( !testTimeStepOutput(tStep) ) {
        return;
    }

    FILE *stream = this->giveOutputStream(tStep);
    fprintf(stream, "# DofMan DataFile Version 1.1\n");
    fprintf( stream, "Output for time %f\n", tStep->giveTargetTime() );

    DofManager *dm;
    double x, y, z, displacement;
    int idof, idm, ndofs;
    Dof *dof;

    Domain *d  = emodel->giveDomain(1);
    int ndm = d->giveNumberOfDofManagers();

    for ( idm = 1; idm <= ndm; idm++ ) {
        dm = d->giveDofManager(idm);
        x = dm->giveCoordinate(1);
        y = dm->giveCoordinate(2);
        z = dm->giveCoordinate(3);
        fprintf(stream, "%d %g %g %g ", dm->giveNumber(), x, y, z);
        ndofs = dm->giveNumberOfDofs();
        for ( idof = 1; idof <= ndofs; idof++ ) {
            dof = dm->giveDof(idof);
            displacement = dof->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
            fprintf(stream, " %g", displacement);
        }

        fprintf(stream, "\n");
    }

    fclose(stream);
}

FILE *
DofManExportModule :: giveOutputStream(TimeStep *tStep)
{
    FILE *answer;

    std::string fileName = this->giveOutputBaseFileName(tStep) + ".dm";
    if ( ( answer = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR2("DofManExportModule::giveOutputStream: failed to open file %s", fileName.c_str());
    }

    return answer;
}
} // end namespace
