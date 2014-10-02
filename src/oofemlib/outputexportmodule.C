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

#include <vector>
#include <cstdio>

#include "outputexportmodule.h"
#include "engngm.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "timestep.h"
#include "classfactory.h"

namespace oofem {
REGISTER_ExportModule(OutputExportModule)

OutputExportModule :: OutputExportModule(int n, EngngModel *e) : ExportModule(n, e)
{
}

IRResultType
OutputExportModule :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    nodeSets.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, nodeSets, _IFT_OutputExportModule_nodeSets);
    elementSets.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, elementSets, _IFT_OutputExportModule_elementSets);

    return ExportModule :: initializeFrom(ir);
}

void
OutputExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }

    if ( domainMask.containsOnlyZeroes() && domainMask.giveSize() > 0 ) {
        return;
    }

    FILE *file = emodel->giveOutputStream();

    fprintf(file, "\n==============================================================");
    fprintf(file, "\nOutput for time % .8e ", tStep->giveTargetTime() * emodel->giveVariableScale(VST_Time) );
    fprintf(file, "\n==============================================================\n");

    for ( int idomain = 1; idomain <= emodel->giveNumberOfDomains(); idomain++ ) {
        if ( domainMask.giveSize() > 0 && domainMask.at(idomain) == 0 ) {
            continue;
        }
        fprintf(file, "Output for domain %3d\n", idomain );

        Domain *domain = emodel->giveDomain(idomain);
        this->doDofManOutput(file, domain, tStep);
        this->doElementOutput(file, domain, tStep);
        ///@todo This module should handle printing reaction forces, we need more information from the domain/engineering model though.
        //emodel->printReactionForces(tStep, idomain);
    }
}

void
OutputExportModule :: doDofManOutput(FILE *file, Domain *domain, TimeStep *tStep)
{
    int ndofman = domain->giveNumberOfDofManagers();

    fprintf(file, "\n\nDofManager output:\n------------------\n");

    int setNum = nodeSets.isEmpty() ? 0 : nodeSets.at(domain->giveNumber()) ;

    if ( nodeSets.isEmpty() ) {
        for ( int i = 1; i <= ndofman; i++ ) {
            DofManager *dman = domain->giveDofManager(i);

            if ( dman->giveParallelMode() == DofManager_null ) {
                continue;
            }

            dman->printOutputAt(file, tStep);
        }
    } else {
        Set *set = domain->giveSet(setNum);
        const IntArray &nodes = set->giveNodeList();

        for ( int i = 1; i <= ndofman; i++ ) {
            DofManager *dman = domain->giveDofManager(i);

            if ( domain->giveDofManager(i)->giveParallelMode() == DofManager_null ) {
                continue;
            }

            if ( nodes.containsSorted(i) ) {
                dman->printOutputAt(file, tStep);
            }
        }
    }

    fprintf(file, "\n\n");
}

void
OutputExportModule :: doElementOutput(FILE *file, Domain *domain, TimeStep *tStep)
{
    int nelem = domain->giveNumberOfElements();

    fprintf(file, "\n\nElement output:\n---------------\n");

    int setNum = elementSets.isEmpty() ? 0 : elementSets.at(domain->giveNumber()) ;

    if ( setNum == 0 ) {
        for ( int i = 1; i <= nelem; i++ ) {
            Element *element = domain->giveElement(i);

            if ( element->giveParallelMode() == Element_remote ) {
                continue;
            }

            element->printOutputAt(file, tStep);
        }
    } else {
        Set *set = domain->giveSet(setNum);
        IntArray elements = set->giveElementList();
        elements.sort();
        for ( int i = 1; i <= nelem; i++ ) {
            Element *element = domain->giveElement(i);

            if ( element->giveParallelMode() == Element_remote ) {
                continue;
            }

            if ( elements.containsSorted(i) ) {
                element->printOutputAt(file, tStep);
            }
        }
    }

    fprintf(file, "\n\n");
}

} // end namespace oofem
