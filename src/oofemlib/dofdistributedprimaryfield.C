/* $Header: /home/cvs/bp/oofem/oofemlib/src/primaryfield.C,v 1.2.4.1 2004/04/05 15:19:43 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

#include "dofdistributedprimaryfield.h"
#include "spatiallocalizer.h"
#include "dofmanager.h"
#include "dof.h"
#include "element.h"
#include "timestep.h"

namespace oofem {
DofDistributedPrimaryField :: DofDistributedPrimaryField(EngngModel *a, int idomain,
                                                         FieldType ft, EquationID ut, int nHist) :
    PrimaryField(a, idomain, ft, ut, nHist)
{ }

DofDistributedPrimaryField :: ~DofDistributedPrimaryField()
{ }

double
DofDistributedPrimaryField :: giveUnknownValue(Dof *dof, ValueModeType mode, TimeStep *atTime)
{
    return dof->giveUnknown(this->ut, mode, atTime);
}


FloatArray *
DofDistributedPrimaryField :: giveSolutionVector(TimeStep *atTime)
{
    _error("giveSolutionVector: not supported");
    return NULL;
}



void
DofDistributedPrimaryField :: advanceSolution(TimeStep *atTime)
{ }


contextIOResultType
DofDistributedPrimaryField :: saveContext(DataStream *stream, ContextMode mode)
{
    // all the job is done by dofs alone
    return CIO_OK;
}

contextIOResultType
DofDistributedPrimaryField :: restoreContext(DataStream *stream, ContextMode mode)
{
    return CIO_OK;
}
} // end namespace oofem
