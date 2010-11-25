/* $Header: /home/cvs/bp/oofem/oofemlib/src/structuralmaterial.C,v 1.19.4.1 2004/04/05 15:19:44 bp Exp $ */
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


//   file STRUCTURALMATERIAL.CC

#include "fluiddynamicmaterial.h"
#include "domain.h"
#include "verbose.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"

#include "mathfem.h"

#include "engngm.h"
#include "fieldmanager.h"

#ifndef __MAKEDEPEND
 #include <stdlib.h>
 #include <math.h>
#endif
#include "contextioerr.h"

namespace oofem {
void
FluidDynamicMaterial :: updateInternalState(const FloatArray &vec, GaussPoint *gp, TimeStep *)
{ }


FluidDynamicMaterialStatus :: FluidDynamicMaterialStatus(int n, Domain *d, GaussPoint *g) :
    MaterialStatus(n, d, g), deviatoricStressVector()
{ }

void
FluidDynamicMaterialStatus :: printOutputAt(FILE *File, TimeStep *tNow)
// Prints the strains and stresses on the data file.
{
    int i, n;

    fprintf(File, "\n deviatoric stresses");
    n = deviatoricStressVector.giveSize();
    for ( i = 1; i <= n; i++ ) {
        fprintf( File, " % .4e", deviatoricStressVector.at(i) );
    }

    fprintf(File, "\n");
}

void
FluidDynamicMaterialStatus :: updateYourself(TimeStep *tStep)
// Performs end-of-step updates.
{
    MaterialStatus :: updateYourself(tStep);
}


void
FluidDynamicMaterialStatus :: initTempStatus()
//
// initialize record at the begining of new load step
//
{
    MaterialStatus :: initTempStatus();
}


contextIOResultType
FluidDynamicMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full ms context (saves state variables, that completely describe
// current state)
// saving the data in  TDictionary is left to material (yield crit. level).
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        _error("saveContex : can't write into NULL stream");
    }

    if ( ( iores = MaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = deviatoricStressVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
FluidDynamicMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full material context (saves state variables, that completely describe
// current state)
//
{
    // FloatArray *s;
    contextIOResultType iores;
    if ( stream == NULL ) {
        _error("saveContex : can't write into NULL stream");
    }

    if ( ( iores = MaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = deviatoricStressVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
} // end namespace oofem
