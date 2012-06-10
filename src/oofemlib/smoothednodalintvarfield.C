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

#include "smoothednodalintvarfield.h"
#include "spatiallocalizer.h"
#include "usrdefsub.h"
#include "element.h"

namespace oofem {
SmoothedNodalInternalVariableField :: SmoothedNodalInternalVariableField(InternalStateType ist, FieldType ft, NodalRecoveryModel::NodalRecoveryModelType st, Domain *d) : Field(ft)
{
    this->istType = ist;
    this->stype = st;
    this->smoother = CreateUsrDefNodalRecoveryModel(this->stype, d);
    // setup smoother (whole domain recovery)
    this->smoother->setRecoveryMode(0, IntArray());
    this->domain = d;
}

SmoothedNodalInternalVariableField :: ~SmoothedNodalInternalVariableField()
{
    if ( smoother ) {
        delete smoother;
    }
}

int
SmoothedNodalInternalVariableField :: evaluateAt(FloatArray &answer, FloatArray &coords, ValueModeType mode, TimeStep *atTime)
{
    int i, region = 1, result = 0; // assume ok
    FloatArray lc, n;
    const FloatArray *nodalValue;

    this->smoother->recoverValues(istType, atTime);
    // request element containing target point
    Element* elem = this->domain->giveSpatialLocalizer()->giveElementContainingPoint(coords);
    if (elem) { // ok element containing target point found
        FEInterpolation* interp = elem->giveInterpolation();
        if (interp) {
            // map target point to element local coordinates
            if (interp->global2local(lc, coords, FEIElementGeometryWrapper(elem))) {
                // evaluate interpolation functions at target point
                interp->evalN(n, lc, FEIElementGeometryWrapper(elem));
                // loop over element nodes
                for (i=1; i<=n.giveSize(); i++) {
                    // request nodal value
                    this->smoother->giveNodalVector(nodalValue, elem->giveDofManagerNumber(i), region);
                    // multiply nodal value by value of corresponding shape function and add this to answer
                    answer.add(n.at(i), *nodalValue);
                }
            } else { // mapping from global to local coordinates failed
                result = 1; // failed
            }
        } else {  // element without interpolation
            result = 1; // failed
        }
    } else { // no element containing given point found
        result = 1; // failed
    }
    return result;
}

int
SmoothedNodalInternalVariableField::evaluateAt(FloatArray &answer, DofManager* dman, ValueModeType mode, TimeStep *atTime)
{
    int region = 1;
    if (dman->hasCoordinates()) {
        const FloatArray* val;
        int result = this->smoother->giveNodalVector(val, dman->giveNumber(), region);
        answer=*val;
        return (result == 1);

    } else {
        return 1; // failed -> dman without coordinates
    }
}

contextIOResultType
SmoothedNodalInternalVariableField :: saveContext(DataStream *stream, ContextMode mode)
{
    // int i, type_id = SmoothedNodalInternalVariableFieldClass;
    // contextIOResultType iores;
    // write class header
    // if (!stream->write(&type_id,1)) return CIO_IOERR;

    return CIO_OK;
}

contextIOResultType
SmoothedNodalInternalVariableField :: restoreContext(DataStream *stream, ContextMode mode)
{
    // int i, class_id;
    // contextIOResultType iores;
    // read class header
    // if (!stream->read(&class_id,1)) return CIO_IOERR;
    // if (class_id != SmoothedNodalInternalVariableField) return CIO_BADVERSION;

    return CIO_OK;
}
} // end namespace oofem
