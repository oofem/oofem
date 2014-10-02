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

#include "smoothednodalintvarfield.h"
#include "spatiallocalizer.h"
#include "classfactory.h"
#include "element.h"
#include "feinterpol.h"
#include "dofmanager.h"
#include "set.h"

namespace oofem {
SmoothedNodalInternalVariableField :: SmoothedNodalInternalVariableField(InternalStateType ist, FieldType ft, NodalRecoveryModel :: NodalRecoveryModelType st, Domain *d) : Field(ft)
{
    this->istType = ist;
    this->stype = st;
    this->smoother = classFactory.createNodalRecoveryModel(this->stype, d);
    this->domain = d;
}

SmoothedNodalInternalVariableField :: ~SmoothedNodalInternalVariableField()
{
    if ( smoother ) {
        delete smoother;
    }
}

int
SmoothedNodalInternalVariableField :: evaluateAt(FloatArray &answer, FloatArray &coords, ValueModeType mode, TimeStep *tStep)
{
    int result = 0; // assume ok
    FloatArray lc, n;
    const FloatArray *nodalValue;

    // use whole domain recovery
    // create a new set containing all elements
    Set elemSet(0, this->domain);
    elemSet.addAllElements();
    this->smoother->recoverValues(elemSet, istType, tStep);
    // request element containing target point
    Element *elem = this->domain->giveSpatialLocalizer()->giveElementContainingPoint(coords);
    if ( elem ) { // ok element containing target point found
        FEInterpolation *interp = elem->giveInterpolation();
        if ( interp ) {
            // map target point to element local coordinates
            if ( interp->global2local( lc, coords, FEIElementGeometryWrapper(elem) ) ) {
                // evaluate interpolation functions at target point
                interp->evalN( n, lc, FEIElementGeometryWrapper(elem) );
                // loop over element nodes
                for ( int i = 1; i <= n.giveSize(); i++ ) {
                    // request nodal value
                    this->smoother->giveNodalVector( nodalValue, elem->giveDofManagerNumber(i) );
                    // multiply nodal value by value of corresponding shape function and add this to answer
                    answer.add(n.at(i), * nodalValue);
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
SmoothedNodalInternalVariableField :: evaluateAt(FloatArray &answer, DofManager *dman, ValueModeType mode, TimeStep *tStep)
{
    if ( dman->hasCoordinates() ) {
        const FloatArray *val;
        int result = this->smoother->giveNodalVector( val, dman->giveNumber() );
        answer = * val;
        return ( result == 1 );
    } else {
        return 1; // failed -> dman without coordinates
    }
}

contextIOResultType
SmoothedNodalInternalVariableField :: saveContext(DataStream &stream, ContextMode mode)
{
    return CIO_OK;
}

contextIOResultType
SmoothedNodalInternalVariableField :: restoreContext(DataStream &stream, ContextMode mode)
{
    return CIO_OK;
}
} // end namespace oofem
