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

#include "load.h"
#include "verbose.h"
#include "timestep.h"
#include "function.h"
#include "dynamicinputrecord.h"

namespace oofem {
Load :: Load(int i, Domain *aDomain) :
    GeneralBoundaryCondition(i, aDomain), componentArray(), dofExcludeMask()
{
    timeFunction = 0;
}


const FloatArray &
Load :: giveComponentArray() const
{
    return componentArray;
}


void
Load :: computeValues(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, const IntArray &dofids, ValueModeType mode)
{
    ///@note Backwards compatibility with input files that don't specify dofs.
#if 1
    if ( this->dofs.giveSize() == 0 ) {
        this->computeValueAt(answer, tStep, coords, mode);
        return;
    }
#endif

    FloatArray loaded_dofs;
    this->computeValueAt(loaded_dofs, tStep, coords, mode);

    answer.resize(dofids.giveSize());
    answer.zero();
    for ( int i = 0; i < dofids.giveSize(); ++i ) {
        int index = this->dofs.findFirstIndexOf(dofids[i]);
        if ( index > 0 ) {
            answer[i] = loaded_dofs.at(index);
        }
    }
}


void
Load :: computeComponentArrayAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
// Returns an array, the load induced at tStep by the receiver.
{
    double factor;

    factor = this->giveTimeFunction()->evaluate(tStep, mode);
    answer = componentArray;
    answer.times(factor);

    if ( !isImposed(tStep) ) {
        answer.zero();
    }
}


IRResultType
Load :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

#  ifdef VERBOSE
    // VERBOSE_PRINT1 ("Instanciating load ",number)
#  endif

    IR_GIVE_FIELD(ir, componentArray, _IFT_Load_components);

    int size = componentArray.giveSize();
    dofExcludeMask.resize(size);
    dofExcludeMask.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, dofExcludeMask, _IFT_Load_dofexcludemask);
    if ( dofExcludeMask.giveSize() != size ) {
        OOFEM_WARNING("dofExcludeMask and componentArray size mismatch");
        return IRRT_BAD_FORMAT;
    } else {
        for ( int i = 1; i <= size; i++ ) {
            if ( dofExcludeMask.at(i) ) {
                componentArray.at(i) = 0.0;
            }
        }
    }

    return GeneralBoundaryCondition :: initializeFrom(ir);
}


void Load :: giveInputRecord(DynamicInputRecord &input)
{
    GeneralBoundaryCondition :: giveInputRecord(input);
    input.setField(this->componentArray, _IFT_Load_components);
    if ( !this->dofExcludeMask.containsOnlyZeroes() ) {
        input.setField(this->dofExcludeMask, _IFT_Load_dofexcludemask);
    }
}


int
Load :: isDofExcluded(int indx)
{
    if ( ( indx > 0 ) && ( indx <= dofExcludeMask.giveSize() ) ) {
        return dofExcludeMask.at(indx);
    } else {
        OOFEM_ERROR("dof index out of range");
    }

    return 0;
}
} // end namespace oofem
