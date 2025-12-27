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

#include "usrdefboundaryload.h"
#include "function.h"
#include "floatarray.h"
#include "timestep.h"
#include "dynamicinputrecord.h"
#include "valuemodetype.h"
#include "domain.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"

namespace oofem {
REGISTER_BoundaryCondition(UsrDefBoundaryLoad);

UsrDefBoundaryLoad :: UsrDefBoundaryLoad(int i, Domain * d) : BoundaryLoad(i, d) {
}

void
UsrDefBoundaryLoad :: computeComponentArrayAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    OOFEM_ERROR("not supported");
}


void
UsrDefBoundaryLoad :: computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode)
{
    // Evaluates the value at specific integration point
    Function* f = domain->giveFunction(this->intensityFunction);
    std :: map< std :: string, FunctionArgument > args = {
        {"x", coords}    
    };
    if ( mode == VM_Total ) {
        args["t"] = FunctionArgument(tStep->giveTargetTime());
        f->evaluate(answer, args);
        return;
    } else if (mode == VM_TotalIntrinsic) {
        args["t"] = FunctionArgument(tStep->giveIntrinsicTime());
        f->evaluate(answer, args);
        return;
    } else if (mode == VM_Incremental) {
        FloatArray answerPrev;
        args["t"] = FunctionArgument(tStep->giveTargetTime());
        f->evaluate(answer, args);
        args["t"] = FunctionArgument(tStep->giveTargetTime()-tStep->giveTimeIncrement());
        f->evaluate(answerPrev, args);
        answer.subtract(answerPrev);  
        return;
    } else {
        OOFEM_ERROR("unknown mode");
    }
}


void
UsrDefBoundaryLoad :: initializeFrom(InputRecord &ir)
{
    BoundaryLoad :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, intensityFunction, _IFT_UsrDefBoundaryLoad_intensityfunction);
    int mgt = (int) SurfaceLoadBGT;
    IR_GIVE_OPTIONAL_FIELD(ir, mgt, _IFT_UsrDefBoundaryLoad_GeomType);
    myGeomType = (bcGeomType) mgt;
    IR_GIVE_OPTIONAL_FIELD(ir, approxOrder, _IFT_UsrDefBoundaryLoad_approxorder);
}


void
UsrDefBoundaryLoad :: giveInputRecord(DynamicInputRecord &input)
{
    BoundaryLoad :: giveInputRecord(input);
    input.setField(this->intensityFunction, _IFT_UsrDefBoundaryLoad_intensityfunction);
    input.setField((int)this->myGeomType, _IFT_UsrDefBoundaryLoad_GeomType);
    input.setField(this->approxOrder, _IFT_UsrDefBoundaryLoad_approxorder);
}


void
UsrDefBoundaryLoad :: saveContext(DataStream &stream, ContextMode mode)
{
    BoundaryLoad :: saveContext(stream, mode);

    if ( mode & CM_Definition ) {
        if ( !stream.write(intensityFunction) ) {
          THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(myGeomType) ) {
          THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(approxOrder) ) {
          THROW_CIOERR(CIO_IOERR);
        }
    }
}


void
UsrDefBoundaryLoad :: restoreContext(DataStream &stream, ContextMode mode)
{
    BoundaryLoad :: restoreContext(stream, mode);

    if ( mode & CM_Definition ) {
        int _val;
        if ( !stream.read(intensityFunction) ) {
          THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(_val) ) {
          THROW_CIOERR(CIO_IOERR);
        }
        myGeomType = (bcGeomType) _val;

        if ( !stream.read(approxOrder) ) {
          THROW_CIOERR(CIO_IOERR);
        }
    }
}

} // end namespace oofem
