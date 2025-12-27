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

#ifndef expressionfield_h
#define expressionfield_h

#include "field.h"
#include "floatarray.h"
#include "intarray.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "dofmanager.h"
#include "error.h"
#include <iostream>


///@name Input fields for PythonField
//@{
#define _IFT_ExpressionField_Name "exprfield"
//@}

namespace oofem {
/**
 * Custom user supplied expression for field.
 * The expression is evaluated at given point (x,y,z) and time (t).
 * @code{.py}
 * evaluateAt(FloatArray &answer, const FloatArray &coords, ValueModeType mode, TimeStep *tStep) # returns int 
 * @endcode
 * 
 */
class OOFEM_EXPORT ExpressionField : public Field
{
private:
    /// Name of python module containing evaluating function function
    std::string expression;

public:
//     Constructor.
    ExpressionField(void);
   
    int evaluateAt(FloatArray &answer, const FloatArray &coords,
                           ValueModeType mode, TimeStep *tStep) override;
                           
                           
    int evaluateAt(FloatArray &answer, DofManager *dman,
                           ValueModeType mode, TimeStep *tStep) override;
                           
    void saveContext(DataStream &stream) override { }
    void restoreContext(DataStream &stream) override { }
    
    const char *giveClassName() const override { return "ExpressionField"; }

    void initializeFrom(InputRecord &ir) override;
};

} // end namespace oofem
#endif // expressionfield_h
