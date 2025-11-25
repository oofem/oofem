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

#ifndef valuemodetype_h
#define valuemodetype_h

#include "meta_enum.hpp"

namespace oofem {
//
// following mode determines the mode of particular unknown
// which can be requested on DOF.
// particular DOF contain for example displacement type unknown,
// but we can request total value, increment of value or velocity or acceleration of
// this unknown. This has been done mainly in  order to  improve runtime checking
// ad Dof level.
//
// see also isUnknownModeIncrementalMode() function (cltypes.C)
// when adding new ValueModeType mode.
/*
  Note:
  VM_Total:   total value evaluated at the end of solution (time) step
  VM_TotalIntrinsic: total value evaluated at intrinsic time 
*/

/**
 * Type representing the mode of UnknownType or CharType, or similar types.
 * Afore mentioned types usually describes the physical meaning of
 * value and ValueModeType provides the further necessary classification. For example "DisplacementVector"
 * value can be further classified to be total displacement (TotalMode) or  velocity of
 * displacement (VelocityMode) an so on.
 */
meta_enum(ValueModeType,int,
    VM_Unknown=0,
    VM_Total=1,
    VM_Velocity=2,
    VM_Acceleration=3,
    VM_Incremental=4,
    VM_RhsTotal=5,
    VM_RhsIncremental=6,
    VM_RhsInitial=7,
    VM_Intermediate=8,
    VM_TotalIntrinsic=9,
    VM_Residual=99
);

constexpr auto __ValueModeTypeToString=ValueModeType_value_to_string;
} // end namespace oofem
#endif // valuemodetype_h
