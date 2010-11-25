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

//
// FILE: internalstatevaluetype.h
//

#ifndef internalstatevaluetype_h
#define internalstatevaluetype_h

namespace oofem {
/// enum determining the type of internal variable
enum InternalStateValueType {
    ISVT_UNDEFINED,
    ISVT_SCALAR,
    ISVT_VECTOR,
    ISVT_TENSOR_S3, // symmetric 3x3 tensor
    ISVT_TENSOR_S3E, // symmetric 3x3 tensor, packed with off diagonal components multiplied by 2
                     // (engineering strain vector, for example)
    ISVT_TENSOR_G // general tensor
};
} // end namespace oofem
#endif // internalstatevaluetype_h

