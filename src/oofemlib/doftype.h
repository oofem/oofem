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
// FILE: doftype.h
//

#ifndef doftype_h
#define doftype_h

#include "enumitem.h"

namespace oofem {

/* Dof Type, determines the type of DOF created
 */
#define dofType_DEF \
    ENUM_ITEM_WITH_VALUE(DT_master, 0) \
    ENUM_ITEM_WITH_VALUE(DT_simpleSlave, 1) \
    ENUM_ITEM_WITH_VALUE(DT_slave, 2)

enum dofType {
    dofType_DEF
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h


const char *__dofTypeToString(dofType _value);

} // end namespace oofem
#endif // doftype_h
