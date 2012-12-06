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

#ifndef dofiditemh
#define dofiditemh

#include "enumitem.h"

namespace oofem {
/**
 * Mask defining the physical meaning of particular DOF in node.
 * Mask arrays are also used in elements, where these arrays
 * are determining required DOFs needed by element and which are then
 * requested on particular nodes. They are used in function
 * Node::giveLocationArray which returns equations numbers
 * corresponding to selected dofs.
 * @deprecated{Use DofIDItem}
 */
typedef char DofID;

#define DofIDItem_DEF \
    ENUM_ITEM_WITH_VALUE(Undef, 0) /* Error value */ \
    ENUM_ITEM_WITH_VALUE(D_u, 1) /* u-displacement (in direction of x-axis) */ \
    ENUM_ITEM_WITH_VALUE(D_v, 2) /* v-displacement (in direction of y-axis) */ \
    ENUM_ITEM_WITH_VALUE(D_w, 3) /* w-displacement (in direction of z-axis) */ \
    ENUM_ITEM_WITH_VALUE(R_u, 4) /* Rotation around x-axis (right hand rule assumed) */ \
    ENUM_ITEM_WITH_VALUE(R_v, 5) /* Rotation around y-axis */ \
    ENUM_ITEM_WITH_VALUE(R_w, 6) /* Rotation around z-axis */ \
  \
    ENUM_ITEM_WITH_VALUE(V_u, 7) /* u-velocity (in direction of x-axis) */ \
    ENUM_ITEM_WITH_VALUE(V_v, 8) /* v-velocity (in direction of y-axis) */ \
    ENUM_ITEM_WITH_VALUE(V_w, 9) /* w-velocity (in direction of z-axis) */ \
  \
    ENUM_ITEM_WITH_VALUE(T_f, 10) /* Temperature field */ \
    ENUM_ITEM_WITH_VALUE(P_f, 11) /* Pressure field */ \
    ENUM_ITEM_WITH_VALUE(G_0, 12) /* DOF for gradient formulation no. 0 */ \
    ENUM_ITEM_WITH_VALUE(G_1, 13) /* DOF for gradient formulation no. 1 */ \
    ENUM_ITEM_WITH_VALUE(C_1, 14) /* Mass concentration of the first constituent */ \
    ENUM_ITEM_WITH_VALUE(X_1, 15) /* Start of xfemManager xfemdof pool */ \
    ENUM_ITEM_WITH_VALUE(X_N, 30) /* End of xfemManager xfemdof pool */ \
    ENUM_ITEM_WITH_VALUE(w_u, 31) /* u-component of change in director field (in direction of x-axis) */ \
	ENUM_ITEM_WITH_VALUE(w_v, 32) /* v-component of change in director field (in direction of y-axis) */ \
	ENUM_ITEM_WITH_VALUE(w_w, 33) /* w-component of change in director field (in direction of z-axis) */ \
	ENUM_ITEM_WITH_VALUE(gam, 34) /* inhomogenous thickness strain in direction of the directorfield m */ \
/**
 * Type representing particular dof type. Values of this type describe the physical meaning of
 * available DOFs.
 * @note{The implementation of Node::computeGNTransformation rely on D_u, D_v and D_w (R_u, R_v, R_w) order.
 * Do not change their order and do not insert any values between these values.}
 */
enum DofIDItem {
    DofIDItem_DEF
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h

const char *__DofIDItemToString(DofIDItem _value);

} // end namespace oofem
#endif // dofiditem_h
