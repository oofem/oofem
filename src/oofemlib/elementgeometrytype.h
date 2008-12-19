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
// FILE: elementgeometrytype.h
//

#ifndef elementgeometrytype_h
#define elementgeometrytype_h

#include "enumitem.h"

/**
 * Enumerative type used to classify element geometry
 * Poosible values are:
 * EGT_line_1 - line elements with two nodes  1-------2
 * EGT_line_2 - line element with three nodes 1---2---3
 * EGT_triangle_1 - triangle element with three nodes
 * EGT_triangle_2 - triangle element with 6 nodes
 *                     3
 *                  6     5
 *               1     4     2
 *
 * EGT_quad_1 - quadrialateral with 4 nodes
 * EGT_tetra_1 - tetrahedron with 4 nodes
 * EGT_hexa_1  - hexahedron with 8 nodes
 * EGT_hexa_2  - hexahedron with 20 nodes
 */
#define Element_Geometry_Type_DEF \
    ENUM_ITEM(EGT_line_1) /* line elements with two nodes  1-------2 */   \
    ENUM_ITEM(EGT_line_2) /* line element with three nodes 1---2---3 */   \
    ENUM_ITEM(EGT_triangle_1) /* triangle element with three nodes */ \
    ENUM_ITEM(EGT_triangle_2) /* triangle element with 6 nodes */ \
    ENUM_ITEM(EGT_quad_1)   /* quadrialateral with 4 nodes */   \
    ENUM_ITEM(EGT_tetra_1)  /* tetrahedron with 4 nodes */   \
    ENUM_ITEM(EGT_hexa_1)   /* hexahedron with 8 nodes */   \
    ENUM_ITEM(EGT_hexa_2)   /* hexahedron with 20 nodes */   \
    ENUM_ITEM(EGT_unknown)  /* unknown element geometry type */

enum Element_Geometry_Type {
    Element_Geometry_Type_DEF
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h


const char *__Element_Geometry_TypeToString(Element_Geometry_Type _value);


#endif // elementgeometrytype_h

