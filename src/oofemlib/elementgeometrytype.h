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

#ifndef elementgeometrytype_h
#define elementgeometrytype_h

#include "enumitem.h"

namespace oofem {
#define Element_Geometry_Type_DEF \
    ENUM_ITEM(EGT_point) /* point in space */ \
    ENUM_ITEM(EGT_line_1) /* line elements with two nodes  1-------2 */ \
    ENUM_ITEM(EGT_line_2) /* line element with three nodes 1---3---2 */ \
    ENUM_ITEM(EGT_triangle_1) /* triangle element with three nodes */ \
    ENUM_ITEM(EGT_triangle_2) /* triangle element with 6 nodes */ \
    ENUM_ITEM(EGT_quad_1)   /* quadrialateral with 4 nodes */ \
    ENUM_ITEM(EGT_quad_1_interface)   /* same as quad_1 but with numbering according to interface element convention */ \
    ENUM_ITEM(EGT_quad_21_interface)   /* quadrialateral with 6 nodes (quadratic-linear) but with numbering according to interface element convention */ \
    ENUM_ITEM(EGT_quad_2)   /* quadratic quadrialateral with 8 nodes */ \
    ENUM_ITEM(EGT_quad9_2)       /* quadratic quadrialateral with 9 nodes */ \
    ENUM_ITEM(EGT_tetra_1)  /* tetrahedron with 4 nodes */ \
    ENUM_ITEM(EGT_tetra_2)  /* tetrahedron with 10 nodes */ \
    ENUM_ITEM(EGT_hexa_1)   /* hexahedron with 8 nodes */ \
    ENUM_ITEM(EGT_hexa_2)   /* hexahedron with 20 nodes */ \
    ENUM_ITEM(EGT_hexa_27)   /* hexahedron with 27 nodes */ \
    ENUM_ITEM(EGT_wedge_1)  /* wedge with 6 nodes */ \
    ENUM_ITEM(EGT_wedge_2)  /* wedge with 15 nodes */ \
    ENUM_ITEM(EGT_Composite) /* Composite geometry, vtk export supported by individual elements */ \
    ENUM_ITEM(EGT_unknown)  /* unknown element geometry type */

/**
 * Enumerative type used to classify element geometry
 * Possible values are:
 * EGT_point  - point in space
 * EGT_line_1 - line elements with two nodes  1-------2
 * EGT_line_2 - line element with three nodes 1---3---2
 * EGT_triangle_1 - triangle element with three nodes
 * EGT_triangle_2 - triangle element with 6 nodes
 *                     3
 *                  6     5
 *               1     4     2
 *
 * EGT_quad_1 - quadrialateral with 4 nodes
 * EGT_quad_2 - quadrialateral with 8 nodes
 * EGT_tetra_1 - tetrahedron with 4 nodes
 * EGT_tetra_2 - tetrahedron with 10 nodes
 * EGT_hexa_1  - hexahedron with 8 nodes
 * EGT_hexa_2  - hexahedron with 20 nodes
 * EGT_hexa_27  - hexahedron with 27 nodes
 */
enum Element_Geometry_Type {
    Element_Geometry_Type_DEF
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h


const char *__Element_Geometry_TypeToString(Element_Geometry_Type _value);
} // end namespace oofem
#endif // elementgeometrytype_h
