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

#ifndef bcgeomtype_h
#define bcgeomtype_h

namespace oofem {
/// Type representing the geometric character of loading.
enum bcGeomType {
    UnknownBGT,     ///< Unknown type.
    NodalLoadBGT,   ///< Concentrated nodal load.
    BodyLoadBGT,    ///< Distributed body load.
    EdgeLoadBGT,    ///< Distributed edge load.
    SurfaceLoadBGT, ///< Distributed surface load.
    PointLoadBGT,   ///< Concentrated point load (placed anywhere).
    GravityPressureBGT, ///<Pressure due to distributed body load.
};
} // end namespace oofem
#endif // bcgeomtype_h
