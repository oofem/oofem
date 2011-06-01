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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef dofmantransftype_h
#define dofmantransftype_h

namespace oofem {
/**
 * Enumerative type, used to specify type of transformation required from dofManager (node).
 * Then global vector @f$ f_g @f$ can be obtained by following operation
 * @f$ f_g = T\cdot f_n@f$, where @f$T@f$ is transformation matrix and @f$f_n@f$ is vector expressed in
 * nodal coordinate system.
 */
enum DofManTransfType {
    _toGlobalCS, ///< Transformation from global c.s in node to node-dependent coordinate system.
    _toNodalCS, ///< Transformation from node-dependent coordinate system to global coordinate system in node.
};
} // end namespace oofem
#endif // dofmantransftype_h
