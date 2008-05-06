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
// FILE: dofmantransftype.h
//

#ifndef dofmantransftype_h
#define dofmantransftype_h

/**
 * Enumerative type, used to specify type of trasformation required from dofManager (node).
 * The _toGlobalCS value requires transformation from node-depenedent coordinate system
 * to gbal coordinte system in node to be assebled. Then global vector fg can be obtained by
 * followwing operation fg = T fn, where T is transformation matrix and fn is vector expressed in
 * nodal coordinate system).
 * The _toNodalCS value represent transformation from global c.s in node to node-dependent
 * coordinate system.
 */
enum DofManTransfType {
    _toGlobalCS,
    _toNodalCS
};

#endif // dofmantransftype_h
