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

#ifndef nmstatus_h
#define nmstatus_h

namespace oofem {
/**
 * Mask defining NumMetod Status; which can be asked after
 * finishing computation by Numerical Method.
 * this mask should report some situation.
 */
typedef unsigned long NM_Status;

#define NM_None         0
#define NM_Success      ( 1L << 1 ) ///< Numerical method exited with success.
#define NM_NoSuccess    ( 1L << 2 ) ///< Numerical method failed to solve problem.
#define NM_KeepTangent  ( 1L << 3 ) ///< Don't assemble new tangent, but use previous.
#define NM_ForceRestart ( 1L << 4 )
} // end namespace oofem
#endif // nmstatus_h
