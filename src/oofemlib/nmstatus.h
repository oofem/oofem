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
// FILE: nmstatus.h
//

#ifndef nmstatus_h
#define nmstatus_h

/* mask defing NumMetod Status; which can be asked after
 * finishing computation by Numerical Method.
 * this mask should report some sitiuation, for
 * exaple:
 * None,
 * Success, NoSuccess  -> at now inform about succes,
 * if no succes, currently NoSuccess is not reported
 * but NoSuccess causes exit();
 *
 * some more detailed messages can be incorporated,
 * like:
 *
 * KeepTangent -> used by some non-linear solvers telling
 *              don't assemble new tangent, but use previous.
 *
 * note: every mask flag begins with NM_ to avoid possible multiple
 * definition
 */

typedef unsigned long NM_Status;
/* Mask selecting status */

#define NM_None         0
#define NM_Success      ( 1L << 1 )
#define NM_NoSuccess    ( 1L << 2 )
#define NM_KeepTangent  ( 1L << 3 )
#define NM_ForceRestart ( 1L << 4 )

#endif // nmstatus_h
