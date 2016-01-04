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

#ifndef linsystsolvertype_h
#define linsystsolvertype_h

namespace oofem {
/**
 * The values of this type should be related not to specific solvers,
 * but more to specific packages that provide linear solver interface
 * (possibly with many solver types) and are represented by a class
 * derived from SparseLinearSystemNM.
 * The selection of particular solver from package should be done using keywords,
 * related to particular package.
 */
enum LinSystSolverType {
    ST_Direct = 0,
    ST_IML    = 1,
    ST_Spooles= 2,
    ST_Petsc  = 3,
    ST_DSS    = 4,
    ST_Feti   = 5,
    ST_MKLPardiso = 6
};
} // end namespace oofem
#endif // linsystsolvertype_h
