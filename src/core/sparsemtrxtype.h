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

#ifndef sparsematrixtype_h
#define sparsematrixtype_h

namespace oofem {
/**
 * Enumerative type used to identify the sparse matrix type.
 */
enum SparseMtrxType {
    SMT_Skyline,       ///< Symmetric skyline.
    SMT_SkylineU,      ///< Unsymmetric skyline.
    SMT_CompCol,       ///< Compressed column.
    SMT_DynCompCol,    ///< Dynamically growing compressed column.
    SMT_SymCompCol,    ///< Symmetric compressed column.
    SMT_DynCompRow,    ///< Dynamically growing compressed row.
    SMT_SpoolesMtrx,   ///< Spooles sparse mtrx representation.
    SMT_PetscMtrx,     ///< PETSc library mtrx representation.
    SMT_DSS_sym_LDL,   ///< Richard Vondracek's sparse direct solver.
    SMT_DSS_sym_LL,    ///< Richard Vondracek's sparse direct solver.
    SMT_DSS_unsym_LU   ///< Richard Vondracek's sparse direct solver.
};
} // end namespace oofem
#endif // sparsematrixtype_h
