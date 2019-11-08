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

#include "sm/Materials/LatticeMaterials/latticestructuralmaterial.h"
#include "domain.h"
#include "verbose.h"
#include "sm/Materials/structuralms.h"
#include "sm/Elements/structuralelement.h"
#include "sm/Elements/nlstructuralelement.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "floatmatrixf.h"
#include "floatarrayf.h"
#include "mathfem.h"
#include "engngm.h"
#include "fieldmanager.h"
#include "dynamicinputrecord.h"

namespace oofem {
LatticeStructuralMaterial :: LatticeStructuralMaterial(int n, Domain *d) : StructuralMaterial(n, d) { }


bool
LatticeStructuralMaterial :: hasMaterialModeCapability(MaterialMode mode) const
//
// returns whether receiver supports given mode
//
{
    return mode == _3dLattice || mode == _2dLattice || mode == _1dLattice;
}


double
LatticeStructuralMaterial :: giveLatticeStress1d(double strain, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("1dLattice mode not supported");
}

FloatArrayF< 3 >
LatticeStructuralMaterial :: giveLatticeStress2d(const FloatArrayF< 3 > &strain, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("2dLattice mode not supported");
}

FloatArrayF< 6 >
LatticeStructuralMaterial :: giveLatticeStress3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("3dLattice mode not supported");
}

FloatMatrixF< 1, 1 >
LatticeStructuralMaterial :: give1dLatticeStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
//
// return material stiffness matrix for 1dlattice
//
{
    OOFEM_ERROR("No general implementation provided");
}

FloatMatrixF< 3, 3 >
LatticeStructuralMaterial :: give2dLatticeStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
//
// return material stiffness matrix for 2dlattice
//
{
    OOFEM_ERROR("No general implementation provided");
}

FloatMatrixF< 6, 6 >
LatticeStructuralMaterial :: give3dLatticeStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
//
// return material stiffness matrix for 2dlattice
//
{
    OOFEM_ERROR("No general implementation provided");
}
} // end namespace oofem
