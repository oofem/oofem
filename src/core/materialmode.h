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

#ifndef materialmode_h
#define materialmode_h


#include "meta_enum.hpp"

namespace oofem {
meta_enum(MaterialMode,int,
    _Unknown,
    _3dMat,
    _3dMatGrad,  /* 3d model with gradient of internal variable */
    _1dMatGrad, /* 1d model with gradient of internal variable */
    _PlaneStressGrad,  /* plane stress with gradient of internal variable */
    _PlaneStrainGrad,  /* plane strain with gradient of internal variable */
    _PlaneStress,
    _PlaneStrain,
    _2dPlate,
    _2dPlateSubSoil, /* Subsoil mode for plates */
    _1dMat,
    _2dBeam,
    _3dBeam,
    _3dShell,
    _3dShellRot,
    _3dDegeneratedShell,
    _3dBeamSubSoil, /* Subsoil model for beams */

    _PlateLayer,
    _2dBeamLayer,
    _PlaneStressRot, /* Plane stress with rotation around z */

    _Fiber,
    _3dMicroplane,
    _3dInterface,
    _2dInterface,
    _1dInterface,

    _1dHeat, /* 1d heat or 1d mass*/
    _1dHeMo, /* 1d heat and mass (one component, transfer */
    _2dHeat, /* 2d heat or 2d mass */
    _2dHeMo, /* 2d heat and mass (one component, transfer */
    _3dHeat, /* 3d heat or 3d mass */
    _3dHeMo, /* 3d heat and mass (one component, transfer */

    _2dFlow,
    _2dAxiFlow,
    _3dFlow,

    _2dUP,
    _3dUP,
    _2dUPV,
    _3dUPV,

    _1dLattice,
    _2dLattice,
    _3dLattice,
    _2dMTLattice,
    _3dMTLattice,
    _Warping
)

// constexpr auto __MaterialModeToString=MaterialMode_value_to_string;
#define __MaterialModeToString(v) std::string(MaterialMode_value_to_string(v)).c_str()
} // end namespace oofem
#endif // materialmode_h
