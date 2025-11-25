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

#ifndef matesponsemode_h
#define matesponsemode_h

#include "meta_enum.hpp"

namespace oofem {
meta_enum(MatResponseMode,int,
    TangentStiffness=0,
    SecantStiffness=1,
    ElasticStiffness=2,
    Stress=3,
    Conductivity=4,  /* element level conductivity matrix */
    Conductivity_ww=5, /* material level conductivity submatrix */
    Conductivity_hh=6, /* material level conductivity submatrix */
    Conductivity_hw=7, /* material level conductivity submatrix */
    Conductivity_wh=8, /* material level conductivity submatrix */
    Capacity=9,
    Capacity_ww=10, /* material level capacity submatrix */
    Capacity_hh=11, /* material level capacity submatrix */
    Capacity_hw=12, /* material level capacity submatrix */
    Capacity_wh=13, /* material level capacity submatrix */
    IntSource=14,
    IntSource_ww=15, /* material level internal source submatrix - water source */
    IntSource_hh=16, /*  - heat source */
    IntSource_hw=17, /*  - heat source dependency on water content change */
    IntSource_wh=18, /*  - water source dependency on temperature change */
    Permeability=19,
    FluidMassBalancePressureContribution=20,
    BiotConstant=21,
    CompressibilityCoefficient=22,
    FluidViscosity=23,
    Flux=24,
    DSigmaDT=25,
    ElasticBulkModulus=26,
    ElasticBulkModulusInverse=27,
    MRM_ScalarOne=28,
    DeviatoricStiffness=29,
    DeviatoricStress=30
);

// constexpr auto __MatResponseModeToString=MatResponseMode_value_to_string;
#define __MatResponseModeToString(v) std::string(MatResponseMode_value_to_string(v)).c_str()
} // end namespace oofem
#endif // matesponsemode_h
