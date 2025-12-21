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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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






#include "interpolationcatalogue.h"

namespace oofem {

FEI1dLin InterpolationCatalogueType::fei1dlin_x = FEI1dLin(1);
FEI1dLin InterpolationCatalogueType::fei1dlin_y = FEI1dLin(2);
FEI1dLin InterpolationCatalogueType::fei1dlin_z = FEI1dLin(3);
#ifdef __MPM_MODULE
ConstantInterpolation InterpolationCatalogueType::feiconst = ConstantInterpolation();
LinearInterpolation InterpolationCatalogueType::feilin = LinearInterpolation();
QuadraticInterpolation InterpolationCatalogueType::feiquad = QuadraticInterpolation();
#endif

const FEInterpolation* InterpolationCatalogueType :: getInterpolationByName (std::string name) {
    if ( name == "fei1dlin_x" ) {
        return &this->fei1dlin_x;
    } else if (name == "fei1dlin_y") {
        return &fei1dlin_y;
    } else if (name == "fei1dlin_z") {
        return &fei1dlin_z;
#ifdef __MPM_MODULE
    } else if (name == "feiconst") {
        return &feiconst;
    } else if (name == "feilin") {
        return &feilin;
    } else if (name == "feiquad") {
        return &feiquad;
#endif
    } else {
        OOFEM_ERROR("Interpolation %s not found in catalogue", name.c_str());
    }
}


InterpolationCatalogueType interpolationCatalogue;

} // end namespace oofem
