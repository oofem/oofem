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

#ifndef interpolationcatalogue_h
#define interpolationcatalogue_h

#include "fei1dlin.h"
#ifdef __MPM_MODULE
#include "prototype2.h"
#endif

namespace oofem {

/**
 * This class represent a repository of interpolation class instances. As the interpolation classes can be shared between elements, it makes a sence to create them only once (here) and reuse them in elements.
 * There is a global instance of this class (created in interpolationcatalogue.C) , which is used by elements to obtain interpolation classes, named interpolationCatalogue.
 * 
 */
class InterpolationCatalogueType {
public:
    static  FEI1dLin fei1dlin_x ;
    static  FEI1dLin fei1dlin_y ;
    static  FEI1dLin fei1dlin_z ;
#ifdef __MPM_MODULE
    static  ConstantInterpolation feiconst;
    static  LinearInterpolation feilin;
    static  QuadraticInterpolation feiquad;
#endif
    const FEInterpolation* getInterpolationByName (std::string name);

};

extern InterpolationCatalogueType interpolationCatalogue;

} // end namespace oofem
#endif // interpolationcatalogue_h
