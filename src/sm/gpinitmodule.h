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

#ifndef gpinitmodule_h
#define gpinitmodule_h

#include "initmodule.h"

#define _IFT_GPInitModule_Name "gpinitmodule"

namespace oofem {
/**
 * Represents GP (Gauss point) initialization module.
 * This module reads certain internal variables
 * of all Gauss points from a file.
 * In this way, one can specify e.g. initial damage
 * and initial stresses computed by another model.
 */
class GPInitModule : public InitModule
{
public:
    /// Creates empty GPInitModule.
    GPInitModule(int n, EngngModel * e);

    /// Destructor
    virtual ~GPInitModule();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void doInit();
    virtual const char *giveClassName() const { return "GPInitModule"; }
};
} // namespace oofem
#endif
