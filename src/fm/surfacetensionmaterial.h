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

#ifndef surfacetensionmaterial_h
#define surfacetensionmaterial_h

#include "material.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "matconst.h"
#include "structuralelement.h"
#include "matstatus.h"

///@name Input fields for SurfaceTensionMaterial
//@{
#define _IFT_SurfaceTensionMaterial_isotropic "isotropic"
//@}

namespace oofem {
class GaussPoint;

/**
 * Currently unused for isotropic surface tension.
 */
class SurfaceTensionMaterialStatus : public MaterialStatus
{
public:
    /// Constructor
    SurfaceTensionMaterialStatus(int n, Domain *d, GaussPoint *g) : MaterialStatus(n, d, g) { }
    /// Destructor
    virtual ~SurfaceTensionMaterialStatus() { }

    virtual const char *giveClassName() const { return "SurfaceTensionMaterialStatus"; }
    virtual classType giveClassID() const { return SurfaceTensionMaterialStatusClass; }
};

/**
 * Material model for surface tension.
 * Currently the model only stores the surface tension energy.
 * @deprecated This material will be removed.
 */
class SurfaceTensionMaterial : public Material
{
public:
    /// Constructor
    SurfaceTensionMaterial(int n, Domain *d) : Material(n, d) { }
    /// Destructor.
    virtual ~SurfaceTensionMaterial() { }

    /**
     * Initializes the material.
     * - g Isotropic surface tension (gamma) (optional, default 0.0)
     */
    virtual IRResultType initializeFrom(InputRecord *ir)
    {
        const char *__proc = "initializeFrom";
        IRResultType result;
        double value = 0.0;
        IR_GIVE_OPTIONAL_FIELD(ir, value, IFT_SurfaceTensionMaterial_isotropic, "g"); // Macro
        this->propertyDictionary->add('g', value);
        return result;
    }

    virtual const char *giveClassName() const { return "SurfaceTensionMaterial"; }
    virtual classType giveClassID() const { return SurfaceTensionMaterialClass; }
};
} // end namespace oofem
#endif // surfacetensionmaterial_h
