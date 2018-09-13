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

#ifndef anisolinearelasticmaterial_h
#define anisolinearelasticmaterial_h

#include "linearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matconst.h"
#include "element.h"

///@name Input fields for AnisotropicLinearElasticMaterial
//@{
#define _IFT_AnisotropicLinearElasticMaterial_Name "anisole"
/// Stiffness coefficients arranged by rows from the diagonal to the right (21 values)
#define _IFT_AnisotropicLinearElasticMaterial_stiff "stiff"
/// Thermal expansion, 6 components in strain-Voigt order,
#define _IFT_AnisotropicLinearElasticMaterial_talpha "talpha"
//@}

namespace oofem {
class GaussPoint;

/**
 * This class implements a general anisotropic linear elastic material in a finite
 * element problem.
 */
class AnisotropicLinearElasticMaterial : public LinearElasticMaterial
{
public:
    AnisotropicLinearElasticMaterial(int n, Domain *d) : LinearElasticMaterial(n, d) {}
    virtual ~AnisotropicLinearElasticMaterial() {}

    IRResultType initializeFrom(InputRecord *ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    const char *giveInputRecordName() const override { return _IFT_AnisotropicLinearElasticMaterial_Name; }
    const char *giveClassName() const override { return "AnisotropicLinearElasticMaterial"; }
};
} // end namespace oofem
#endif // anisolinearelasticmaterial_h
