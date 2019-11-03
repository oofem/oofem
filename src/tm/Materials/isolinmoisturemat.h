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

#ifndef isolinmoisturemat_h
#define isolinmoisturemat_h

#include "tm/Materials/isomoisturemat.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for IsotropicLinMoistureTransferMaterial
//@{
#define _IFT_IsotropicLinMoistureTransferMaterial_Name "isolinmoisturemat"
#define _IFT_IsotropicLinMoistureTransferMaterial_perm "perm" ///< Moisture permeability
#define _IFT_IsotropicLinMoistureTransferMaterial_capa "capa" ///< Moisture capacity
//@}

namespace oofem {
/**
 * This class implements a isotropic moisture tranport material. A material
 * is an attribute of a domain. It is usually also attribute of many elements.
 */
class IsotropicLinMoistureTransferMaterial : public IsotropicMoistureTransferMaterial
{
protected:
    double moistureCapacity = 0.;
    double permeability = 0.;

public:
    IsotropicLinMoistureTransferMaterial(int n, Domain * d) : IsotropicMoistureTransferMaterial(n, d) { }

    void initializeFrom(InputRecord &ir) override;
    double givePermeability(GaussPoint *gp, TimeStep *tStep) const override;
    double giveMoistureCapacity(GaussPoint *gp, TimeStep *tStep) const override;

    const char *giveInputRecordName() const override { return _IFT_IsotropicLinMoistureTransferMaterial_Name; }
    const char *giveClassName() const override { return "IsotropicLinMoistureTransferMaterial"; }
};
} // end namespace oofem
#endif // isolinmoisturemat_h
