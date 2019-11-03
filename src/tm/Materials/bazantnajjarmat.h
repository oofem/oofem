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

#ifndef bazantnajjarmat_h
#define bazantnajjarmat_h

#include "tm/Materials/isomoisturemat.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for BazantNajjarMoistureTransferMaterial
//@{
#define _IFT_BazantNajjarMoistureTransferMaterial_Name "bazantnajjarmoisturemat"
#define _IFT_BazantNajjarMoistureTransferMaterial_c1 "c1"
#define _IFT_BazantNajjarMoistureTransferMaterial_n "n"
#define _IFT_BazantNajjarMoistureTransferMaterial_alpha0 "alpha0"
#define _IFT_BazantNajjarMoistureTransferMaterial_hc "hc"
#define _IFT_BazantNajjarMoistureTransferMaterial_capa "capa"
//@}

namespace oofem {
/**
 * This class implements a isotropic moisture tranport material. A material
 * is an attribute of a domain. It is usually also attribute of many elements.
 */
class BazantNajjarMoistureTransferMaterial : public IsotropicMoistureTransferMaterial
{
protected:
    /// sorption isotherm derivative [kg/m^3]
    double moistureCapacity = 0.;

    /// maximal permeability [kg/ m s]
    double C1 = 0.;
    /// exponent in nonlinear permeability function [-]
    double n = 0.;
    /// fraction minimal/maximal permeability [-]
    double alpha0 = 0.;
    /// nonlinear threshold [-]
    double hC = 0.;

public:
    BazantNajjarMoistureTransferMaterial(int n, Domain * d) : IsotropicMoistureTransferMaterial(n, d) { }

    /// evaluates permeability according to Bazant - Najjar function for diffusivity
    double givePermeability(GaussPoint *gp, TimeStep *tStep) const override;
    /// evaluates slope of the sorption isotherm
    double giveMoistureCapacity(GaussPoint *gp, TimeStep *tStep) const override;

    const char *giveInputRecordName() const override { return _IFT_BazantNajjarMoistureTransferMaterial_Name; }
    const char *giveClassName() const override { return "BazantNajjarMoistureTransferMaterial"; }

    void initializeFrom(InputRecord &ir) override;

    double giveHumidity(GaussPoint *gp, ValueModeType mode) const override;
};
} // end namespace oofem
#endif // bazantnajjarmat_h
