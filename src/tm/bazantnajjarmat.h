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

#include "isomoisturemat.h"
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
    double moistureCapacity;

    /// maximal permeability [kg/ m s]
    double C1;
    /// exponent in nonlinear permeability function [-]
    double n;
    /// fraction minimal/maximal permeability [-]
    double alpha0;
    /// nonlinear threshold [-]
    double hC;

public:
    BazantNajjarMoistureTransferMaterial(int n, Domain * d) : IsotropicMoistureTransferMaterial(n, d) { }
    virtual ~BazantNajjarMoistureTransferMaterial() { }

    /// evaluates permeability according to Bazant - Najjar function for diffusivity
    virtual double givePermeability(GaussPoint *gp, TimeStep *tStep);
    /// evaluates slope of the sorption isotherm
    virtual double giveMoistureCapacity(GaussPoint *gp, TimeStep *tStep);

    virtual const char *giveInputRecordName() const { return _IFT_BazantNajjarMoistureTransferMaterial_Name; }
    virtual const char *giveClassName() const { return "BazantNajjarMoistureTransferMaterial"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual double giveHumidity(GaussPoint *gp, ValueModeType mode);
};
} // end namespace oofem
#endif // bazantnajjarmat_h
