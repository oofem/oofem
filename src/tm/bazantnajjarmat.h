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
    BazantNajjarMoistureTransferMaterial(int n, Domain *d) : IsotropicMoistureTransferMaterial(n, d) { }
    virtual ~BazantNajjarMoistureTransferMaterial() { }

    /// evaluates permeability according to Bazant - Najjar function for diffusivity
    virtual double givePermeability(GaussPoint *gp, TimeStep *atTime);
    /// evaluates slope of the sorption isotherm
    virtual double giveMoistureCapacity(GaussPoint *gp, TimeStep *atTime);

    virtual const char *giveClassName() const { return "BazantNajjarMoistureTransferMaterial"; }
    virtual classType giveClassID() const { return BazantNajjarMoistureTransferMaterialClass; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual double giveHumidity(GaussPoint *gp, ValueModeType mode);

    /*
     * virtual double give(int aProperty, GaussPoint *gp);
     */

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new TransportMaterialStatus(1, domain, gp);  }
};
} // end namespace oofem
#endif // bazantnajjarmat_h
