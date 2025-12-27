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

#ifndef isoheatmat_h
#define isoheatmat_h

#include "tm/Materials/transportmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "scalarfunction.h"
#include "floatmatrixf.h"

///@name Input fields for IsotropicHeatTransferMaterial
//@{
#define _IFT_IsotropicHeatTransferMaterial_Name "isoheat"
#define _IFT_IsotropicHeatTransferMaterial_k "k" ///< Conductivity
#define _IFT_IsotropicHeatTransferMaterial_c "c" ///< Specific heat
#define _IFT_IsotropicHeatTransferMaterial_maturityT0 "maturityt0" ///< Baseline for maturity method
#define _IFT_IsotropicHeatTransferMaterial_d "td"
//@}

namespace oofem {

/**
 * This class implements an isotropic linear heat  material. A material
 * is an attribute of a domain. It is usually also attribute of many elements.
 */
class IsotropicHeatTransferMaterial : public TransportMaterial
{
protected:
    ScalarFunction conductivity; ///< Conductivity (k in input file).
    ScalarFunction capacity;     ///< Capacity (c in input file).
    ScalarFunction density;      ///< Density (td in input file).
    double maturityT0 = 0.;           ///< Baseline for maturity mathod

public:
    IsotropicHeatTransferMaterial(int n, Domain * d);

    FloatArrayF<3> computeFlux3D(const FloatArrayF<3> &grad, double field, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    virtual double giveIsotropicConductivity(GaussPoint *gp, TimeStep *tStep) const;

    double giveCharacteristicValue(MatResponseMode mode,
                                   GaussPoint *gp,
                                   TimeStep *tStep) const override;

    virtual double giveMaturityT0() const { return maturityT0; }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    const char *giveInputRecordName() const override { return _IFT_IsotropicHeatTransferMaterial_Name; }
    const char *giveClassName() const override { return "IsotropicHeatTransferMaterial"; }

    void initializeFrom(InputRecord &ir) override;

    double giveProperty(int aProperty, GaussPoint *gp, TimeStep *tStep) const;
    double giveTemperature(GaussPoint *gp) const;
};

} // end namespace oofem
#endif // isoheatmat_h
