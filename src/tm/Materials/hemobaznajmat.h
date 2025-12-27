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

#ifndef hemobaznajmat_h
#define hemobaznajmat_h

#include "tm/Materials/transportmaterial.h"

///@name Input fields for HeMoBazNajMaterial
//@{
#define _IFT_HeMoBazNajMaterial_Name "hemobaznajmat"
#define _IFT_HeMoBazNajMaterial_c1 "c1"
#define _IFT_HeMoBazNajMaterial_n "n"
#define _IFT_HeMoBazNajMaterial_alpha0 "alpha0"
#define _IFT_HeMoBazNajMaterial_hc "hc"
#define _IFT_HeMoBazNajMaterial_capa "capa"
#define _IFT_HeMoBazNajMaterial_k "k" ///< Conductivity
#define _IFT_HeMoBazNajMaterial_c "c" ///< Specific heat
//@}

namespace oofem {
/**
 */
class HeMoBazNajMaterial : public TransportMaterial
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

    double heatConductivity = 0.; ///< Conductivity (k in input file).
    double heatCapacity = 0.;     ///< Capacity (c in input file).

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    HeMoBazNajMaterial(int n, Domain *d) : TransportMaterial(n, d) { }

    std::pair<FloatArrayF<3>, FloatArrayF<3>> computeHeMoFlux3D(const FloatArrayF<3> &grad_t, const FloatArrayF<3> &grad_w, double t, double h, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    double giveCharacteristicValue(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override;

    void initializeFrom(InputRecord &ir) override;

    double give(int aProperty, GaussPoint *gp) const override;

    bool hasMaterialModeCapability(MaterialMode mode) const override;
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    const char *giveInputRecordName() const override { return _IFT_HeMoBazNajMaterial_Name; }
    const char *giveClassName() const override { return "HeMoBazNajMaterial"; }

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<HeMoTransportMaterialStatus>(gp); }

protected:
    double computeCapacityCoeff(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const;

    double perm_mm(double h, double T) const;
    double perm_mh(double h, double T) const;
    double perm_hm(double h, double T) const;
    double perm_hh(double h, double T) const;
};
} // end namespace oofem
#endif // hemobaznajmat_h
