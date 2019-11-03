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

#ifndef hemotkmat_h
#define hemotkmat_h

#include "tm/Materials/transportmaterial.h"

///@name Input fields for HeMoTKMaterial
//@{
#define _IFT_HeMoTKMaterial_Name "hemotk"
#define _IFT_HeMoTKMaterial_a_0 "a_0"
#define _IFT_HeMoTKMaterial_nn "nn"
#define _IFT_HeMoTKMaterial_phi_c "phi_c"
#define _IFT_HeMoTKMaterial_delta_wet "delta_wet"
#define _IFT_HeMoTKMaterial_w_h "w_h"
#define _IFT_HeMoTKMaterial_n "n"
#define _IFT_HeMoTKMaterial_a "a"
#define _IFT_HeMoTKMaterial_latent "latent"
#define _IFT_HeMoTKMaterial_c "c"
#define _IFT_HeMoTKMaterial_rho "rho"
#define _IFT_HeMoTKMaterial_chi_eff "chi_eff"
#define _IFT_HeMoTKMaterial_por "por"
#define _IFT_HeMoTKMaterial_rho_gws "rho_gws"
//@}

namespace oofem {
/**
 * This class implements a coupled heat and mass transfer material model. It computes conductivity and capacity
 * matrices for coupled heat and moisture transfer.
 * Assumptions: Water vapor is the only driving mechanism; relative humidity is from range 0.2 - 0.98 (I and II region).
 * @author Tomas Krejci
 * Source: T.K. Doctoral Thesis; Bazant and Najjar, 1972; Pedersen, 1990.
 */
class HeMoTKMaterial : public TransportMaterial
{
protected:
    double a_0 = 0.;       ///< Constant (obtained from experiments) [Bazant and Najjar, 1972]
    double nn = 0.;        ///< Constant-exponent (obtained from experiments) [Bazant and Najjar, 1972]
    double phi_c = 0.;     ///< Constant-relative humidity  (obtained from experiments) [Bazant and Najjar, 1972]
    double delta_wet = 0.; ///< Constant-water vapor permeability (obtained from experiments) [Bazant and Najjar, 1972]

    double w_h = 0.;       ///< Constant water content (obtained from experiments) [Pedersen, 1990]
    double n = 0.;         ///< Constant-exponent (obtained from experiments) [Pedersen, 1990]
    double a = 0.;         ///< Constant (obtained from experiments) [Pedersen, 1990]

    double latent = 0.;    ///< Latent heat of evaporation.
    double c = 0.;         ///< Thermal capacity.
    double rho = 0.;       ///< Volume density.
    double chi_eff = 0.;   ///< Effective thermal conductivity.

    double por = 0.;       ///< Porosity.
    double rho_gws = 0.;   ///< Saturation volume density.

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    HeMoTKMaterial(int n, Domain * d) : TransportMaterial(n, d) { }

    std::pair<FloatArrayF<3>, FloatArrayF<3>> computeHeMoFlux3D(const FloatArrayF<3> &grad_t, const FloatArrayF<3> &grad_w, double t, double h, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    double giveCharacteristicValue(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override;

    void initializeFrom(InputRecord &ir) override;

    double give(int aProperty, GaussPoint *gp) const override;

    bool hasMaterialModeCapability(MaterialMode mode) const override;
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    // identification
    const char *giveInputRecordName() const override { return _IFT_HeMoTKMaterial_Name; }
    const char *giveClassName() const override { return "HeMoTKMaterial"; }

    double sorption_isotherm(double phi) const;
    double inverse_sorption_isotherm(double w) const;
    double give_dphi_dw(double w) const;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new HeMoTransportMaterialStatus(gp); }

protected:
    double computeCapacityCoeff(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const;
    /**
     * Returns positive value of humidity, use VM_Velocity for previous (equilibrated) value
     */
    double giveHumidity(GaussPoint *gp, ValueModeType mode) const override;

    double get_latent(double w, double t) const;
    double get_ceff(double w, double t) const;
    double get_chi(double w, double t) const;

    double perm_wt(double w, double t) const;
    double perm_ww(double w, double t) const;
    double give_delta_gw(double phi) const;
    double give_dpgw_dt(double t, double phi) const;

    double get_b(double w, double t) const;
    double get_sat(double w, double t) const;
    double give_p_gws(double t) const;
};
} // end namespace oofem
#endif // hemotkmat_h
