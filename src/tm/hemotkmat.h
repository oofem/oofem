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

#ifndef hemotkmat_h
#define hemotkmat_h

#include "transportmaterial.h"

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
    double a_0;       ///< Constant (obtained from experiments) [Bazant and Najjar, 1972]
    double nn;        ///< Constant-exponent (obtained from experiments) [Bazant and Najjar, 1972]
    double phi_c;     ///< Constant-relative humidity  (obtained from experiments) [Bazant and Najjar, 1972]
    double delta_wet; ///< Constant-water vapor permeability (obtained from experiments) [Bazant and Najjar, 1972]

    double w_h;       ///< Constant water content (obtained from experiments) [Pedersen, 1990]
    double n;         ///< Constant-exponent (obtained from experiments) [Pedersen, 1990]
    double a;         ///< Constant (obtained from experiments) [Pedersen, 1990]

    double latent;    ///< Latent heat of evaporation.
    double c;         ///< Thermal capacity.
    double rho;       ///< Volume density.
    double chi_eff;   ///< Effective thermal conductivity.

    double por;       ///< Porosity.
    double rho_gws;   ///< Saturation volume density.

public:

    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    HeMoTKMaterial(int n, Domain *d) : TransportMaterial(n, d) { }
    /// Destructor.
    virtual ~HeMoTKMaterial() { }

    virtual void giveCharacteristicMatrix(FloatMatrix &answer,
                                          MatResponseForm form,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *atTime);

    virtual double giveCharacteristicValue(MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime);

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual double give(int aProperty, GaussPoint *gp);

    virtual int hasMaterialModeCapability(MaterialMode mode);

    // identification
    virtual const char *giveClassName() const { return "HeMoTKMaterial"; }
    virtual classType giveClassID() const { return HeMoTKMaterialClass; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new TransportMaterialStatus(1, domain, gp); }

    double sorption_isotherm(double phi);
    double inverse_sorption_isotherm(double w);
    double give_dphi_dw(double w);

protected:
    void computeConductivityMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void matcond1d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime);
    void matcond2d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime);
    void matcond3d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime);

    double computeCapacityCoeff(MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    /**
     * Returns positive value of humidity, use VM_Velocity for previous (equilibrated) value
     */
    double giveHumidity(GaussPoint *gp, ValueModeType mode);

    double get_latent(double w, double t);
    double get_ceff(double w, double t);
    double get_chi(double w, double t);

    double perm_wt(double w, double t);
    double perm_ww(double w, double t);
    double give_delta_gw(double phi);
    double give_dpgw_dt(double t, double phi);

    double get_b(double w, double t);
    double get_sat(double w, double t);
    double give_p_gws(double t);

    // post-processing, poi export
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
};
} // end namespace oofem
#endif // hemotkmat_h
