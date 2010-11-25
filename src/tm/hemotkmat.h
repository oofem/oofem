/* $Header: */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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
#include "flotarry.h"
#include "flotmtrx.h"

#include "matconst.h"

namespace oofem {
class GaussPoint;

/**
 * This class implements a coupled heat and mass transfer material model. It computes condictivity and capacity
 * matrices for coupled heat and moisture transfer.
 * Autor: Tomas Krejci
 * Source: T.K. Doctoral Thesis; Bazant and Najjar, 1972; Pedersen, 1990;
 * Assumptions:      water vapor is the only driving mechanism; relative humidity is from range 0.2 - 0.98 (I and II region);
 */
class HeMoTKMaterial : public TransportMaterial
{
protected:
    double a_0;       //constant (obtained from experiments) a_0 [Bazant and Najjar, 1972]
    double nn;        //constant-exponent (obtained from experiments) n [Bazant and Najjar, 1972]
    double phi_c;     //constant-relative humidity  (obtained from experiments) phi_c [Bazant and Najjar, 1972]
    double delta_wet; //constant-water vapor permeability (obtained from experiments) delta_wet [Bazant and Najjar, 1972]

    double w_h;       //constant water content (obtained from experiments) w_h [Pedersen, 1990]
    double n;         //constant-exponent (obtained from experiments) n [Pedersen, 1990]
    double a;         //constant (obtained from experiments) A [Pedersen, 1990]

    double latent;    //latent heat of evaporation
    double c;         //thermal capacity
    double rho;       //volume density
    double chi_eff;   //effective thermal conductivity

    double por;       //porosity
    double rho_gws;   //saturation volume density

public:

    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n material number
     * @param d domain to which new material will belong
     */
    HeMoTKMaterial(int n, Domain *d) : TransportMaterial(n, d) { }
    /// Destructor.
    ~HeMoTKMaterial()                { }

    /**
     * Computes the characteristic matrix of receiver in given integration point, respecting its history.
     * The algorithm should use temporary or equlibrium  history variables stored in integration point status
     * to compute and return required result.
     * @param answer contains result
     * @param form material response form
     * @param mode  material response mode
     * @param gp integration point
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    virtual void  giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime);

    /**
     * Computes the characteristic value of receiver in given integration point, respecting its history.
     * The algorithm should use temporary or equlibrium  history variables stored in integration point status
     * to compute and return required result.
     * @param mode material response mode
     * @param gp integration point
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    virtual double  giveCharacteristicValue(MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *atTime);


    /**
     * Returns true if stiffness matrix of receiver is symmetric
     * Default implementation returns true.
     */
    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode);

    /**
     * Initializes receiver acording to object description stored in input record.
     */
    IRResultType initializeFrom(InputRecord *ir);

    // non-standard - returns time independent material constant
    double   give(int, GaussPoint *);

    /**
     * Request material extension.
     * @param ext material extension tested
     * @return nonzero if implemented
     */
    virtual int testMaterialExtension(MaterialExtension ext) { return ( ( ext == Material_TransportCapability ) ? 1 : 0 ); }
    int hasMaterialModeCapability(MaterialMode mode);
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "HeMoTKMaterial"; }
    /// Returns classType id of receiver.
    classType giveClassID()         const { return HeMoTKMaterialClass; }

    /**
     * Creates new copy of associated status
     * and inserts it into given integration point.
     * @param gp Integration point where newly created status will be stored.
     * @return reference to new status.
     */
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new TransportMaterialStatus(1, domain, gp); }

    double sorption_isotherm(double phi);
    double inverse_sorption_isotherm(double w);
    double give_dphi_dw(double w);

protected:


    //void matcond (matrix &d,long ncomp,long ipp);
    void computeConductivityMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void matcond1d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime);
    void matcond2d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime);
    void matcond3d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime);

    //void matcap (double &c,long ipp);
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

#ifdef __OOFEG
#endif
};
} // end namespace oofem
#endif // hemotkmat_h
