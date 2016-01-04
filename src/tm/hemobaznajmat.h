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

#ifndef hemobaznajmat_h
#define hemobaznajmat_h

#include "transportmaterial.h"

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
    double moistureCapacity;

    /// maximal permeability [kg/ m s]
    double C1;
    /// exponent in nonlinear permeability function [-]
    double n;
    /// fraction minimal/maximal permeability [-]
    double alpha0;
    /// nonlinear threshold [-]
    double hC;

    double heatConductivity; ///< Conductivity (k in input file).
    double heatCapacity;     ///< Capacity (c in input file).


public:

    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    HeMoBazNajMaterial(int n, Domain *d) : TransportMaterial(n, d) { }
    /// Destructor.
    virtual ~HeMoBazNajMaterial() { }

    virtual void giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep);

    virtual void giveCharacteristicMatrix(FloatMatrix &answer,
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
    virtual const char *giveInputRecordName() const { return _IFT_HeMoBazNajMaterial_Name; }
    virtual const char *giveClassName() const { return "HeMoBazNajMaterial"; }

protected:
    void computeConductivityMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void matcond1d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime);
    void matcond2d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime);
    void matcond3d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime);

    double computeCapacityCoeff(MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);

    double perm_mm(double h, double T);
    double perm_mh(double h, double T);
    double perm_hm(double h, double T);
    double perm_hh(double h, double T);


    // post-processing, poi export
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime);
};
} // end namespace oofem
#endif // hemobaznajmat_h
