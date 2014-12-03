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

#ifndef latticetransmat_h
#define latticetransmat_h

#include "transportmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "cltypes.h"
#include "matstatus.h"

///@name Input fields for Lattice2d_mt
//@{
#define _IFT_LatticeTransportMaterial_Name "latticetransmat"
#define _IFT_LatticeTransportMaterial_vis "vis"
#define _IFT_LatticeTransportMaterial_k "k"
#define _IFT_LatticeTransportMaterial_thetas "thetas"
#define _IFT_LatticeTransportMaterial_thetar "thetar"
#define _IFT_LatticeTransportMaterial_contype "contype"
#define _IFT_LatticeTransportMaterial_m "m"
#define _IFT_LatticeTransportMaterial_a "a"
#define _IFT_LatticeTransportMaterial_thetam "thetam"
#define _IFT_LatticeTransportMaterial_paev "paev"
#define _IFT_LatticeTransportMaterial_ctor "ctor"
#define _IFT_LatticeTransportMaterial_clim "clim"
#define _IFT_LatticeTransportMaterial_c "c"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to LatticeTransportMaterial.
 */
class LatticeTransportMaterialStatus : public TransportMaterialStatus
{
protected:
    /// Liquid mass in element
    double mass;
    double oldPressure;

public:
    ///Constructor
    LatticeTransportMaterialStatus(int n, Domain * d, GaussPoint * g);
    /// Destructor
    virtual ~LatticeTransportMaterialStatus() { }

    void printOutputAt(FILE *, TimeStep *);

    /// Returns pressure
    double givePressure() { return field.at(1); }

    double giveOldPressure() { return oldPressure; }

    /// Sets the mass
    void setMass(double input) { this->mass = input; }

    /// Returns mass
    double giveMass() { return this->mass; }

    virtual void updateYourself(TimeStep *tStep);

    virtual void initTempStatus();

    virtual const char *giveClassName() const { return "LatticeTransportMaterialStatus"; }
};

/**
 * This class implements a transport constitutive model for saturated and unsaturated porous materials for lattice elements.
 */
class LatticeTransportMaterial : public TransportMaterial
{
protected:
    ///Viscosity of fluid
    double viscosity;

    ///Parameter of van Genuchten law
    double paramM;

    ///Parameter of van Genuchten law
    double paramA;

    ///Intrinsic permeability of porous material
    double permeability;

    ///Porosity of porous material
    double porosity;

    ///Density of fluid
    double density;

    ///Type of conductivity and capcity laws.
    int conType;

    ///Type of conductivity and capcity laws.
    int capacity;

    ///Relative saturated water content
    double thetaS;

    ///Residual water content
    double thetaR;

    ///modified water content
    double thetaM;

    ///crack tortuosity
    double crackTortuosity;

    ///crack limit
    double crackLimit;


    ///suction air entry value
    double suctionAirEntry;

    /// Material mode for convenient access.
    MaterialMode matMode;

public:
    LatticeTransportMaterial(int n, Domain * d) : TransportMaterial(n, d) { }

    virtual ~LatticeTransportMaterial() { }

    virtual void giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep);

    virtual void  giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *tStep) {; }

    virtual double  giveCharacteristicValue(MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *tStep);

    /**
     * Computes the conductivity.
     * @param suction Capillary stress
     * @param gp Integration point.
     * @param tStep Time step.
     */
    double computeConductivity(double suction, GaussPoint *gp, TimeStep *tStep);

    /**
     * Computes the capacity.
     * @param suction Capillary stress
     * @param gp Integration point.
     */
    double computeCapacity(double suction, GaussPoint *gp);

    /**
     * Computes the mass.
     * @param stateVector Capillary stress
     * @param gp Integration point.
     */

    double computeMass(FloatArray &stateVector, GaussPoint *gp);

    virtual const char *giveInputRecordName() const { return _IFT_LatticeTransportMaterial_Name; }
    virtual const char *giveClassName() const { return "LatticeTransportMaterial"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual double give(int, GaussPoint *gp);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
};
} // end namespace oofem
#endif
