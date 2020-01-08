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

#include "tm/Materials/transportmaterial.h"
#include "dictionary.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
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
    double mass = 0.;
    double oldPressure = 0.;

public:
    ///Constructor
    LatticeTransportMaterialStatus(GaussPoint *g);

    void printOutputAt(FILE *, TimeStep *) const override;

    /// Returns pressure
    double givePressure() const { return field; }

    double giveOldPressure() const { return oldPressure; }

    /// Sets the mass
    void setMass(double input) { this->mass = input; }

    /// Returns mass
    double giveMass() const { return this->mass; }

    void updateYourself(TimeStep *tStep) override;

    void initTempStatus() override;

    const char *giveClassName() const override { return "LatticeTransportMaterialStatus"; }
};

/**
 * This class implements a transport constitutive model for saturated and unsaturated porous materials for lattice elements.
 */
class LatticeTransportMaterial : public TransportMaterial
{
protected:
    /// Viscosity of fluid
    double viscosity = 0.;

    /// Parameter of van Genuchten law
    double paramM = 0.;

    /// Parameter of van Genuchten law
    double paramA = 0.;

    /// Intrinsic permeability of porous material
    double permeability = 0.;

    /// Porosity of porous material
    double porosity = 0.;

    /// Density of fluid
    double density = 0.;

    /// Type of conductivity and capcity laws.
    int conType = 0.;

    /// Type of conductivity and capcity laws.
    int capacity = 0.;

    /// Relative saturated water content
    double thetaS = 0.;

    /// Residual water content
    double thetaR = 0.;

    /// Modified water content
    double thetaM = 0.;

    /// Crack tortuosity
    double crackTortuosity = 0.;

    /// Crack limit
    double crackLimit = 0.;

    /// Suction air entry value
    double suctionAirEntry = 0.;


public:
    LatticeTransportMaterial(int n, Domain *d) : TransportMaterial(n, d) { }

    void initializeFrom(InputRecord &ir) override;

    FloatArrayF< 3 >computeFlux3D(const FloatArrayF< 3 > &grad, double field, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF< 3, 3 >computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override { return {}; }

    double  giveCharacteristicValue(MatResponseMode mode,
                                    GaussPoint *gp,
                                    TimeStep *tStep) const override;

    /**
     * Computes the conductivity.
     * @param suction Capillary stress
     * @param gp Integration point.
     * @param tStep Time step.
     */
    double computeConductivity(double suction, GaussPoint *gp, TimeStep *tStep) const;

    /**
     * Computes the capacity.
     * @param suction Capillary stress
     * @param gp Integration point.
     */
    double computeCapacity(double suction, GaussPoint *gp) const;

    /**
     * Computes the mass.
     * @param stateVector Capillary stress
     * @param gp Integration point.
     */

    double computeMass(FloatArray &stateVector, GaussPoint *gp) const;

    const char *giveInputRecordName() const override { return _IFT_LatticeTransportMaterial_Name; }
    const char *giveClassName() const override { return "LatticeTransportMaterial"; }

    double give(int, GaussPoint *gp) const override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;
};
} // end namespace oofem
#endif
