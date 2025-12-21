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

#ifndef latticeslip_h
#define latticeslip_h

#include "latticelinearelastic.h"
#include "latticematstatus.h"

///@name Input fields for LatticeSlip
//@{
#define _IFT_LatticeSlip_Name "latticeslip"
#define _IFT_LatticeSlip_talpha "talpha"
#define _IFT_LatticeSlip_e "e"
#define _IFT_LatticeSlip_a1 "a1"
#define _IFT_LatticeSlip_a2 "a2"
#define _IFT_LatticeSlip_t0 "t0"
//@}

namespace oofem {
class LatticeSlipStatus : public LatticeMaterialStatus
{
protected:


public:
    LatticeSlipStatus(GaussPoint *g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "LatticeSlipStatus"; }

    void initTempStatus() override;

    void updateYourself(TimeStep *) override;

    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;
};



/**
 * This class implements a slip model for interface elements between lattice and beam elements.
 */

class LatticeSlip : public LatticeLinearElastic
{
protected:

    ///Normal modulus
    double eNormal = 0.;

    ///Ratio of shear and normal modulus
    double alphaOne = 0.;

    ///Ratio of torsion and normal modulus
    double alphaTwo = 0.;

    ///Strength for slip component
    double tauZero = 0.;

    /// coefficient variation of the Gaussian distribution
    double coefficientOfVariation = 0.;

    /// flag which chooses between no distribution (0) and Gaussian distribution (1)
    double localRandomType = 0.;

    double cAlpha = 0.;

    double tAlphaMean = 0.;

public:

    /// Constructor
    LatticeSlip(int n, Domain *d);


    const char *giveInputRecordName() const override { return _IFT_LatticeSlip_Name; }
    const char *giveClassName() const override { return "LatticeSlip"; }

    void initializeFrom(InputRecord &ir) override;

    FloatArrayF< 6 >giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const override;


    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }


    double computeDeltaDissipation(GaussPoint *gp, TimeStep *atTime) const;

    bool hasMaterialModeCapability(MaterialMode mode) const override;


    Interface *giveInterface(InterfaceType) override;

    FloatArrayF< 6 >giveLatticeStress3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;


protected:

    int giveIPValue(FloatArray &answer,
                    GaussPoint *gp,
                    InternalStateType type,
                    TimeStep *atTime) override;
};
} // end namespace oofem

#endif
