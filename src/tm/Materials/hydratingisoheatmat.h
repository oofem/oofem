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

#ifndef hydratingisoheatmat_h
#define hydratingisoheatmat_h

#include "tm/Materials/isoheatmat.h"
#include "tm/Materials/hydram.h"

///@name Input fields for HydratingIsoHeatMaterial
//@{
#define _IFT_HydratingIsoHeatMaterial_Name "hisoheat"
#define _IFT_HydratingIsoHeatMaterial_hydration "hydration"
#define _IFT_HydratingIsoHeatMaterial_mix "mix"
#define _IFT_HydratingIsoHeatMaterial_noHeat "noheat"
#define _IFT_HydratingIsoHeatMaterial_noLHS "nolhs"
//@}

namespace oofem {
/**
 * Isotropic material for heat with hydration.
 */
class HydratingTransportMaterialStatus : public TransportMaterialStatus, public HydrationModelStatusInterface
{
public:
    HydratingTransportMaterialStatus(GaussPoint * g) : TransportMaterialStatus(g), HydrationModelStatusInterface() { }

    Interface *giveInterface(InterfaceType t) override;
    const char *giveClassName() const override { return "HydratingTransportMaterialStatus"; }

    void updateYourself(TimeStep *tStep) override {
        HydrationModelStatusInterface :: updateYourself(tStep);
        TransportMaterialStatus :: updateYourself(tStep);
    }
    void printOutputAt(FILE *file, TimeStep *tStep) const override;
};

/**
 * This class implements a isotropic linear heat material in a finite element problem.
 * A material is an attribute of a domain. It is usually also attribute of many elements.
 * Isotropic Linear Heat Material with interface to the Hydration Model
 */
class HydratingIsoHeatMaterial : public IsotropicHeatTransferMaterial, public HydrationModelInterface
{
protected:
    bool hydration = false, hydrationHeat = false, hydrationLHS = false;

public:
    HydratingIsoHeatMaterial(int n, Domain * d) : IsotropicHeatTransferMaterial(n, d), HydrationModelInterface() { }

    void setMixture(MixtureType mix);

    /// Return true if hydration heat source is present.
    bool hasInternalSource() const override;
    void computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const override;
    void updateInternalState(const FloatArray &state, GaussPoint *gp, TimeStep *tStep) override;

    double giveCharacteristicValue(MatResponseMode mode,
                                   GaussPoint *gp,
                                   TimeStep *tStep) const override;

    void saveIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp) override;
    void restoreIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp) override;

    // identification and auxiliary functions
    const char *giveInputRecordName() const override { return _IFT_HydratingIsoHeatMaterial_Name; }
    const char *giveClassName() const override { return "HydratingIsoHeatMaterial"; }

    void initializeFrom(InputRecord &ir) override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    MaterialStatus *CreateStatus(GaussPoint *gp) const override;
};
} // end namespace oofem
#endif // hydratingisoheatmat_h
