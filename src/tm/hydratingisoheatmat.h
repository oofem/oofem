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

#ifndef hydratingisoheatmat_h
#define hydratingisoheatmat_h

#include "isoheatmat.h"
#include "../sm/hydram.h"

///@name Input fields for HydratingIsoHeatMaterial
//@{
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
    HydratingTransportMaterialStatus(int n, Domain *d, GaussPoint *g) : TransportMaterialStatus(n, d, g), HydrationModelStatusInterface() { }
    virtual ~HydratingTransportMaterialStatus() { }

    virtual Interface *giveInterface(InterfaceType t);
    virtual const char *giveClassName() const { return "HydratingTransportMaterialStatus"; }
    virtual classType giveClassID() const { return HydratingTransportMaterialStatusClass; }

    virtual void updateYourself(TimeStep *atTime) {
        HydrationModelStatusInterface :: updateYourself(atTime);
        TransportMaterialStatus :: updateYourself(atTime);
    }
    virtual void printOutputAt(FILE *file, TimeStep *atTime);
};

/**
 * This class implements a isotropic linear heat material in a finite element problem.
 * A material is an attribute of a domain. It is usually also attribute of many elements.
 * Isotropic Linear Heat Material with interface to the Hydration Model
 */
class HydratingIsoHeatMaterial : public IsotropicHeatTransferMaterial, public HydrationModelInterface
{
protected:
    int hydration, hydrationHeat, hydrationLHS;

public:
    HydratingIsoHeatMaterial(int n, Domain *d) : IsotropicHeatTransferMaterial(n, d), HydrationModelInterface() { }
    virtual ~HydratingIsoHeatMaterial() { }

    void setMixture(MixtureType mix);

    /// Return true if hydration heat source is present.
    virtual int hasInternalSource();
    virtual void computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *atTime, ValueModeType mode);
    virtual void updateInternalState(const FloatArray &state, GaussPoint *gp, TimeStep *tStep);

    virtual double giveCharacteristicValue(MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime);

    virtual contextIOResultType saveIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp);
    virtual contextIOResultType restoreIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp);

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "HydratingIsoHeatMaterial"; }
    virtual classType giveClassID() const { return HydratingIsoHeatMaterialClass; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    // post-processing
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);

protected:
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
};
} // end namespace oofem
#endif // hydratingisoheatmat_h
