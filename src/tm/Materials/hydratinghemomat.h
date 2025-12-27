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

#ifndef hydratinghemomat_h
#define hydratinghemomat_h

#include "tm/Materials/hemotkmat.h"
#include "tm/Materials/hydram.h"

///@name Input fields for HydratingHeMoMaterial
//@{
#define _IFT_HydratingHeMoMaterial_Name "hhemotk"
#define _IFT_HydratingHeMoMaterial_hydration "hydration"
#define _IFT_HydratingHeMoMaterial_mix "mix"
#define _IFT_HydratingHeMoMaterial_noHeat "noheat"
#define _IFT_HydratingHeMoMaterial_noLHS "nolhs"
//@}

namespace oofem {
/**
 * Heat and moisture transport material with hydration.
 */
class HydratingHeMoMaterial : public HeMoTKMaterial, public HydrationModelInterface
{
protected:
    bool hydration = false, hydrationHeat = false, hydrationLHS = false, teplotaOut = false;

public:
    HydratingHeMoMaterial(int n, Domain * d) : HeMoTKMaterial(n, d), HydrationModelInterface() { }

    void setMixture(MixtureType mix);

    bool hasInternalSource() const override;
    void computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const override;
    void updateInternalState(const FloatArray &state, GaussPoint *gp, TimeStep *tStep) override;

    double giveCharacteristicValue(MatResponseMode mode,
                                   GaussPoint *gp,
                                   TimeStep *tStep) const override;

    // saves current context(state) into stream
    void saveIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp) override;
    void restoreIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp) override;

    // identification and auxiliary functions
    const char *giveInputRecordName() const override { return _IFT_HydratingHeMoMaterial_Name; }
    const char *giveClassName() const override { return "HydratingHeMoMaterial"; }

    void initializeFrom(InputRecord &ir) override;

    // post-processing
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;
};
} // end namespace oofem
#endif // hydratinghemomat_h
