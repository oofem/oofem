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

#ifndef winklerpasternak_h
#define winklerpasternak_h

#include "sm/Materials/structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matconst.h"
#include "matstatus.h"
#include "stressstrainprincmode.h"

///@name Input fields for WinklerPasternakMaterial
//@{
#define _IFT_WinklerPasternakMaterial_Name "winklerpasternak"
#define _IFT_WinklerPasternakMaterial_C1 "c1"
#define _IFT_WinklerPasternakMaterial_C2 "c2"
#define _IFT_WinklerPasternakMaterial_C2X "c2x"
#define _IFT_WinklerPasternakMaterial_C2Y "c2y"
//@}

namespace oofem {

class GaussPoint;
/**
 * Implementation of 2D winkler-pasternak model for plate (and potentiaonnaly beam) subsoil model.
 */
class WinklerPasternakMaterial : public StructuralMaterial
{
protected:
    /// C1 constant, defined as $\int_0^hE_{oed}(z)\left\(d\Psi(z)\over dz\right\)^2\ dz$
    double c1 = 0.;
    /// C2 constants in x and y directions, defined as $\int_0^hG_{x,y}(z)Psi^2(z)\ dz$
    double c2x = 0., c2y = 0.;

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    WinklerPasternakMaterial(int n, Domain * d);

    bool hasMaterialModeCapability(MaterialMode mode) const override;
    const char *giveClassName() const override { return "WinklerPasternakMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_WinklerPasternakMaterial_Name; }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    FloatArrayF<3> giveRealStressVector_2dPlateSubSoil(const FloatArrayF<3> &reducedE, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> give2dPlateSubSoilStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;
    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;
};
} // end namespace oofem
#endif // winklerpasternak_h
