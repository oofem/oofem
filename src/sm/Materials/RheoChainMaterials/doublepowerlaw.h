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

#ifndef doublepowerlaw_h
#define doublepowerlaw_h

#include "maxwellChM.h"

///@name Input fields for DoublePowerLawMaterial
//@{
#define _IFT_DoublePowerLawMaterial_Name "doublepowerlaw"
#define _IFT_DoublePowerLawMaterial_e28 "e28"
#define _IFT_DoublePowerLawMaterial_fi1 "fi1"
#define _IFT_DoublePowerLawMaterial_m "m"
#define _IFT_DoublePowerLawMaterial_n "n"
#define _IFT_DoublePowerLawMaterial_alpha "alpha"
//@}

namespace oofem {
/**
 * This class implements a rheologic double power law material model.
 */
class DoublePowerLawMaterial : public MaxwellChainMaterial
{
protected:
    double E28 = 0.; ///< Young modulus at age of 28 days [MPa].
    double fi1 = 0.; ///< Basic creep coefficient.
    double m = 0., n = 0.;
    double alpha = 0.;

public:
    DoublePowerLawMaterial(int n, Domain *d) : MaxwellChainMaterial(n, d) { }

    const char *giveClassName() const override { return "DoublePowerLawMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_DoublePowerLawMaterial_Name; }
    void initializeFrom(InputRecord &ir) override;

    double computeCreepFunction(double t, double t_prime, GaussPoint *gp, TimeStep *tStep) const override;
};
} // end namespace oofem
#endif // doublepowerlaw_h
