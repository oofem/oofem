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
    double E28; ///< Young modulus at age of 28 days [MPa].
    double fi1; ///< Basic creep coefficient.
    double m, n;
    double alpha;

public:
    DoublePowerLawMaterial(int n, Domain *d) : MaxwellChainMaterial(n, d) { }
    virtual ~DoublePowerLawMaterial() { }

    virtual const char *giveClassName() const { return "DoublePowerLawMaterial"; }
    virtual const char *giveInputRecordName() const { return _IFT_DoublePowerLawMaterial_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual double computeCreepFunction(double t, double t_prime);

protected:
};
} // end namespace oofem
#endif // doublepowerlaw_h
