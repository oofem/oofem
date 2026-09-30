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
 *               Copyright (C) 1993 - 2026   Borek Patzak
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

#ifndef latticedamageshell_h
#define latticedamageshell_h

#include "latticedamage.h"

///@name Input fields for LatticeDamageShell
//@{
#define _IFT_LatticeDamageShell_Name "latticedamageshell"
//@}

namespace oofem {

/**
 * LatticeDamage for shell cross-sections. Differs from the base only in the damage
 * driver: the out-of-plane (transverse) shear component (2) is excluded from the
 * equivalent strain, so transverse shear cannot initiate cracking.
 */
class LatticeDamageShell : public LatticeDamage
{
public:
    LatticeDamageShell(int n, Domain *d) : LatticeDamage(n, d) { }

    const char *giveInputRecordName() const override { return _IFT_LatticeDamageShell_Name; }
    const char *giveClassName() const override { return "LatticeDamageShell"; }

    double computeEquivalentStrain(const FloatArrayF< 6 > &strain, GaussPoint *gp) const override;
};

} // end namespace oofem
#endif
