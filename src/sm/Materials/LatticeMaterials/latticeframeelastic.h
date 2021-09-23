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
 *               Copyright (C) 1993 - 2019   Borek Patzak
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

#ifndef latticeframeelastic_h
#define latticeframeelastic_h

#include "latticestructuralmaterial.h"
#include "cltypes.h"
#include "randommaterialext.h"
#include "strainvector.h"
#include "stressvector.h"
#include "latticematstatus.h"

///@name Input fields for LatticeFrameElastic
//@{
#define _IFT_LatticeFrameElastic_Name "latticeframeelastic"
#define _IFT_LatticeFrameElastic_talpha "talpha"
#define _IFT_LatticeFrameElastic_e "e"
#define _IFT_LatticeFrameElastic_n "n"
//@}

namespace oofem {
/**
 * This class implements a local random linear elastic model for lattice elements.
 */
class LatticeFrameElastic : public LatticeStructuralMaterial
{
protected:
    ///Normal modulus
    double e;

    ///Ratio of shear and normal modulus
    double nu;

public:
    LatticeFrameElastic(int n, Domain *d) : LatticeStructuralMaterial(n, d) { };


    FloatArrayF< 6 >giveThermalDilatationVector(GaussPoint *gp,  TimeStep *tStep) const override;

    const char *giveInputRecordName() const override { return _IFT_LatticeFrameElastic_Name; }

    const char *giveClassName() const override { return "LatticeFrameElastic"; }

    void initializeFrom(InputRecord &ir) override;


    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    FloatArrayF< 6 >giveFrameForces3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) override;

    FloatMatrixF< 6, 6 >give3dFrameStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;


    bool hasMaterialModeCapability(MaterialMode mode) const override;

    Interface *giveInterface(InterfaceType) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    MaterialStatus *giveStatus(GaussPoint *gp) const override;

protected:
};
} // end namespace oofem

#endif
