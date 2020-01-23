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

#ifndef latticelinearelastic_h
#define latticelinearelastic_h

#include "latticestructuralmaterial.h"
#include "cltypes.h"
#include "randommaterialext.h"
#include "strainvector.h"
#include "stressvector.h"
#include "latticematstatus.h"

///@name Input fields for LatticeLinearElastic
//@{
#define _IFT_LatticeLinearElastic_Name "latticelinearelastic"
#define _IFT_LatticeLinearElastic_talpha "talpha"
#define _IFT_LatticeLinearElastic_e "e"
#define _IFT_LatticeLinearElastic_n "n"
#define _IFT_LatticeLinearElastic_a1 "a1"
#define _IFT_LatticeLinearElastic_a2 "a2"
#define _IFT_LatticeLinearElastic_localrandomtype "randomtype"
#define _IFT_LatticeLinearElastic_cov "cov"
#define _IFT_LatticeLinearElastic_calpha "calpha"
//@}

namespace oofem {
/**
 * This class implements a local random linear elastic model for lattice elements.
 */
class LatticeLinearElastic : public LatticeStructuralMaterial, public RandomMaterialExtensionInterface
{
protected:
    ///Normal modulus
    double eNormalMean = 0.;

    ///Ratio of shear and normal modulus
    double alphaOne = 0.;

    ///Ratio of torsion and normal modulus
    double alphaTwo = 0.;

    /// coefficient variation of the Gaussian distribution
    double coefficientOfVariation = 0.;

    /// flag which chooses between no distribution (0) and Gaussian distribution (1)
    double localRandomType = 0.;

    /// parameter which allows to prescribed thermal displacement
    double cAlpha = 0.;

public:
    LatticeLinearElastic(int n, Domain *d) : LatticeStructuralMaterial(n, d), RandomMaterialExtensionInterface() { };


    LatticeLinearElastic(int n, Domain *d, double eNormalMean, double alphaOne, double alphaTwo);

    const char *giveInputRecordName() const override { return _IFT_LatticeLinearElastic_Name; }

    const char *giveClassName() const override { return "LatticeLinearElastic"; }

    void initializeFrom(InputRecord &ir) override;

    FloatArrayF< 6 >giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    FloatArrayF< 6 >giveLatticeStress3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) override;

    FloatMatrixF< 6, 6 >give3dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF< 3, 3 >give2dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;


    bool hasMaterialModeCapability(MaterialMode mode) const override;


    Interface *giveInterface(InterfaceType) override;

    virtual void giveRandomParameters(FloatArray &param);

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    MaterialStatus *giveStatus(GaussPoint *gp) const override;

    double  give(int aProperty, GaussPoint *gp) const override;

protected:
};
} // end namespace oofem

#endif
