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
#define _IFT_LatticeLinearElastic_a3 "a3"
#define _IFT_LatticeLinearElastic_tcrit "tcrit"
#define _IFT_LatticeLinearElastic_nu "nu"
#define _IFT_LatticeLinearElastic_em "em"
#define _IFT_LatticeLinearElastic_bio "bio"
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

    ///Ratio of bending and normal modulus
    double alphaTwo = 0.;

  ///Ratio of torsion and normal modulus
  double alphaThree = 0.;

    /// coefficient variation of the Gaussian distribution
    double coefficientOfVariation = 0.;

    /// flag which chooses between no distribution (0) and Gaussian distribution (1)
    double localRandomType = 0.;


    /// parameter which allows to prescribed thermal displacement
    double cAlpha = 0.;

    /// Poisson's ratio; triggers the Fahy/Griffiths-Mustoe shortcut for
    /// the spring ratios when supplied with `em`.
    double nu = 0.;
    bool nuWasGiven = false;

    /// Macroscopic Young's modulus (continuum E). Defaults to eNormalMean
    /// if `em` is not given; decoupled from eNormalMean when both are set.
    double emMacro = 0.;

    /// Biot coefficient — fraction of the pore (fluid) pressure that enters the
    /// effective normal stress. 0 = no poromechanical coupling, 1 = full pore
    /// pressure. Shared by the derived damage/plasticity-damage materials.
    double biotCoefficient = 0.;

public:
    LatticeLinearElastic(int n, Domain *d) : LatticeStructuralMaterial(n, d), RandomMaterialExtensionInterface() { };


    LatticeLinearElastic(int n, Domain *d, double eNormalMean, double alphaOne, double alphaTwo, double alphaThree);

    const char *giveInputRecordName() const override { return _IFT_LatticeLinearElastic_Name; }

    const char *giveClassName() const override { return "LatticeLinearElastic"; }

    void initializeFrom(const std::shared_ptr<InputRecord> &ir) override;

    FloatArrayF< 6 >giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    FloatArrayF< 6 >giveLatticeStress3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) override;

    FloatMatrixF< 6, 6 >give3dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF< 3, 3 >give2dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;


    bool hasMaterialModeCapability(MaterialMode mode) const override;


    Interface *giveInterface(InterfaceType) override;

    virtual void giveRandomParameters(FloatArray &param);

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

    MaterialStatus *giveStatus(GaussPoint *gp) const override;

    double  give(int aProperty, GaussPoint *gp) const override;

protected:
    /// Pore (fluid) pressure acting on this Gauss point's lattice element.
    /// In a coupled StaggeredProblem run it is read *live* from the dual
    /// transport element (`couplingNumbers` + the transport slave problem);
    /// when run standalone it falls back to the element's static `givePressures`
    /// field. The caller multiplies by the Biot coefficient before adding it to
    /// the normal stress. Returns 0 if no coupling/pressure applies.
    double giveCouplingPressure(GaussPoint *gp, TimeStep *tStep);

    /// Biot coefficient used to scale the pore pressure in the stress. Base
    /// implementation returns the constant `biotCoefficient`; LatticeDamage
    /// overrides it for the damage-evolving form (`btype 1`).
    virtual double computeBiot(double omega, double kappa, double le) const { return biotCoefficient; }
};
} // end namespace oofem

#endif
