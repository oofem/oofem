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

#ifndef latticeframesteelplastic_h
#define latticeframesteelplastic_h

#include "latticestructuralmaterial.h"
#include "cltypes.h"
#include "randommaterialext.h"
#include "strainvector.h"
#include "stressvector.h"
#include "latticematstatus.h"

///@name Input fields for LatticeFrameSteelPlastic
//@{
#define _IFT_LatticeFrameSteelPlastic_Name "LatticeFrameSteelPlastic"
#define _IFT_LatticeFrameSteelPlastic_talpha "talpha"
#define _IFT_LatticeFrameSteelPlastic_e "e"
#define _IFT_LatticeFrameSteelPlastic_n "n"
#define _IFT_LatticeFrameSteelPlastic_nx0 "nx0"
#define _IFT_LatticeFrameSteelPlastic_mx0 "mx0"
#define _IFT_LatticeFrameSteelPlastic_my0 "my0"
#define _IFT_LatticeFrameSteelPlastic_mz0 "mz0"
#define _IFT_LatticeFrameSteelPlastic_tol "tol"
#define _IFT_LatticeFrameSteelPlastic_iter "iter"
#define _IFT_LatticeFrameSteelPlastic_sub "sub"
#define _IFT_LatticeFrameSteelPlastic_plastic "plastic"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to LatticeFrameSteelPlastic.
 * @authors: Gumaa Abdelrhim, Peter Grassl
 */

class LatticeFrameSteelPlasticStatus : public LatticeMaterialStatus
{
public:

    enum state_flag_values {
        LatticeFrameSteelPlastic_Elastic,
        LatticeFrameSteelPlastic_Unloading,
        LatticeFrameSteelPlastic_Plastic,
    };


    enum LatticeFrameSteelPlastic_ReturnResult {
        RR_NotConverged,
        RR_Converged
    };


protected:

    int tempReturnResult = LatticeFrameSteelPlasticStatus::RR_NotConverged;



public:

    /// Constructor
    LatticeFrameSteelPlasticStatus(int n, Domain *d, GaussPoint *g);


    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "LatticeFrameSteelPlasticStatus"; }

    void letTempReturnResultBe(const int result) { tempReturnResult = result; }

    int giveTempReturnResult() const { return tempReturnResult; }
};


/**
 * This class implements a local random linear elastic model for lattice elements.
 */
class LatticeFrameSteelPlastic : public LatticeStructuralMaterial
{
protected:

    ///Normal modulus
    double e;

    ///Ratio of shear and normal modulus
    double nu;

    ///maximum axial force in x-axis x-axis nx0
    double nx0;

    ///maximum  bending moment about x-axis mx0
    double mx0;

    ///maximum  bending moment about x-axis my0
    double my0;

    ///maximum  bending moment about x-axis mz0
    double mz0;

    /// yield tolerance
    double yieldTol;

    /// maximum number of iterations for stress return
    double newtonIter;

    ///number Of SubIncrements
    double numberOfSubIncrements;

    ///plastic flag
    double plasticFlag;

    enum LatticeFrameSteelPlastic_ReturnResult { RR_NotConverged, RR_Converged };
    //   mutable LatticeFrameSteelPlastic_ReturnResult returnResult = RR_NotConverged; /// FIXME: This must be removed. Not thread safe. Shouldn't be stored at all.

    double initialYieldStress = 0.;

    //

public:
    LatticeFrameSteelPlastic(int n, Domain *d) : LatticeStructuralMaterial(n, d) { };

    FloatArrayF< 4 >computeFVector(const FloatArrayF< 4 > &sigma, GaussPoint *gp, TimeStep *tStep) const;

    FloatMatrixF< 4, 4 >computeDMMatrix(const FloatArrayF< 4 > &sigma, GaussPoint *gp, TimeStep *tStep) const;

    FloatArrayF< 6 >giveThermalDilatationVector(GaussPoint *gp,  TimeStep *tStep) const override;

    FloatArrayF< 6 >giveReducedLatticeStrain(GaussPoint *gp, TimeStep *tStep) const;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime) override;

    virtual FloatArrayF< 6 >giveReducedStrain(GaussPoint *gp, TimeStep *tStep) const;

    FloatArrayF< 6 >performPlasticityReturn(GaussPoint *gp, const FloatArrayF< 6 > &reducedStrain, TimeStep *tStep) const;

    void performRegularReturn(FloatArrayF< 4 > &stress, double yieldValue, GaussPoint *gp, TimeStep *tStep) const;

    double computeYieldValue(const FloatArrayF< 4 > &sigma, GaussPoint *gp, TimeStep *tStep) const;

    FloatMatrixF< 5, 5 >computeJacobian(const FloatArrayF< 4 > &sigma, const double deltaLambda, GaussPoint *gp, TimeStep *tStep) const;

    FloatArrayF< 6 >giveFrameForces3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) override;

    const char *giveInputRecordName() const override { return _IFT_LatticeFrameSteelPlastic_Name; }

    const char *giveClassName() const override { return "LatticeFrameSteelPlastic"; }

    void initializeFrom(InputRecord &ir) override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    Interface *giveInterface(InterfaceType) override;

    FloatMatrixF< 6, 6 >give3dFrameStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    bool hasMaterialModeCapability(MaterialMode mode) const override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    MaterialStatus *giveStatus(GaussPoint *gp) const override;

protected:
};
} // end namespace oofem

#endif
