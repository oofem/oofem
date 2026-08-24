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

#ifndef latticeframesteelplastic_h
#define latticeframesteelplastic_h

#include "latticeframeelastic.h"
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
#define _IFT_LatticeFrameSteelPlastic_h "h"
#define _IFT_LatticeFrameSteelPlastic_hlength "hlength"
#define _IFT_LatticeFrameSteelPlastic_htype "htype"
#define _IFT_LatticeFrameSteelPlastic_h_eps "h_eps"
#define _IFT_LatticeFrameSteelPlastic_h_function_eps "h(eps)"
#define _IFT_LatticeFrameSteelPlastic_capmode "capmode"
#define _IFT_LatticeFrameSteelPlastic_sigmay "sigmay"

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

    double kappa = 0.;

    double tempKappa = 0.;

public:

    /// Constructor
    LatticeFrameSteelPlasticStatus(int n, Domain *d, GaussPoint *g);


    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "LatticeFrameSteelPlasticStatus"; }

    void letTempReturnResultBe(const int result) { tempReturnResult = result; }

    int giveTempReturnResult() const { return tempReturnResult; }

    void updateYourself(TimeStep *atTime) override;

    void initTempStatus() override;

    double giveKappa() const { return kappa; }
    void   setTempKappa(double newKappa) { tempKappa = newKappa; }
};


/**
 * This class implements a local random linear elastic model for lattice elements.
 */
class LatticeFrameSteelPlastic : public LatticeFrameElastic
{
protected:


    ///capacity mode to allow for shell mode 
    int capacityMode = 0;

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

    ///maximum number Of SubIncrements
    double numberOfSubIncrements;

    double hardeningLength;

    /// type of hardening function
    int hType;

    /// user-defined hardening (yield stress - kappa)
    FloatArray h_eps, h_function_eps;

    /// Hardening modulus.
    double H = 0.;

  

  
    enum LatticeFrameSteelPlastic_ReturnResult { RR_NotConverged, RR_Converged };

    double initialYieldStress = 0.;

    double sigmay = 0.0;
  
    bool autoCapacitiesFromSigmay = false;

public:
    LatticeFrameSteelPlastic(int n, Domain *d) : LatticeFrameElastic(n, d) { };

    FloatArrayF< 5 >computeFVector(const FloatArrayF< 4 > &sigma, const double kappa, GaussPoint *gp, TimeStep *tStep) const;
    FloatArrayF< 5 >computeMVector(const FloatArrayF< 4 > &sigma, const double kappa, GaussPoint *gp, TimeStep *tStep) const;

    FloatMatrixF< 5, 5 >computeDMMatrix(const FloatArrayF< 4 > &sigma, const double kappa, GaussPoint *gp, TimeStep *tStep) const;

    double computeDHardeningDKappa(const double kappa) const;

    FloatArrayF< 6 >giveThermalDilatationVector(GaussPoint *gp,  TimeStep *tStep) const override;

    FloatArrayF< 6 >giveReducedLatticeStrain(GaussPoint *gp, TimeStep *tStep) const;

  void giveEffectiveCapacities(double &nx0Eff, double &mx0Eff, double &my0Eff, double &mz0Eff, GaussPoint *gp, TimeStep *tStep) const;
  
    virtual FloatArrayF< 6 >giveReducedStrain(GaussPoint *gp, TimeStep *tStep) const;

    FloatArrayF< 6 >performPlasticityReturn(GaussPoint *gp, const FloatArrayF< 6 > &reducedStrain, TimeStep *tStep) const;

    double performRegularReturn(FloatArrayF< 4 > &stress, double yieldValue, GaussPoint *gp, TimeStep *tStep) const;

    double computeYieldValue(const FloatArrayF< 4 > &sigma, const double kappa, GaussPoint *gp, TimeStep *tStep) const;

    FloatMatrixF< 6, 6 >computeJacobian(const FloatArrayF< 4 > &sigma, const double deltaLambda, double kappa, GaussPoint *gp, TimeStep *tStep) const;

    FloatArrayF< 6 >giveFrameForces3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) override;

    void giveBaseCapacitiesFromSection(double &nx0Base, double &mx0Base, double &my0Base, double &mz0Base, GaussPoint *gp) const;
  
    const char *giveInputRecordName() const override { return _IFT_LatticeFrameSteelPlastic_Name; }

    const char *giveClassName() const override { return "LatticeFrameSteelPlastic"; }

    void initializeFrom(const std::shared_ptr<InputRecord> &ir) override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    Interface *giveInterface(InterfaceType) override;

    bool hasMaterialModeCapability(MaterialMode mode) const override;

    double computeHardening(const double kappa) const;

    std::unique_ptr< MaterialStatus > CreateStatus(GaussPoint *gp) const override;

    MaterialStatus *giveStatus(GaussPoint *gp) const override;

protected:
};
} // end namespace oofem

#endif
