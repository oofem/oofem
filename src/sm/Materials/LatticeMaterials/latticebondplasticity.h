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

#ifndef latticebondplasticity_h
#define latticebondplasticity_h

#include "latticelinearelastic.h"
#include "latticematstatus.h"

///@name Input fields for LatticePlasticityDamage
//@{
#define _IFT_LatticeBondPlasticity_Name "latticebondplast"
#define _IFT_LatticeBondPlasticity_tol "tol"
#define _IFT_LatticeBondPlasticity_iter "iter"
#define _IFT_LatticeBondPlasticity_sub "sub"
#define _IFT_LatticeBondPlasticity_fc "fc"
#define _IFT_LatticeBondPlasticity_angle1 "angle1"
#define _IFT_LatticeBondPlasticity_ef "ef"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to LatticeBondPlasticity.
 */
class LatticeBondPlasticityStatus : public LatticeMaterialStatus
{
protected:

    double kappaP = 0.;
    double tempKappaP = 0.;

    int surfaceValue = 0.;

public:

    /// Constructor
    LatticeBondPlasticityStatus(int n, Domain *d, GaussPoint *g);

    /// Returns the last equilibrated scalar measure of the largest strain level
    double giveKappaP() const { return kappaP; }
    /// Returns the temp. scalar measure of the largest strain level
    double giveTempKappaP() const { return tempKappaP; }
    /// Sets the temp scalar measure of the largest strain level to given value
    void setTempKappaP(double newKappa) { tempKappaP = newKappa; }

    void setSurfaceValue(int val)
    {
        surfaceValue = val;
    }

    int giveSurfaceValue()
    {
        return surfaceValue;
    }


    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "LatticeBondPlasticityStatus"; }

    void initTempStatus() override;

    void updateYourself(TimeStep *) override;


    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;
};

class LatticeBondPlasticity : public LatticeLinearElastic
{
protected:

    enum LatticeBondPlasticity_SurfaceType { ST_Unknown, ST_Vertex, ST_Shear, ST_Compression };

    enum LatticeBondPlasticity_ReturnResult { RR_NotConverged, RR_Converged, RR_Elastic };

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime) override;

    /// compressive strength
    double fc = 0.;

    /// frictional angle of the yield surface
    double frictionAngleOne = 0.;
    double frictionAngleTwo = 0.;

    double flowAngle = 0.;

    /// yield tolerance
    double yieldTol = 0.;

    /// maximum number of iterations for stress return
    int newtonIter = 0;

    /// maximum number of subincrements
    int numberOfSubIncrements = 0;

    //parameter in hardening law
    double ef = 0.;

public:

    /// Constructor
    LatticeBondPlasticity(int n, Domain *d);

    const char *giveInputRecordName() const override { return _IFT_LatticeBondPlasticity_Name; }
    const char *giveClassName() const override { return "LatticeBondPlasticity"; }

    double computeHardening(double kappa) const;

    double computeDHardeningDKappa(double kappa) const;

    double computeParamA(const double kappa) const;
    double computeDParamADKappa(const double kappa) const;

    double computeShift(const double kappa) const;
    double computeDShiftDKappa(const double kappa) const;


    void initializeFrom(InputRecord &ir) override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }



    bool checkForVertexCase(FloatArrayF< 3 > &stress, GaussPoint *gp) const;

    void performVertexReturn(FloatArrayF< 3 > &stress, GaussPoint *gp) const;

    double computeTransition(const double kappa, GaussPoint *gp) const;

    bool hasMaterialModeCapability(MaterialMode mode) const override;


    /* Compute the A matrix used for closest point return
     */
    FloatMatrixF< 3, 3 >computeAMatrix(const FloatArrayF< 3 > &sigma,
                                       const double tempKappa,
                                       const double deltaLambda,
                                       int transitionFlag,
                                       GaussPoint *gp) const;

    FloatMatrixF< 4, 4 >computeJacobian(const FloatArrayF< 3 > &sigma,
                                        const double tempKappa,
                                        const double deltaLambda,
                                        int transitionFlag,
                                        GaussPoint *gp) const;

    FloatArrayF< 3 >computeFVector(const FloatArrayF< 3 > &sigma,
                                   const double kappa,
                                   int transitionFlag,
                                   GaussPoint *gp) const;
    FloatArrayF< 3 >computeMVector(const FloatArrayF< 3 > &sigma,
                                   const double kappa,
                                   int transitionFlag,
                                   GaussPoint *gp) const;
    FloatMatrixF< 3, 3 >computeDMMatrix(const FloatArrayF< 3 > &sigma,
                                        const double deltaLambda,
                                        int transitionFlag,
                                        GaussPoint *gp) const;

    FloatArrayF< 6 >giveLatticeStress3d(const FloatArrayF< 6 > &jump, GaussPoint *gp, TimeStep *tStep) override;


    FloatArrayF< 3 >performPlasticityReturn(GaussPoint *gp, const FloatArrayF< 3 > &totalStrain, TimeStep *) const;

    double performRegularReturn(FloatArrayF< 3 > &stress,
                                LatticeBondPlasticity_ReturnResult &returnResult,
                                LatticeBondPlasticity_SurfaceType &surfaceType,
                                double yieldValue,
                                int transitionFlag,
                                GaussPoint *gp) const;

    double computeYieldValue(const FloatArrayF< 3 > &sigma,
                             const double tempKappa,
                             int transitionFlag,
                             GaussPoint *gp) const;


    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

    virtual FloatArrayF< 6 >giveReducedStrain(GaussPoint *gp, TimeStep *tStep) const;
};
} // end namespace oofem
#endif
