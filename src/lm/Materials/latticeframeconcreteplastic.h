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

#ifndef latticeframeconcreteplastic_h
#define latticeframeconcreteplastic_h


#include "latticeframeelastic.h"
#include "cltypes.h"
#include "randommaterialext.h"
#include "strainvector.h"
#include "stressvector.h"
#include "latticematstatus.h"
///@name Input fields for LatticeFrameConcretePlastic
//@{
#define _IFT_LatticeFrameConcretePlastic_Name "LatticeFrameConcretePlastic"
#define _IFT_LatticeFrameConcretePlastic_talpha "talpha"
#define _IFT_LatticeFrameConcretePlastic_e "e"
#define _IFT_LatticeFrameConcretePlastic_n "n"
#define _IFT_LatticeFrameConcretePlastic_nx0 "nx0"
#define _IFT_LatticeFrameConcretePlastic_vy0 "vy0"
#define _IFT_LatticeFrameConcretePlastic_vz0 "vz0"
#define _IFT_LatticeFrameConcretePlastic_mx0 "mx0"
#define _IFT_LatticeFrameConcretePlastic_my0 "my0"
#define _IFT_LatticeFrameConcretePlastic_mz0 "mz0"
#define _IFT_LatticeFrameConcretePlastic_nx01 "nx01"
#define _IFT_LatticeFrameConcretePlastic_vy01 "vy01"
#define _IFT_LatticeFrameConcretePlastic_vz01 "vz01"
#define _IFT_LatticeFrameConcretePlastic_mx01 "mx01"
#define _IFT_LatticeFrameConcretePlastic_my01 "my01"
#define _IFT_LatticeFrameConcretePlastic_mz01 "mz01"
#define _IFT_LatticeFrameConcretePlastic_tol "tol"
#define _IFT_LatticeFrameConcretePlastic_iter "iter"
#define _IFT_LatticeFrameConcretePlastic_sub "sub"
#define _IFT_LatticeFrameConcretePlastic_plastic "plastic"
#define _IFT_LatticeFrameConcretePlastic_wu "wu"
#define _IFT_LatticeFrameConcretePlastic_wf "wf"
#define _IFT_LatticeFrameConcretePlastic_wftwo "wftwo"
#define _IFT_LatticeFrameConcretePlastic_qzero "qzero"
#define _IFT_LatticeFrameConcretePlastic_capmode "capmode"
#define _IFT_LatticeFrameConcretePlastic_sigmay "sigmay"

//@}

namespace oofem {
/**
 * This class implements associated Material Status to LatticeFrameConcretePlastic.
 * @authors: Gumaa Abdelrhim, Peter Grassl
 */

class LatticeFrameConcretePlasticStatus : public LatticeMaterialStatus
{
public:

    enum state_flag_values {
        LatticeFrameConcretePlastic_Elastic,
        LatticeFrameConcretePlastic_Unloading,
        LatticeFrameConcretePlastic_Plastic,
    };
    enum LatticeFrameConcretePlastic_ReturnResult {
        RR_NotConverged,
        RR_Converged,
    };

protected:

    double kappaD = 0.;
    double tempKappaD = 0.;
    double damage = 0.;
    double tempDamage = 0.;


    int tempReturnResult = LatticeFrameConcretePlasticStatus::RR_NotConverged;

public:

    /// Constructor
    LatticeFrameConcretePlasticStatus(int n, Domain *d, GaussPoint *g);


    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "LatticeFrameConcretePlasticStatus"; }

    void letTempReturnResultBe(const int result) { tempReturnResult = result; }

    int giveTempReturnResult() const { return tempReturnResult; }

    double giveKappaD() const { return this->kappaD; }

    double giveTempKappaD() const { return this->tempKappaD; }

    void   setTempKappaD(double newKappa) { this->tempKappaD = newKappa; }

    double giveDamage() const { return this->damage; }

    double giveTempDamage() const { return this->tempDamage; }

    void   setTempDamage(double newDamage) { this->tempDamage = newDamage; }

    void initTempStatus() override;

    void updateYourself(TimeStep *) override;

    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;
};


/**
 * This class implements a concrete plasticity model for concrete for frame elements.
 */
class LatticeFrameConcretePlastic : public LatticeFrameElastic

{
protected:


    struct CapacitySet {
        double nx0  = 0.0;
        double vy0  = 0.0;
        double vz0  = 0.0;
        double mx0  = 0.0;
        double my0  = 0.0;
        double mz0  = 0.0;

        double nx01 = 0.0;
        double vy01 = 0.0;
        double vz01 = 0.0;
        double mx01 = 0.0;
        double my01 = 0.0;
        double mz01 = 0.0;
    };

   
    int capacityMode = 0;
  
    ///maximum axial force in x-axis x-axis nx0
    double nx0;

    double vy0;

    double vz0;

    ///maximum  bending moment about x-axis mx0
    double mx0;

    ///maximum  bending moment about x-axis my0
    double my0;

    ///maximum  bending moment about x-axis mz0
    double mz0;

    ///maximum axial force in x-axis x-axis nx01
    double nx01;

    double vy01;

    double vz01;

    ///maximum  bending moment about x-axis mx01
    double mx01;

    ///maximum  bending moment about x-axis my01
    double my01;

    ///maximum  bending moment about x-axis mz01
    double mz01;

    /// yield tolerance
    double yieldTol;

    /// maximum number of iterations for stress return
    double newtonIter;

    ///number Of SubIncrements
    double numberOfSubIncrements;

    ///plastic flag
    double plasticFlag;

    ///wu
    double wu;

    ///wf
    double wf;

    ///wftwo
    double wftwo;

    ///qzero
    double qzero;

    double sigmay = 0.0;

    bool autoCapacitiesFromSigmay = false;

    enum LatticeFrameConcretePlastic_ReturnResult { RR_NotConverged, RR_Converged, RR_Unknown, RR_known };

    double initialYieldStress = 0.;

public:
    LatticeFrameConcretePlastic(int n, Domain *d) : LatticeFrameElastic(n, d) { };

    FloatArrayF< 6 >computeFVector(const FloatArrayF< 6 > &stress, const FloatArrayF< 6 > &k, GaussPoint *gp, TimeStep *tStep) const;

    FloatMatrixF< 6, 6 >computeDMMatrix(const FloatArrayF< 6 > &stress,  const FloatArrayF< 6 > &k, GaussPoint *gp, TimeStep *tStep) const;

    FloatArrayF< 6 >giveThermalDilatationVector(GaussPoint *gp,  TimeStep *tStep) const override;

    FloatArrayF< 6 >giveLatticeStrain(GaussPoint *gp, TimeStep *tStep) const;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime) override;
    virtual FloatArrayF< 6 >giveStrain(GaussPoint *gp, TimeStep *tStep) const;

    FloatArrayF< 6 >performPlasticityReturn(GaussPoint *gp, const FloatArrayF< 6 > &strain, TimeStep *tStep) const;

    void performRegularReturn(FloatArrayF< 6 > &stress, LatticeFrameConcretePlastic_ReturnResult,  const FloatArrayF< 6 > &k, double yieldValue, GaussPoint *gp, TimeStep *tStep) const;

    void giveSwitches(IntArray &answer, int k);

    double computeYieldValue(const FloatArrayF< 6 > &stress,  const FloatArrayF< 6 > &k, GaussPoint *gp, TimeStep *tStep) const;

    FloatMatrixF< 7, 7 >computeJacobian(const FloatArrayF< 6 > &stress, const FloatArrayF< 6 > &k, const double deltaLambda, GaussPoint *gp, TimeStep *tStep) const;

    IntArray checkStatus(const FloatArrayF< 6 > &stress, IntArray &k, GaussPoint *gp, TimeStep *tStep) const;

    FloatArrayF< 6 >giveFrameForces3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) override;

    const char *giveInputRecordName() const override { return _IFT_LatticeFrameConcretePlastic_Name; }

    const char *giveClassName() const override { return "LatticeFrameConcretePlastic"; }

    void initializeFrom(const std::shared_ptr<InputRecord> &ir) override;

    CapacitySet giveEffectiveCapacities(GaussPoint *gp, TimeStep *tStep) const;
  
    CapacitySet giveBaseCapacitiesFromSection(GaussPoint *gp) const;
  
    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    Interface *giveInterface(InterfaceType) override;

    FloatMatrixF< 6, 6 >give3dFrameStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    bool hasMaterialModeCapability(MaterialMode mode) const override;

    std::unique_ptr< MaterialStatus > CreateStatus(GaussPoint *gp) const override;

    MaterialStatus *giveStatus(GaussPoint *gp) const override;

protected:
    const FloatArray checkStatus(FloatArrayF< 6 >f, GaussPoint *pPoint, TimeStep *pStep) const;
    IntArray checkStatus(const FloatArrayF< 6 > &stress, IntArray &k) const;
    const FloatArrayF< 6 >checkTransition(const FloatArrayF< 6 > &stress, GaussPoint *gp, TimeStep *tStep) const;
};
} // end namespace oofem

#endif
