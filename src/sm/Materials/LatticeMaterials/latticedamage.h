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

#ifndef latticedamage_h
#define latticedamage_h

#include "latticelinearelastic.h"
#include "latticematstatus.h"
#include "randommaterialext.h"

///@name Input fields for LatticeDamage
//@{
#define _IFT_LatticeDamage_Name "latticedamage"
#define _IFT_LatticeDamage_softeningType "stype"
#define _IFT_LatticeDamage_wf "wf"
#define _IFT_LatticeDamage_wfOne "wf1"
#define _IFT_LatticeDamage_e0Mean "e0"
#define _IFT_LatticeDamage_e0OneMean "e01"
#define _IFT_LatticeDamage_coh "coh"
#define _IFT_LatticeDamage_ec "ec"
#define _IFT_LatticeDamage_bio "bio"

#define _IFT_LatticeDamage_btype "btype"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to LatticeDamage.
 */
class LatticeDamageStatus : public LatticeMaterialStatus
{
protected:


    /// scalar measure of the largest strain level ever reached in material
    double kappa = 0.;

    /// non-equilibrated scalar measure of the largest strain level
    double tempKappa = 0.;

    /// scalar measure of the largest strain level ever reached in material
    double equivStrain = 0.;

    /// non-equilibrated scalar measure of the largest strain level
    double tempEquivStrain = 0.;

    /// damage level of material
    double damage = 0.;

    /// non-equilibrated damage level of material
    double tempDamage = 0.;

    /// random material parameter stored in status, since each gp has a differnet value.
    double e0 = 0.;

    /// computed biot coefficient
    double biot = 0.;

public:
    LatticeDamageStatus(GaussPoint *g);

    /// Returns the last equilibrated scalar measure of the largest strain level
    double giveKappa() const { return kappa; }
    /// Returns the temp. scalar measure of the largest strain level
    double giveTempKappa() const { return tempKappa; }

    /// Sets the temp scalar measure of the largest strain level to given value
    void   setTempKappa(double newKappa) { tempKappa = newKappa; }

    /// Returns the last equilibrated scalar measure of the largest strain level
    double giveEquivalentStrain() const { return equivStrain; }
    /// Returns the temp. scalar measure of the largest strain level
    double giveTempEquivalentStrain() const { return tempEquivStrain; }
    /// Sets the temp scalar measure of the largest strain level to given value
    void   setTempEquivalentStrain(double newEquivStrain) { tempEquivStrain = newEquivStrain; }


    /// Returns the last equilibrated damage level
    double giveDamage() const { return damage; }
    /// Returns the temp. damage level
    double giveTempDamage() const { return tempDamage; }
    /// Sets the temp damage level to given value
    void   setTempDamage(double newDamage) { tempDamage = newDamage; }

    void printOutputAt(FILE *file, TimeStep *tStep) const override;


    const char *giveClassName() const override { return "LatticeDamageStatus"; }

    void initTempStatus() override;

    void updateYourself(TimeStep *) override;

    ///Set random e0
    void setE0(double val) { e0 = val; }

    void setBiotCoefficientInStatus(double variable) { biot = variable; }

    void saveContext(DataStream &stream, ContextMode mode) override;

    void restoreContext(DataStream &stream, ContextMode mode) override;
};


/**
 * This class implements a local random damage model for quasi-brittle materials for lattice (1D, 2D and 3D) elements.
 */
class LatticeDamage : public LatticeLinearElastic
{
protected:


    /// max effective strain at peak
    double e0Mean = 0.;
    double e0OneMean = 0.;

    ///tensile strength
    double ftMean = 0.;
    double ftOneMean = 0.;

    /**parameter which determines the typ of the softeningFunction
     * 1 = linear softening
     * 2 = bilinear softening
     * 3 = exponential softening
     **/
    int softeningType = 0.;

    /// determines the softening -> corresponds to crack opening when tension stress vanishes
    double wf = 0., wfOne = 0.;

    double ultimateFactor = 0.;

    //parameter for the elliptic equivalent strain function
    double coh = 0.;

    //parameter for the elliptic equivalent strain function
    double ec = 0.;

    /// flag which chooses between no distribution (0) and Gaussian distribution (1)
    double localRandomType = 0.;

    /// Biot's coefficient
    double biotCoefficient = 0.;

    /// Parameter specifying how the biot coefficient changes with the crack opening
    int biotType = 0;


public:
    LatticeDamage(int n, Domain *d);

    const char *giveInputRecordName() const override { return _IFT_LatticeDamage_Name; }
    const char *giveClassName() const override { return "LatticeDamage"; }

    bool hasAnalyticalTangentStiffness() const override { return true; }

    void initializeFrom(InputRecord &ir) override;


    FloatArrayF< 6 >giveLatticeStress3d(const FloatArrayF< 6 > &jump, GaussPoint *gp, TimeStep *tStep) override;

    FloatMatrixF< 6, 6 >give3dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF< 3, 3 >give2dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;


    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    bool hasMaterialModeCapability(MaterialMode mode) const override;

    void performDamageEvaluation(GaussPoint *gp, FloatArrayF< 6 > &reducedStrain) const;

    virtual double computeEquivalentStrain(const FloatArrayF< 6 > &strain, GaussPoint *gp) const;

    virtual double computeBiot(double omega, double kappa, double le) const;

    virtual double computeDamageParam(double kappa, GaussPoint *gp) const;
    ///Compute increment of dissipation for post-processing reasons
    double computeDeltaDissipation2d(double omega, const FloatArrayF< 3 > &reducedStrain, GaussPoint *gp, TimeStep *atTime) const;

    double computeDeltaDissipation3d(double omega, const FloatArrayF< 6 > &reducedStrain, GaussPoint *gp, TimeStep *atTime) const;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

    double give(int aProperty, GaussPoint *gp) const override;


protected:
    double computeReferenceGf(GaussPoint *gp) const;
    double computeIntervals(double testDissipation, double referenceGf) const;

    int giveIPValue(FloatArray &answer,
                    GaussPoint *gp,
                    InternalStateType type,
                    TimeStep *atTime) override;
};
} // end namespace oofem

#endif
