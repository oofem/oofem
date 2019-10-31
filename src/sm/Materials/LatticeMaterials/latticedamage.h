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
    double kappa;

    /// non-equilibrated scalar measure of the largest strain level
    double tempKappa;

    /// scalar measure of the largest strain level ever reached in material
    double equivStrain;

    /// non-equilibrated scalar measure of the largest strain level
    double tempEquivStrain;

    /// damage level of material
    double damage;

    /// non-equilibrated damage level of material
    double tempDamage;

    /// random material parameter stored in status, since each gp has a differnet value.
    double e0;

    /// computed biot coefficient
    double biot;

public:

    /// Constructor
    LatticeDamageStatus(GaussPoint *g);
    /// Destructor
    ~LatticeDamageStatus() {}


    /// Returns the last equilibrated scalar measure of the largest strain level
    double giveKappa() { return kappa; }
    /// Returns the temp. scalar measure of the largest strain level
    double giveTempKappa() { return tempKappa; }
    /// Sets the temp scalar measure of the largest strain level to given value
    void   setTempKappa(double newKappa) { tempKappa = newKappa; }

    /// Returns the last equilibrated scalar measure of the largest strain level
    double giveEquivalentStrain() { return equivStrain; }
    /// Returns the temp. scalar measure of the largest strain level
    double giveTempEquivalentStrain() { return tempEquivStrain; }
    /// Sets the temp scalar measure of the largest strain level to given value
    void   setTempEquivalentStrain(double newEquivStrain) { tempEquivStrain = newEquivStrain; }


    /// Returns the last equilibrated damage level
    double giveDamage() { return damage; }
    /// Returns the temp. damage level
    double giveTempDamage() { return tempDamage; }
    /// Sets the temp damage level to given value
    void   setTempDamage(double newDamage) { tempDamage = newDamage; }


    void printOutputAt(FILE *file, TimeStep *tStep);


    const char *giveClassName() const override { return "LatticeDamageStatus"; }

    void initTempStatus();

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
    double e0Mean;

    double e0OneMean;

    /**parameter which determines the typ of the softeningFunction
     * 1 = linear softening
     * 2 = bilinear softening
     * 3 = exponential softening
     **/
    int softeningType;

    /// determines the softening -> corresponds to crack opening when tension stress vanishes
    double wf, wfOne;

    double ultimateFactor;

    //parameter for the elliptic equivalent strain function
    double coh;

    //parameter for the elliptic equivalent strain function
    double ec;

    /// flag which chooses between no distribution (0) and Gaussian distribution (1)
    double localRandomType;

    double biotCoefficient;

    /// Parameter specifying how the biot coefficient changes with the crack opening
    int biotType;


public:

    /// Constructor
    LatticeDamage(int n, Domain *d);

    /// Destructor
    virtual ~LatticeDamage();

    const char *giveInputRecordName() const override { return _IFT_LatticeDamage_Name; }
    const char *giveClassName() const override { return "LatticeDamage"; }

    virtual bool hasAnalyticalTangentStiffness() const { return true; }

    IRResultType initializeFrom(InputRecord *ir) override;

    virtual FloatArrayF< 3 >giveLatticeStress2d(const FloatArrayF< 3 > &jump, GaussPoint *gp, TimeStep *tStep) override;

    virtual FloatArrayF< 6 >giveLatticeStress3d(const FloatArrayF< 6 > &jump, GaussPoint *gp, TimeStep *tStep) override;

    virtual FloatMatrixF< 3, 3 >give2dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    virtual FloatMatrixF< 6, 6 >give3dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }



    int hasMaterialModeCapability(MaterialMode mode);


    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime);

    virtual double computeBiot(double omega, double kappa, double le);

    virtual void computeDamageParam(double &omega, double kappa, GaussPoint *gp);
    ///Compute increment of dissipation for post-processing reasons
    double computeDeltaDissipation2d(double omega, FloatArray &reducedStrain, GaussPoint *gp, TimeStep *atTime);

    double computeDeltaDissipation3d(double omega, FloatArray &reducedStrain, GaussPoint *gp, TimeStep *atTime);

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    double give(int aProperty, GaussPoint *gp) const;


protected:

    int giveIPValue(FloatArray &answer,
                    GaussPoint *gp,
                    InternalStateType type,
                    TimeStep *atTime) override;
};
} // end namespace oofem

#endif
