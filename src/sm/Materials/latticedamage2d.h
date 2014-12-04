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
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef latticedamage2d_h
#define latticedamage2d_h

#include "material.h"
#include "Materials/linearelasticmaterial.h"
#include "latticematstatus.h"
#include "cltypes.h"
#include "../sm/Materials/structuralmaterial.h"
#include "randommaterialext.h"
#include "strainvector.h"
#include "stressvector.h"

///@name Input fields for LatticeDamage2d
//@{
#define _IFT_LatticeDamage2d_Name "latticedamage2d"
#define _IFT_LatticeDamage2d_eNormal "e"
#define _IFT_LatticeDamage2d_alphaOne "a1"
#define _IFT_LatticeDamage2d_alphaTwo "a2"
#define _IFT_LatticeDamage2d_softeningType "stype"
#define _IFT_LatticeDamage2d_wf "wf"
#define _IFT_LatticeDamage2d_wfOne "wf1"
#define _IFT_LatticeDamage2d_localrandomtype "randomtype"
#define _IFT_LatticeDamage2d_coefficientOfVariation "cov"
#define _IFT_LatticeDamage2d_e0Mean "e0"
#define _IFT_LatticeDamage2d_e0OneMean "e01"
#define _IFT_LatticeDamage2d_coh "coh"
#define _IFT_LatticeDamage2d_ec "ec"
#define _IFT_LatticeDamage2d_calpha "calpha"
#define _IFT_LatticeDamage2d_bio "bio"
#define _IFT_LatticeDamage2d_btype "btype"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to LatticeDamage2d.
 */
class LatticeDamage2dStatus : public LatticeMaterialStatus, public RandomMaterialStatusExtensionInterface
{
protected:
    /// Scalar measure of the largest strain level ever reached in material.
    double kappa;
    /// Non-equilibrated scalar measure of the largest strain level.
    double tempKappa;

    /// Scalar measure of the strain.
    double equivStrain;
    /// Non-equilibrated scalar measure of the strain.
    double tempEquivStrain;

    /// Damage level of material.
    double damage;
    /// Non-equilibrated damage level of material.
    double tempDamage;

    /// Dissipation.
    double dissipation;
    /// Non-equilibrated dissipation..
    double tempDissipation;

    /// Increment of dissipation.
    double deltaDissipation;
    /// Non-equilibrated increment of dissipation.
    double tempDeltaDissipation;

    /// Characteristic length.
    double le;

    /// Random material parameter stored in status, since each gp has a different value.
    double e0;

    /// Set biot coefficient
    double biot;

    /**
     * The crack_flag indicates if the gp is cracked:
     * crack_flag = 0 gp is uncracked
     * crack_flag = 1 gp is cracked and damage grows
     * crack_flag = 2 gp is cracked and damage does not grow
     */
    int crack_flag;

    /// Non-equilibrated temp flag.
    int temp_crack_flag;

    /// Crack width.
    double crackWidth;
    
    /// Old crack width
    double oldCrackWidth;

    /// Non-equilibrated crack width.
    double tempCrackWidth;

    /// equilibrated normal stress
    double normalStress;

    /// old normal stress
    double oldNormalStress;

    /// nonequilibrated normal stress
    double tempNormalStress;

    /// Reduced strain.
    FloatArray reducedStrain;

    /// Non-equilibrated reduced strain.
    FloatArray tempReducedStrain;

    int updateFlag;

public:

    /// Constructor
    LatticeDamage2dStatus(int n, Domain * d, GaussPoint * g);
    /// Destructor
    virtual ~LatticeDamage2dStatus() { }


    /// Returns the last equilibrated scalar measure of the largest strain level
    double giveKappa() { return kappa; }
    /// Returns the temp. scalar measure of the largest strain level
    double giveTempKappa() { return tempKappa; }
    /// Sets the temp scalar measure of the largest strain level to given value
    void setTempKappa(double newKappa) { tempKappa = newKappa; }

    /// Returns the last equilibrated scalar measure of the largest strain level
    double giveEquivalentStrain() { return equivStrain; }
    /// Returns the temp. scalar measure of the largest strain level
    double giveTempEquivalentStrain() { return tempEquivStrain; }
    /// Sets the temp scalar measure of the largest strain level to given value
    void setTempEquivalentStrain(double newEquivStrain) { tempEquivStrain = newEquivStrain; }

    /// Returns the last equilibrated dissipation
    double giveDissipation() { return dissipation; }
    /// Returns the temp. dissipation
    double giveTempDissipation() { return tempDissipation; }
    /// Sets the temp dissipation
    void setTempDissipation(double newDiss) { tempDissipation = newDiss; }

    /// Returns the last equilibrated increment of dissipation
    double giveDeltaDissipation() { return deltaDissipation; }
    /// Returns the temp. increment dissipation
    double giveTempDeltaDissipation() { return tempDeltaDissipation; }
    /// Sets the temp. increment dissipation
    void setTempDeltaDissipation(double newDiss) { tempDeltaDissipation = newDiss; }


    /// Returns the last equilibrated damage level
    double giveDamage() { return damage; }
    /// Returns the temp. damage level
    double giveTempDamage() { return tempDamage; }
    /// Sets the temp damage level to given value
    void setTempDamage(double newDamage) { tempDamage = newDamage; }

    /// Gives the temp value of plastic strain.
    const FloatArray &giveTempReducedStrain() const { return tempReducedStrain; }
    /// Gives the old equilibrated value of plastic strain.
    const FloatArray &giveReducedStrain() const { return reducedStrain; }

    /**
     * Assign the temp value of plastic strain.
     * @param v New temp value of plastic strain.
     */
    void letTempReducedStrainBe(FloatArray v)  { tempReducedStrain = std :: move(v); }

    /// Prints the receiver state to given stream
    void printOutputAt(FILE *file, TimeStep *tStep);

    /// Returns characteristic length stored in receiver
    double giveLe() { return le; }

    /// Sets characteristic length to given value
    void setLe(double ls) { le = ls; }

    ///Sets the temp_crack_flag
    void setTempCrackFlag(int val);

    /// Returns the crack_flag
    int giveCrackFlag();

    /// Set random e0
    void setE0(double val) { e0 = val; }

    /// Sets the temp crack width
    void setTempCrackWidth(double val);

    /// Gives the last equilibrated crack width
    double giveCrackWidth() { return this->crackWidth; }


    virtual double giveNormalStress() { return this->normalStress; }

    virtual double giveOldNormalStress() { return this->oldNormalStress; }

    virtual int hasBeenUpdated() { return this->updateFlag; }


    /// Sets the temp normalStress
    void setTempNormalStress(double val);

    /// Sets the old normalStress
    void setOldNormalStress(double val);

    // definition
    virtual const char *giveClassName() const { return "LatticeDamage2dStatus"; }

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *);

    virtual Interface *giveInterface(InterfaceType);

    virtual void setVariableInStatus(double variable);

    void setBiotCoefficientInStatus(double variable);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
};


/**
 * This class implements a local random isotropic damage model for concrete in tension for 2D lattice elements.
 */
class LatticeDamage2d : public StructuralMaterial, public RandomMaterialExtensionInterface
{
protected:

    /// Normal modulus
    double eNormal;
    /// Shear modulus
    double eShear;
    /// Torsion modulus
    double eTorsion;
    /// Ratio of shear and normal modulus
    double alphaOne;
    /// Ratio of torsion and normal modulus
    double alphaTwo;

    /// Mean effective strain at peak
    double e0Mean;

    /// Mean effective strain at sigma1
    double e0OneMean;

    ///coefficient used for modelling eigendisplacements
    double cAlpha;

    /**
     * Parameter which determines the typ of the softeningFunction
     * 1 = linear softening
     * 2 = bilinear softening
     * 3 = exponential softening
     */
    int softeningType;

    /// Determines the softening -> corresponds to threshold of crack opening (not strain)
    double wf, wfOne;

    /// Parameter for the elliptic equivalent strain function
    double ec;

    /// Parameter setting ratio of shear and tensile strength
    double coh;

    /// Parameter controlling the amount of fluid pressure added to the mechanical stress (Biot's coefficient)
    double biotCoefficient;

    /// Parameter specifying how the biot coefficient changes with the crack opening
    int biotType;

    /// Coefficient variation of the Gaussian distribution
    double coefficientOfVariation;

    /// Flag which chooses between no distribution (0) and Gaussian distribution (1)
    double localRandomType;

public:

    /// Constructor
    LatticeDamage2d(int n, Domain * d);
    /// Destructor
    virtual ~LatticeDamage2d();

    virtual const char *giveInputRecordName() const { return _IFT_LatticeDamage2d_Name; }
    virtual const char *giveClassName() const { return "LatticeDamage2d"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual void  giveStiffnessMatrix(FloatMatrix &answer,
                                      MatResponseMode mode,
                                      GaussPoint *gp,
                                      TimeStep *tStep);

    virtual void computeStressIndependentStrainVector(FloatArray &answer,
                                                      GaussPoint *gp,
                                                      TimeStep *tStep,
                                                      ValueModeType mode);

    virtual void giveSecantStiffnessMatrix(FloatMatrix &answer,
                                   GaussPoint *gp,
                                   TimeStep *tStep);


    virtual void giveElasticStiffnessMatrix(FloatMatrix &answer,
                                    GaussPoint *gp,
                                    TimeStep *tStep);

    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);

    virtual void computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp);

    virtual double computeBiot(double omega, double kappa, double le);

    virtual Interface *giveInterface(InterfaceType);


    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *,
                                      const FloatArray &, TimeStep *);

    virtual void giveRandomParameters(FloatArray &param);


    ///Compute increment of dissipation for post-processing reasons
    double computeDeltaDissipation(double omega, FloatArray &reducedStrain, GaussPoint *gp, TimeStep *tStep);

    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp,  TimeStep *tStep);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual MaterialStatus *giveStatus(GaussPoint *gp) const;

    virtual double give(int aProperty, GaussPoint *gp);


protected:

    virtual int giveIPValue(FloatArray &answer,
                            GaussPoint *gp,
                            InternalStateType type,
                            TimeStep *tStep);
};
} // end namespace oofem
#endif
