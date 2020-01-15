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

#ifndef rcm2_h
#define rcm2_h

#include "material.h"
#include "sm/Materials/linearelasticmaterial.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "intarray.h"

///@name Input fields for RCM2Material
//@{
#define _IFT_RCM2Material_gf "gf"
#define _IFT_RCM2Material_ft "ft"
//@}

namespace oofem {
// material contant's keys for give()
#define pscm_Ee 300
#define pscm_Et 301
#define pscm_Gf 302
#define pscm_Beta 303
#define pscm_G  304
#define pscm_Ft 305

// crack statuses
#define pscm_NONE 0
#define pscm_OPEN 1
#define pscm_SOFTENING 2
#define pscm_RELOADING 3
#define pscm_UNLOADING 4
#define pscm_CLOSED 5

#define pscm_NEW_CRACK    20
#define pscm_NEW_FULLY_OPEN_CRACK    21
#define pscm_REOPEN_CRACK 22

#define rcm_SMALL_STRAIN 1.e-7
#define rcm2_BIGNUMBER 1.e8

/**
 * This class implements associated Material Status to SmearedCrackingMaterail.
 */
class RCM2MaterialStatus : public StructuralMaterialStatus
{
protected:
    /// One value from (pscm_NONE, pscm_OPEN, pscm_SOFTENING, pscm_RELOADING, pscm_UNLOADING, pscm_CLOSED
    IntArray crackStatuses, tempCrackStatuses;
    /// Max crack strain reached.
    FloatArray maxCrackStrains, tempMaxCrackStrains;
    /// Components of crack strain vector.
    FloatArray crackStrainVector, oldCrackStrainVector;
    /// Storing direction of cracks in columwise format.
    FloatMatrix crackDirs, tempCrackDirs;

    FloatArray charLengths;

    // FloatArray minEffStrainsForFullyOpenCrack;
    // FloatArray *crackStrainVector, *crackStrainIncrementVector;
    // charLengths and minEffStrainsForFullyOpenCrack are not temp,
    // because they are set only if a new crack is created,
    // if a new crack is created this is recognized based on value in
    // tempCrackStatuses->at(i), so we need not create temp version
    // of theese variables.
    //
    FloatArray principalStrain, oldPrincipalStrain;
    FloatArray principalStress, oldPrincipalStress;
    IntArray crackMap;

public:
    RCM2MaterialStatus(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    const FloatArray &getPrincipalStrainVector() const { return principalStrain; }
    const FloatArray &getPrincipalStressVector() const { return principalStress; }
    const FloatArray &givePrevPrincStrainVector() const { return oldPrincipalStrain; }
    const FloatArray &givePrevPrincStressVector() const { return oldPrincipalStress; }
    void letPrincipalStrainVectorBe(FloatArray pv) { principalStrain = std :: move(pv); }
    void letPrincipalStressVectorBe(FloatArray pv) { principalStress = std :: move(pv); }

    const IntArray & giveCrackMap() const { return crackMap; }
    void letCrackMapBe(IntArray map) { crackMap = std :: move(map); }
    virtual int isCrackActive(int i) const;
    virtual int giveNumberOfActiveCracks() const;
    virtual int giveNumberOfTempActiveCracks() const;
    int giveTempAlreadyCrack() const { return this->giveNumberOfTempActiveCracks(); }

    //double giveMinCrackStrainsForFullyOpenCrack (int icrack) {return minEffStrainsForFullyOpenCrack.at(icrack);}
    //void   setMinCrackStrainsForFullyOpenCrack (int icrack, double val) {minEffStrainsForFullyOpenCrack.at(icrack) = val;}

    const FloatMatrix &giveTempCrackDirs() { return tempCrackDirs; }
    void letTempCrackDirsBe(FloatMatrix a) { tempCrackDirs = std :: move(a); }
    //const FloatArray & giveTempMaxCrackStrain() { return tempMaxCrackStrains;}
    double giveTempMaxCrackStrain(int icrack) { return tempMaxCrackStrains.at(icrack); }
    void setTempMaxCrackStrain(int icrack, double val) { tempMaxCrackStrains.at(icrack) = val; }
    const IntArray &giveTempCrackStatus() { return tempCrackStatuses; }
    int giveTempCrackStatus(int icrack) const { return tempCrackStatuses.at(icrack); }
    void setTempCrackStatus(int icrack, int val) { tempCrackStatuses.at(icrack) = val; }

    const FloatArray &giveCrackStrainVector() const { return crackStrainVector; }
    double giveCrackStrain(int icrack) const { return crackStrainVector.at(icrack); }
    const FloatArray &giveOldCrackStrainVector() { return oldCrackStrainVector; }
    void letCrackStrainVectorBe(FloatArray a) { crackStrainVector = std :: move(a); }
    void letOldCrackStrainVectorBe(FloatArray a) { oldCrackStrainVector = std :: move(a); }

    double giveCharLength(int icrack) const {
        if ( icrack ) {
            return charLengths.at(icrack);
        } else {
            return 0.0;
        }
    }
    void setCharLength(int icrack, double val) { charLengths.at(icrack) = val; }

    // query for non-tem variables (usefull for postprocessing)
    const FloatMatrix &giveCrackDirs() { return crackDirs; }
    const IntArray &giveCrackStatus() { return crackStatuses; }
    int giveAlreadyCrack() const { return this->giveNumberOfActiveCracks(); }

    const char *giveClassName() const override { return "RCM2MaterialStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};

/**
 * This class implements a Rotating Crack Model for fracture in smeared fashion
 * ( only material stiffness modification is required, no changes in
 * mesh topology) coupled with plastic behaviour.
 * In this class we follow a paper written by R. de Borst
 * "Smeared Cracking, plasticity, creep - Unified approach "
 * which allows cracking phenomena to be easily incorporated into
 * existing code for plastic materials.
 *
 * A model is based on Fracture Energy criterion, with (linear softening
 * phenomena and with shear retention factor), controlled un-re loading.
 */
class RCM2Material : public StructuralMaterial
{
protected:
    LinearElasticMaterial *linearElasticMaterial = nullptr;
    double Gf = 0., Ft = 0.;
    //double beta = 0.;

public:
    RCM2Material(int n, Domain * d);

    bool hasMaterialModeCapability(MaterialMode mode) const override;

    const char *giveClassName() const override { return "RCM2Material"; }

    void initializeFrom(InputRecord &ir) override;

    double give(int aProperty, GaussPoint *gp) const override;

    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                              const FloatArray &reducedStrain, TimeStep *tStep) override;

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<RCM2Material*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<4> giveRealStressVector_PlaneStrain(const FloatArrayF<4> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<RCM2Material*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<3> giveRealStressVector_PlaneStress(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<RCM2Material*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<1> giveRealStressVector_1d(const FloatArrayF<1> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<RCM2Material*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<2> giveRealStressVector_2dBeamLayer(const FloatArrayF<2> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<RCM2Material*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<5> giveRealStressVector_PlateLayer(const FloatArrayF<5> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<RCM2Material*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new RCM2MaterialStatus(gp); }

    FloatArrayF<6> giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const override
    { return linearElasticMaterial->giveThermalDilatationVector(gp, tStep); }

protected:

    virtual void checkForNewActiveCracks(IntArray &answer, GaussPoint *gp, const FloatArray &,
                                         const FloatArray &, FloatArray &, const FloatArray &);
    virtual void updateCrackStatus(GaussPoint *gp, const FloatArray &crackStrain);

    virtual void checkIfClosedCracks(GaussPoint *gp, FloatArray &crackStrainVector, IntArray &);
    virtual int checkSizeLimit(GaussPoint *gp, double) { return 0; }
    virtual double giveNormalCrackingStress(GaussPoint *gp, double eps_cr, int i) = 0;
    virtual double giveMinCrackStrainsForFullyOpenCrack(GaussPoint *gp, int i) = 0;
    virtual double computeStrength(GaussPoint *gp, double) = 0;
    virtual void updateStatusForNewCrack(GaussPoint *, int, double);
    virtual double giveCharacteristicElementLength(GaussPoint *gp, const FloatArray &crackPlaneNormal);
    virtual double giveCrackingModulus(MatResponseMode rMode, GaussPoint *gp,
                                       double effStrain, int i) { return 1.e20; }

    virtual void giveMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode,
                                             GaussPoint *gp,
                                             TimeStep *tStep);

    void giveCrackedStiffnessMatrix(FloatMatrix &answer,
                                    MatResponseMode rMode,
                                    GaussPoint *gp,
                                    TimeStep *tStep);
    /*
     * void  computeTrialStressIncrement (FloatArray& answer, GaussPoint *gp,
     * const FloatArray& strainIncrement, TimeStep* tStep);
     */
    virtual void giveEffectiveMaterialStiffnessMatrix(FloatMatrix &answer,
                                                      MatResponseMode rMode,
                                                      GaussPoint *gp, TimeStep *tStep);

    void giveRealPrincipalStressVector3d(FloatArray &answer, GaussPoint *,
                                         FloatArray &, FloatMatrix &, TimeStep *);
    void giveNormalElasticStiffnessMatrix(FloatMatrix &answer,
                                          bool reduce, MatResponseMode,
                                          GaussPoint *, TimeStep *tStep,
                                          const FloatMatrix &);
    void updateActiveCrackMap(GaussPoint *gp, const IntArray *activatedCracks = NULL);
    // Give3dMaterialStiffnessMatrix should return 3d material stiffness matrix
    // taking into account possible failure or fracture of material
    double giveResidualStrength() { return 0.01 * this->Ft; }


    FloatMatrixF<3,3> givePlaneStressStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<4,4> givePlaneStrainStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<1,1> give1dStressStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<2,2> give2dBeamLayerStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<5,5> givePlateLayerStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;
};
} // end namespace oofem
#endif // rcm2_h
