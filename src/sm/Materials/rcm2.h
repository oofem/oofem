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
#include "Materials/linearelasticmaterial.h"
#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
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
    RCM2MaterialStatus(int n, Domain * d, GaussPoint * g);
    virtual ~RCM2MaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

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

    // definition
    virtual const char *giveClassName() const { return "RCM2MaterialStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    // saves current context(state) into stream
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
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
    LinearElasticMaterial *linearElasticMaterial;
    double Gf, Ft;
    //double beta;

public:
    RCM2Material(int n, Domain * d);
    virtual ~RCM2Material();

    // identification and auxiliary functions
    virtual int hasNonLinearBehaviour() { return 1; }
    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual const char *giveClassName() const { return "RCM2Material"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual double give(int aProperty, GaussPoint *gp);

    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);

    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_2dBeamLayer(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_PlateLayer(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }


    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new RCM2MaterialStatus(1, domain, gp); }

protected:

    virtual void checkForNewActiveCracks(IntArray &answer, GaussPoint *gp, const FloatArray &,
                                         const FloatArray &, FloatArray &, const FloatArray &);
    virtual void updateCrackStatus(GaussPoint *gp, const FloatArray &crackStrain);

    virtual void checkIfClosedCracks(GaussPoint *gp, FloatArray &crackStrainVector, IntArray &);
    virtual int  checkSizeLimit(GaussPoint *gp, double) { return 0; }
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
    FloatMatrix *GiveCrackTransformationMtrx(GaussPoint *gp, int i);
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


    virtual void givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode mmode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);
    virtual void givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseMode mmode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);
    virtual void give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode mmode,
                                       GaussPoint *gp,
                                       TimeStep *tStep);
    virtual void give2dBeamLayerStiffMtrx(FloatMatrix &answer, MatResponseMode mmode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);
    virtual void givePlateLayerStiffMtrx(FloatMatrix &answer, MatResponseMode mmode,
                                         GaussPoint *gp,
                                         TimeStep *tStep);
};
} // end namespace oofem
#endif // rcm2_h
