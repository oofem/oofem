/* $Header: /home/cvs/bp/oofem/sm/src/rcm2.h,v 1.5 2003/04/06 14:08:31 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   ****************************************************
//   *** CLASS ROTATING SMEARED CRACK MODEL  ************
//   ****************************************************

#ifndef rcm2_h
#define rcm2_h

#include "material.h"
#include "linearelasticmaterial.h"
#include "structuralmaterial.h"
#include "structuralms.h"

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

class GaussPoint;


//class PlasticSmearedCrackingMaterialStatus : public PerfectlyPlasticMaterialStatus
class RCM2MaterialStatus : public StructuralMaterialStatus
{
    /*
     * This class implements associated Material Status to SmearedCrackingMaterail.
     * It is atribute of matStatusDictionary at every GaussPoint, for which this material
     * is active.
     * DESCRIPTION:
     * Idea used there is that we have variables
     * describing:
     * 1) state at previous equilibrium state (variables without temp)
     * 2) state during searching new equilibrium (variables with temp)
     * when we start search new state from previous equilibrium one we copy
     * non-tem variables into temp ones. And after we reach new equilibrium
     * (now decribed by temp variables) we copy tem-var into non-tepm ones
     * (see function updateYourself).
     *
     * variables description:
     *
     * alreadyCrack - non zero if material already cracked during whole process.
     * crackStatuses- one value from (pscm_NONE, pscm_OPEN , pscm_SOFTENING, pscm_RELOADING,
     * pscm_UNLOADING, pscm_CLOSED).
     * maxCrackStrains - max  Crack Strain reached;
     * maxTotalStrains - max  Total strain reached (means elastic+cracking strain).
     * reachedSofteningStress - reched softening stress
     * crackDirs - current crack direction in gp.
     * minEffStrainsForFullyOpenCrack - for each crack (if exist) treshold value for
     * fully open crack. (Elastic+cracking).
     * charLengths - characteristic Lengthsfor each active crack.
     * crackStrainVector, crackStrainIncrementVector - components of crack strain vector in gp.
     *
     * TASK:
     *
     */

protected:

    IntArray crackStatuses, tempCrackStatuses;
    FloatArray maxCrackStrains, tempMaxCrackStrains;
    FloatArray crackStrainVector, oldCrackStrainVector;
    FloatMatrix crackDirs, tempCrackDirs;
    // floatMatrix storing direction of cracks in columwise format
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
    RCM2MaterialStatus(int n, Domain *d, GaussPoint *g);
    ~RCM2MaterialStatus();

    void   printOutputAt(FILE *file, TimeStep *tStep);

    // saves current context(state) into stream
    // contextIOResultType    saveContext (FILE* stream, void *obj = NULL);
    // contextIOResultType    restoreContext(FILE* stream, void *obj = NULL);
    void getPrincipalStrainVector(FloatArray &answer) const { answer =  principalStrain; }
    void getPrincipalStressVector(FloatArray &answer) const { answer =  principalStress; }
    void givePrevPrincStrainVector(FloatArray &answer) const { answer = oldPrincipalStrain; }
    void givePrevPrincStressVector(FloatArray &answer) const { answer = oldPrincipalStress; }
    void        letPrincipalStrainVectorBe(const FloatArray &pv) {
        principalStrain = pv;
    }
    void        letPrincipalStressVectorBe(const FloatArray &pv) {
        principalStress = pv;
    }

    void giveCrackMap(IntArray &answer) const { answer = crackMap; }
    void letCrackMapBe(IntArray &map) { crackMap = map; }
    virtual int          isCrackActive(int i) const;
    virtual int          giveNumberOfActiveCracks() const;
    virtual int    giveNumberOfTempActiveCracks() const;
    int    giveTempAlreadyCrack() const { return this->giveNumberOfTempActiveCracks(); }

    //double giveMinCrackStrainsForFullyOpenCrack (int icrack) {return minEffStrainsForFullyOpenCrack.at(icrack);}
    //void   setMinCrackStrainsForFullyOpenCrack (int icrack, double val) {minEffStrainsForFullyOpenCrack.at(icrack) = val;}

    void giveTempCrackDirs(FloatMatrix &answer) { answer =  tempCrackDirs; }
    void letTempCrackDirsBe(const FloatMatrix &a) { tempCrackDirs = a; }
    //void giveTempMaxCrackStrain(FloatArray& answer) {answer = tempMaxCrackStrains;}
    double giveTempMaxCrackStrain(int icrack) { return tempMaxCrackStrains.at(icrack); }
    void   setTempMaxCrackStrain(int icrack, double val) { tempMaxCrackStrains.at(icrack) = val; }
    void giveTempCrackStatus(IntArray &answer) { answer = tempCrackStatuses; }
    int  giveTempCrackStatus(int icrack) const { return tempCrackStatuses.at(icrack); }
    void setTempCrackStatus(int icrack, int val) { tempCrackStatuses.at(icrack) = val; }

    void giveCrackStrainVector(FloatArray &answer) { answer = crackStrainVector; }
    double giveCrackStrain(int icrack) const { return crackStrainVector.at(icrack); }
    void giveOldCrackStrainVector(FloatArray &answer) { answer = oldCrackStrainVector; }
    void letCrackStrainVectorBe(const FloatArray &a)
    { crackStrainVector = a; }
    void letOldCrackStrainVectorBe(const FloatArray &a)
    { oldCrackStrainVector = a; }



    double giveCharLength(int icrack) const { if ( icrack ) { return charLengths.at(icrack); } else { return 0.0; } }
    void   setCharLength(int icrack, double val) { charLengths.at(icrack) = val; }

    // query for non-tem variables (usefull for postprocessing)
    void giveCrackDirs(FloatMatrix &answer) { answer = crackDirs; }
    void giveCrackStatus(IntArray &answer) { answer = crackStatuses; }
    int  giveAlreadyCrack() const { return this->giveNumberOfActiveCracks(); }


    // definition
    const char *giveClassName() const { return "RCM2MaterialStatus"; }
    classType             giveClassID() const { return RCMMaterialStatusClass; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *); // update after new equilibrium state reached

    // saves current context(state) into stream
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};



class RCM2Material : public StructuralMaterial
{
    /*
     *
     * DESCRIPTION
     * This class implements a Rotating Crack Model for fracture in smeared fashion
     * ( only material stiffness modification is required, no changes in
     * mesh topology) coupled with plastic behaviour.
     * In this class we follow a paper written by R. de Borst
     * "Smeared Crcking, plasticity, creep - Unified approach "
     * which allows cracking phenomena to be easily incorporated into
     * existing code for plastic materials.
     *
     * A model is based on Fracture Energy criterion, with (linear softening
     * phenomena and with shear retention factor), controlled un-re loading.
     *
     * TASK
     * - Returning standard material stiffness and flexibility marices for 3d-case.
     * according to current state determined by using data stored
     * in Gausspoint.
     * - Returning a material property (method 'give'). Only for non-standard elements.
     * - Returning real stress state vector(tensor) at gauss point for 3d - case.
     * - Storing & restoring Material Status sored in gp matStatusDictionary.
     */

protected:
    LinearElasticMaterial *linearElasticMaterial;
    double Gf, Ft;
    // double beta;

public:

    RCM2Material(int n, Domain *d);
    ~RCM2Material();

    // identification and auxiliary functions
    int hasNonLinearBehaviour()   { return 1; }
    int hasMaterialModeCapability(MaterialMode mode);
    const char *giveClassName() const { return "RCM2Material"; }
    classType giveClassID()         const { return RCMMaterialClass; }

    // contextIOResultType    saveContext (FILE* stream, void *obj = NULL);
    // contextIOResultType    restoreContext(FILE* stream, void *obj = NULL);
    IRResultType initializeFrom(InputRecord *ir);
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    // non-standard - returns time independent material constant
    double   give(int, GaussPoint*);

    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix & answer,
                                               MatResponseForm, MatResponseMode,
                                               GaussPoint * gp,
                                               TimeStep * atTime);

    void giveRealStressVector(FloatArray & answer,  MatResponseForm, GaussPoint *,
                              const FloatArray &, TimeStep *);


#ifdef __OOFEG
#endif

    /**
     * Returns the integration point corresponding value in Reduced form.
     * @param answer contain corresponding ip value, zero sized if not available
     * @param aGaussPoint integration point
     * @param type determines the type of internal variable
     * @param type determines the type of internal variable
     * @returns nonzero if ok, zero if var not supported
     */
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    /**
     * Returns the mask of reduced indexes of Internal Variable component .
     * @param answer mask of Full VectorSize, with components beeing the indexes to reduced form vectors.
     * @param type determines the internal variable requested (physical meaning)
     * @returns nonzero if ok or error is generated for unknown mat mode.
     */
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);


    /**
     * Returns the type of internal variable (scalar, vector, tensor,...).
     * @param type determines the type of internal variable
     * @returns type of internal variable
     */
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    /**
     * Returns the corresponding integration point  value size in Reduced form.
     * @param type determines the type of internal variable
     * @returns var size, zero if var not supported
     */
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);

    MaterialStatus *CreateStatus(GaussPoint *gp) const { return new RCM2MaterialStatus(1, domain, gp); }

protected:

    // two functions used to initialize and updating temporary variables in
    // gp's status. These variables are used to controll process, when
    // we try to find equlibrium state.
    virtual void initTempStatus(GaussPoint *gp);
    virtual void checkForNewActiveCracks(IntArray &answer, GaussPoint *gp, const FloatArray &,
                                         const FloatArray &, FloatArray &, const FloatArray &);
    virtual void updateCrackStatus(GaussPoint *gp, const FloatArray &crackStrain);

    virtual void checkIfClosedCracks(GaussPoint *gp, FloatArray &crackStrainVector, IntArray &);
    virtual int  checkSizeLimit(GaussPoint *gp, double) { return 0; }
    virtual double giveNormalCrackingStress(GaussPoint *gp, double eps_cr, int i) = 0;
    virtual double giveMinCrackStrainsForFullyOpenCrack(GaussPoint *gp, int i) = 0;
    virtual double computeStrength(GaussPoint *, double) = 0;
    virtual void   updateStatusForNewCrack(GaussPoint *, int, double);
    virtual double giveCharacteristicElementLenght(GaussPoint *, const FloatArray &);
    virtual double giveCrackingModulus(MatResponseMode rMode, GaussPoint *gp,
                                       double effStrain, int i) { return 1.e20; }

    virtual void giveMaterialStiffnessMatrix(FloatMatrix & answer, MatResponseForm, MatResponseMode,
                                             GaussPoint * gp,
                                             TimeStep * atTime);

    void    giveCrackedStiffnessMatrix(FloatMatrix &answer,
                                       MatResponseMode rMode,
                                       GaussPoint *gp,
                                       TimeStep *atTime);
    /*
     * void  computeTrialStressIncrement (FloatArray& answer, GaussPoint *gp,
     * const FloatArray& strainIncrement, TimeStep* atTime);
     */

    // FloatMatrix* GiveNMatrix (GaussPoint* gp);
    FloatMatrix *GiveCrackTransformationMtrx(GaussPoint *gp, int i);
    virtual void giveEffectiveMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form,
                                                      MatResponseMode rMode,
                                                      GaussPoint *gp, TimeStep *atTime);

    void giveRealPrincipalStressVector3d(FloatArray &answer, GaussPoint *,
                                         FloatArray &, FloatMatrix &, TimeStep *);
    void giveNormalElasticStiffnessMatrix(FloatMatrix & answer,
                                          MatResponseForm, MatResponseMode,
                                          GaussPoint *, TimeStep * atTime,
                                          const FloatMatrix &);
    void         updateActiveCrackMap(GaussPoint *gp, const IntArray *activatedCracks = NULL);
    // Give3dMaterialStiffnessMatrix should return 3d material stiffness matrix
    // taking into account possible failure or fracture of material
    double giveResidualStrength() { return 0.01 * this->Ft; }


    void givePlaneStressStiffMtrx(FloatMatrix & answer, MatResponseForm form, MatResponseMode,
                                  GaussPoint * gp,
                                  TimeStep * atTime);
    void givePlaneStrainStiffMtrx(FloatMatrix & answer, MatResponseForm form, MatResponseMode,
                                  GaussPoint * gp,
                                  TimeStep * atTime);
    void give1dStressStiffMtrx(FloatMatrix & answer, MatResponseForm form, MatResponseMode,
                               GaussPoint * gp,
                               TimeStep * atTime);
    void give2dBeamLayerStiffMtrx(FloatMatrix & answer, MatResponseForm form, MatResponseMode,
                                  GaussPoint * gp,
                                  TimeStep * atTime);
    void give2dPlateLayerStiffMtrx(FloatMatrix & answer, MatResponseForm form, MatResponseMode,
                                   GaussPoint * gp,
                                   TimeStep * atTime);
    void give3dShellLayerStiffMtrx(FloatMatrix & answer, MatResponseForm form, MatResponseMode,
                                   GaussPoint * gp,
                                   TimeStep * atTime);
};

} // end namespace oofem
#endif // rcm2_h
