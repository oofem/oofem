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
 *               Copyright (C) 1993 - 2016 Borek Patzak
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

#ifndef fcm_h
#define fcm_h

#include "material.h"
#include "isolinearelasticmaterial.h"
#include "structuralmaterial.h"
#include "structuralms.h"
#include "intarray.h"

///@name Input fields for FCMMaterial
//@{
#define _IFT_FCM_nAllowedCracks "ncracks"
#define _IFT_FCM_crackSpacing "crackspacing"
#define _IFT_FCM_multipleCrackShear "multiplecrackshear"
#define _IFT_FCM_ecsm "ecsm"
//@}

namespace oofem {
// crack statuses
#define pscm_NONE 0
#define pscm_JUST_INIT 1
#define pscm_SOFTENING 2
#define pscm_UNLO_RELO 3
#define pscm_CLOSED 4

#define fcm_SMALL_STRAIN 1.e-12
#define fcm_BIGNUMBER 1.e6
#define fcm_TOLERANCE 1.e-6


/**
 * This class implements associated Material Status to FCMMaterial (fixed crack material).
 */
class FCMMaterialStatus : public StructuralMaterialStatus
{
protected:
    /// crack statuses (none, just initialized, softenin, unlo-relo, closed)
    IntArray crackStatuses, tempCrackStatuses;
    /// Max. crack strain reached in the entire previous history
    FloatArray maxCrackStrains, tempMaxCrackStrains;
    /// Components of crack strain vector (normal as well as shear).
    FloatArray crackStrainVector, tempCrackStrainVector;
    /// Storing direction of cracks (crack normals) in columwise format.
    FloatMatrix crackDirs;
    /// Characteristic lengths computed from the crack orientation and element geometry
    FloatArray charLengths;

    /// transformation matrix converting stress from global to local coordinate system
    FloatMatrix transMatrix_G2Lstress;
    /// transformation matrix converting strain from global to local coordinate system
    FloatMatrix transMatrix_G2Lstrain;
    /// transformation matrix converting stress from local to global coordinate system
    FloatMatrix transMatrix_L2Gstress;
    /// transformation matrix converting strain from local to global coordinate system
    FloatMatrix transMatrix_L2Gstrain;

    /// number of maximum possible cracks (optional parameter)
    int nMaxCracks;

public:
    FCMMaterialStatus(int n, Domain *d, GaussPoint *g);
    virtual ~FCMMaterialStatus();
    /// prints the output into the output file
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    /// returns number of cracks from the previous time step (equilibrated value)
    virtual int giveNumberOfCracks() const;
    /// returns temporary number of cracks
    virtual int giveNumberOfTempCracks() const;

    /// returns vector with maximum cracking strains (max 3 components)
    const FloatArray &giveMaxCrackStrainVector() { return maxCrackStrains; }
    /// returns maximum crack strain for the i-th crack (equilibrated value)
    double giveMaxCrackStrain(int icrack) { return maxCrackStrains.at(icrack); }
    /// sets value of the maximum crack strain for the i-th crack (equilibrated value)
    void setMaxCrackStrain(int icrack, double val) { maxCrackStrains.at(icrack) = val; }

    /// returns maximum crack strain for the i-th crack (temporary value)
    double giveTempMaxCrackStrain(int icrack) { return tempMaxCrackStrains.at(icrack); }
    /// sets value of the maximum crack strain for the i-th crack (temporary value)
    void setTempMaxCrackStrain(int icrack, double val) { tempMaxCrackStrains.at(icrack) = val; }

    /// returns vector of temporary crack statuses
    const IntArray &giveTempCrackStatus() { return tempCrackStatuses; }
    /// returns temporary value of status associated with i-th crack direction
    int giveTempCrackStatus(int icrack) const { return tempCrackStatuses.at(icrack); }
    /// sets temporary value of status for of the i-th crack
    void setTempCrackStatus(int icrack, int val) { tempCrackStatuses.at(icrack) = val; }
    /// return equilibrated value of status associated with i-th crack direction
    int giveCrackStatus(int icrack) const { return crackStatuses.at(icrack); }

    /// return equilibrated crack strain vector (max 6 components)
    const FloatArray &giveCrackStrainVector() const { return crackStrainVector; }
    /// return temporary crack strain vector (max 6 components)
    const FloatArray &giveTempCrackStrainVector() { return tempCrackStrainVector; }
    /// returns i-th component of the crack strain vector (equilibrated)
    double giveCrackStrain(int icrack) const { return crackStrainVector.at(icrack); }
    /// returns i-th component of the crack strain vector (temporary)
    double giveTempCrackStrain(int icrack) const { return tempCrackStrainVector.at(icrack); }
    /// sets temporary vector of cracking strains (max 6 components)
    void setTempCrackStrainVector(FloatArray a) { tempCrackStrainVector = std :: move(a); }
    /// sets temporary value of i-th cracking strain (max 6 components)
    void setTempCrackStrain(int icrack, double val) { tempCrackStrainVector.at(icrack) = val; }
    /// sets equilibrated vector of cracking strains (max 6 components)
    void setCrackStrainVector(FloatArray a) { crackStrainVector = std :: move(a); }
    /// sets transformation matrix for stress transformation from global to local coordinate system
    void setG2LStressVectorTransformationMtrx(FloatMatrix t) { transMatrix_G2Lstress = std :: move(t); }
    /// sets transformation matrix for strain transformation from global to local coordinate system
    void setG2LStrainVectorTransformationMtrx(FloatMatrix s) { transMatrix_G2Lstrain = std :: move(s); }
    /// sets transformation matrix for stress transformation from local to global coordinate system
    void setL2GStressVectorTransformationMtrx(FloatMatrix t) { transMatrix_L2Gstress = std :: move(t); }
    /// sets transformation matrix for stress transformation from global to local coordinate system
    void setL2GStrainVectorTransformationMtrx(FloatMatrix s) { transMatrix_L2Gstrain = std :: move(s); }

    /// returns transformation matrix for stress transformation from global to local coordinate system
    const FloatMatrix &giveG2LStressVectorTransformationMtrx() { return transMatrix_G2Lstress; }
    /// sets transformation matrix for strain transformation from global to local coordinate system
    const FloatMatrix &giveG2LStrainVectorTransformationMtrx() { return transMatrix_G2Lstrain; }
    /// sets transformation matrix for stress transformation from local to global coordinate system
    const FloatMatrix &giveL2GStressVectorTransformationMtrx() { return transMatrix_L2Gstress; }
    /// sets transformation matrix for stress transformation from global to local coordinate system
    const FloatMatrix &giveL2GStrainVectorTransformationMtrx() { return transMatrix_L2Gstrain; }

    /// returns characteristic length associated with i-th crack direction
    double giveCharLength(int icrack) const {
        if ( icrack ) {
            return charLengths.at(icrack);
        } else {
            return 0.0;
        }
    }
    /// sets characteristic length for i-th crack
    void setCharLength(int icrack, double val) { charLengths.at(icrack) = val; }
    /// returns crack directions
    const FloatMatrix &giveCrackDirs() { return crackDirs; }
    /// returns crack statuses
    const IntArray &giveCrackStatus() { return crackStatuses; }
    /// sets matrix with crack directions (normal vectors)
    void setCrackDirs(FloatMatrix a) { crackDirs = std :: move(a); }
    /// returns maximum number of cracks associated with current mode
    virtual int giveMaxNumberOfCracks(GaussPoint *gp);

    // definition
    virtual const char *giveClassName() const { return "FCMMaterialStatus"; }

    /// initializes temporary status
    virtual void initTempStatus();
    /// replaces equilibrated values with temporary values
    virtual void updateYourself(TimeStep *tStep);

    /// saves current context(state) into stream
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    /// restores context(state) from stream
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
};

/**
 * This class implements a Fixed Crack Model for fracture (after initiation the crack directions cannot rotate).
 * After the onset of cracking, the cracks can still transfer tractions (i.e. cohesive crack approach) in
 * both normal and shear directions. Cracks can develop only in mutually perpendicular directions.
 * In elastic state this model is isotropic.
 * This is class is purely abstract, it can be used only in the derived classes (e.g. ConcreteFCM)
 */
class FCMMaterial : public StructuralMaterial
{
protected:
    IsotropicLinearElasticMaterial *linearElasticMaterial;

public:
    FCMMaterial(int n, Domain *d);
    virtual ~FCMMaterial();

    // identification and auxiliary functions
    virtual int hasNonLinearBehaviour() { return 1; }
    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual const char *giveClassName() const { return "FCMMaterial"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual double give(int aProperty, GaussPoint *gp);

    IsotropicLinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);


    virtual void givePlaneStressStiffMtrx(FloatMatrix &answer,
                                          MatResponseMode mmode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);

    virtual void givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                          MatResponseMode mmode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);

    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void initializeCrack(GaussPoint *gp, FloatMatrix &base, int nCrack);

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

    /// uses temporary cracking strain and characteristic length to obtain the crack opening
    virtual double computeNormalCrackOpening(GaussPoint *gp, int i);
    /// uses maximum equilibrated cracking strain and characteristic length to obtain the maximum reached crack opening
    virtual double computeMaxNormalCrackOpening(GaussPoint *gp, int i);

    /// computes total shear slip on a given crack plane (i = 1, 2, 3); the slip is computed from the temporary cracking strain
    virtual double computeShearSlipOnCrack(GaussPoint *gp, int i);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new FCMMaterialStatus(1, domain, gp); }

protected:

    /// allowed number of cracks (user-defined)
    int nAllowedCracks;

    /// Method used for evaluation of characteristic element size
    ElementCharSizeMethod ecsMethod;

    /// if true = takes shear compliance of all cracks, false = only dominant crack contribution, default value is false
    bool multipleCrackShear;

    /// comutes tensile strength
    virtual double giveTensileStrength(GaussPoint *gp) = 0;

    /// checks possible snap-back
    virtual void checkSnapBack(GaussPoint *gp, int crack) = 0;

    /// updates crack statuses
    virtual void updateCrackStatus(GaussPoint *gp);

    /// computes normal stress associated with i-th crack direction
    virtual double giveNormalCrackingStress(GaussPoint *gp, double eps_cr, int i) = 0;

    /// returns characteristic element length in given direction
    virtual double giveCharacteristicElementLength(GaussPoint *gp, const FloatArray &crackPlaneNormal);

    /// returns stiffness in the normal direction of the i-th crack
    virtual double giveCrackingModulus(MatResponseMode rMode, GaussPoint *gp, int i) = 0;

    /// returns Geff which is necessary in the global stiffness matrix
    virtual double computeEffectiveShearModulus(GaussPoint *gp, int i) = 0;

    /// shear modulus for a given shear direction (4, 5, 6)
    virtual double computeTotalD2Modulus(GaussPoint *gp, int i);

    /// shear modulus for a given crack plane (1, 2, 3)
    virtual double computeD2ModulusForCrack(GaussPoint *gp, int icrack) = 0;

    /// computes the maximum value of the shear stress; if the shear stress exceeds this value, it is cropped
    virtual double maxShearStress(GaussPoint *gp, int i) = 0;

    /// returns true for closed or no cracks associated to given shear direction (i = 4, 5, 6)
    virtual bool isIntactForShear(GaussPoint *gp, int i);

    /// returns true for closed or no crack (i = 1, 2, 3)
    virtual bool isIntact(GaussPoint *gp, int icrack);

    /// checks if the globalStress does not exceed strength in the direction of newBase for n-th crack
    virtual bool checkStrengthCriterion(FloatMatrix &newBase, const FloatArray &globalStress, GaussPoint *gp, TimeStep *tStep, int nCrack);

    /// compares trial stress with strength. Returns true if the strength is exceeded. Function oveloaded in the nonlocal approach for the fiber reinforced composites
    virtual bool isStrengthExceeded(const FloatMatrix &base, GaussPoint *gp, TimeStep *tStep, int iCrack, double trialStress);

    /// function calculating ratio used to split shear slips on two crack planes
    virtual double computeShearStiffnessRedistributionFactor(GaussPoint *gp, int ithCrackPlane, int jthCrackDirection);

    /// value of crack spacing (allows to "have" more parallel cracks in one direction if the element size exceeds user-defined or computed crack spacing).
    double crackSpacing;

    /// returns either user-provided value of crack spacing or a value computed from composition
    virtual double giveCrackSpacing(void);

    /// returns number of fictiotious parallel cracks in the direction of i-th crack
    virtual double giveNumberOfCracksInDirection(GaussPoint *gp, int iCrack);

    /// returns number of cracks for given shear direction (i = 4, 5, 6) which is treated as the maximum of the two associated normal directions
    virtual double giveNumberOfCracksForShearDirection(GaussPoint *gp, int i);


    virtual void giveMaterialStiffnessMatrix(FloatMatrix & answer, MatResponseMode,
                                             GaussPoint * gp,
                                             TimeStep * tStep);

    /// returns local stiffness matrix of the crack
    virtual void giveLocalCrackedStiffnessMatrix(FloatMatrix &answer,
                                                 MatResponseMode rMode,
                                                 GaussPoint *gp, TimeStep *tStep);

    /// returns overall Young's modulus
    virtual double computeOverallElasticStiffness(void) { return linearElasticMaterial->giveYoungsModulus(); }

    /// returns overall shear modulus
    virtual double computeOverallElasticShearModulus(void) { return linearElasticMaterial->giveShearModulus(); }
};
} // end namespace oofem
#endif // fcm_h
