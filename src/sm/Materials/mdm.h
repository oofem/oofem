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

#ifndef mdm_h
#define mdm_h

/**
 * Select the mapping algorithm. The IDM_USE_MMAShapeFunctProjection does not work, since
 * this mapper does not preserve the max. property of damage and equivalent strain.
 */
#define MDM_USE_MMAClosestIPTransfer
//#define MDM_USE_MMAShapeFunctProjection
//#define MDM_USE_MMALeastSquareProjection


// Allows to select different mapping algorithms
#define MDM_MAPPING_DEBUG 1

#include "microplanematerial.h"
#include "structuralnonlocalmaterialext.h"
#include "matconst.h"
#include "../sm/Materials/structuralms.h"
#include "materialmapperinterface.h"
#include "mmaclosestiptransfer.h"

#ifdef MDM_MAPPING_DEBUG
 #include "mmashapefunctprojection.h"
 #include "mmaleastsquareprojection.h"

#else

 #ifdef MDM_USE_MMAShapeFunctProjection
  #include "mmashapefunctprojection.h"
 #endif
 #ifdef MDM_USE_MMALeastSquareProjection
  #include "mmaleastsquareprojection.h"
 #endif
#endif

#include "mmashapefunctprojection.h"
#include "mmaleastsquareprojection.h"

///@name Input fields for MDM
//@{
#define _IFT_MDM_Name "mdm"
#define _IFT_MDM_talpha "talpha"
#define _IFT_MDM_parmd "parmd"
#define _IFT_MDM_nonloc "nonloc"
#define _IFT_MDM_r "r"
#define _IFT_MDM_efp "efp"
#define _IFT_MDM_ep "ep"
#define _IFT_MDM_gf "gf"
#define _IFT_MDM_ft "ft"
#define _IFT_MDM_formulation "formulation"
#define _IFT_MDM_mode "mode"
#define _IFT_MDM_mapper "mapper"
#define _IFT_MDM_sourceRegionSet "sourceregset"
//@}

namespace oofem {
class Microplane;

/**
 * Material status class MDMStatus associated to MDM matarial
 */
class MDMStatus : public StructuralMaterialStatus, public StructuralNonlocalMaterialStatusExtensionInterface
{
protected:
    /// Damage values on individual microplanes.
    FloatArray Psi, PsiTemp;
    /// Macroscopic damage tensor.
    FloatMatrix DamageTensor, DamageTensorTemp, localDamageTensor;
    /// Principal damage directions.
    FloatArray tempDamageTensorEigenValues, damageTensorEigenValues;
    FloatMatrix tempDamageTensorEigenVectors, damageTensorEigenVectors;

public:
    MDMStatus(int n, int nsd, int nmplanes, Domain * d, GaussPoint * g);
    virtual ~MDMStatus();

    void setTempDamageTensorEigenVals(FloatArray src) { tempDamageTensorEigenValues = std :: move(src); }
    void setTempDamageTensorEigenVec(FloatMatrix src) { tempDamageTensorEigenVectors = std :: move(src); }
    const FloatArray &giveTempDamageTensorEigenVals() { return tempDamageTensorEigenValues; }
    const FloatArray &giveDamageTensorEigenVals() { return damageTensorEigenValues; }
    const FloatMatrix &giveTempDamageTensorEigenVec() { return tempDamageTensorEigenVectors; }
    const FloatMatrix &giveDamageTensorEigenVec() { return damageTensorEigenVectors; }

    double giveMicroplaneTempDamage(int m) { return PsiTemp.at(m); }
    double giveMicroplaneDamage(int m) { return Psi.at(m); }
    void setMicroplaneTempDamage(int m, double val) { PsiTemp.at(m) = val; }
    const FloatArray & giveMicroplaneDamageValues() { return Psi; }
    void setMicroplaneTempDamageValues(FloatArray src) { PsiTemp = std :: move(src); }

    const FloatMatrix &giveTempDamageTensor() { return DamageTensorTemp; }
    const FloatMatrix &giveDamageTensor() { return DamageTensor; }
    void setTempDamageTensor(FloatMatrix src) { DamageTensorTemp = std :: move(src); }
    void setLocalDamageTensorForAverage(FloatMatrix src) { localDamageTensor = std :: move(src); }
    const FloatMatrix &giveLocalDamageTensorForAverage() { return localDamageTensor; }
    const FloatMatrix *giveLocalDamageTensorForAveragePtr() { return & localDamageTensor; }


    // definition
    virtual const char *giveClassName() const { return "MDMStatus"; }

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    /**
     * Interface requesting service.  In the case of nonlocal constitutive models,
     * the use of multiple inheritance is assumed. Typically, the class representing nonlocal
     * constitutive model status is derived both from class representing local status and from class
     * NonlocalMaterialStatusExtension or from one of its derived classes
     * (which declare services and variables corresponding to specific analysis type).
     * @return In both cases, this function returns pointer to this object, obtained by
     * returning adress of component or using pointer conversion from receiver to base class
     * NonlocalMaterialStatusExtension. If no nonlocal extension exists, NULL pointer is returned.
     */
    virtual Interface *giveInterface(InterfaceType it);
};


/**
 * Implementation of microplane damage material (According to M.Jirasek).
 */
class MDM : public MicroplaneMaterial, public StructuralNonlocalMaterialExtensionInterface, public MaterialModelMapperInterface
{
public:
    enum MDMFormulatrionType { COMPLIANCE_DAMAGE, STIFFNESS_DAMAGE };
    enum MDMModeType { mdm_3d, mdm_2d };

protected:
    int ndc; ///< Number of damage components.
    int nsd; ///< Number of spatial dimensions.
    int type_dam;
    int type_soft;

    /// Parameter controlling the elastic limit.
    double mdm_Ep;
    /// Prescribed value of ef-ep.
    double mdm_Efp;
    double ParMd; ///< (m/E*Ep)
    double tempDillatCoeff; ///< Temperature dilatation coeff.

    /// Fracture energy (necessary to determine Ep and Efp if not given).
    double Gf;
    /// Macroscopic tensile strength (necessary to determine Ep and Efp if not given).
    double Ft;

    MDMFormulatrionType formulation;
    MDMModeType mdmMode;

    /// Reference to bulk (undamaged) material.
    StructuralMaterial *linearElasticMaterial;

    /// Flag indicating local or nonlocal mode.
    int nonlocal;
    /// Interaction radius, related to the nonlocal characteristic length of material.
    double R;

    ///cached source element set used to map internal variables (adaptivity), created on demand
    Set *sourceElemSet;

#ifdef MDM_MAPPING_DEBUG
    /// Mapper used to map internal variables in adaptivity.
    static MMAShapeFunctProjection mapperSFT;
    static MMALeastSquareProjection mapperLST;

    // determines the mapping algorithm to be used
    enum MDMMapperType { mdm_cpt=0, mdm_sft=1, mdm_lst=2 };
    MDMMapperType mapperType;
#else

 #ifdef MDM_USE_MMAClosestIPTransfer
    /// Mapper used to map internal variables in adaptivity.
    static MMAClosestIPTransfer mapper;
 #endif
 #ifdef MDM_USE_MMAShapeFunctProjection
    /// Mapper used to map internal variables in adaptivity.
    static MMAShapeFunctProjection mapper;
 #endif
 #ifdef MDM_USE_MMALeastSquareProjection
    /// Mapper used to map internal variables in adaptivity.
    static MMALeastSquareProjection mapper;
 #endif
#endif

    /// Mapper used to map stresses in adaptivity.
    static MMAClosestIPTransfer mapper2;

public:
    /**
     * Constructor. Creates Microplane Material belonging to domain d, with number n.
     * @param n Material number.
     * @param d Domain to which newly created material belongs.
     */
    MDM(int n, Domain * d) : MicroplaneMaterial(n, d), StructuralNonlocalMaterialExtensionInterface(d), MaterialModelMapperInterface()
    {
        linearElasticMaterial = NULL;
        nonlocal = 0;
        type_dam = 0;
        type_soft = 0;
        mdm_Ep = mdm_Efp = -1.0;
        sourceElemSet = NULL;
    }
    /// Destructor.
    virtual ~MDM() {
        if ( linearElasticMaterial ) {
            delete linearElasticMaterial;
        }
        if ( sourceElemSet ) {
            delete sourceElemSet;
        }
    }

    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);


    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    // identification and auxiliary functions
    virtual const char *giveInputRecordName() const { return _IFT_MDM_Name; }
    virtual const char *giveClassName() const { return "MDM"; }

    virtual void giveRealMicroplaneStressVector(FloatArray &answer, Microplane *mplane,
                                                const FloatArray &strain, TimeStep *tStep) { }

    virtual void initializeData(int numberOfMicroplanes);

    virtual void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep);

    virtual double computeWeightFunction(const FloatArray &src, const FloatArray &coord);

    virtual int hasBoundedSupport() { return 1; }
    /**
     * Determines the width (radius) of limited support of weighting function.
     */
    virtual void giveSupportRadius(double &radius) { radius = this->R; }

    virtual Interface *giveInterface(InterfaceType it);

    virtual int MMI_map(GaussPoint *gp, Domain *oldd, TimeStep *tStep);
    virtual int MMI_update(GaussPoint *gp, TimeStep *tStep, FloatArray *estrain = NULL);
    virtual int MMI_finish(TimeStep *tStep);

    virtual int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip);
    virtual int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip);
    virtual int estimatePackSize(DataStream &buff, GaussPoint *ip);
    virtual double predictRelativeComputationalCost(GaussPoint *gp);
    virtual double predictRelativeRedistributionCost(GaussPoint *gp) { return 1.0; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
    /// Returns reference to undamaged (bulk) material.
    StructuralMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual MaterialStatus *CreateMicroplaneStatus(GaussPoint *gp)
    { return NULL; }
    void computeDamageTensor(FloatMatrix &tempDamageTensor, const FloatArray &totalStrain,
                             GaussPoint *gp, TimeStep *tStep);
    void computeLocalDamageTensor(FloatMatrix &tempDamageTensor, const FloatArray &totalStrain,
                                  GaussPoint *gp, TimeStep *tStep);
    double computeDamageOnPlane(GaussPoint *gp, Microplane *mplane, const FloatArray &strain);
    void computePDC(FloatMatrix &tempDamageTensor, FloatArray &tempDamageTensorEigenVals,
                    FloatMatrix &tempDamageTensorEigenVec);
    void transformStrainToPDC(FloatArray &answer, FloatArray &strain,
                              FloatMatrix &t, GaussPoint *gp);
    void applyDamageTranformation(FloatArray &strainPDC, const FloatArray &tempDamageTensorEigenVals);
    void transformStressFromPDC(FloatArray &answer, const FloatArray &stressPDC, const FloatMatrix &t, GaussPoint *gp);
    void computeEffectiveStress(FloatArray &stressPDC, const FloatArray &strainPDC,
                                GaussPoint *gp, TimeStep *tStep);
    void giveMaterialStiffnessMatrix(FloatMatrix &answer,
                                     MatResponseMode mode, GaussPoint *gp,
                                     TimeStep *tStep);
    void applyDamageToStiffness(FloatMatrix &d, GaussPoint *gp);
    void transformStiffnessfromPDC(FloatMatrix &de, const FloatMatrix &t);

    virtual void givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode mmode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);

    virtual void givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                          MatResponseMode mmode, GaussPoint *gp,
                                          TimeStep *tStep);

    void rotateTensor4(FloatMatrix &Dlocal, const FloatMatrix &t);
    void formTransformationMatrix(FloatMatrix &answer, const FloatMatrix &t, int n);

    //double giveParameterEfp(const FloatArray &reducedStrain, GaussPoint *gp);
    void giveRawMDMParameters(double &Efp, double &Ep, const FloatArray &reducedStrain, GaussPoint *gp);
};
} // end namespace oofem
#endif // mdm_h
