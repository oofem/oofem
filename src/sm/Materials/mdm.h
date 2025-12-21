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
#include "sm/Materials/isolinearelasticmaterial.h"
#include "sm/Materials/structuralms.h"
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
    MDMStatus(GaussPoint * g, int nsd, int nmplanes);

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

    const char *giveClassName() const override { return "MDMStatus"; }

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

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
    Interface *giveInterface(InterfaceType it) override;
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
    int ndc = 0; ///< Number of damage components.
    int nsd = 0; ///< Number of spatial dimensions.
    int type_dam = 0;
    int type_soft = 0;

    /// Parameter controlling the elastic limit.
    mutable double mdm_Ep = 0.;
    /// Prescribed value of ef-ep.
    mutable double mdm_Efp = 0.; /// THREAD UNSAFE. These are modified when evaluating the material. Why? This must be avoided at all cost.
    double ParMd = 0.; ///< (m/E*Ep)
    double tempDillatCoeff = 0.; ///< Temperature dilatation coeff.

    /// Fracture energy (necessary to determine Ep and Efp if not given).
    double Gf = 0.;
    /// Macroscopic tensile strength (necessary to determine Ep and Efp if not given).
    double Ft = 0.;

    MDMFormulatrionType formulation = COMPLIANCE_DAMAGE;
    MDMModeType mdmMode = mdm_3d;

    /// Reference to bulk (undamaged) material.
    IsotropicLinearElasticMaterial linearElasticMaterial;

    /// Flag indicating local or nonlocal mode.
    int nonlocal = 0;
    /// Interaction radius, related to the nonlocal characteristic length of material.
    double R = 0.;

    ///cached source element set used to map internal variables (adaptivity), created on demand
    std::unique_ptr<Set> sourceElemSet;

#ifdef MDM_MAPPING_DEBUG
    /// Mapper used to map internal variables in adaptivity.
    static MMAShapeFunctProjection mapperSFT;
    static MMALeastSquareProjection mapperLST;

    // determines the mapping algorithm to be used
    enum MDMMapperType { mdm_cpt=0, mdm_sft=1, mdm_lst=2 };
    MDMMapperType mapperType = mdm_cpt;
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
    MDM(int n, Domain * d) : MicroplaneMaterial(n, d), 
        StructuralNonlocalMaterialExtensionInterface(d),
        MaterialModelMapperInterface(),
        linearElasticMaterial(n, d)
    {
        mdm_Ep = mdm_Efp = -1.0;
    }

    bool hasMaterialModeCapability(MaterialMode mode) const override;

    void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                              const FloatArray &reducedStrain, TimeStep *tStep) const override;

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<MDM*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<4> giveRealStressVector_PlaneStrain(const FloatArrayF<4> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<MDM*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<3> giveRealStressVector_PlaneStress(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<MDM*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }

    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    FloatArrayF<6> giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const override;

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    const char *giveInputRecordName() const override { return _IFT_MDM_Name; }
    const char *giveClassName() const override { return "MDM"; }

    void initializeData(int numberOfMicroplanes) override;

    void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep) const override;

  double computeWeightFunction(const double cl, const FloatArray &src, const FloatArray &coord) const override;

    int hasBoundedSupport() const override { return 1; }
    /**
     * Determines the width (radius) of limited support of weighting function.
     */
    virtual void giveSupportRadius(double &radius) { radius = this->R; }

    Interface *giveInterface(InterfaceType it) override;

    int MMI_map(GaussPoint *gp, Domain *oldd, TimeStep *tStep) override;
    int MMI_update(GaussPoint *gp, TimeStep *tStep, FloatArray *estrain = nullptr) override;
    int MMI_finish(TimeStep *tStep) override;

    int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override;
    int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override;
    int estimatePackSize(DataStream &buff, GaussPoint *ip) override;
    double predictRelativeComputationalCost(GaussPoint *gp) override;
    double predictRelativeRedistributionCost(GaussPoint *gp) override { return 1.0; }

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

protected:
    void computeDamageTensor(FloatMatrix &tempDamageTensor, const FloatArray &totalStrain,
                             GaussPoint *gp, TimeStep *tStep) const;
    void computeLocalDamageTensor(FloatMatrix &tempDamageTensor, const FloatArray &totalStrain,
                                  GaussPoint *gp, TimeStep *tStep) const;
    double computeDamageOnPlane(GaussPoint *gp, int mnumber, const FloatArray &strain) const;
    void computePDC(FloatMatrix &tempDamageTensor, FloatArray &tempDamageTensorEigenVals,
                    FloatMatrix &tempDamageTensorEigenVec) const;
    void transformStrainToPDC(FloatArray &answer, FloatArray &strain,
                              FloatMatrix &t, GaussPoint *gp) const;
    void applyDamageTranformation(FloatArray &strainPDC, const FloatArray &tempDamageTensorEigenVals) const;
    void transformStressFromPDC(FloatArray &answer, const FloatArray &stressPDC, const FloatMatrix &t, GaussPoint *gp) const;
    void computeEffectiveStress(FloatArray &stressPDC, const FloatArray &strainPDC,
                                GaussPoint *gp, TimeStep *tStep) const;
    void giveMaterialStiffnessMatrix(FloatMatrix &answer,
                                     MatResponseMode mode, GaussPoint *gp,
                                     TimeStep *tStep) const;
    void applyDamageToStiffness(FloatMatrix &d, GaussPoint *gp) const;
    void transformStiffnessfromPDC(FloatMatrix &de, const FloatMatrix &t) const;

    FloatMatrixF<3,3> givePlaneStressStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<4,4> givePlaneStrainStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;

    void rotateTensor4(FloatMatrix &Dlocal, const FloatMatrix &t) const;
    void formTransformationMatrix(FloatMatrix &answer, const FloatMatrix &t, int n) const;

    //double giveParameterEfp(const FloatArray &reducedStrain, GaussPoint *gp);
    void giveRawMDMParameters(double &Efp, double &Ep, const FloatArray &reducedStrain, GaussPoint *gp) const;
};
} // end namespace oofem
#endif // mdm_h
