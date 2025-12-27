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

#ifndef idm1_h
#define idm1_h

/**
 * Select the mapping algorithm. The IDM_USE_MMAShapeFunctProjection does not work, since
 * this mapper does not preserve the max. property of damage and equivalent strain.
 */
//#define IDM_USE_MMAClosestIPTransfer
#define IDM_USE_MMAContainingElementProjection
//#define IDM_USE_MMAShapeFunctProjection
//#define IDM_USE_MMALeastSquareProjection

/*
 * Selects the use of mapped strain or projected strain from element.
 */
#define IDM_USE_MAPPEDSTRAIN

#include "material.h"
#include "sm/Materials/linearelasticmaterial.h"
#include "sm/Materials/isodamagemodel.h"
#include "sm/Materials/structuralms.h"
#include "randommaterialext.h"
#include "materialmapperinterface.h"

#ifdef IDM_USE_MMAClosestIPTransfer
 #include "mmaclosestiptransfer.h"
#endif

#ifdef IDM_USE_MMAContainingElementProjection
 #include "mmacontainingelementprojection.h"
#endif

#ifdef IDM_USE_MMAShapeFunctProjection
 #include "mmashapefunctprojection.h"
#endif

#ifdef IDM_USE_MMALeastSquareProjection
 #include "mmaleastsquareprojection.h"
#endif

///@name Input fields for IsotropicDamageMaterial1
//@{
#define _IFT_IsotropicDamageMaterial1_Name "idm1"
#define _IFT_IsotropicDamageMaterial1_e0 "e0"
#define _IFT_IsotropicDamageMaterial1_ef "ef"
#define _IFT_IsotropicDamageMaterial1_wf "wf"
#define _IFT_IsotropicDamageMaterial1_equivstraintype "equivstraintype"
#define _IFT_IsotropicDamageMaterial1_damageLaw "damlaw"
#define _IFT_IsotropicDamageMaterial1_k "k"
#define _IFT_IsotropicDamageMaterial1_md "md"
#define _IFT_IsotropicDamageMaterial1_ecsm "ecsm"
#define _IFT_IsotropicDamageMaterial1_At "at"
#define _IFT_IsotropicDamageMaterial1_Bt "bt"
#define _IFT_IsotropicDamageMaterial1_ft "ft"
#define _IFT_IsotropicDamageMaterial1_wkwf "wkwf"
#define _IFT_IsotropicDamageMaterial1_e1ef "e1ef"
#define _IFT_IsotropicDamageMaterial1_skft "skft"
#define _IFT_IsotropicDamageMaterial1_s1 "s1"
#define _IFT_IsotropicDamageMaterial1_sk "sk"
#define _IFT_IsotropicDamageMaterial1_wk "wk"
#define _IFT_IsotropicDamageMaterial1_e1 "e1"
#define _IFT_IsotropicDamageMaterial1_ek "ek"
#define _IFT_IsotropicDamageMaterial1_gf "gf"
#define _IFT_IsotropicDamageMaterial1_gft "gft"
#define _IFT_IsotropicDamageMaterial1_ep "ep"
#define _IFT_IsotropicDamageMaterial1_e2 "e2"
#define _IFT_IsotropicDamageMaterial1_nd "nd"
#define _IFT_IsotropicDamageMaterial1_checkSnapBack "checksnapback"
#define _IFT_IsotropicDamageMaterial1_n "griff_n"
#define _IFT_IsotropicDamageMaterial1_c1 "c1"
#define _IFT_IsotropicDamageMaterial1_c2 "c2"
#define _IFT_IsotropicDamageMaterial1_alphaps "alphaps"
#define _IFT_IsotropicDamageMaterial1_h "h"
#define _IFT_IsotropicDamageMaterial1_w_k "w_k"
#define _IFT_IsotropicDamageMaterial1_w_r "w_r"
#define _IFT_IsotropicDamageMaterial1_w_f "w_f"
#define _IFT_IsotropicDamageMaterial1_f_k "f_k"
#define _IFT_IsotropicDamageMaterial1_f_r "f_r"
//@}

namespace oofem {
#define IDM1_ITERATION_LIMIT 1.e-9

/**
 * This class implements associated Material Status to IsotropicDamageMaterial1.
 * Stores the characteristic length of the element.
 */
class IsotropicDamageMaterial1Status : public IsotropicDamageMaterialStatus, public RandomMaterialStatusExtensionInterface
{
public:
    /// Constructor
    IsotropicDamageMaterial1Status(GaussPoint *g);

    const char *giveClassName() const override { return "IsotropicDamageMaterial1Status"; }

    Interface *giveInterface(InterfaceType it) override;
};

/**
 * This class implements a simple local isotropic damage model for concrete in tension.
 * A model is based on isotropic damage concept, assuming that damage evolution law
 * is postulated in explicit form, relation damage parameter (omega) to scalar measure
 * of the largest strain level ever reached in material (kappa).
 */
class IsotropicDamageMaterial1 : public IsotropicDamageMaterial,
    public RandomMaterialExtensionInterface,
    public MaterialModelMapperInterface
{
protected:
    /// Equivalent strain at stress peak (or a similar parameter).
    double e0 = 0.;
    /// Determines ductility -> corresponds to fracturing strain.
    double ef = 0.;
    /// Determines ductility -> corresponds to crack opening in the cohesive crack model.
    double wf = 0.;

    /**
     * Determines the softening -> corresponds to the initial fracture energy. For a linear law, it is the area
     * under the stress/strain curve. For an exponential law, it is the area bounded by the elastic range
     * and a tangent to the softening part of the curve at the peak stress. For a bilinear law,
     * gf corresponds to area bounded by elasticity and the first linear softening line projected to zero stress.
     */
    double gf = 0.;

    /// Determines the softening for the bilinear law -> corresponds to the total fracture energy.
    double gft = 0.;

    /// Determines the softening for the bilinear law -> corresponds to the strain at the knee point.
    double ek = 0.;

    /// Determines the softening for the bilinear law -> corresponds to the crack opening at the knee point.
    double wk = 0.;

    /// Determines the softening for the bilinear law -> corresponds to the stress at the knee point.
    double sk = 0.;

    /// Parameters used in Hordijk's softening law
    double c1 = 3., c2 = 6.93;  // default value of Hordijk parameter

    /// Parameters used in Trilinear_Cohesive_Crack softening law
    double w_k = 0., w_r = 0., w_f = 0., f_k = 0., f_r = 0.;

    /** Type characterizing the algorithm used to compute equivalent strain measure.
     *  Note that the assigned numbers to enum values have to correspond to values
     *  used in initializeFrom to resolve EquivStrainType. If not, the consistency
     *  between initializeFrom and giveInputRecord methods is lost.
     */
    enum EquivStrainType {
        EST_Mazars=0,
        EST_Rankine_Smooth=1,
        EST_ElasticEnergy=2,
        EST_Mises=3,
        EST_Rankine_Standard=4,
        EST_ElasticEnergyPositiveStress=5,
        EST_ElasticEnergyPositiveStrain=6,
        EST_Griffith=7,
        EST_Unknown = 100
    };
    /// Parameter specifying the definition of equivalent strain.
    EquivStrainType equivStrainType = EST_Unknown;

    /// Parameter used in Mises definition of equivalent strain.
    double k = 0.;

    /// Parameter used in Griffith's criterion
    double griff_n = 8.;

    /// Temporary parameter reading type of softening law, used in other isotropic damage material models.
    int damageLaw = 0;

    /** Type characterizing the formula for the damage law. For example, linear softening can be specified
     *   with fracturing strain or crack opening.
     */
    enum SofteningType { ST_Unknown, ST_Exponential, ST_Linear, ST_Mazars, ST_Smooth, ST_SmoothExtended, ST_Exponential_Cohesive_Crack, ST_Linear_Cohesive_Crack, ST_BiLinear_Cohesive_Crack, ST_Disable_Damage, ST_PowerExponential, ST_DoubleExponential, ST_Hordijk_Cohesive_Crack, ST_ModPowerExponential, ST_Trilinear_Cohesive_Crack };

    /// Parameter specifying the type of softening (damage law).
    SofteningType softType = ST_Unknown;

    /// Parameters used in Mazars damage law.
    double At = 0., Bt = 0.;
    /// Parameter used in "smooth damage law".
    double md = 1.;

    /// Parameters used if softType = 7 (extended smooth damage law)
    double e1 = 0., e2 = 0., s1 = 0., nd = 0.;
    /// Check possible snap back flag
    int checkSnapBack = 0;

    /// auxiliary input variablesfor softType == ST_SmoothExtended
    double ep = 0., ft = 0.;

    /// Parameters used by the model with permanent strain
    double ps_alpha = 0., ps_H = 0.;

    /// Method used for evaluation of characteristic element size
    ElementCharSizeMethod ecsMethod = ECSM_Unknown;

    /// Cached source element set used to map internal variables (adaptivity), created on demand
    Set *sourceElemSet = nullptr;

#ifdef IDM_USE_MMAClosestIPTransfer
    /// Mapper used to map internal variables in adaptivity.
    static MMAClosestIPTransfer mapper;
#endif
#ifdef IDM_USE_MMAContainingElementProjection
    /// Mapper used to map internal variables in adaptivity.
    static MMAContainingElementProjection mapper;
#endif
#ifdef IDM_USE_MMAShapeFunctProjection
    /// Mapper used to map internal variables in adaptivity.
    static MMAShapeFunctProjection mapper;
#endif
#ifdef IDM_USE_MMALeastSquareProjection
    /// Mapper used to map internal variables in adaptivity.
    static MMALeastSquareProjection mapper;
#endif

public:
    /// Constructor
    IsotropicDamageMaterial1(int n, Domain *d);
    /// Destructor
    virtual ~IsotropicDamageMaterial1();

    const char *giveClassName() const override { return "IsotropicDamageMaterial1"; }
    const char *giveInputRecordName() const override { return _IFT_IsotropicDamageMaterial1_Name; }
    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    /**
     * Computes invariants I1 and J2 of the strain tensor
     * from the strain components stored in a vector.
     * @param strainVector Input strain components.
     * @param[out] I1e Output value of strain invariant I1.
     * @param[out] J2e Output value of strain invariant J2.
     */
    static void computeStrainInvariants(const FloatArray &strainVector, double &I1e, double &J2e);

    bool isCrackBandApproachUsed() const { return ( this->softType == ST_Exponential_Cohesive_Crack || this->softType == ST_Linear_Cohesive_Crack || this->softType == ST_BiLinear_Cohesive_Crack || this->softType == ST_Trilinear_Cohesive_Crack || this->gf != 0. ); }
    double computeEquivalentStrain(const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) const override;

    void computeEta(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) const override;
    double computeDamageParam(double kappa, const FloatArray &strain, GaussPoint *gp) const override;
    /**
     * computes the value of damage parameter omega,
     * based on a given value of equivalent strain,
     * using iterations to achieve objectivity,
     * based on the crack band concept (effective element size used)
     * @param[out] omega Contains the resulting damage.
     * @param kappa Equivalent strain measure.
     * @param gp Integration point.
     */
    double computeDamageParamForCohesiveCrack(double kappa, GaussPoint *gp) const;
    /**
     * Returns the value of damage parameter
     * corresponding to a given value
     * of the damage-driving variable kappa,
     * depending on the type of selected damage law,
     * using a simple dependence (no adjustment for element size).
     * @param kappa Equivalent strain measure.
     * @param gp Integration point.
     */
    double damageFunction(double kappa, GaussPoint *gp) const;
    /**
     * Returns the value of compliance parameter
     * corresponding to a given value
     * of the damage-driving variable kappa,
     * depending on the type of selected damage law,
     * using a simple dependence (no adjustment for element size).
     * The compliance parameter gamma is defined as
     * gamma = omega/(1-omega)
     * where omega is the damage.
     * @param kappa Equivalent strain measure.
     * @param gp Integration point.
     */
    /**
     * Returns the value of derivative of damage function
     * wrt damage-driving variable kappa corresponding
     * to a given value of the  kappa, depending on
     * the type of selected damage law.
     * @param kappa Equivalent strain measure.
     * @param gp Integration point.
     */
    double damageFunctionPrime(double kappa, GaussPoint *gp) const override;
    /**
     * Returns the value of compliance parameter
     * corresponding to a given value
     * of the damage-driving variable kappa,
     * depending on the type of selected damage law,
     * using a simple dependence (no adjustment for element size).
     * The compliance parameter gamma is defined as
     * gamma = omega/(1-omega)
     * where omega is the damage.
     * @param kappa Equivalent strain measure.
     * @param gp Integration point.
     */
    double complianceFunction(double kappa, GaussPoint *gp) const;

    double evaluatePermanentStrain(double kappa, double omega) const override;

    Interface *giveInterface(InterfaceType it) override;

    int MMI_map(GaussPoint *gp, Domain *oldd, TimeStep *tStep) override;
    int MMI_update(GaussPoint *gp, TimeStep *tStep, FloatArray *estrain = nullptr) override;
    int MMI_finish(TimeStep *tStep) override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;
    MaterialStatus *giveStatus(GaussPoint *gp) const override;

    double give(int aProperty, GaussPoint *gp) const override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    void restoreContext(DataStream &stream, ContextMode mode) override;
    void saveContext(DataStream &stream, ContextMode mode) override;
protected:
    /**
     * Performs initialization, when damage first appear. The characteristic length is
     * computed from the direction of largest positive principal strain and stored
     * in corresponding status.
     * @param kappa Scalar measure of strain level.
     * @param totalStrainVector Current total strain vector.
     * @param gp Integration point.
     */
    void initDamaged(double kappa, FloatArray &totalStrainVector, GaussPoint *gp) const override;
};
} // end namespace oofem
#endif // idm1_h
