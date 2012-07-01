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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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
//#define IDM_USE_MAPPEDSTRAIN

#include "material.h"
#include "linearelasticmaterial.h"
#include "isodamagemodel.h"
#include "structuralms.h"
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
    IsotropicDamageMaterial1Status(int n, Domain *d, GaussPoint *g);
    /// Destructor
    virtual ~IsotropicDamageMaterial1Status() { }

    // definition
    virtual const char *giveClassName() const { return "IsotropicDamageMaterial1Status"; }
    virtual classType giveClassID() const { return IsotropicDamageMaterialStatusClass; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual Interface *giveInterface(InterfaceType it);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
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
    double e0;
    /// Determines ductility -> corresponds to fracturing strain.
    double ef;
    /// Determines ductility -> corresponds to crack opening in the cohesive crack model.
    double wf;

    /**
     * Determines the softening -> corresponds to the initial fracture energy. For a linear law, it is the area
     * under the stress/strain curve. For an exponential law, it is the area bounded by the elastic range
     * and a tangent to the softening part of the curve at the peak stress. For a bilinear law,
     * gf corresponds to area bounded by elasticity and the first linear softening line projected to zero stress.
     */
    double gf;

    /// Determines the softening for the bilinear law -> corresponds to the total fracture energy.
    double gft;

    /// Determines the softening for the bilinear law -> corresponds to the strain at the knee point.
    double ek;

    /// Type characterizing the algorithm used to compute equivalent strain measure.
    enum EquivStrainType { EST_Unknown, EST_Mazars, EST_Rankine_Smooth, EST_Rankine_Standard, EST_ElasticEnergy, EST_ElasticEnergyPositiveStress, EST_ElasticEnergyPositiveStrain, EST_Mises };
    /// Parameter specifying the definition of equivalent strain.
    EquivStrainType equivStrainType;

    /// Parameter used in Mises definition of equivalent strain.
    double k;

    /// Remporary parameter reading type of softening law, used in other isotropic damage material models.
    int damageLaw;

    /** Type characterizing the formula for the damage law. For example, linear softening can be specified
     *   with fracturing strain or crack opening.
     */
    enum SofteningType { ST_Unknown, ST_Exponential, ST_Linear, ST_Mazars, ST_Smooth, ST_SmoothExtended, ST_Exponential_Cohesive_Crack, ST_Linear_Cohesive_Crack, ST_BiLinear_Cohesive_Crack, ST_Disable_Damage };

    /// Parameter specifying the type of softening (damage law).
    SofteningType softType;

    /// Parameters used in Mazars damage law.
    double At, Bt;
    /// Parameter used in "smooth damage law".
    double md;

    /// Parameters used if softType = 7 (extended smooth damage law)
    double e1, e2, s1, nd;
    /// Check possible snap back flag
    int checkSnapBack;

    /// Method used for evaluation of characteristic element size
    ElementCharSizeMethod ecsMethod;

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

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "IsotropicDamageMaterial1"; }
    virtual classType giveClassID() const { return IsotropicDamageMaterial1Class; }
    virtual const char *giveInputRecordName() const { return "idm1"; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    /**
     * Computes invariants I1 and J2 of the strain tensor
     * from the strain components stored in a vector.
     * @param strainVector Input strain components.
     * @param[out] I1e Output value of strain invariant I1.
     * @param[out] J2e Output value of strain invariant J2.
     */
    static void computeStrainInvariants(const FloatArray &strainVector, double &I1e, double &J2e);

    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp);
    /**
     * computes the value of damage parameter omega,
     * based on a given value of equivalent strain,
     * using iterations to achieve objectivity,
     * based on the crack band concept (effective element size used)
     * @param[out] omega Contains the resulting damage.
     * @param kappa Equivalent strain measure.
     * @param gp Integration point.
     */
    void computeDamageParamForCohesiveCrack(double &omega, double kappa, GaussPoint *gp);
    /**
     * Returns the value of damage parameter
     * corresponding to a given value
     * of the damage-driving variable kappa,
     * depending on the type of selected damage law,
     * using a simple dependence (no adjustment for element size).
     * @param kappa Equivalent strain measure.
     * @param gp Integration point.
     */
    double damageFunction(double kappa, GaussPoint *gp);
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
    double complianceFunction(double kappa, GaussPoint *gp);

    virtual Interface *giveInterface(InterfaceType it);

    virtual int MMI_map(GaussPoint *gp, Domain *oldd, TimeStep *tStep);
    virtual int MMI_update(GaussPoint *gp, TimeStep *tStep, FloatArray *estrain = NULL);
    virtual int MMI_finish(TimeStep *tStep);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
    virtual MaterialStatus *giveStatus(GaussPoint *gp) const;

    virtual double give(int aProperty, GaussPoint *gp);

protected:
    /**
     * Performs initialization, when damage first appear. The characteristic length is
     * computed from the direction of largest positive principal strain and stored
     * in corresponding status.
     * @param kappa Scalar measure of strain level.
     * @param totalStrainVector Current total strain vector.
     * @param gp Integration point.
     */
    virtual void initDamaged(double kappa, FloatArray &totalStrainVector, GaussPoint *gp);
};
} // end namespace oofem
#endif // idm1_h

