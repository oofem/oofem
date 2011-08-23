/* $Header: /home/cvs/bp/oofem/sm/src/idm1.h,v 1.9 2003/04/06 14:08:30 bp Exp $ */
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


//   *****************************************************************************
//   *** CLASS SIMPLE ISOTROPIC DAMAGE MODEL FOR CONCRETE IN TENSION  ************
//   *****************************************************************************

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

/**
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
 */
class IsotropicDamageMaterial1Status : public IsotropicDamageMaterialStatus, public RandomMaterialStatusExtensionInterface
{
protected:

    /**
     * characteristic element length, computed when damage initialized from direction of max. positive
     * principal strain. Fixed during further loading
     */
    double le;


public:

    /// Constructor
    IsotropicDamageMaterial1Status(int n, Domain *d, GaussPoint *g);
    /// Destructor
    ~IsotropicDamageMaterial1Status() { }

    // void   printOutputAt (FILE *file, TimeStep* tStep) ;

    /// Returns characteristic length stored in receiver
    double giveLe()  { return le; }
    /// Sets characteristic length to given value
    void   setLe(double ls) { le = ls; }

    // definition
    const char *giveClassName() const { return "IsotropicDamageMaterial1Status"; }
    classType             giveClassID() const { return IsotropicDamageMaterialStatusClass; }

    /**
     * Initializes the temporary internal variables, describing the current state according to
     * previously reached equilibrium internal variables.
     */
    virtual void initTempStatus();
    /**
     * Update equilibrium history variables according to temp-variables.
     * Invoked, after new equilibrium state has been reached.
     */
    virtual void updateYourself(TimeStep *); // update after new equilibrium state reached

    /** Interface requesting service */
    virtual Interface *giveInterface(InterfaceType);

    // saves current context(state) into stream
    /**
     * Stores context of receiver into given stream.
     * Le attribute is stored. Corresponding parent method invoked.
     * @param stream stream where to write data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     */
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores context of receiver from given stream.
     * Le attribute is restored. Corresponding parent method invoked.
     * @param stream stream where to read data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     */
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};

/**
 * This class implements a simple local isotropic damage model for concrete in tension.
 * A model is based on isotropic damage concept, assuming that damage evolution law
 * is postulated in explicit form, relatin damage parameter (omega) to scalar measure
 * of the largest strain level ever reached in material (kappa).
 */
class IsotropicDamageMaterial1 : public IsotropicDamageMaterial, public RandomMaterialExtensionInterface
    , public MaterialModelMapperInterface
{
    /*
     *
     * DESCRIPTION
     * This class implements a simple local isotropic damage model for concrete in tension
     *
     * A model is based on isotropic damage concept, assuming that damage evolution law
     * is postulated in explicit form, relating damage parameter (omega) to scalar measure
     * of the largest strain level ever reached in material (kappa).
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

    /// equivalent strain at stress peak (or a similar parameter)
    double e0;
    /// determines ductility -> corresponds to fracturing strain
    double ef;
    /// determines ductility -> corresponds to crack opening in the cohesive crack model
    double wf;
    
    /**Determines the softening -> corresponds to the initial fracture energy. For a linear law, it is the area 
    *  under the stress/strain curve. For an exponential law, it is the area bounded by the elastic range 
    *  and a tangent to the softening part of the curve at the peak stress. For a bilinear law, 
    *  Gf corresponds to area bounded by elasticity and the first linear softening line projected to zero stress
    */
    double gf;
    
    /// Determines the softening for the bilinear law -> corresponds to the total fracture energy
    double gft;
    
    /// Determines the softening for the bilinear law -> corresponds to the strain at the knee point
    double ek;
    
    
    /// type characterizing the algorithm used to compute equivalent strain measure
    enum EquivStrainType { EST_Unknown, EST_Mazars, EST_Rankine, EST_ElasticEnergy, EST_Mises };
    /// parameter specifying the definition of equivalent strain
    EquivStrainType equivStrainType;

    /// parameter used in Mises definition of equivalent strain
    double k;

    /// temporary parameter reading type of softening law, used in other isotropic damage material models
    int damageLaw;
    
    /// type characterizing the formula for the damage law
    enum SofteningType { ST_Unknown, ST_Exponential, ST_Linear, ST_Mazars, ST_Smooth, ST_SmoothExtended, ST_Exponential_Cohesive_Crack, ST_Linear_Cohesive_Crack, ST_BiLinear_Cohesive_Crack };
    /// parameter specifying the type of softening (damage law)
    SofteningType softType;

    /// parameters used in Mazars damage law
    double At, Bt;
    /// parameter used in "smooth damage law"
    double md;

#ifdef IDM_USE_MMAClosestIPTransfer
    /// Mapper used to map internal variables in adaptivity
    static MMAClosestIPTransfer mapper;
#endif
#ifdef IDM_USE_MMAContainingElementProjection
    /// Mapper used to map internal variables in adaptivity
    static MMAContainingElementProjection mapper;
#endif
#ifdef IDM_USE_MMAShapeFunctProjection
    /// Mapper used to map internal variables in adaptivity
    static MMAShapeFunctProjection mapper;
#endif
#ifdef IDM_USE_MMALeastSquareProjection
    /// Mapper used to map internal variables in adaptivity
    static MMALeastSquareProjection mapper;
#endif

public:

    /// Constructor
    IsotropicDamageMaterial1(int n, Domain *d);
    /// Destructor
    ~IsotropicDamageMaterial1();

    // identification and auxiliary functions
    const char *giveClassName() const { return "IsotropicDamageMaterial1"; }
    classType giveClassID()         const { return IsotropicDamageMaterial1Class; }
    /// Returns input record name of the receiver.
    const char *giveInputRecordName() const { return "idm1"; }

    /**
     * Initializes receiver acording to object description stored in input record..
     * The density of material is read into property dictionary (keyword 'd')
     */
    IRResultType initializeFrom(InputRecord *ir);
    /** Setups the input record string of receiver
     * @param str string to be filled by input record
     * @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    /**
     * Computes invariants I1 and J2 of the strain tensor
     * from the strain components stored in a vector.
     * @param strainVector input strain components
     * @param I1e output value of strain invariant I1
     * @param J2e output value of strain invariant J2
     */
    void computeStrainInvariants(FloatArray *strainVector, double *I1e, double *J2e);
    /**
     * Computes the equivalent strain measure from given strain vector (full form).
     * @param kappa return param, comtaining the corresponding equivalent strain
     * @param strain total strain vector in full form
     * @param gp integration point
     * @param atTime timeStep
     */
    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime);
    /**
     * computes the value of damage parameter omega,
     * based on a given value of equivalent strain
     * @param omega contains the resulting damage
     * @param kappa equivalent strain measure
     */
    virtual void computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp);
    /**
     * computes the value of damage parameter omega,
     * based on a given value of equivalent strain,
     * using iterations to achieve objectivity,
     * based on the crack band concept (effective element size used)
     * @param omega contains the resulting damage
     * @param kappa equivalent strain measure
     */
    void computeDamageParamForCohesiveCrack(double &omega, double kappa, GaussPoint *gp);
    /**
     * Returns the value of damage parameter
     * corresponding to a given value
     * of the damage-driving variable kappa,
     * depending on the type of selected damage law,
     * using a simple dependence (no adjustment for element size).
     * @param kappa equivalent strain measure
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
     * @param kappa equivalent strain measure
     */
    double complianceFunction(double kappa, GaussPoint *gp);

    /** Interface requesting service */
    virtual Interface *giveInterface(InterfaceType);

    /**
     * @name The interface required by MaterialModelMapperInterface
     */
    //@{
    /** Maps the required internal state variables from
     * old mesh oldd to given ip. The result is stored in gp status.
     * @param gp Integration point belonging to new domain which values will be mapped
     * @param oldd old mesh reference
     * @param tStep time step
     * @return nonzero if o.k.
     */
    virtual int MMI_map(GaussPoint *gp, Domain *oldd, TimeStep *tStep);
    /** Updates the required internal state variables from previously mapped values.
     * The result is stored in gp status. This map and update splitting is necessary,
     * for example for nonlocal models tahe local quantity to be averaged must be mapped in all ips
     * and then update can happen, because it may depend on nonlocal variable, which is computed
     * from local values.
     * @param gp Integration point belonging to new domain which values will be mapped
     * @param tStep time step
     * @return nonzero if o.k.
     */
    virtual int MMI_update(GaussPoint *gp, TimeStep *tStep, FloatArray *estrain = NULL);
    /**
     * Finishes the mapping for given time step. Used to perform cleanup.
     * Typically some mappers reguire to compute some global mesh data related to
     * current step, which are valid for example to all IPs - so they are computed only once for
     * all IPs, stored and they need to be dealocated. These mappers are typically class variables,
     * but their finish is invoked by all members.
     * @return nonzero if ok.
     */
    virtual int MMI_finish(TimeStep *tStep);

    //@}

    /// Creates corresponding status
    MaterialStatus *CreateStatus(GaussPoint *gp) const;
    /**
     * Returns material status of receiver in given integration point.
     * If status does not exist yet, it is created using CreateStatus  member function.
     * @param gp Returns reference to material status belonging to integration
     * point gp.
     * @return material status associated with given integration point.
     */
    MaterialStatus *giveStatus(GaussPoint *gp) const;
    /**
     * Returns the value of material property 'aProperty'. Property must be identified
     * by unique int id. Intgeration point also passed to allow for materials with spatially
     * varying properties
     * @param aProperty id of peroperty requested
     * @param gp integration point,
     * @return property value
     */
    double   give(int aProperty, GaussPoint *gp);

protected:
    /**
     * Perfoms initialization, when damage first appear. The Le characteristic length is
     * computed from the direction of largest positive principal strain and stored
     * in corresponding status.
     * @param kappa scalar measure of strain level
     * @param totalStrainVector current total strain vector
     * @param gp integration point
     */
    void initDamaged(double kappa, FloatArray &totalStrainVector, GaussPoint *gp);
};
} // end namespace oofem
#endif // idm1_h

