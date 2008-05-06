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
#define IDM_USE_MMAClosestIPTransfer
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



#include "materialmapperinterface.h"

#ifdef IDM_USE_MMAClosestIPTransfer
#include "mmaclosestiptransfer.h"
#endif

#ifdef IDM_USE_MMAShapeFunctProjection
#include "mmashapefunctprojection.h"
#endif

#ifdef IDM_USE_MMALeastSquareProjection
#include "mmaleastsquareprojection.h"
#endif


// material contant's keys for give()
#define IDM1_ITERATION_LIMIT 1.e-8

/// Type characterizing the algorithm used to compute equivalent strain measure
enum EquivStrainType { EST_Mazars, EST_Rankine, EST_ElasticEnergy };

/**
 * This class implements associated Material Status to IsotropicDamageMaterial1.
 */
class IsotropicDamageMaterial1Status : public IsotropicDamageMaterialStatus
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
class IsotropicDamageMaterial1 : public IsotropicDamageMaterial
    , public MaterialModelMapperInterface
{
    /*
     *
     * DESCRIPTION
     * This class implements a simple local isotropic damage model for concrete in tension
     *
     * A model is based on isotropic damage concept, assuming that damage evolution law
     * is postulated in explicit form, relatin damage parameter (omega) to scalar measure
     * of the largest strain level ever reached in material (kappa).
     *
     * TASK
     * - Returning standard material stiffness and flexibility marices for 3d-case.
     * according to current state determined by using data stored
     * in Gausspoint.
     * - Returning a material property (method 'give'). Only for non-standard elements.
     * - Returning real stress state vector(tensor) at gauss point for 3d - case.
     * - Storing & restoring Material Staus sored in gp matStatusDictionary.
     */

protected:

    /// max effective strain at peak
    double e0;
    /// determines the softening -> corresponds to crack opening (not strain) when tension stress vanishes
    double ef;
    /// Algorithm used to compute equiv. strain
    EquivStrainType equivStrainType;

#ifdef IDM_USE_MMAClosestIPTransfer
    /// Mapper used to map internal variables in adaptivity
    static MMAClosestIPTransfer mapper;
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
     * Computes the equivalent strain measure from given strain vector (full form).
     * @param kappa return param, comtaining the corresponding equivalent strain
     * @param strain total strain vector in full form
     * @param gp integration point
     * @param atTime timeStep
     */
    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime);
    /**
     * computes the value of damage parameter omega, based on given value of equivalent strain
     * @param omega contains result
     * @param kappa equivalent strain measure
     */
    virtual void computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp);
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
    MaterialStatus *CreateStatus(GaussPoint *gp) const { return new IsotropicDamageMaterial1Status(1, IsotropicDamageMaterial1 :: domain, gp); }

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


#endif // idm1_h

