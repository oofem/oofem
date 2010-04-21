/* $Header: /home/cvs/bp/oofem/sm/src/idmnl1.h,v 1.9.4.1 2004/04/05 15:19:47 bp Exp $ */
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

//   ******************************************************************************
//   *** CLASS NONLOCAL ISOTROPIC DAMAGE MODEL FOR CONCRETE IN TENSION ************
//   ******************************************************************************

#ifndef idmnl1_h
#define idmnl1_h

#include "idm1.h"
#include "structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"


#include "sparsemtrx.h"
#include "dynalist.h"

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#include "conTable.h"
#endif

namespace oofem {

/* Flag indicating a special averaging on elements:
 * 0=no special averaging
 * 1=boundary layer method averages the strains over the whole finite element, the influence radius is close to zero
 */

class GaussPoint;

/**
 * This class implements associated Material Status to IDNLMaterial.
 */
class IDNLMaterialStatus : public IsotropicDamageMaterial1Status, public StructuralNonlocalMaterialStatusExtensionInterface
{
protected:

    /// Equivalent strain for avaraging
    double localEquivalentStrainForAverage;

    /*
     * // Variables used to track loading/reloading
     * public:
     * enum LastStateType {LST_elastic, LST_loading, LST_unloading};
     * LastStateType lst;
     */

public:

    /// constructor
    IDNLMaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    ~IDNLMaterialStatus();

    /// Prints the receiver state to given stream
    void   printOutputAt(FILE *file, TimeStep *tStep);

    /// Returns the local  equivalent strain to be averaged
    double giveLocalEquivalentStrainForAverage()     { return localEquivalentStrainForAverage; }
    /// Sets the localEquivalentStrainForAverage to given value
    void   setLocalEquivalentStrainForAverage(double ls) { localEquivalentStrainForAverage = ls; }

    // definition
    const char *giveClassName() const { return "IDNLMaterialStatus"; }
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
     * Corresponding parent method invoked.
     * @param stream stream where to write data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     */
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores context of receiver from given stream.
     * Corresponding parent method invoked.
     * @param stream stream where to read data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     */
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Interface requesting service.
     * In the case of nonlocal constitutive models,
     * the use of multiple inheritance is assumed. Typically, the class representing nonlocal
     * constitutive model status is derived both from class representing local status and from class
     * NonlocalMaterialStatusExtensionInterface or from one of its derived classes
     * (which declare services and variables corresponding to specific analysis type).
     * @return In both cases, this function returns pointer to this object, obtained by
     * returning adress of component or using pointer conversion from receiver to base class
     * NonlocalMaterialStatusExtensionInterface.
     */
    virtual Interface *giveInterface(InterfaceType);
};


/**
 * This class implements a Nonlocal Isotropic Damage Model for Concrete in Tension
 * Model based on nonlocal averaging of equivalent strain.
 */
class IDNLMaterial : public IsotropicDamageMaterial1, public StructuralNonlocalMaterialExtensionInterface,
    public NonlocalMaterialStiffnessInterface
{
protected:
    /// Interaction radius, related to the nonlocal characteristic length of material.
    double R;

public:

    /// Constructor
    IDNLMaterial(int n, Domain *d);
    /// Destructor
    ~IDNLMaterial();

    // identification and auxiliary functions
    const char *giveClassName() const { return "IDNLMaterial"; }
    classType giveClassID()         const { return IsotropicDamageMaterial1Class; }
    /// Returns input record name of the receiver.
    const char *giveInputRecordName() const { return "idmnl1"; }

    /// Initializes the receiver from given record
    IRResultType initializeFrom(InputRecord *ir);
    /** Setups the input record string of receiver
     * @param str string to be filled by input record
     * @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    /** Interface requesting service */
    virtual Interface *giveInterface(InterfaceType);


    /**
     * Computes the equivalent nonlocal strain measure from given strain vector (full form).
     * @param kappa return param, comtaining the corresponding equivalent strain
     * @param strain total strain vector in full form
     * @param gp integration point
     * @param atTime time step
     */
    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime);
    /**
     * Computes the equivalent local strain measure from given strain vector (full form).
     * @param kappa return param, comtaining the corresponding equivalent strain
     * @param strain total strain vector in full form
     * @param gp integration point
     * @param atTime time step
     */
    void computeLocalEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime)
    { IsotropicDamageMaterial1 :: computeEquivalentStrain(kappa, strain, gp, atTime); }

    /**
     * Implements the service updating local variables in given integration points,
     * which take part in nonlocal average process. Actually, no update is necessary,
     * because the value used for nonlocal averaging is strain vector used for nonlocal secant stiffness
     * computation. It is therefore necessary only to store local strain in corresponding status.
     * This service is declared at StructuralNonlocalMaterial level.
     * @param equivalentStrain equivalent strain vector in given integration point.
     * @param gp integration point to update.
     * @param atTime solution step indicating time of update.
     */
    virtual void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *atTime);

    /**
     * Computes the value of nonlocal weight function in given point.
     * @param src coordinates of source point.
     * @param coord coordinates of point, where nonlocal weight function is evaluated.
     * @return value of weight function.
     */
    virtual double computeWeightFunction(const FloatArray &src, const FloatArray &coord);
    /**
     * Determines, whether receiver has bounded weighting function (limited support)
     * @return true if weighting function bounded, zero otherwise
     */
    virtual int hasBoundedSupport() { return 1; }
    /**
     * Determines the width (radius) of limited support of weighting function
     */
    virtual void giveSupportRadius(double &radius) { radius = this->R; }

#ifdef __OOFEG
    ///Plots the sparse structure of stiffness contribution.
    virtual void NonlocalMaterialStiffnessInterface_showSparseMtrxStructure(GaussPoint *gp, oofegGraphicContext &gc, TimeStep *atTime);
#endif
    /**
     * Computes the damage parameter from given equivalent strain in given integration point.
     */
    void computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *g);

    /**@name Services required by NonlocalMaterialStiffnessInterface and related ones to support Nonlocal Stiffness*/
    //@{
    /// compute and add IP contributions to destination matrix
    virtual void NonlocalMaterialStiffnessInterface_addIPContribution(SparseMtrx &dest, const UnknownNumberingScheme &s, 
                                                                      GaussPoint *gp, TimeStep *atTime);
    /**
     * Returns integration list of receiver. Contains localIntegrationRecord structures, containing
     * references to integration points and their weights that influence to nonlocal average in
     * receiver's associated integration point.
     */
    virtual dynaList< localIntegrationRecord > *NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint *gp);
    /**
     * Computes the "local" part of nonlocal stiffness contribution assembled for given integration point.
     * @param gp source integration point
     * @param loc local code numbers
     * @param s determines the equation numbering scheme
     * @param lcontrib "local" contribution
     * @return nonzero if local point contributes (loading) or zero if not (unloading in elastic range, elastic)
     */
    int     giveLocalNonlocalStiffnessContribution(GaussPoint *gp, IntArray &loc, const UnknownNumberingScheme &s, 
                                                   FloatArray &lcontrib, TimeStep *atTime);
    /**
     * Computes the "remote" part of nonlocal stiffness contribution assembled for given integration point.
     * @param gp remote integration point
     * @param loc remote element code numbers
     * @param s determines the equation numbering scheme
     * @param rcontrib "remote" contribution
     */
    void     giveRemoteNonlocalStiffnessContribution(GaussPoint *gp, IntArray &rloc, const UnknownNumberingScheme &s, 
                                                     FloatArray &rcontrib, TimeStep *atTime);
    /**
     * Computes elastic stiffness for normal stress components
     * @param answer result of size (3,3)
     * @param mode determines the MatResponseMode
     * @param gp integration point
     * @param atTime time step
     */
    void giveNormalElasticStiffnessMatrix(FloatMatrix &answer,
                                          MatResponseMode rMode,
                                          GaussPoint *gp, TimeStep *atTime);

    //@}

#ifdef __PARALLEL_MODE
    /**
     * Updates domain before nonloc average (using updateDomainBeforeNonlocAverage service)
     * to ensure, that the localStrainVectorForAverage variable is correctly updated and
     * pack this localStrainVectorForAverage into given buffer.
     * @see Material::packUnknowns for description.
     * @param buff communication buffer
     * @param stepN solution step
     * @param ip integration point
     */
    int packUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip);
    /**
     * Unpack localStrainVectorForAverage value from given buffer.
     * @see Material::unpackAndUpdateUnknowns service.
     * @param buff communication buffer
     * @param stepN solution step.
     * @param ip integration point
     */
    int unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip);
    /**
     * Estimates the necessary pack size to hold all packed data of receiver.
     */
    int estimatePackSize(CommunicationBuffer &buff, GaussPoint *ip);
    /**
     *  Returns the weight representing relative computational cost of receiver
     *  The reference material model  is linear isotropic material - its weight is set to 1.0
     *  The other material models should compare to this reference model.
     */
    virtual double predictRelativeComputationalCost(GaussPoint *gp);
    /**
     * Returns the relative redistribution cost of the receiver
     */
    virtual double predictRelativeRedistributionCost(GaussPoint *gp) { return 1.0; }
#endif

    /// Creates the corresponding material status
    MaterialStatus *CreateStatus(GaussPoint *gp) const { return new IDNLMaterialStatus(1, IsotropicDamageMaterial1 :: domain, gp); }

protected:
    virtual void initDamaged(double kappa, FloatArray &totalStrainVector, GaussPoint *gp) { }
};

} // end namespace oofem
#endif // idmnl1_h
