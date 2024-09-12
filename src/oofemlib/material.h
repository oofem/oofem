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

#ifndef material_h
#define material_h

#include "femcmpnn.h"

#include "matconst.h"
#include "matstatus.h"
#include "materialmode.h"
#include "timestep.h"
#include "internalstatetype.h"
#include "internalstatevaluetype.h"
#include "matresponsemode.h"
#include "dictionary.h"
#include "chartype.h"

///@name Input fields for Material
//@{
#define _IFT_Material_density "d"
#define _IFT_Material_castingtime "castingtime"
#define _IFT_Material_preCastingTimeMat "precastingtimemat"
//@}

namespace oofem {
#define STRAIN_STEPS 10.0

class GaussPoint;
class Dictionary;
class FloatArray;
class FloatMatrix;
class Element;
class ProcessCommunicator;

/**
 * Abstract base class for all material models. Declares the basic common interface
 * to all material models. Derived classes should expand this interface, because they are
 * assumed to be  base classes for analysis specific tasks (for example mechanical or
 * thermal analysis).
 *
 * Instance of integration point class is assumed to be implicit argument to
 * all method, depending on internal state in point of consideration.
 * To provide opportunity for storing arbitrary material model related history variables
 * in integration points, associated material status class is introduced.
 * Each new material model class should be declared together with its associated status class
 * (derived from MaterialStatus class). This status can be seen as simple container,
 * storing necessary history variables and providing some access and modification methods.
 * Each integration point can contain material status. Material model should create
 * unique copy of its associated status in each integration point.
 * Because integration point is parameter of all messages to material model
 * class, material model therefore can easily access  all history variables it needs.
 *
 * The attribute 'propertyDictionary' contains all the properties of a material
 * like its Young modulus, its mass density or Poisson ratio.
 *
 * Its task is to indicate whether there required material mode is valid for receiver
 * (method hasMaterialModeCapability). Note: for some material models and linear materials
 * there need not exist support for assembling material char matrix at material level,
 * all is handled properly at crossSection level (_2dBeam mode, 3dShellMode, ...).
 * But this function must indicate whether mode is valid or not for real stress computation.
 *
 * @see MaterialStatus class
 * @see GaussPoint class
 */
class OOFEM_EXPORT Material : public FEMComponent
{
protected:
    /**
     * Property dictionary.
     * Can be used to store constant material parameters, which
     * are same for all integration points.
     * Note: Try to avoid using the dictionary because of a very slow access. Use rather separate variables to
     * store material parameters.
     */
    Dictionary propertyDictionary;

    /**
     * Casting time. For solution time less than casting time the material
     * is assumed to have no stiffness etc. This attribute is declared here,
     * but support for this functionality must be incorporated by particular
     * material model
     */
    double castingTime;
    
    /// Material existing before casting time - optional parameter, zero by default
    int preCastingTimeMat;
    

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    Material(int n, Domain *d);
    /// Destructor.
    virtual ~Material() = default;

    /**
     * Returns true if stiffness matrix of receiver is symmetric
     * Default implementation returns true.
     */
    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const { return true; }
    /**
     * @brief Returns characteristic matrix of the receiver
     * 
     */
    virtual void giveCharacteristicMatrix(FloatMatrix &answer, MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const {}
    /**
     * @brief Returns characteristic vector of the receiver
     * 
     */
    virtual void giveCharacteristicVector(FloatArray &answer, FloatArray& flux, MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const {}
    /**
     * @brief Returns characteristic value of the receiver
     * 
     */
    virtual double giveCharacteristicValue(MatResponseMode type, GaussPoint* gp, TimeStep *tStep) const ;
    /**
     * Returns the value of material property 'aProperty'. Property must be identified
     * by unique int id. Integration point also passed to allow for materials with spatially
     * varying properties
     * @param aProperty ID of property requested.
     * @param gp Integration point,
     * @return Property value.
     */
    virtual double give(int aProperty, GaussPoint *gp) const;
    /**
     * Returns true if 'aProperty' exists on material.
     * @param aProperty ID of property requested.
     * @param gp Integration point.
     * @return True if 'aProperty' exists.
     */
    virtual bool hasProperty(int aProperty, GaussPoint *gp) const;
    /**
     * Modify 'aProperty', which already exists on material. Intended for evolving material properties.
     * @param aProperty ID of a property requested.
     * @param value Assigned value.
     * @param gp Integration point.
     */
    virtual void modifyProperty(int aProperty, double value, GaussPoint *gp);
    /**
     * @return Casting time of the receiver.
     */
    double giveCastingTime() const { return this->castingTime; }
    /**
     * @param tStep Time step to check activity for.
     * @return True if material is activated for given solution step.
     */
    virtual bool isActivated(TimeStep *tStep) const {
        if ( tStep ) {
            return ( tStep->giveTargetTime() >= this->castingTime );
        } else {
            return true;
        }
    }

    // identification and auxiliary functions
    /**
     * Tests if material supports material mode.
     * @param mode Required material mode.
     * @return Nonzero if supported, zero otherwise.
     */
    virtual bool hasMaterialModeCapability(MaterialMode mode) const;

    /**
     * Tests if material supports casting time
     * @return Nonzero if supported, zero otherwise.
     */
    virtual bool hasCastingTimeSupport() const;

    ///@name Access functions for internal states. Usually overloaded by new material models.
    //@{
    /**
     * Sets the value of a certain variable at a given integration point to the given value.
     * @param value Contains the value(s) to be set (in reduced form).
     * @param gp Integration point.
     * @param type Determines the type of internal variable.
     * @param type Determines the type of internal variable.
     * @returns Nonzero if ok, zero if var not supported.
     */
    virtual int setIPValue(const FloatArray &value, GaussPoint *gp, InternalStateType type)
    { return 0; }
    /**
     * Returns the integration point corresponding value in Reduced form.
     * @param answer Contain corresponding ip value, zero sized if not available.
     * @param gp Integration point to which the value refers.
     * @param type Determines the type of internal variable.
     * @param tStep Determines the time step.
     * @returns Nonzero if the assignment can be done, zero if this type of variable is not supported.
     */
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    //@}

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    void printYourself() override;

    /**
     * Stores integration point state to output stream.
     * @param stream Output stream.
     * @param mode Determines amount of info required in stream (state, definition, ...).
     * @param gp integration point.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    virtual void saveIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp);
    /**
     * Reads integration point state to output stream.
     * @param stream Output stream.
     * @param mode Determines amount of info required in stream (state, definition, ...).
     * @param gp integration point.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    virtual void restoreIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp);

    /**
     * Allows programmer to test some internal data, before computation begins.
     * For example, one may use this function, to ensure that element has material with
     * required capabilities is assigned to element. This must be done after all
     * mesh components are instanciated.
     * @return Nonzero if receiver is consistent.
     */
    int checkConsistency() override;
    /**
     * Restores consistency of the status, i.e., computes or corrects
     * the values of certain status variables such that the state is admissible.
     * For instance, if the initial values of some internal variables
     * are read from a file, other internal variables are adjusted accordingly.
     */
    virtual void restoreConsistency(GaussPoint *gp) { }
    /**
     * Optional function to call specific procedures when initializing a material.
     * For example, multiscale simulations need to create master and slave material statuses on specific integration points before the computation.
     * @param element Pointer to element.
     * @return Zero on error.
     */
    virtual int initMaterial(Element *element);
    /**
     * Returns material status of receiver in given integration point.
     * If status does not exist yet, it is created using CreateStatus member function.
     * @param gp Returns reference to material status belonging to integration
     * point gp.
     * @return Material status associated with given integration point.
     */
    virtual MaterialStatus *giveStatus(GaussPoint *gp) const;
    /*
     * In the case of nonlocal constitutive models,
     * the use of multiple inheritance is assumed. Typically, the class representing nonlocal
     * constitutive model is derived both from class representing local model and from class
     * NonlocalMaterialExtension or from one of its derived classes
     * (which declare services and variables corresponding to specific analysis type).
     * Or alternatively, material model can introduce stronger form of this relation,
     * it can possess object of NonlocalMaterialExtension.
     * @return In both cases, this function returns pointer to this object, obtained by
     * returning address of component or using pointer conversion from receiver to base class
     * NonlocalMaterialExtension. If no nonlocal extension exists, NULL pointer is returned.
     */
    //virtual NonlocalMaterialExtension* giveNonlocalMaterialExtensionPtr () {return NULL;}

    /**
     * Pack all necessary data of integration point (according to element parallel_mode)
     * into given communication buffer. The nature of packed data is material model dependent.
     * Typically, for material of "local" response (response depends only on integration point local state)
     * no data are exchanged. For "nonlocal" constitutive models the send/receive of local values which
     * undergo averaging is performed between local and corresponding remote elements.
     * @param buff Communication buffer.
     * @param tStep Solution step.
     * @param ip Integration point.
     */
    virtual int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) { return 1; }
    /**
     * Unpack and updates all necessary data of given integration point (according to element parallel_mode)
     * into given communication buffer.
     * @see packUnknowns service.
     * @param buff Communication buffer.
     * @param tStep Solution step.
     * @param ip Integration point.
     */
    virtual int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) { return 1; }
    /**
     * Estimates the necessary pack size to hold all packed data of receiver.
     */
    virtual int estimatePackSize(DataStream &buff, GaussPoint *ip) { return 0; }
    /**
     * Returns the weight representing relative computational cost of receiver
     * The reference material model is linear isotropic material - its weight is set to 1.0
     * The other material models should compare to this reference model.
     */
    virtual double predictRelativeComputationalCost(GaussPoint *gp) { return 1.0; }
    /**
     * Returns the relative redistribution cost of the receiver
     */
    virtual double predictRelativeRedistributionCost(GaussPoint *gp) { return 1.0; }

    /**
     * Creates new copy of associated status and inserts it into given integration point.
     * @param gp Integration point where newly created status will be stored.
     * @return Reference to new status.
     */
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const
    { return nullptr; }

    /**
     * Initializes temporary variables stored in integration point status
     * at the beginning of new time step.
     * Temporary history variables (they describe state of material during
     * solution of time step) are initialized according to history variables, which
     * describe state corresponding to previous equilibrium solution.
     * Default implementation simply extracts status from integration point and
     * calls its initTempStatus method.
     */
    virtual void initTempStatus(GaussPoint *gp) const;

    /**
     * Stores receiver state to output stream.
     * @param stream Output stream.
     * @param mode Determines amount of info required in stream (state, definition, ...).
     * @exception throws an ContextIOERR exception if error encountered.
     */
    void saveContext(DataStream &stream, ContextMode mode) override;
    /**
     * Restores the receiver state previously written in stream.
     * @see saveContext
     * @param stream Input stream.
     * @param mode Determines amount of info available in stream (state, definition, ...).
     * @exception throws an ContextIOERR exception if error encountered.
     */
    void restoreContext(DataStream &stream, ContextMode mode) override;
};
} // end namespace oofem
#endif // material_h
