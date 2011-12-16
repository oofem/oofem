/* $Header: /home/cvs/bp/oofem/sm/src/slipmodelcebfip90.h,v 1.1.2.1 2004/04/05 15:19:47 bp Exp $ */
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

//   **********************************************************************
//   *** CLASS SIMPLE INTERFACE MATERIAL FOR INTERFACE ELEMENTS   ************
//   **********************************************************************

#ifndef simpleinterfacemat_h
#define simpleinterfacemat_h

#include "material.h"
#include "linearelasticmaterial.h"
#include "structuralmaterial.h"
#include "structuralms.h"

namespace oofem {
// material contant's keys for give()
class GaussPoint;

/**
 * This class implements associated Material Status to IsoInterfaceDamageMaterial.
 * It is atribute of matStatusDictionary at every GaussPoint, for which this material
 * is active.
 */
class SimpleInterfaceMaterialStatus : public StructuralMaterialStatus
{
    /*
     * This class implements associated Material Status to SimpleInterfaceMaterial.
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
     */

protected:
    FloatArray shearStressShift, tempShearStressShift;

public:
    /// Constructor
    SimpleInterfaceMaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    ~SimpleInterfaceMaterialStatus();

    /// Prints the receiver state to stream
    void   printOutputAt(FILE *file, TimeStep *tStep);

    // definition
    const char *giveClassName() const { return "SimpleInterfaceMaterialStatus"; }
    classType             giveClassID() const { return SimpleInterfaceMaterialStatusClass; }

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

    FloatArray giveShearStressShift();
    FloatArray giveTempShearStressShift();
    void setTempShearStressShift(FloatArray newShearStressShift) { tempShearStressShift = newShearStressShift; };

    // saves current context(state) into stream
    /**
     * Stores context of receiver into given stream.
     * Only non-temp internal history variables are stored.
     * @param stream stream where to write data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     */
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores context of receiver from given stream.
     * @param stream stream where to read data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     */
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};


/**
 * Base class representing general isotropic damage model.
 * It is based on isotropic damage concept, assuming that damage evolution law
 * is postulated in explicit form, relatin damage parameter (omega) to scalar measure
 * of the largest strain level ever reached in material (kappa).
 */
class SimpleInterfaceMaterial : public StructuralMaterial
{
protected:
    double kn;
    double stiffCoeff;
    double frictCoeff;

public:
    /// Constructor
    SimpleInterfaceMaterial(int n, Domain *d);
    /// Destructor
    ~SimpleInterfaceMaterial();

    /// Returns nonzero indicating that receiver is nonlinear
    int hasNonLinearBehaviour()   { return 1; }
    /**
     * Tests, if material supports material mode.
     * @param mode required material mode
     * @return nonzero if supported, zero otherwise
     */
    int hasMaterialModeCapability(MaterialMode mode);
    const char *giveClassName() const { return "SimpleInterfaceMaterial"; }
    classType giveClassID()         const { return SimpleInterfaceMaterialClass; }

    /**
     * Computes full 3d material stiffness matrix at given integration point, time, respecting load history
     * in integration point.
     * @param answer computed results
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseForm from,
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *atTime);

    /**
     * Computes the real stress vector for given strain increment and integration point.
     * Temp vars are updated accordingly
     * @param answer contains result
     * @param form material response form
     * @param gp integration point
     * @param reducedStrain strain  vector in reduced form
     * @param tStep current time step (most models are able to respond only when atTime is current time step)
     */
    void giveRealStressVector(FloatArray & answer,  MatResponseForm, GaussPoint *,
                              const FloatArray &, TimeStep *);

    virtual void  giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime);

    virtual int giveStressStrainComponentIndOf(MatResponseForm, MaterialMode mmode, int);
    virtual void giveStressStrainMask(IntArray & answer, MatResponseForm, MaterialMode mmode) const;
    virtual int giveSizeOfReducedStressStrainVector(MaterialMode);
    void giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *,
                                         const FloatArray &charVector3d);
    void giveFullCharacteristicVector(FloatArray &answer,  GaussPoint *,
                                      const FloatArray &);

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
    /**
     * Returns a vector of coefficients of thermal dilatation in direction
     * of each material principal (local) axis.
     * @param answer vector of thermal dilatation coefficients
     * @param gp integration point
     * @param tStep time step (most models are able to respond only when atTime is current time step)
     */
    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *, TimeStep *);

    /**
     * Instanciates the receiver from input record
     */
    IRResultType initializeFrom(InputRecord *ir);
    /** Setups the input record string of receiver
     *  @param str string to be filled by input record
     *  @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    /// Creates corresponding material status
    MaterialStatus *CreateStatus(GaussPoint *gp) const { return new SimpleInterfaceMaterialStatus(1, domain, gp); }
};
} // end namespace oofem
#endif // simpleinterfacemat_h
