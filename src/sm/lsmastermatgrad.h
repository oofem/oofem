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

//   *******************************
//   *** CLASS Mises material
//   *******************************

#ifndef lsmastermatgrad_h
#define lsmastermatgrad_h

#include "lsmastermat.h"
#include "structuralmaterial.h"
#include "structuralms.h"
#include "linearelasticmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {
class GaussPoint;
class Domain;

class LsMasterMatGrad : public LsMasterMat
{
    /*
     * This class implements an isotropic elastoplastic material
     * with Mises yield condition, associated flow rule
     * and linear isotropic hardening.
     *
     * It differs from other similar materials (such as J2Mat)
     * by implementation - here we use the radial return, which
     * is the most efficient algorithm for this specific model.
     * Also, an extension to large strain will be available.
     *
     */

protected:
    /// reference to the basic elastic material
    
public:
    LsMasterMatGrad(int n, Domain *d);
    ~LsMasterMatGrad();
    
    /// specifies whether a given material mode is supported by this model
    int hasMaterialModeCapability(MaterialMode mode);
    IRResultType initializeFrom(InputRecord *ir);

    void giveCharacteristicMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime);
    void give3dKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void give3dGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void giveInternalLength(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    const char *giveClassName() const { return "LsMasterMatGrad"; }
    classType giveClassID()         const { return LsMasterMatClass; }

     MaterialStatus *CreateStatus(GaussPoint *gp) const;
      /// evaluates the material stiffness matrix
    void give3dMaterialStiffnessMatrix(FloatMatrix & answer,
                                       MatResponseForm, MatResponseMode,
                                       GaussPoint * gp,
                                       TimeStep * atTime);

    /// evaluates the stress
    void giveRealStressVector(FloatArray & answer,  MatResponseForm, GaussPoint *,
                              const FloatArray &, TimeStep *);

   
protected:


     

    //   virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);

    /**
     * Returns the mask of reduced indexes of Internal Variable component .
     * @param answer mask of Full VectorSize, with components beeing the indexes to reduced form vectors.
     * @param type determines the internal variable requested (physical meaning)
     * @returns nonzero if ok or error is generated for unknown mat mode.
     */
    //    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);

    /**
     * Returns the type of internal variable (scalar, vector, tensor,...).
     * @param type determines the type of internal variable
     * @returns type of internal variable
     */
    //   virtual InternalStateValueType giveIPValueType(InternalStateType type);

    /**
     * Returns the corresponding integration point  value size in Reduced form.
     * @param type determines the type of internal variable
     * @returns var size, zero if var not supported
     */
    //   virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);
};

//=============================================================================


class LsMasterMatGradStatus : public LsMasterMatStatus
{
protected:

public:
  LsMasterMatGradStatus(int n, Domain *d, GaussPoint *g, int s);
    ~LsMasterMatGradStatus();


     void printOutputAt(FILE *file, TimeStep *tStep);

    /// initializes the temporary status
    virtual void initTempStatus();

    /// updates the state after a new equilibrium state has been reached
    virtual void updateYourself(TimeStep *);

    /// saves the current context(state) into a stream
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    /// restores the state from a stream
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /// identifies this class by its name
    const char *giveClassName() const { return "LsMasterMatGradStatus"; }

    /// identifies this class by its ID number
    classType             giveClassID() const
    { return LsMasterMatStatusClass; }
};
} // end namespace oofem
#endif // misesmat_h
