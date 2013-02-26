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
//   *** CLASS Large-Strain Master Material
//   *******************************

#ifndef lsmastermat_h
#define lsmastermat_h

#include "structuralmaterial.h"
#include "structuralms.h"
#include "linearelasticmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {
class GaussPoint;
class Domain;

class LsMasterMat : public StructuralMaterial
{
    /*
     * This class implements a large strain master material.
     * stress and stiffness are computed from small strain(slaveMat) material model 
     * using a strain tensor from the Seth-Hill strain tensors familly(depends on parameter m,
     * m = 0 logarithmic strain,m = 1 Green-Lagrange strain ...)
     * then stress and stiffness are transformed 2.PK stress and appropriate stiffness 
     */

protected:
    /// reference to the basic elastic material
    LinearElasticMaterial *linearElasticMaterial;

    /// 'slave' material model number
    int slaveMat;
    /// specifies the strain tensor
    double m;

 
public:
    LsMasterMat(int n, Domain *d);
    ~LsMasterMat();
    
    /// specifies whether a given material mode is supported by this model
    int hasMaterialModeCapability(MaterialMode mode);

    /// reads the model parameters from the input file
    IRResultType initializeFrom(InputRecord *ir);

    // identification and auxiliary functions
    int hasNonLinearBehaviour()   { return 1; }
    const char *giveClassName() const { return "LsMasterMat"; }
    classType giveClassID()     const { return LsMasterMatClass; }

    /// returns a reference to the basic elastic material
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    //    virtual int         giveSizeOfFullHardeningVarsVector();
    //    virtual int         giveSizeOfReducedHardeningVarsVector(GaussPoint *);
    /// confirms that the stiffness matrix is symmetric
    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    /// creates a new material status  corresponding to this class
    MaterialStatus *CreateStatus(GaussPoint *gp) const;

      /// evaluates the material stiffness matrix
    void give3dMaterialStiffnessMatrix(FloatMatrix & answer,
                                       MatResponseForm, MatResponseMode,
                                       GaussPoint * gp,
                                       TimeStep * atTime);

    /// evaluates the stress
    void giveRealStressVector(FloatArray & answer,  MatResponseForm, GaussPoint *,
                              const FloatArray &, TimeStep *);

    /// transformation matrix
    void constructTransformationMatrix(FloatMatrix F, GaussPoint *gp);

    

protected:


     

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
};

//=============================================================================


class LsMasterMatStatus : public StructuralMaterialStatus
{
protected:
  FloatMatrix Pmatrix,TLmatrix, transformationMatrix;
  int slaveMat;

public:
  LsMasterMatStatus(int n, Domain *d, GaussPoint *g, int s);
    ~LsMasterMatStatus();


    void givePmatrix(FloatMatrix &answer)
    { answer = Pmatrix; }
    void giveTLmatrix(FloatMatrix &answer)
    { answer = TLmatrix; }
    void giveTransformationMatrix(FloatMatrix &answer)
    { answer = transformationMatrix;}

    void setPmatrix(FloatMatrix values) { Pmatrix = values; }
    void setTLmatrix(FloatMatrix values) { TLmatrix = values; }
    void setTransformationMatrix(FloatMatrix values){transformationMatrix = values;}
    /// prints the output variables into the *.out file
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
    const char *giveClassName() const { return "LsMasterMatStatus"; }

    /// identifies this class by its ID number
    classType             giveClassID() const
    { return LsMasterMatStatusClass; }
};
} // end namespace oofem
#endif // misesmat_h
