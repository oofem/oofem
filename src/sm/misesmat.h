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

#ifndef misesmat_h
#define misesmat_h

#include "structuralmaterial.h"
#include "structuralms.h"
#include "linearelasticmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"

class GaussPoint;
class Domain;

class MisesMat : public StructuralMaterial
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
    LinearElasticMaterial *linearElasticMaterial;

    /// elastic shear modulus
    double G;

    /// elastic bulk modulus
    double K;

    /// hardening modulus
    double H;   
 
    /// initial (uniaxial) yield stress 
    double sig0; 

public:
    MisesMat(int n, Domain *d);
    ~MisesMat();

    /// specifies whether a given material mode is supported by this model
    int hasMaterialModeCapability(MaterialMode mode);

    /// reads the model parameters from the input file
    IRResultType initializeFrom(InputRecord *ir);

    // identification and auxiliary functions
    int hasNonLinearBehaviour()   { return 1; }
    const char *giveClassName() const { return "MisesMat"; }
    classType giveClassID()         const { return MisesMatClass; }

    /// returns a reference to the basic elastic material
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    //    virtual int         giveSizeOfFullHardeningVarsVector();
    //    virtual int         giveSizeOfReducedHardeningVarsVector(GaussPoint *);
    /// confirms that the stiffness matrix is symmetric
    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return true; }

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

    /**
     * Returns the integration point corresponding value in Reduced form.
     * @param answer contain corresponding ip value, zero sized if not available
     * @param aGaussPoint integration point
     * @param type determines the type of internal variable
     * @param type determines the type of internal variable
     * @returns nonzero if ok, zero if var not supported
     */
    //    int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);

    /**
     * Returns the mask of reduced indexes of Internal Variable component .
     * @param answer mask of Full VectorSize, with components beeing the indexes to reduced form vectors.
     * @param type determines the internal variable requested (physical meaning)
     * @returns nonzero if ok or error is generated for unknown mat mode.
     */
    //    int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);
};

//=============================================================================


class MisesMatStatus : public StructuralMaterialStatus
{
protected:
    /// plastic strain vector
    FloatArray plasticStrainVector;
    FloatArray tempPlasticStrainVector;

    /// cumulative plastic strain (hardening variable)
    double kappa;
    double tempKappa;

public:
    MisesMatStatus(int n, Domain *d, GaussPoint *g);
    ~MisesMatStatus();

    void givePlasticStrainVector(FloatArray& answer)
    {answer = plasticStrainVector;}

    void giveCumulativePlasticStrain(double& answer)
    {answer = kappa;}

    void letTempPlasticStrainVectorBe(FloatArray values)
    {tempPlasticStrainVector = values;}

    void setTempCumulativePlasticStrain(double value)
    {tempKappa = value;}

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
    const char *giveClassName() const { return "MisesMatStatus"; }

    /// identifies this class by its ID number
    classType             giveClassID() const
    { return MisesMatStatusClass; }
};

#endif // misesmat_h
