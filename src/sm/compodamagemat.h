/* v 1.8 2009/04/27 vs  */
/*
 *
 *****    *****   ******  ******  ***   ***
 **   **  **   **  **      **      ** *** **
 **   **  **   **  ****    ****    **  *  **
 **   **  **   **  **      **      **     **
 **   **  **   **  **      **      **     **
 *****    *****   **      ******  **     **
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 2009   Vit Smilauer
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

#ifndef compodamagemat_h
#define compodamagemat_h

#include "material.h"
#include "linearelasticmaterial.h"
#include "structuralmaterial.h"
#include "structuralms.h"
#include "cltypes.h"

namespace oofem {

/**
 * Material damage model for transversely orthotropic material. Fixed cracks are induced in principal material coordinates, always perpendicular to material axis. Six cracking modes are implemented in three orthogonal directions - three for tension/compression failure and three for shear failure. Material orientation can be specified on each element with "lmcs" keyword, otherwise global material orientation is assumed.
 *
 * Evolution of cracks is determined separately in compression and tension. Linear softening is assumed (can be easily changed to exponential but requires then inner iterations)
 *
 * For derivation of the model see the book Bazant and Planas, Fracture and Size Effect in Concrete and Other Quasibrittle Materials, pp.236 or article Bazant and Oh: Crack band theory for fracture of concrete, Materials and Structures, 1983.
 *
 * The model is aimed for 3D problems but extension for 1D truss works. In this particular case, only the first array component is used in all involved variables.
 *
 * Prefix temp* refers to unequilibrated values, e.g. tempOmega[] is a temporal damage array
 */


//class for manipulating values at GP
class CompoDamageMatStatus : public StructuralMaterialStatus
{
public:
    /// Constructor
    CompoDamageMatStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    ~CompoDamageMatStatus();

    /// Prints the receiver state to stream
    void   printOutputAt(FILE *file, TimeStep *tStep);

    /// Initializes temp variables - called at the beginning of each time increment (not iteration), resets tempStressVector, tempStrainVector
    void initTempStatus();

    /// Update after new equilibrium state reached
    void updateYourself(TimeStep *);

    //tempVal are values used during iteration, Val are equilibrated values, updated after last iteration in previous time step

    /// Highest strain ever reached at IP. Can be unequilibrated from last iterations [6 tension, 6 compression]
    FloatArray tempKappa;

    /// Highest strain ever reached in all previous equilibrated steps [6 tension, 6 compression]
    FloatArray kappa;

    /// Highest damage ever reached at IP. Can be unequilibrated from last iterations [6 for tension and compression]
    FloatArray tempOmega;

    /// Checks whether snapback occured at IP.
    IntArray hasSnapBack;

    /// Highest damage ever reached in all previous equilibrated steps at IP [6 for tension and compression]
    FloatArray omega;

    ///Iteration in the time step
    int Iteration;

    /// Stress at which damage starts. For uniaxial loading is equal to given maximum stress in the input. The stress is linearly interpolated between increments at IP [6 tension, 6 compression]
    FloatArray initDamageStress;
    /// Strain when damage is initiated at IP. In literature denoted eps_0 [6 tension, 6 compression]
    FloatArray strainAtMaxStress;
    /// Maximum strain when stress becomes zero due to complete damage (omega = 1) at IP. Determined from fracture energy and characteristic length. Derived for mode I only on the uniaxial loading case. In literature denoted often eps_f [6 tension, 6 compression]
    FloatArray maxStrainAtZeroStress;

    /// only for printing purposes in CompoDamageMatStatus
    FloatArray tempStressMLCS;

    /// Characteristic element length at IP in three perpendicular planes aligned with material orientation
    FloatArray elemCharLength;

    /**
     * Stores context of receiver into given stream.
     * Only non-temp internal history variables are stored.
     * @param stream stream where to write data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     */
    contextIOResultType   saveContext(DataStream *stream, ContextMode mode, void *obj);
    /**
     * Restores context of receiver from given stream.
     * @param stream stream where to read data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     */
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj);

    const char *giveClassName() const { return "CompoDamageMatStatus"; }
    classType             giveClassID() const { return CompoDamageMatStatusClass; }
private:
};




//implementation of composite damage model
class CompoDamageMat : public StructuralMaterial
{
public:
    /// Constructor
    CompoDamageMat(int n, Domain *d);
    /// Destructor
    ~CompoDamageMat();
    const char *giveClassName() const { return "CompositeDamageMaterial"; }
    classType             giveClassID() const { return CompositeDamageMaterialClass; }
    /// Returns input record name of the receiver.
    const char *giveInputRecordName() const { return "compodamagemat"; }

    /**
     * Initializes receiver acording to object description stored in input record..
     * The density of material is read into property dictionary (keyword 'd')
     */
    IRResultType initializeFrom(InputRecord *ir);
    /** Setups the input record string of receiver
     * @param str string to be filled by input record
     * @param keyword print record keyword (default true)
     */
    int giveInputRecordString(std :: string &str, bool keyword = true);


    /// updates nonTemp variables when equilibrium reached
    void updateYourself(TimeStep *); // update after new equilibrium state reached
    /// Creates corresponding material status
    MaterialStatus *CreateStatus(GaussPoint *gp) const { return new CompoDamageMatStatus(1, domain, gp); }

    /**
     * Computes full 3d material secant stiffness matrix at given integration point, time, respecting load history in integration point. The matrix is rotated to global coordinate system
     * @param answer computed results
     * @param form material response form (no effect)
     * @param mode material response mode (no effect)
     * @param gp integration point
     * @param atTime time step
     */
    void give3dMaterialStiffnessMatrix(FloatMatrix & answer,
                                       MatResponseForm, MatResponseMode,
                                       GaussPoint * gp,
                                       TimeStep * atTime);
    /**
     * Computes the real stress vector for given strain increment and integration point based on reduced strain. Contains constitutive equation of damage model. Called in each iteration. Temp vars are updated accordingly
     * @param answer contains result
     * @param form material response form (unused)
     * @param gp integration point
     * @param reducedStrain strain vector in reduced form
     * @param tStep current time step
     */
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
    int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    /**
     * Returns the type of internal variable (scalar, vector, tensor,...).
     * @param type determines the type of internal variable
     * @returns type of internal variable
     */
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    /**
     * Returns the mask of reduced indexes of Internal Variable component .
     * @param answer mask of Full VectorSize, with components beeing the indexes to reduced form vectors.
     * @param type determines the internal variable requested (physical meaning)
     * @returns nonzero if ok or error is generated for unknown mat mode.
     */
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);
    /**
     * Returns the corresponding integration point  value size in Reduced form.
     * @param type determines the type of internal variable
     * @returns var size, zero if var not supported
     */
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);

    ///Optional parameter determining after how many iterations within the time step the damage is calculated. This is important for stress evaluation which is unequilibrated in the beginning. Variables strainAtMaxStress, initDamageStress, maxStrainAtZeroStress are evaluated afterIter.
    int afterIter;

protected:
    /**
     * Returns 3D material stiffness matrix [6x6] in unrotated form. The matrix is reduced by by omega variables
     * @param answer full symmetric matrix
     * @param gp integration point
     */
    void giveUnrotated3dMaterialStiffnessMatrix(FloatMatrix &answer, GaussPoint *gp);
    /**
     * Returns [6x6] rotation matrix in global coordinate system. The matrix relates local c.s. to global c.s. Local c.s. can be specified with 'lmcs' flag defined on element.
     * @param gp integration point
     * @returns pointer to matrix
     */
    int giveMatStiffRotationMatrix(FloatMatrix &answer, GaussPoint *gp);

    /// Six stress components of tension components read from the input file
    FloatArray inputTension;
    /// Six stress components of compression components read from the input file
    FloatArray inputCompression;

    /// Stress components which are allowed for snap back [6 tension, 6 compression]
    IntArray allowSnapBack;

    /**
     * Fills array elemCharLength with characteristic length related to three perpendicular planes. The planes are of the same orientation as material.
     * @param status pointer to integration point's status
     * @param gp integration point
     * @param elementCs material orientation matrix
     */
    void giveCharLength(CompoDamageMatStatus *status, GaussPoint *gp, FloatMatrix &elementCs);
    /**
     * Computes characteristic length for fixed planes of material orientation.
     * @param charLenModes returns six lengths
     * @param gp integration point
     */
    void giveCharLengthForModes(FloatArray &charLenModes, GaussPoint *gp);
    /**
    * Check that elemnt is small or Gf large enough to prevent snap-back
    * @param gp integration point
    */
    void checkSnapBack(GaussPoint *gp, MaterialMode mMode);
};

} // end namespace oofem
#endif // compodamagemat_h
