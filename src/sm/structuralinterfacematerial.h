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
#ifndef structuralinterfacematerial_h
#define structuralinterfacematerial_h

#include "material.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matconst.h"
#include "matstatus.h"

///@name Input fields for StructuralInterfaceMaterial
//@{

#define _IFT_StructuralInterfaceMaterial_useNumericalTangent "use_num_tangent"

//@}

namespace oofem {
class GaussPoint;

/**
 * Abstract base class for all "structural" interface models. It declares common  services provided
 * by all interface material models. The implementation of these services is partly left on derived classes,
 * which will implement constitutive model dependent part.
 * Some general purpose services are implemented on this level. For details, how to store
 * material model related history variables in integration points, see base class @ref Material documentation.
 *
 * Its capabilities can be examined using hasMaterialModeCapability  service.
 * It is generally assumed, that results obtained from constitutive model services are according to
 * valid material mode. This mode is determined from integration point, which is compulsory parameter of all material
 * services.
 * Full and reduced formats of stress/strain vectors are introduced for convenience.
 * Methods for converting vectors between full and reduced format are provided.
 *
 * Structural material introduces following basic stress/strain modes
 * - 3d state - all components of general stress/strain vector are generally nonzero.
 *   General 3d jump vector has following components {j_x, j_y, j_z}
 *
 * Generally speaking, there are following major tasks, covered by declared services.
 * - Computing engineering/first PK traction vector at integration point for given jump/discontinuity and updating its
 *   state (still temporary state, after overall equilibrium is reached).
 * - Updating its state (final state), when equilibrium has been reached.
 * - Storing/restoring its context to stream.
 * - Returning its material properties.
 *
 * @author Jim Brouzoulis
 */
class StructuralInterfaceMaterial : public Material
{
public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    StructuralInterfaceMaterial(int n, Domain *d);
    /// Destructor.
    virtual ~StructuralInterfaceMaterial() { }


    /**
     * Computes the first Piola-Kirchoff traction vector for given total jump/gap and integration point.
     * The total gap is computed from the displacement field at the given time step.
     * The service should use previously reached equilibrium history variables. Also
     * it should update temporary history variables in status according to newly reached state.
     * The temporary history variables are moved into equilibrium ones after global structure
     * equilibrium has been reached by iteration process.
     * @param answer Contains result.
     * @param gp Integration point.
     * @param reducedF Deformation gradient in in reduced form.
     * @param tStep Current time step (most models are able to respond only when tStep is current time step).
     */
    virtual void giveFirstPKTraction_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump,
                                        const FloatArray &reducedF, TimeStep *tStep)
    { _error("giveFirstPKTraction_1d: not implemented "); }
    virtual void giveFirstPKTraction_2d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump,
                                        const FloatArray &reducedF, TimeStep *tStep)
    { _error("giveFirstPKTraction_2d: not implemented "); }
    virtual void giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump,
                                        const FloatMatrix &F, TimeStep *tStep)
    { _error("giveFirstPKTraction_3d: not implemented "); }

    virtual void giveEngTraction_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep);
    virtual void giveEngTraction_2d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep);
    virtual void giveEngTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep);

    /**
     * Gives the tangent: @f$ \frac{\partial T}{\partial j} @f$.
     * Where T is the first PK traction and j is the spatial jump between the two sides, x(+) - x(-)
     * @param answer The computed tangent from the last evaluated first-PK-stress.
     * @param rMode Material mode.
     * @param gp Gauss point.
     * @param tStep Time step.
     */
    virtual void give1dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
    { _error("give1dStiffnessMatrix_dTdj: not implemented "); }
    virtual void give2dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
    { _error("give2dStiffnessMatrix_dTdj: not implemented "); }
    virtual void give3dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    virtual void give1dStiffnessMatrix_Eng(FloatMatrix &answer,  MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void give2dStiffnessMatrix_Eng(FloatMatrix &answer,  MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void give3dStiffnessMatrix_Eng(FloatMatrix &answer,  MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    // Numerical stiffness (intended to work regardless of dimension)
    virtual void giveStiffnessMatrix_dTdj_Num(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);


    /**
     * Tells if the model has implemented analytical tangent stiffness.
     * If not, the tangent must be computed numerically.
     */
    virtual bool hasAnalyticalTangentStiffness() const = 0;

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "StructuralInterfaceMaterial"; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);


    //virtual int setIPValue(const FloatArray &value, GaussPoint *gp, InternalStateType type);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    bool useNumericalTangent; ///@todo make private

protected:
    /**
     * Transforms traction vector into another coordinate system.
     * @param answer Transformed traction vector
     * @param base Transformation matrix. The columns in the matrix corresponds to the base vectors of the new
     * coordinate system to which we do transformation. These vectors must
     * be expressed in the same coordinate system as source strainVector.
     * @param strainVector 3d traction.
     * @param transpose Determines if we transpose matrix before transforming.
     */
    static void transformTractionTo(FloatArray &answer, const FloatMatrix &base,
                                    const FloatArray &strainVector, bool transpose = false);

    /**
     * Computes jump vector transformation matrix from standard vector transformation matrix.
     * @param answer Transformation matrix for strain vector.
     * @param base A (3,3) matrix, where on each column are stored unit direction vectors of
     * local coordinate axes to which we do transformation.
     * @param transpose Determines if we transpose matrix before transforming.
     */
    static void giveJumpTranformationMtrx(FloatMatrix &answer, const FloatMatrix &base,
                                          bool transpose = false);
};
} // end namespace oofem
#endif // StructuralInterfaceMaterial_h
