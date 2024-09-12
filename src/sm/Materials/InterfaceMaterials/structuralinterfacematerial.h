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
#include "floatarrayf.h"

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
    StructuralInterfaceMaterial(int n, Domain * d);

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
    virtual double giveFirstPKTraction_1d(double jump, double reducedF, GaussPoint *gp, TimeStep *tStep) const;
    virtual FloatArrayF<2> giveFirstPKTraction_2d(const FloatArrayF<2> &jump, const FloatMatrixF<2,2> &reducedF, GaussPoint *gp, TimeStep *tStep) const;
    virtual FloatArrayF<3> giveFirstPKTraction_3d(const FloatArrayF<3> &jump, const FloatMatrixF<3,3> &F, GaussPoint *gp, TimeStep *tStep) const
    { OOFEM_ERROR("not implemented ");  }

    virtual double giveEngTraction_1d(double jump, GaussPoint *gp, TimeStep *tStep) const;
    virtual FloatArrayF<2> giveEngTraction_2d(const FloatArrayF<2> &jump, GaussPoint *gp, TimeStep *tStep) const;
    virtual FloatArrayF<3> giveEngTraction_3d(const FloatArrayF<3> &jump, GaussPoint *gp, TimeStep *tStep) const;

    /**
     * Gives the tangent: @f$ \frac{\partial T}{\partial j} @f$.
     * Where T is the first PK traction and j is the spatial jump between the two sides, x(+) - x(-)
     * @param answer The computed tangent from the last evaluated first-PK-stress.
     * @param rMode Material mode.
     * @param gp Gauss point.
     * @param tStep Time step.
     */
    virtual FloatMatrixF<1,1> give1dStiffnessMatrix_dTdj(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const;
    virtual FloatMatrixF<2,2> give2dStiffnessMatrix_dTdj(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const;
    virtual FloatMatrixF<3,3> give3dStiffnessMatrix_dTdj(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const;

    virtual FloatMatrixF<1,1> give1dStiffnessMatrix_Eng(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const;
    virtual FloatMatrixF<2,2> give2dStiffnessMatrix_Eng(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const;
    virtual FloatMatrixF<3,3> give3dStiffnessMatrix_Eng(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const;

    FloatMatrixF<1,1> give1dStiffnessMatrix_dTdj_Num(GaussPoint *gp, TimeStep *tStep) const;
    FloatMatrixF<2,2> give2dStiffnessMatrix_dTdj_Num(GaussPoint *gp, TimeStep *tStep) const;
    FloatMatrixF<3,3> give3dStiffnessMatrix_dTdj_Num(GaussPoint *gp, TimeStep *tStep) const;

    FloatMatrixF<1,1> give1dStiffnessMatrix_Eng_Num(GaussPoint *gp, TimeStep *tStep) const;
    FloatMatrixF<2,2> give2dStiffnessMatrix_Eng_Num(GaussPoint *gp, TimeStep *tStep) const;
    FloatMatrixF<3,3> give3dStiffnessMatrix_Eng_Num(GaussPoint *gp, TimeStep *tStep) const;

    /**
     * Tells if the model has implemented analytical tangent stiffness.
     * If not, the tangent must be computed numerically.
     */
    virtual bool hasAnalyticalTangentStiffness() const = 0;

    // identification and auxiliary functions
    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    virtual FloatArray giveInterfaceStrength() { return {0}; }

    //virtual int setIPValue(const FloatArray &value, GaussPoint *gp, InternalStateType type);
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    bool useNumericalTangent; ///@todo make private
};
} // end namespace oofem
#endif // structuralinterfacematerial_h
