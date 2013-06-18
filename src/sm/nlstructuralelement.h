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

#ifndef nlstructuralelement_h
#define nlstructuralelement_h

#include "structuralelement.h"
#include "domain.h"
#include "floatmatrix.h"

///@name Input fields for NLStructuralElement
//@{
#define _IFT_NLStructuralElement_nlgeoflag "nlgeo"
//@}

namespace oofem {
class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * Abstract base class for "structural" finite elements with geometrical nonlinearities.
 * The base StructuralElement supports nonlinearities at constitutive level.
 * The new implementation of relevant services is extended to incorporate general geometric
 * nonlinearity. The services for computing nonlinear parts of geometrical equations for
 * particular strain components are provided. The updated services for stiffness,
 * strain and internal forces computations are formulated.
 *
 * The activation of non-linear effects can be generally controlled on element level
 * using nlGeometry, allowing elements to be used also in linear computations.
 *
 * @note
 * Methods ComputeStressVector and computeStrainVector as they are implemented here
 * assume Total Lagrange approach.
 *
 * Tasks :
 * - Defining itself
 * - Calculating its contribution to the problem:
 *   - calculating its mass matrix M, its stiffness matrix K, its load vector
 *     f, its location array ;
 *   - calculating its contribution to the LHS and RHS of the linear system,
 *     using Static, Newmark,etc, formula. These contributions are usually
 *     combinations of M,K,f.
 * - Performing end-of-step operations:
 *   - calculating the strains and stresses at its Gauss points ;
 *   - printing its output in the data file and updating itself ;
 */
class NLStructuralElement : public StructuralElement
{
protected:
    /// Flag indicating if geometrical nonlinearities apply.
    int nlGeometry;

public:
    /**
     * Constructor. Creates element with given number, belonging to given domain.
     * @param n Element number.
     * @param d Domain to which new material will belong.
     */
    NLStructuralElement(int n, Domain *d);
    /// Destructor.
    virtual ~NLStructuralElement() { }

    /**
     * Return nlgeometry mode
     * 0 - small stran mode
     * 1 - G-L strain mode
     * 2 - deformation gradient mode
    */
    int giveGeometryMode(){return nlGeometry;}


    void computeSecondPKStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);



    /**
     * Computes the stiffness matrix of receiver.
     * The response is evaluated using @f$ \int (B_1+B_2(r))^{\mathrm{T}} D(B_1+B_2(r))\;\mathrm{d}v @f$, where
     * @f$ B_2 @f$ is nonlinear contribution evaluated using computeNLBMatrixAt service for each strain component
     * (@f$ B_2(i) = \Delta r^{\mathrm{T}} A(i) @f$).
     * Necessary transformations and reduced integration are taken into account.
     *
     * @param answer Computed stiffness matrix (symmetric).
     * @param rMode Response mode.
     * @param tStep Time step.
     */
    virtual void computeStiffnessMatrix(FloatMatrix &answer,
                                        MatResponseMode rMode, TimeStep *tStep);

    virtual void OLDcomputeStiffnessMatrix(FloatMatrix &answer,
                                        MatResponseMode rMode, TimeStep *tStep);

    /**
     * Computes numerically stiffness matrix of receiver. The response is evaluated
     * using @f$ \int (B_1+B_2(r))^{\mathrm{T}}D(B_1+B_2(r))\;\mathrm{d}v @f$, where
     * @f$ B_2 @f$ is nonlinear contribution evaluated using computeNLBMatrixAt service for each strain component
     * (@f$ B_2(i) = \Delta r^{\mathrm{T}} A(i) @f$).
     * Numerical integration procedure uses integrationRulesArray
     * for numerical integration. This implementation regards element integration rules as that they represent sub-cells
     * so that the integration is performed over all subcells for all terms.
     * For higher numerical performance, only one half of stiffness matrix is computed and answer is then symmetrized.
     * Therefore, if element matrix will be generally nonsymmetric, one must specialize this method.
     * Finally, the result is transformed into global coordinate system (or nodal coordinate system, if it is defined).
     * @param answer Computed stiffness matrix (symmetric).
     * @param rMode Response mode.
     * @param tStep Time step.
     */
    void computeStiffnessMatrix_withIRulesAsSubcells(FloatMatrix &answer,
                                                     MatResponseMode rMode, TimeStep *tStep);


    // stress equivalent vector (vector of internal forces) - for nonLinear Analysis.
    /**
     * Evaluates nodal representation of real internal forces.
     * The response is evaluated using @f$ F = \int (B+B_2)^{\mathrm{T}}\sigma\;\mathrm{d}V @f$ formula, where
     * @f$ B @f$ is linear strain-displacement contribution, @f$ B_2 @f$ is nonlinear contribution evaluated using
     * computeNLBMatrixAt service for each strain component (@f$ B_2(i) = \Delta r^{\mathrm{T}} A(i) @f$).
     * Necessary transformations are taken into account.
     *
     * @param answer Equivalent nodal forces vector.
     * @param tStep Time step
     * @param useUpdatedGpRecord If equal to zero, the stresses in integration points are computed (slow but safe).
     */
    virtual void OLDgiveInternalForcesVector(FloatArray &answer,
                                          TimeStep *tStep, int useUpdatedGpRecord = 0);


    virtual void giveInternalForcesVector(FloatArray &answer,
                                          TimeStep *tStep, int useUpdatedGpRecord = 0);

    /**
     * Evaluates nodal representation of real internal forces.
     * The response is evaluated using @f$ F = \int (B+B_2)^{\mathrm{T}}\sigma\;\mathrm{d}V @f$ formula, where
     * @f$ B @f$ is linear strain-displacement contribution, @f$ B_2 @f$ is nonlinear contribution evaluated using
     * computeNLBMatrixAt service for each strain component (@f$ B_2(i) = \Delta r^{\mathrm{T}} A(i) @f$).
     * Numerical integration procedure uses integrationRulesArray
     * for numerical integration. This implementation regards element integration rules as that they represent sub-cells
     * so that the integration is performed over all subcells for all terms.
     * Necessary transformations are taken into account.
     *
     * @param answer Equivalent nodal forces vector.
     * @param tStep Time step.
     * @param useUpdatedGpRecord If equal to zero, the stresses in integration points are computed (slow but safe).
     */
    void giveInternalForcesVector_withIRulesAsSubcells(FloatArray &answer,
                                                       TimeStep *tStep, int useUpdatedGpRecord = 0);

    /**
     * Compute strain vector of receiver evaluated at given integration point at time
     * step stepN from element displacement vector.
     * Total (green) strains are computed using following scheme
     * @f$ \varepsilon_G = Br+B_2(r) @f$, where @f$ B_2 @f$ is obtained using
     * computeNLBMatrixAt service for each strain component (@f$ B_2(i) = \Delta r^{\mathrm{T}} A(i) @f$)
     * and @f$ B @f$ is usual linear strain-displacement matrix.
     * Necessary transformation are taken into account.
     *
     * @param answer Element strain vector.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    // should be small strain only
    void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    void OLDcomputeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);


    void computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);
    void computeGreenLagrangeStrainVector(FloatArray &answer, FloatArray &F, MaterialMode matMode);

    //void trans_dSdE_2_dPdF(FloatMatrix &answer, const FloatMatrix &dSdE, const FloatArray &F, const FloatArray &S );

    // data management
    virtual IRResultType initializeFrom(InputRecord *ir);

    // definition
    virtual const char *giveClassName() const { return "NLStructuralElement"; }
    virtual classType giveClassID() const { return NLStructuralElementClass; }

protected:
    /**
     * Computes the nonlinear part of strain-displacement (geometrical) equation
     * related to i-th component of strain vector.
     * @param answer Returned nonlinear strain vector component contribution.
     * @param gp Integration point.
     * @param i Determines the component of strain vector for which contribution is assembled.
     * @see computeStrainVector
     */
    virtual void computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *gp, int i) {
        answer.resize(0, 0);
        return;
    }

    void dyadicProductBelow(FloatMatrix &answer, FloatArray &A, FloatArray &B);
    int giveVoigtIndexSym(int ind1, int ind2);

    void computeGLBMatrixAt(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep); 
    int checkConsistency();
    /**
     * Computes a matrix which, multiplied by the column matrix of nodal displacements,
     * gives the displacement gradient stored by columns.
     * The components of this matrix are derivatives of the shape functions,
     * but they are arranged in a somewhat different way from the usual B matrix.
     * @param gp Integration point.
     * @param answer BF matrix at this point.
     */ //@todo rename to computeBHmatrix as it is more appropriate /JB
    virtual void computeBFmatrixAt(GaussPoint *gp, FloatMatrix &answer) {
        OOFEM_ERROR("NLStructuralElement::computeBFmatrixAt : method not implemented for this element");
        return;
    }

    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer) {
        OOFEM_ERROR("NLStructuralElement::computeBHmatrixAt : method not implemented for this element");
        return;
    }
    friend class GradDpElement;
};
} // end namespace oofem
#endif // nlstructuralelement_h
