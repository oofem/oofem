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

#ifndef nlstructuralelement_h
#define nlstructuralelement_h

#include "../sm/Elements/structuralelement.h"

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
 * The base class StructuralElement supports nonlinearities at constitutive level.
 * The new implementation of relevant services is extended to incorporate general geometric
 * nonlinearity. The services for computing nonlinear parts of geometrical equations for
 * particular strain components are provided. The updated services for stiffness,
 * strain and internal forces computations are formulated.
 *
 * The activation of non-linear effects can be generally controlled on element level
 * using nlGeometry, allowing elements to be used also in linear computations.
 * nlGeometry = 0 - Small strain theory based on the engineering stress and strain
 * nlGeometry = 1 - Finite deformation theory formulated in terms of the deformation gradient F
 *                  and first Piola-Kirchoff stress P in the virtual work
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
 *
 * @author Jim Brouzoulis (among others)
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
    NLStructuralElement(int n, Domain * d);
    /// Destructor.
    virtual ~NLStructuralElement() { }

    /**
     * Returns the geometry mode describing the formulation used in the internal work
     * 0 - Engineering (small deformation) stress-strain mode
     * 1 - First Piola-Kirchhoff - Deformation gradient mode, P is defined as FS
     * 2 - Second Piola-Kirchhoff - Green-Lagrange strain mode with deformation gradient as input (deprecated and not supported)
     */
    int giveGeometryMode() { return nlGeometry; }

    /**
     * Computes the first Piola-Kirchhoff stress tensor on Voigt format. This method will
     * be called if nlGeo = 1 and mode = TL. This method computes the deformation gradient F and passes
     * it on to the crossection which then asks for the stress from the material.
     * @note P is related to S through F*S.
     *
     * @param answer Computed stress vector in Voigt form.
     * @param gp Gauss point at which the stress is evaluated.
     * @param tStep Time step.
     */
    void computeFirstPKStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    /**
     * Computes the Cauchy stress tensor on Voigt format. This method will
     * be called if nlGeo = 1 and mode = UL. This method computes the deformation gradient F and passes
     * it on to the crossection which then asks for the stress from the material.
     *
     * @param answer Computed stress vector in Voigt form.
     * @param gp Gauss point at which the stress is evaluated.
     * @param tStep Time step.
     */
    void computeCauchyStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    /**
     * Computes the stiffness matrix of receiver.
     * The response is evaluated using @f$ \int B_{\mathrm{H}}^{\mathrm{T}} D B_{\mathrm{H}} \;\mathrm{d}v @f$, where
     * @f$ B_{\mathrm{H}} @f$ is the B-matrix which produces the displacement gradient vector @f$ H_{\mathrm{V}} @f$ when multiplied with
     * the solution vector a.
     * Reduced integration are taken into account.
     *
     * @param answer Computed stiffness matrix.
     * @param rMode Response mode.
     * @param tStep Time step.
     */
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);


    /**
     * Computes the initial stiffness matrix of receiver. This method is used only if mode = UL
     * The response is evaluated using @f$ \int B ({\mathrm{\sigma}}\otimes \delta )B_{\mathrm{H}} \;\mathrm{d}v @f$, where
     * @f$ B @f$ is the classical B-matrix, but computed wrt updated node position
     *
     * @param answer Computed initial stiffness matrix.
     * @param tStep Time step.
     */
    virtual void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep);

    /**
     * Computes the stiffness matrix of receiver.
     * The response is evaluated using @f$ \int B_{\mathrm{H}}^{\mathrm{T}} D B_{\mathrm{H}} \;\mathrm{d}v @f$, where
     * @f$ B_{\mathrm{H}} @f$ is the B-matrix which produces the displacement gradient vector @f$ H_{\mathrm{V}} @f$ when multiplied with
     * the solution vector a.
     * @note Reduced intergration is not taken into account.
     * The integration procedure uses an integrationRulesArray for numerical integration. Each integration rule is
     * considered to represent a separate sub-cell/element. Typically this would be used when integration of the element
     * domain needs special treatment, e.g. when using the XFEM.
     *
     * @param answer Computed stiffness matrix.
     * @param rMode Response mode.
     * @param tStep Time step.
     */
    void computeStiffnessMatrix_withIRulesAsSubcells(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    /**
     * Evaluates nodal representation of real internal forces.
     * Necessary transformations are taken into account. @todo what is meant?
     *
     * @param answer Equivalent nodal forces vector.
     * @param tStep Time step
     * @param useUpdatedGpRecord If equal to zero, the stresses in integration points are computed (slow but safe).
     */

    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);

    /**
     * Evaluates nodal representation of real internal forces.
     *
     * Numerical integration procedure uses integrationRulesArray
     * for numerical integration. The integration procedure uses an integrationRulesArray for numerical integration.
     * Each integration rule is considered to represent a separate sub-cell/element. Typically this would be used when
     * integration of the element domain needs special treatment, e.g. when using the XFEM.
     *
     * @param answer Equivalent nodal forces vector.
     * @param tStep Time step.
     * @param useUpdatedGpRecord If equal to zero, the stresses in the integration points are computed (slow but safe).
     */
    void giveInternalForcesVector_withIRulesAsSubcells(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);

    /**
     * Computes the deformation gradient in Voigt form at integration point ip and at time
     * step tStep. Computes the displacement gradient and adds an identitiy tensor.
     *
     * @param answer Deformation gradient vector
     * @param gp Gauss point.
     * @param tStep Time step.
     */
    virtual void computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    /**
      * Computes the current volume of element
      */
    double computeCurrentVolume(TimeStep *tStep);

    // data management
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    // definition
    virtual const char *giveClassName() const { return "NLStructuralElement"; }

protected:
    int checkConsistency();
    /**
     * Computes a matrix which, multiplied by the column matrix of nodal displacements,
     * gives the displacement gradient stored by columns.
     * The components of this matrix are derivatives of the shape functions,
     * but they are arranged in a somewhat different way from the usual B matrix.
     * @param gp Integration point.
     * @param answer BF matrix at this point.
     */

    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer) {
        OOFEM_ERROR("method not implemented for this element");
        return;
    }
    friend class GradDpElement;
    friend class PhaseFieldElement;
    friend class XfemStructuralElementInterface;
};
} // end namespace oofem
#endif // nlstructuralelement_h
