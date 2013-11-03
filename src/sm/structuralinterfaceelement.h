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

#ifndef structuralinterfaceelement_h
#define structuralinterfaceelement_h

#include "element.h"
#include "femcmpnn.h"
#include "domain.h"
#include "floatmatrix.h"
#include "loadtimefunction.h"
#include "matresponsemode.h"
#include "valuemodetype.h"
#include "integrationdomain.h"
#include "dofmantransftype.h"
#include "structuralinterfacecrosssection.h"

namespace oofem {

class TimeStep;
class Node;
class StructuralInterfaceMaterial;
class GaussPoint;
class FloatArray;
class IntArray;
class FEInterpolation;


/**
 * Abstract base class for all structural interface elements. It declares a common interface available
 * to all derived elements. The implementation of these services is partly left to the derived classes,
 * some general services are implemented here (but they can be overload by more efficient
 * element implementations).
 * The general implementation provided here is intended for both linear and nonlinear computations.
 *
 * Tasks:
 * - Obtaining its basic data, the element reads in the data file the number
 *   of these objects, then obtains the data from the domain (methods
 *   'giveNode', 'giveMaterial',etc).
 * - Calculating its contribution to the problem :
 *   Calculating its mass matrix M, its stiffness matrix K, its load vector
 *   f, its location array.
 *   Calculating its contribution to the LHS and RHS of the linear system,
 *   using Static,Newmark,etc, formula. These contributions are usually
 *   combinations of M,K,f.
 * - Performing end-of-step operations;
 *   Calculating the strains and stresses at its Gauss points.
 * - Printing its output in the data file and updating itself.
 */
class StructuralInterfaceElement : public Element
{
protected:
    /// Initial displacement vector, describes the initial nodal displacements when element has been casted.
    FloatArray *initialDisplacements;
    FEInterpolation *interpolation;
    virtual FEInterpolation *giveInterpolation() const { return interpolation; };

public:
    /**
     * Constructor. Creates structural element with given number, belonging to given domain.
     * @param n Element number.
     * @param d Domain to which new material will belong.
     */
    StructuralInterfaceElement(int n, Domain *d);
    /// Destructor.
    virtual ~StructuralInterfaceElement();

    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType, TimeStep *tStep);
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep);

    virtual void giveDefaultDofManDofIDMask(int inode, IntArray &answer) const { this->giveDofManDofIDMask(inode, EID_MomentumBalance, answer); }
    virtual void giveDefaultInternalDofManDofIDMask(int inode, IntArray &answer) const { this->giveInternalDofManDofIDMask(inode, EID_MomentumBalance, answer); }

    /**
     * Computes the stiffness/tangent matrix of receiver. Default implementation computes element stiffness using
     * @f$ K=\int_{\Gamma} N^{\mathrm{T}} D N \mathrm{d}V @f$ formulae, where @f$ N @f$ is the element geometric matrix such 
     * that @f$ j = N a @f$ and @f$ D @f$ is the stiffness matrix of the interface material.
     * Numerical integration procedure uses integrationRulesArray for numerical integration. 
     *
     * The geometrical matrix is obtained using computeNmatrixAt service and the constitutive matrix is obtained using
     * computeConstitutiveMatrixAt service.
     *
     * For higher numerical performance, only one half of stiffness matrix is computed and answer is then symmetrized.
     * Therefore, if element matrix will be generally nonsymmetric, one must specialize this method.
     * Finally, the result is transformed into global coordinate system (or nodal coordinate system, if it is defined).
     *
     * @param answer Computed stiffness matrix (symmetric).
     * @param rMode Response mode.
     * @param tStep Time step.
     */
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    // stress equivalent vector in nodes (vector of internal forces)
    // - mainly for nonLinear Analysis.
    /**
     * Returns equivalent nodal forces vectors. Useful for nonlinear analysis.
     * Default implementation computes result as @f$ F=\int_v B^{\mathrm{T}} \sigma \mathrm{d}V @f$, where @f$ \sigma @f$ is the
     * real element stress vector obtained using computeStressVector service (if useUpdatedGpRecord=0) or
     * (if useUpdatedGpRecord=1) from integration point status.
     * The geometric matrix is obtained using computeBmatrixAt service.
     * Integration is performed using default integration rule, which should produce always valid results,
     * assuming that strains used for computation of stresses are valid.
     * @param answer Internal nodal forces vector.
     * @param tStep Time step.
     * @param useUpdatedGpRecord If equal to zero, the stresses in integration points are computed (slow but safe), else if
     * nonzero the stresses are taken directly from integration point status (should be derived from StructuralMaterialStatus)
     * (fast, but engineering model must ensure valid status data in each integration point).
     */
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);
    virtual void computeTraction(FloatArray &traction, IntegrationPoint *ip, FloatArray &jump, TimeStep *tStep);
    virtual void computeSpatialJump(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);


    //@}

    // Overloaded methods.
    virtual void updateInternalState(TimeStep *tStep);
    virtual void updateYourself(TimeStep *tStep);
    virtual int checkConsistency();
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual const char *giveClassName() const { return "StructuralInterfaceElement"; };
    //virtual classType giveClassID() const { return StructuralInterfaceElementClass; }

    StructuralInterfaceCrossSection *giveInterfaceCrossSection();
    virtual double computeAreaAround(GaussPoint *gp) = 0;

protected:



    /**
     * Computes constitutive matrix of receiver. Default implementation uses element cross section
     * giveCharMaterialStiffnessMatrix service.
     * @param answer Constitutive matrix.
     * @param rMode Material response mode of answer.
     * @param gp Integration point for which constitutive matrix is computed.
     * @param tStep Time step.
     */
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer,
                                             MatResponseMode rMode, GaussPoint *gp,
                                             TimeStep *tStep);


    /**
     * Computes volume related to integration point on local surface.
     * @param gp Surface integration point.
     * @param iSurf Surface number.
     */
    virtual double computeSurfaceVolumeAround(GaussPoint *gp, int iSurf) { return 0.; }

    /**
     * Computes global coordinates of integration point on  local surface.
     * @param answer Global coordinates.
     * @param gp Surface integration point.
     * @param iSurf Surface number.
     */
    virtual void computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iSurf) { answer.resize(0); }

    // Global to local element c.s transformation for load vector dofs
    /**
     * Returns transformation matrix from global coordinate system to local
     * element coordinate system for element load vector components.
     * If no transformation is necessary, answer is empty matrix (default);
     * @param answer Transformation matrix.
     * @return Nonzero if transformation matrix is not empty matrix, zero otherwise.
     */
    virtual int computeLoadGToLRotationMtrx(FloatMatrix &answer) {
        answer.beEmptyMtrx();
        return 0;
    }

    // Local edge (LE-local Edge c.s) or surface (LS-local surface c.s) c.s
    // to element local c.s for load vector dofs
    /**
     * Returns transformation matrix from local edge c.s  to element local coordinate system
     * of load vector components. Necessary, because integration must be done in local coordinate
     * system of entity (edge or surface).
     * If no transformation is necessary, answer is empty matrix (default);
     * @param answer Computed rotation matrix.
     * @param iEdge Edge number.
     * @param gp Integration point (point, where transformation is computed, useful for curved edges).
     * @return Nonzero if transformation matrix is not empty matrix, zero otherwise.
     */
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp) {
        answer.beEmptyMtrx();
        return 0;
    }
    /**
     * Returns transformation matrix from local surface c.s  to element local coordinate system
     * of load vector components. Necessary, because integration must be done in local coordinate
     * system of entity (edge or surface).
     * If no transformation is necessary, answer is empty matrix (default);
     * @param answer Computed rotation matrix.
     * @param iSurf Surface number.
     * @param gp Integration point (point, where transformation is computed, useful for curved surfaces)
     * @return Nonzero if transformation matrix is not empty matrix, zero otherwise.
     */
    virtual int computeLoadLSToLRotationMatrix(FloatMatrix &answer, int iSurf, GaussPoint *gp) {
        answer.beEmptyMtrx();
        return 0;
    }
    //@}


    /**
     * Computes the stress vector of receiver at given integration point, at time step stepN.
     * The nature of these stresses depends on the element's type.
     * @param answer Stress vector.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    //virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);

    /**
     * Computes the geometrical matrix of receiver in given integration point.
     * The product of this matrix (assembled at given integration point) and element displacement
     * vector is element strain vector. If lowerIndx and upperIndx parameters are specified,
     * answer is formed only for strains within this interval. This will affects the size of answer.
     *
     * @param gp Integration point for which answer is computed.
     * @param answer Geometric matrix of receiver.

     */
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer){};

    /**
     * Computes modified interpolation matrix (N) for the element which multiplied 
     * with the unknowns vector (u) produces the spatial jump.
     * @param gp Integration point for which answer is assembled.
     * @param answer Interpolation matrix evaluated at gp.
     */
    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer) = 0;

    virtual void computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer) = 0;


    /**
     * Returns maximum approximation order used by receiver.
     * Must be implemented by derived classes
     * @return Order of approximation.
     */
    virtual int giveApproxOrder() { return 0; }
    /**
     * Return desired number of integration points for consistent mass matrix
     * computation, if required.
     * @return Number of integration points for mass matrix.
     */
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 0; }

    virtual int testCrossSectionExtension(CrossSectExtension ext) { return ( ( ext == CS_StructuralInterfaceCapability ) ? 1 : 0 ); }
};
} // end namespace oofem
#endif // structuralinterfaceelement_h
