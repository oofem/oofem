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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#ifndef structuralelementevaluator2_h
#define structuralelementevaluator2_h

#include "elementgeometry.h"
#include "elementevaluator.h"

#include "femcmpnn.h"
#include "domain.h"
#include "floatmatrix.h"
#include "matresponsemode.h"
#include "valuemodetype.h"
#include "integrationdomain.h"
#include "dofmantransftype.h"

///@name Input fields for NLStructuralElement
//@{
#define _IFT_StructuralElementEvaluator2_nlgeoflag "nlgeo"
//@}


namespace oofem {
#define ALL_STRAINS -1

class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatArray;
class IntArray;
class SparseMtrx; // required by addNonlocalStiffnessContributions declaration
class StructuralCrossSection;
class IDNLMaterial;

/**
 * Abstract base class for all "structural" finite elements. It declares common interface provided
 * by all derived elements. The implementation of these services is partly left on derived classes,
 * some general services are implemented here generally (But they can be overload by more efficient
 * element implementation).
 * The general implementation provided here is intended for both linear and nonlinear computations.
 * Both material and geometrical nonlinearities are taken into account. 
 *
 * The basic data of an element are the numbers of its 'numberOfNodes' nodes,
 * stored in 'nodeArray', of its 'material', of its body loads (eg, the dead
 * weight) stored in 'loadArray'. These data are obtained from the domain.
 * The element contains logical reference to Volume object (
 * which stores material,geometrical characteristics, integration points, state of
 * material and so on)
 * The calculated data of an element are its 'massMatrix', its 'stiffnessMatrix'
 * its 'locationArray'. Since the load vector is recalculated at every
 * time step, it is not given the status of attribute.
 *
 * The activation of non-linear effects can be generally controlled on element level
 * using nlGeometry, allowing elements to be used also in linear computations.
 * nlGeometry = 0 - Small strain theory based on the engineering stress and strain
 * nlGeometry = 1 - Finite deformation theory formulated in terms of the deformation gradient F
 *                  and first Piola-Kirchoff stress P in the virtual work
 * @note
 * Methods ComputeStressVector and computeStrainVector as they are implemented here
 * assume Total Lagrange approach.
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
class StructuralElementEvaluator2 : public ElementEvaluator
{
protected:
    /// Initial displacement vector, describes the initial nodal displacements when element has been casted.
    FloatArray *initialDisplacements;
	/// Displacement interpolation number describes the order of displacement interpolation, for pure displacement models din ís equal to one
    int displacementInterpolationNumber;
	 /// Flag indicating if geometrical nonlinearities apply.
    int nlGeometry;


public:
    /**
     * Constructor. Creates structural element evaluator
     * @param din Displacement interpolation number.
     */
    StructuralElementEvaluator2(int din);
	/**
     * Constructor. Creates structural element evaluator with default din equals to one
     */
    StructuralElementEvaluator2();
    /// Destructor.
    virtual ~StructuralElementEvaluator2();

	virtual void giveCharacteristicMatrix(FloatMatrix & answer, CharType, TimeStep * tStep, ElementGeometry* elemeGeometry);
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep, ElementGeometry* elemeGeometry);

//	virtual void giveDefaultDofManDofIDMask(int inode, IntArray &answer, ElementGeometry* elemGeometry) const { elemGeometry->giveDofManDofIDMask(inode, EID_MomentumBalance, answer); }
 //   virtual void giveDefaultInternalDofManDofIDMask(int inode, IntArray &answer, ElementGeometry* elemGeometry) const { elemGeometry->giveInternalDofManDofIDMask(inode, EID_MomentumBalance, answer); }
    /**
     * Computes mass matrix of receiver. Default implementation returns consistent mass matrix and uses
     * numerical integration. Returns result of this->computeConsistentMassMatrix service, transformed into
     * nodal coordinate system.
     * Requires the computeNmatrixAt and giveMassMtrxIntegrationgMask services to be implemented.
     * @param answer Mass matrix.
     * @param tStep Time step.
     */
    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep, ElementGeometry* elemeGeometry);
    /**
     * Computes lumped mass matrix of receiver. Default implementation returns lumped consistent mass matrix.
     * Then returns lumped mass transformed into nodal coordinate system.
     * The lumping procedure zeroes all off-diagonal members and zeroes also all diagonal members
     * corresponding to non-displacement DOFs. Such diagonal matrix is then rescaled, to preserve
     * the element mass.
     * Requires the computeNmatrixAt and giveMassMtrxIntegrationgMask services to be implemented.
     * @param answer Lumped mass matrix.
     * @param tStep Time step.
     */
    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep, ElementGeometry* elemeGeometry);
    /**
     * Computes consistent mass matrix of receiver using numerical integration over element volume.
     * Mass matrix is computed as @f$ M = \int_V N^{\mathrm{T}} \rho N dV @f$, where @f$ N @f$ is displacement approximation matrix.
     * The number of necessary integration points  is determined using this->giveNumberOfIPForMassMtrxIntegration
     * service. Only selected degrees of freedom participate in integration of mass matrix. This is described
     * using dof mass integration mask. This mask is obtained from this->giveMassMtrxIntegrationgMask service.
     * The nonzero mask value at i-th position indicates that i-th element DOF participates in mass matrix
     * computation. The result is in element local coordinate system.
     * @param answer Mass matrix.
     * @param tStep Time step.
     * @param mass Total mass of receiver.
     */
    virtual void computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, ElementGeometry* elemeGeometry, const double *ipDensity = NULL);
    /**
     * Returns mask indicating, which unknowns (their type and ordering is the same as
     * element unknown vector) participate in mass matrix integration.
     * Nonzero value at i-th position
     * indicates that corresponding row in interpolation matrix N will participate in
     * mass matrix integration (typically only displacements are taken into account).
     * @param answer Integration mask, if zero sized, all unknowns participate. This is default.
     */
    virtual void giveMassMtrxIntegrationgMask(IntArray &answer) { answer.resize(0); }
    /**
     * Computes numerically stiffness matrix of receiver. Default implementation computes element stiffness using
     * @f$ K=\int_v B^{\mathrm{T}} D B \mathrm{d}V @f$ formulae, where @f$ B @f$ is element geometric matrix and @f$ D @f$ is material stiffness matrix.
     * No geometrical nonlinearity is taken into account. NUmerical integration procedure uses integrationRulesArray
     * for numerical integration. Support for reduced or selected integration is implemented. The individual integration
     * rules are assumed to correspond to different terms from which the overall matrix is assembled.
     *
     * If numberOfIntegrationRules is equal to
     * 1, the full integration of all coefficients is performed. Otherwise, integration is performed using following rules.
     * Each integration rule can specify start and end strain index of strain vector components for which is valid.
     * It is necessary to ensure that these start and end indexes, dividing geometrical matrix into blocks,
     * are not overlapping and that each strain component is included.
     *
     * Then stiffness matrix is obtained as summation of integrals @f$ I_{ij}=\int_v B^{\mathrm{T}}_i D_{ij} B_j \mathrm{d}V @f$
     * where @f$ B_i @f$ is i-th block of geometrical matrix and @f$ D_{ij} @f$ is corresponding constitutive sub-matrix.
     * The geometrical matrix is obtained using computeBmatrixAt service and the constitutive matrix is obtained using
     * computeConstitutiveMatrixAt service.
     * The @f$ I_{ij} @f$ integral is evaluated using such integration rule, which is valid for i-th or j-th block
     * and has smaller number of integration points.
     *
     * For higher numerical performance, only one half of stiffness matrix is computed and answer is then symmetrized.
     * Therefore, if element matrix will be generally nonsymmetric, one must specialize this method.
     * Finally, the result is transformed into global coordinate system (or nodal coordinate system, if it is defined).
     *
     * @param answer Computed stiffness matrix (symmetric).
     * @param rMode Response mode.
     * @param tStep Time step.
     */
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep, ElementGeometry* elemeGeometry);
    /**
     * @see giveStiffnessMatrix
     */
    void computeStiffnessMatrix_withIRulesAsSubcells(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep, ElementGeometry* elemeGeometry);


	
    virtual void computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode, ElementGeometry* elemeGeometry);
    virtual void computeLocalForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode, ElementGeometry* elemeGeometry);
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
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, ElementGeometry* elemeGeometry, int useUpdatedGpRecord = 0);

    /**
     * @see giveInternalForcesVector
     */
    virtual void giveInternalForcesVector_withIRulesAsSubcells(FloatArray &answer,
                                                               TimeStep *tStep, ElementGeometry* elemeGeometry, int useUpdatedGpRecord = 0);

    /**
     * Compute strain vector of receiver evaluated at given integration point at given time
     * step from element displacement vector.
     * The nature of strain vector depends on the element type.
     * @param answer Requested strain vector.
     * @param gp Integration point where to calculate the strain.
     * @param tStep Time step.
     */
    virtual void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ElementGeometry* elemeGeometry);

	

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep, ElementGeometry* elemeGeometry);
    /**
     * Computes at given time (tStep) the the resulting temperature component array.
     * This is summation of all temperature load  components of  receiver.
     * @param answer Resulting temperature components of receiver.
     * @param tStep Time step.
     * @param gp Integration point.
     * @param mode Determines response mode.
     */
    virtual void computeResultingIPTemperatureAt(FloatArray &answer, TimeStep *tStep, GaussPoint *gp, ValueModeType mode, ElementGeometry* elemeGeometry);
    /**
     * Computes at given time the resulting eigenstrain component array.
     * This is summation of all eigenstrains imposed on the element.
     * @param answer Resulting eigenstrain components of receiver.
     * @param tStep Time step.
     * @param gp Integration point.
     * @param mode Determines response mode.
     */
    virtual void computeResultingIPEigenstrainAt(FloatArray &answer, TimeStep *tStep, GaussPoint *gp, ValueModeType mode, ElementGeometry* elemeGeometry);

    /**@name Methods related to nonlocal models */
    //@{
    /**
     * Updates internal element state (in all integration points of receiver)
     * before nonlocal averaging takes place. Used by so nonlocal materials,
     * because their response in particular point depends not only on state in this point, but
     * depends also on state in point's neighborhood. Nonlocal quantity is computed as nonlocal
     * average of local quantities. Therefore, before updating integration point state depending on
     * nonlocal quantity (or quantities), local quantities in all integration points must be updated
     * in advance. This function updates local quantities of material model using
     * updateBeforeNonlocalAverage member function of structural nonlocal material class.
     * @param tStep Time step.
     */
    virtual void updateBeforeNonlocalAverage(TimeStep *tStep, ElementGeometry* elemGeometry);
    /**
     * Returns the "nonlocal" location array of receiver. This is necessary, when stiffness matrix
     * of nonlocal model is assembled. Since model is nonlocal, the value at given IP depends on
     * other IP (generally belonging to different elements) and as a consequence leads to
     * increase of stiffness matrix profile, to take into account this "remote" dependency.
     * @param locationArray Location arrays from neighboring elements.
     * @param us Unknown numbering scheme.
     */
    virtual void giveNonlocalLocationArray(IntArray &locationArray, const UnknownNumberingScheme &us, ElementGeometry* elemGeometry);
    /**
     * Adds the "nonlocal" contribution to stiffness matrix, to account for nonlocality of
     * material model. Typically, this contribution is obtained by summing up mutual IP contributions.
     */
    virtual void addNonlocalStiffnessContributions(SparseMtrx &dest, const UnknownNumberingScheme &s, TimeStep *tStep, ElementGeometry* elemeGeometry);
    //@}

	/**@name Methods related to geometrical nonlinearity */
    //@{
	/**
     * Returns the geometry mode describing the formulation used in the internal work
     * 0 - Engineering (small deformation) stress-strain mode
     * 1 - First Piola-Kirchhoff - Deformation gradient mode, P is defined as FS
     */
    int giveGeometryMode() { return nlGeometry; } 
	
	/**
     * Computes the deformation gradient in Voigt form at integration point ip and at time
     * step tStep. Computes the displacement gradient and adds an identitiy tensor.
     *
     * @param answer Deformation gradient vector
     * @param gp Gauss point.
     * @param tStep Time step.
     */
    void computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ElementGeometry* elemeGeometry);
	
    /**
     * Computes initial stress matrix for linear stability problem.
     * Default implementation is not provided.
     * Please note, that initial stress matrix depends on normal forces of element,
     * corresponding engineering model must take this into account.
     * @param answer Computed initial stress matrix.
     * @param tStep Time step.
     */
    virtual void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep, ElementGeometry* elemeGeometry)
    {
       // _error("computeInitialStressMatrix: not implemented");
    }
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
    void computeFirstPKStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ElementGeometry* elemeGeometry);
	
    /**
     * Computes the Cauchy stress tensor on Voigt format. This method will
     * be called if nlGeo = 1 and mode = UL. This method computes the deformation gradient F and passes
     * it on to the crossection which then asks for the stress from the material.
     *
     * @param answer Computed stress vector in Voigt form.
     * @param gp Gauss point at which the stress is evaluated.
     * @param tStep Time step.
     */
    void computeCauchyStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ElementGeometry* elemeGeometry);
	//@}
  

    // Overloaded methods.
    virtual int adaptiveUpdate(TimeStep *tStep, ElementGeometry* elemGeometry);
    virtual void updateInternalState(TimeStep *tStep, ElementGeometry* elemGeometry);
    virtual void updateYourself(TimeStep *tStep, ElementGeometry* elemGeometry);
	virtual int checkConsistency(ElementGeometry* elemGeometry);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual const char *giveClassName() const { return "StructuralElementEvaluator"; }

#ifdef __OOFEG
    /**
     * Returns internal state variable (like stress,strain) at node of element in Reduced form,
     * the way how is obtained is dependent on InternalValueType.
     * The value may be local, or smoothed using some recovery technique /
     * returns zero if element is unable to respond to request.
     * @param answer Contains result, zero sized if not supported.
     * @param type Determines the internal variable requested (physical meaning).
     * @param mode Determines the mode of variable (recovered, local, ...).
     * @param node Node number, for which variable is required.
     * @param tStep Time step.
     * @return Nonzero if o.k, zero otherwise.
     */
    virtual int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                        int node, TimeStep *tStep, ElementGeometry* elemeGeometry);
    /// Shows sparse structure
    void showSparseMtrxStructure(CharType mtrx, oofegGraphicContext &gc, TimeStep *tStep, ElementGeometry* elemeGeometry);
    /// Shows extended sparse structure (for example, due to nonlocal interactions for tangent stiffness)
    virtual void showExtendedSparseMtrxStructure(CharType mtrx, oofegGraphicContext &gc, TimeStep *tStep, ElementGeometry* elemeGeometry);

#endif

    // Interface for b.c.s applied by Sets:
    virtual void computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep, ElementGeometry* elemeGeometry);
    virtual void computeBoundaryLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, ElementGeometry* elemeGeometry);

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
                                             TimeStep *tStep, ElementGeometry* elemeGeometry);
    /// Helper function which returns the structural cross-section for the element.
	StructuralCrossSection *giveStructuralCrossSection(ElementGeometry* elemGeometry);

    virtual void createMaterialStatus();

protected:


    /**
     * Computes the load vector due to body load acting on receiver, at given time step.
     * Default implementation computes body load vector numerically as @f$ l=\int_V N^{\mathrm{T}} f \rho\;\mathrm{d}V @f$
     * using default integration rule. Result is transformed to global c.s.
     * @param answer Computed load vector due to body load
     * @param load Body load which contribution is computed.
     * @param tStep Time step.
     * @param mode determines the response mode
     */
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode, ElementGeometry* elemeGeometry);

    /**@name Edge and surface load support services.
     */
    //@{
    /**
     * Computes point load vector contribution of receiver for given load (should has BoundaryLoad Base).
     * @param answer Computed load vector.
     * @param load Point load which contribution is computed.
     * @param tStep Time step.
     * @param mode determines response mode
     */
    virtual void computePointLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode, ElementGeometry* elemeGeometry);
    /**
     * Computes edge load vector contribution of receiver for given load (should has BoundaryLoad Base).
     * Each edge should have unique number assigned to identify it.
     * The default implementation does integration of load vector in local edge space
     * (i.e. one dimensional integration is performed on line). This general implementation requires
     * that element must provide following services:
     * - computeEgdeNMatrixAt - returns interpolation matrix of local edge DOFs in the local edge space.
     * - computeEdgeVolumeAround - returns volume corresponding to integration point at local edge.
     * - giveEdgeDofMapping - returns integer array specifying local dof edge mapping to "global" element dofs.
     *
     * Integration rule is set up automatically, based on element interpolation order and load approximation.
     * Integration points are set-up using standard integration rule services (setUpIntegrationPoints method).
     * Gauss integration rule is used.
     * If derived class overrides this default implementation somehow, the above services
     * must not be implemented.
     *
     * @param answer Computed load vector.
     * @param load Edge load which contribution is computed.
     * @param iEdge Edge number where load applies.
     * @param tStep Time step.
     * @param mode Determines response mode.
     */
    virtual void computeEdgeLoadVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep, ValueModeType mode, ElementGeometry* elemeGeometry);
    /**
     * Computes surface load vector contribution of receiver for given load (should has BoundaryLoad Base).
     * Each surface should have unique number assigned to identify it.
     * The default implementation does integration of load vector in local surface space
     * (i.e. two dimensional integration is performed on triangle or square).
     * This general implementation requires
     * that element must provide following services:
     * - getSurfaceIntegrationRule - returns integration rule for surface for given polynomial order.
     * - computeSurfaceNMatrixAt - returns interpolation matrix of local surface DOFs in the local edge space.
     * - computeSurfaceVolumeAround - returns volume corresponding to integration point of local surface.
     * - giveSurfaceDofMapping - returns integer array specifying local dof surface mapping to "global" element dofs.
     *
     * Integration rule is set up automatically, based on element interpolation order and load approximation.
     * Integration points are set-up using standard integration rule services (setUpIntegrationPoints method).
     * Gauss integration rule is used.
     * If derived class overrides this default implementation somehow, the above services
     * must not be implemented.
     *
     * @param answer Computed load vector.
     * @param load Surface load which contribution is computed.
     * @param iSurf Surface number where load applies.
     * @param tStep Time step.
     * @param mode Determines response mode.
     */
    virtual void computeSurfaceLoadVectorAt(FloatArray &answer, Load *load, int iSurf, TimeStep *tStep, ValueModeType mode, ElementGeometry* elemeGeometry);
    /**
     * Computes Edge interpolation matrix. Interpolation matrix provide way, how to compute
     * local edge unknowns (nonzero element unknowns on edge) at any integration point of edge, based on
     * local edge unknowns in edge nodes.
     * The edge numbering and local edge coordinate system is element dependent.
     * The integration point is specified using one-dimensional iso coordinates.
     * @param answer Interpolation matrix of edge.
     * @param iedge Edge number.
     * @param gp Integration point.
     * @todo Should be a FlotArray instead
     */
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp, ElementGeometry* elemeGeometry) { answer.resize(0, 0); }
    /**
     * Computes surface interpolation matrix. Interpolation matrix provide way, how to compute
     * local surface unknowns (nonzero element unknowns on surface) at any integration point of surface, based on
     * local unknowns in surface nodes.
     * Local coordinate system of surface edge and element surface numbering is element dependent.
     * The integration point is specified using two-dimensional iso coordinates, or using area coordinates
     * for triangular surface.
     * @param answer Interpolation matrix of surface.
     * @param iSurf Surface number.
     * @param gp Integration point.
     */
    virtual void computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *gp, ElementGeometry* elemeGeometry) { answer.resize(0, 0); }
    /**
     * Assembles edge dof mapping mask, which provides mapping between edge local DOFs and "global" element
     * DOFs. Mask can be imagined as local edge code numbers used to localize local edge DOFs to
     * element DOFs.
     * @param answer Edge DOF mapping.
     * @param iEdge Edge number.
     */
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge, ElementGeometry* elemeGeometry) const { answer.resize(0); }
    /**
     * Assembles surface dof mapping mask, which provides mapping between surface local DOFs and "global" element
     * DOFs. Mask can be imagined as local surface code numbers used to localize local DOFs to
     * element DOFs.
     * @param answer Surface DOF mapping.
     * @param iSurf Surface number
     */
    virtual void giveSurfaceDofMapping(IntArray &answer, int iSurf, ElementGeometry* elemeGeometry) const { answer.resize(0); }
    /**
     * Returns integration rule for integration over element surface.
     * @param i order of integrated polynomial
     * @return Best integration rule to integrate polynomial of order i over element surface.
     */
    virtual IntegrationRule *GetSurfaceIntegrationRule(int i, ElementGeometry* elemeGeometry) { return NULL; }
    /**
     * Computes volume related to integration point on local edge.
     * @param gp edge integration point
     * @param iEdge edge number
     */
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge, ElementGeometry* elemeGeometry) { return 0.; }
    /**
     * Computes volume related to integration point on local surface.
     * @param gp Surface integration point.
     * @param iSurf Surface number.
     */
    virtual double computeSurfaceVolumeAround(GaussPoint *gp, int iSurf, ElementGeometry* elemeGeometry) { return 0.; }
    /**
     * Computes global coordinates of integration point on local edge.
     * @param answer Global coordinates.
     * @param gp Edge integration point.
     * @param iEdge Edge number.
     */
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge, ElementGeometry* elemeGeometry) { answer.resize(0); }
    /**
     * Computes global coordinates of integration point on  local surface.
     * @param answer Global coordinates.
     * @param gp Surface integration point.
     * @param iSurf Surface number.
     */
    virtual void computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iSurf, ElementGeometry* elemeGeometry) { answer.resize(0); }

    // Global to local element c.s transformation for load vector dofs
    /**
     * Returns transformation matrix from global coordinate system to local
     * element coordinate system for element load vector components.
     * If no transformation is necessary, answer is empty matrix (default);
     * @param answer Transformation matrix.
     * @return Nonzero if transformation matrix is not empty matrix, zero otherwise.
     */
    virtual int computeLoadGToLRotationMtrx(FloatMatrix &answer) {
		answer.clear();
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
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp, ElementGeometry* elemeGeometry) {
        answer.clear();
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
    virtual int computeLoadLSToLRotationMatrix(FloatMatrix &answer, int iSurf, GaussPoint *gp, ElementGeometry* elemeGeometry) {
        answer.clear();
        return 0;
    }
    //@}

    /**
     * Computes the stress vector of receiver at given integration point, at time step tStep.
     * The nature of these stresses depends on the element's type.
     * @param answer Stress vector.
     * @param strain Strain vector.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep, ElementGeometry* elemeGeometry);

	
    /**
     * Computes the geometrical matrix of receiver in given integration point.
     * The product of this matrix (assembled at given integration point) and element displacement
     * vector is element strain vector. If lowerIndx and upperIndx parameters are specified,
     * answer is formed only for strains within this interval. This will affects the size of answer.
     *
     * @param gp Integration point for which answer is computed.
     * @param answer Geometric matrix of receiver.
     * @param lowerIndx If specified, answer is formed only for strain with index equal and greater than  lowerIndx.
     * This parameter has default value 1 (answer is formed from first strain).
     * @param upperIndx If specified, answer is formed only for strain with index less and equal than  upperIndx.
     * This parameter has default value ALL_STRAINS (answer is formed for all strains).
     */
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemeGeometry,
                                  int lowerIndx = 1, int upperIndx = ALL_STRAINS) = 0;

	/**
     * Computes a matrix which, multiplied by the column matrix of nodal displacements,
     * gives the displacement gradient stored by columns.
     * The components of this matrix are derivatives of the shape functions,
     * but they are arranged in a somewhat different way from the usual B matrix.
     * @param gp Integration point.
     * @param answer BF matrix at this point.
     */

    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemeGeometry) {
        OOFEM_SIMPLE_ERROR("NLStructuralElement::computeBHmatrixAt : method not implemented for this element");
        return;
    }

public:
    /**
     * Computes interpolation matrix for element unknowns.
     * The order and meaning of unknowns is element dependent.
     * @param iLocCoord Local coordinates.
     * @param answer Interpolation matrix evaluated at gp.
     */
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer, ElementGeometry* elemeGeometry);

protected:
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

    /**
     * General service for condensation of stiffness and optionally load vector and mass or initial stress matrices
     * of receiver.
     * @param stiff Stiffness matrix to be condensed. Must be specified.
     * @param mass Mass or initial stress matrix. If parameter is NULL, only stiffness and/or load is condensed.
     * @param load Load vector of receiver. If specified then it is condensed. If no load vector condensation is necessary
     * set parameter to NULL pointer.
     * @param what Integer array. If at i-th position is nonzero, then i-th component is condensed.
     */
    void condense(FloatMatrix *stiff, FloatMatrix *mass, FloatArray *load, IntArray *what);

    friend class IDNLMaterial;
    friend class TrabBoneNL3D;
    friend class MisesMatNl;
    friend class RankineMatNl;
    friend class GradDpElement;
};


class OneDimensionalStructuralElementEvaluator : public StructuralElementEvaluator2
{
	virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemeGeometry,
                                  int lowerIndx = 1, int upperIndx = ALL_STRAINS);

	virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemeGeometry);

	virtual void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const 
	{
		answer.setValues(1, D_u);    
	}
};

class PlaneStressElementEvaluator2 : public StructuralElementEvaluator2
{
	virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemeGeometry,
                                  int lowerIndx = 1, int upperIndx = ALL_STRAINS);

	virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemeGeometry);

	virtual void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const 
	{
		answer.setValues(2, D_u, D_v);    
	}
};

class PlaneStrainElementEvaluator : public StructuralElementEvaluator2
{
	virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemeGeometry,
                                  int lowerIndx = 1, int upperIndx = ALL_STRAINS);

	virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemeGeometry);

	virtual void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const 
	{
		answer.setValues(2, D_u, D_v);    
	}
};

class AxisymmetricStructuralElementEvaluator : public StructuralElementEvaluator2
{
	virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemeGeometry,
                                  int lowerIndx, int upperIndx);

	virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemeGeometry);

	virtual void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const 
	{
		answer.setValues(2, D_u, D_v);    
	}
};



class ThreeDimensionalStructuralElementEvaluator : public StructuralElementEvaluator2
{
	virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemeGeometry,
                                  int lowerIndx = 1, int upperIndx = ALL_STRAINS);

	virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, ElementGeometry* elemeGeometry);

	virtual void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const 
  {
	  answer.setValues(3, D_u, D_v,D_w);
  }
};









} // end namespace oofem
#endif // structuralelementevaluator2_h
