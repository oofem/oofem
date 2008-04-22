/* $Header: /home/cvs/bp/oofem/oofemlib/src/structuralelement.h,v 1.24 2003/04/06 14:08:26 bp Exp $ */
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


//   ********************************
//   *** CLASS STRUCTURAL ELEMENT ***
//   ********************************


#ifndef structuralelement_h
#define structuralelement_h

#include "element.h"
#include "femcmpnn.h"
#include "domain.h"
#include "flotmtrx.h"
#include "cltypes.h"

#define ALL_STRAINS -1

class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;
class SparseMtrx; // required by addNonlocalStiffnessContributions declaration
class IDNLMaterial;
/**
 * Abstract base class for all "structural" finite elements. It declares common interface provided
 * by all derived elements. The implementation of these services is partly left on derived classes,
 * some general services are implemented here generally (But they can be overload by more efficient
 * element implementation).
 * The general implementation provided here is intended for both linear and nonlinear computations.
 * At this level, only material model nonlinearities are taken into account. If particular element type
 * will participate in geometrically - nonlinear computation, it should be derived from derived
 \ref NLStructuralElement class, which provide support for this cases.
 */
class StructuralElement : public Element
{
    /*
     * This abstract class is the most important class of the program. It is the
     * superclass of all classes implementing structural finite elements
     * (bar, shell, etc).
     * An element is an attribute of a domain.
     * DESCRIPTION :
     * The basic data of an element are the numbers of its 'numberOfNodes' nodes,
     * stored in 'nodeArray', of its 'material', of its body loads (eg, the dead
     * weight) stored in 'loadArray'. These data are obtained from the domain.
     * The element contains logical reference to Volume object (
     * which stores material,geometrical characteristics, gaussPoints, state of
     * material and so on)
     * The calculated data of an element are its 'massMatrix', its 'stiffnessMa-
     * trix', its 'locationArray'. Since the load vector is recalculated at every
     * time step, it is not given the status of attribute.
     * TASKS :
     * -defining itself :
     *   .typing itself (methods 'typed' and 'ofType'). When the domain creates
     *    an element, it actually creates a temporary instance of Element, then
     *    asks this element to transform itself into an element of the right type
     *    (PlaneStrain, Truss2D, etc) ;
     *   .obtaining its basic data : the element reads in the data file the num-
     *    ber of these objects, then obtains the data from the domain (methods
     *    'giveNode', 'giveMaterial',etc) ;
     * -calculating its contribution to the problem :
     *   .calculating its mass matrix M, its stiffness matrix K, its load vector
     *    f, its location array ;
     *   .calculating its contribution to the LHS and RHS of the linear system,
     *    using Static,Newmark,etc, formula. These contributions are usually
     *    combinations of M,K,f.
     * -performing end-of-step operations :
     *   .calculating the strains and stresses at its Gauss points ;
     *   .printing its output in the data file and updating itself ;
     */

protected:
    /// Cached transformation matrix of receiver
    FloatMatrix *rotationMatrix;
    /// Flag indicating if tranformation matrix has been already computed
    int rotationMatrixDefined;
    //     FloatMatrix* constitutiveMatrix;

public:
    /**
     * Constructor. Creates structural element with given number, belonging to given domain.
     * @param n element number
     * @param d domain to which new material will belong
     */
    StructuralElement(int, Domain *);           // constructors
    /// Destructor.
    ~StructuralElement();                        // destructor

    // characteristic  matrix
    /**
     * Computes characteristic matrix of receiver of requested type in given Timestep.
     * @param answer requested characteristic matrix
     * @param mtrx   id  of characteristic component requested.
     * @param tStep  time step when answer is computed.
     * @see CharType type.
     * @return If element has no capability to compute requsted type of caharacteristic matrix
     * error function is invoked.
     */
    void giveCharacteristicMatrix(FloatMatrix & answer, CharType, TimeStep *);
    /**
     * Computes characteristic vector of receiver of requested type in given Timestep.
     * @param answer requested characteristic vector
     * @param type   id  of characteristic component requested.
     * @param mode   determines mode of answer.
     * @param tStep  time step when answer is computed.
     * @see CharType.
     * @return If element has no capability to compute requsted type of caharacteristic vector
     * error function is invoked.
     */
    void  giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *);


    /**
     * Computes mass matrix of receiver. Default implementation returns consistent mass matrix and uses
     * numerical integration. Returns result of this->computeConsistentMassMatrix service, transformed into
     * nodal coordinate system.
     * Requires the computeNmatrixAt and giveMassMtrxIntegrationgMask services to be implemented.
     * @param answer mass matrix
     * @param tStep time step
     */
    virtual void          computeMassMatrix(FloatMatrix &answer, TimeStep *);
    /**
     * Computes lumped mass matrix of receiver. Default implementation returns lumped consistent mass matrix.
     * Then returns lumped mass transformed into nodal coordinate system.
     * The lumping procedure zeroes all off-diagonal members and zeroes also all diagonal members
     * corresponding to non-displacement DOFs. Such diagonal matrix is then rescaled, to preserve
     * the element mass.
     * Requires the computeNmatrixAt and giveMassMtrxIntegrationgMask services to be implemented.
     * @param answer mass matrix
     * @param tStep time step
     */
    virtual void          computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *);
    /**
     * Computes consistent mass matrix of receiver using numerical integration over element volume.
     * Mass matrix is computed as \f$M=\int_V N^T \rho N dV\f$, where \f$N\f$ is displacement approximation matrix.
     * The number of necessary integration points  is determined using this->giveNumberOfIPForMassMtrxIntegration
     * service. Only selected degrees of freedom participate in integration of mass matrix. This is described
     * using dof mass integration mask. This mask is obtained from this->giveMassMtrxIntegrationgMask service.
     * The nonzero mask value at i-th position indicates that i-th element DOF participates in mass matrix
     * computation. The result is in element local coordinate system.
     * @param answer mass matrix
     * @param tStep time step
     * @param mass as result contain total mass of receiver.
     */
    virtual void          computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass);
    /**
     * Returns mask indicating, which unknowns (their type and ordering is the same as
     * element unknown vector) participate in mass matrix integration.
     * Nonzero value at i-th position
     * indicates that corresponding row in interpolation matrix N will participate in
     * mass matrix integration (typically only displacements are taken into account).
     * @param answer integration mask, if zero sized, all unknowns participate. This is default.
     * @returns zero sized mask
     */
    virtual void          giveMassMtrxIntegrationgMask(IntArray &answer)
    { answer.resize(0);
      return; }
    /**
     * Computes numerically stiffness matrix of receiver. Default implementation computes element stiffness using
     \f$K=\int_v B^T D B dV\f$ formulae, where \f$B\f$ is element geometric matrix and \f$D\f$ is material stiffness matrix.
     * No geometrical nonlinearity is taken into account. NUmerical integration procedure uses integrationRulesArray
     * for numrical integration. Support for reduced or selected integration is implemented.
     * If numberOfIntegrationRules is equal to
     * 1, the full integration of all coefficients is performed. Otherwise, integration is performed using following rules.
     * Each integration rule can specify start and end strain index of strain vector components for which is valid.
     * It is necessary to ensure that these strat and end indexes, dividing geometrical matrix into blocks,
     * are not overlapping and that each strain component is included.
     * Then stiffness matrix is obtained as summation of integrals \f$I_{ij}=\int_v B^T_i D_{IJ} B_j dV\f$
     * where \f$B_i\f$ is i-th block of geometrical matrix and \f$D_{ij}\f$ is corresponding constitutive sub-matrix.
     * The geometrical matrix is obtained using computeBmatrixAt service and the constitituve matrix is obtained using
     * computeConstitutiveMatrixAt service.
     * The \f$I_{ij}\f$ integral is evaluated using such integration rule, which is valid for i-th or j-th block
     * and has smaller number of integration points.
     * For higher numerical performance, only one half of stiffness matrix is computed and answer is then symmetrized.
     * Therefore, if element matrix will be generally nonsymmetric, one must specialize this method.
     * Finaly, the result is transformed into global coordinate system (or nodal coordinate system, if it is defined).
     * @param answer computed stiffness matrix (symmetric)
     * @param rMode response mode
     * @param tStep time step
     */
    virtual void          computeStiffnessMatrix(FloatMatrix &answer,
                                                 MatResponseMode rMode, TimeStep *tStep);

    /**
     * Computes initial stress matrix for linear stability problem.
     * Default implementation is not provided.
     * Please note, that initial stress matrix depends on normal forces of element,
     * corresponding engineering model must take this into account.
     * @param answer computed initial stress matrix
     * @param tStep time step.
     */
    virtual void          computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep)
    { _error("computeInitialStressMatrix: not implemented");
      return; }

    // load vector
    /**
     * Computes element load vector of receiver induced by non force influences.
     * Parts due to prescribed strains (like temperature) and due to  prescribed displacement are included.
     * (Precisely, result is computePrescribedStrainLoadVectorAt contribution plus
     * contribution of computeBcLoadVectorAt service)
     * @param answer computed element load vector in global c.s.
     * @param stepN time step
     * @param determines response mode
     */
    void                  computeNonForceLoadVector(FloatArray &answer, TimeStep *, ValueModeType mode);
    /**
     * Computes force dependent part of load vector. It is load vector induced by applied force loading.
     * Element body load and eleemnt boundary load (edge and surface load) is included.
     * (precisely result is summation of computeBodyLoadVectorAt, computeEdgeLoadVectorAt and
     * computeSurfaceLoadVectorAt service results contributions)
     * @param answer computed load vector
     * @param tStep time step
     * @param mode determines the response (total, incremental)
     */
    void computeForceLoadVector(FloatArray & answer, TimeStep *, ValueModeType);
    virtual void computeLocalForceLoadVector(FloatArray & answer, TimeStep *, ValueModeType);
    // stress equivalent vector in nodes (vector of internal forces)
    // - mainly for nonLinear Analysis.
    /**
     * Returns equivalent nodal forces vectors. Usefull for nonlinear analysis.
     * Default implementation computes result as \f$F=\int_v B^T \sigma dV\f$, where \f$\sigma\f$ is the
     * real element stress vector obtained using computeStressVector service (if useUpdatedGpRecord=0) or
     * (if useUpdatedGpRecord=1) from integration point status.
     * The geometric matrix is obtained using computeBmatrixAt service.
     * Integration is performed using default integration rule, which should produce always valid results,
     * assuming that strains used for computation of stresses are valid.
     * @param answer equivalent nodal forces vector
     * @param tStep time step
     * @param useUpdatedGpRecord if equal to zero, the stresses in integration points are computed (slow but safe), else if
     * nonzero the stresses are taken directly from integration point status (should be derived from StructuralMaterialStatus)
     * (fast, but engineering model must ensure valid status data in each integration point).
     */
    virtual void giveInternalForcesVector(FloatArray &answer,
                                          TimeStep *, int useUpdatedGpRecord = 0);

    /**
     * Compute strain vector of receiver evaluated at given integration point at time
     * step stepN from element displacement vector.
     * The nature of strain vector depends on the element type.
     * @param answer element strain vector
     * @param gp integration point
     * @param stepN time step
     */
    virtual void   computeStrainVector(FloatArray &answer, GaussPoint *, TimeStep *);

    /**
     * Returns the integration point corresponding value.
     * @param answer contain corresponding ip value, zero sized if not available
     * @param aGaussPoint integration point
     * @param type determines the type of internal variable
     */
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    /**
     * Computes at given time (stepN) the the resulting temperature component array.
     * This is summation of all temperature load  components of  receiver.
     * @param answer resulting temperature components of receiver
     * @param stepN time step
     * @param mode determines response mode.
     */
    virtual void  computeResultingIPTemperatureAt(FloatArray &answer, TimeStep *, GaussPoint *gp, ValueModeType mode);

    /**@name Methods related to nonlocal models */
    //@{
    /**
     * Updates internal element state (in all integration points of receiver)
     * before nonlocal averaging takes place. Used by so nonlocal materials,
     * because their response in particular point depends not only on state in this point, but
     * depends also on state in point's neighbourhood. Nonlocal quatity is computed as nonlocal
     * average of local quantities. Therefore, before updating integration point state depending on
     * nonlocal quantity (or quantities), local quantities in all integration points must be updated
     * in advance. This function updates local quantities of material model using
     * updateBeforeNonlocalAverage member function of structural nonlocal material class.
     */
    virtual void updateBeforeNonlocalAverage(TimeStep *atTime);
    /**
     * Returns the "nonlocal" location array of receiver. This is necessary, when stiffness matrix
     * of nonlocal model is assembled. Since model is nonlocal, the value at given IP depends on
     * other IP (generally belonging to different elements) and as a consequence leads to
     * increase of stifness matrix profile, to take into account this "remote" dependency.
     */
    virtual void giveNonlocalLocationArray(IntArray &locationArray);
    /**
     * Adds the "nonlocal" contribution to stiffness matrix, to account for nonlocality of
     * material model. Typically, this contribution is obtained by summing up mutual IP contributions.
     */
    virtual void addNonlocalStiffnessContributions(SparseMtrx &dest, TimeStep *atTime);
    //@}

    /**
     * Updates the internal state variables stored in all IPs according to
     * already mapped state.
     * @param oldd old mesh reference
     * @param tStep time step
     * @return nonzero if o.k.
     */
    virtual int adaptiveUpdate(TimeStep *tStep);




    // time step termination
    // void                  printOutputAt (FILE *, TimeStep*) ;
    /**
     * Updates element state corresponding to newly reached solution.
     * It computes stress vector in each element integration point (to ensure that data in integration point's
     * statuses are valid).
     * @param tStep finished time step
     */
    void                  updateInternalState(TimeStep *);

    // consistency check
    /**
     * Checks internal data consistency. Default implementation checks if material implements
     * Material_StructuralCapability and if cross section model implements CS_StructuralCapability
     * extension (or interface).
     * @return nonzero if consistency check is o.k.
     */
    virtual int    checkConsistency();

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "StructuralElement"; }
    /// Returns classType id of receiver.
    classType                giveClassID() const
    { return StructuralElementClass; }
    /*      int                   giveNumber ()
     * { return FEMComponent::giveNumber() ;} */

#ifdef __OOFEG
    /**
     * Returns internal state variable (like stress,strain) at node of element in Reduced form,
     * the way how is obtained is dependent on InternalValueType.
     * The value may be local, or smoothed useing some recovery technigue/
     * returns zero if element is unable to respont to request.
     * @param answer contains result, zero sized if not supported
     * @param type determines the internal variable requested (physical meaning)
     * @param mode determines the mode of variable (recovered, local, ...)
     * @param node node number, for which variable is required
     * @param atTime time step
     * @return nonzero if o.k, zero otherwise
     */
    virtual int   giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                          int node, TimeStep *atTime);
    /// Shows sparse structure
    void showSparseMtrxStructure(CharType mtrx, oofegGraphicContext &gc, TimeStep *atTime);
    /// Shows extended sparse structure (for example, due to nonlocal interactios for tangent stiffness)
    virtual void showExtendedSparseMtrxStructure(CharType mtrx, oofegGraphicContext &gc, TimeStep *atTime);

#endif

protected:
    /**
     * Computes constitutive matrix of receiver. Default implementation uses element cross section
     * giveCharMaterialStiffnessMatrix service.
     * @param answer computed answer
     * @param rMode material response mode of answer
     * @param gp integration point for which constitutive matrix is computed
     * @param tStep time step
     */
    virtual void          computeConstitutiveMatrixAt(FloatMatrix &answer,
                                                      MatResponseMode rMode, GaussPoint *,
                                                      TimeStep *tStep);

    /**
     * Computes load vector of receiver due to the prescribed displacements of receiver nodes.
     * Implementation supports the changes of static system (must be also supported by
     * engineering model).
     * @param answer load vector due to precsribed b.c., zero sized answer if load vector is zero.
     * @param tStep time step, vhen load vector is assembled
     * @param mode determines response mode
     */
    void                  computeBcLoadVectorAt(FloatArray &answer, TimeStep *, ValueModeType mode);
    /**
     * Computes the load vector due to body load acting on receiver, at given time step.
     * Default implementation computes body load vector numerically as \f$l=\int_V N^T f \rho dV\f$
     * using default integration rule. Result is transfromed to global c.s.
     * @param answer computed load vector due to body load
     * @param forLoad pointer to Bodyload object, which contribution is computed
     * @param stepN time step
     * @param mode determines the response mode
     */
    virtual void          computeBodyLoadVectorAt(FloatArray &answer, Load *, TimeStep *, ValueModeType mode);

    // BEGIN START edge and surface load support
    /**@name Edge and surface load support services.
     */
    //@{
    /**
     * Computes point load vector contribution of receiver for given load (should has BoundaryLoad Base).
     * @param answer computed load vector
     * @param load edge load object
     * @param tStep time step
     * @param mode determines response mode
     */
    virtual void   computePointLoadVectorAt(FloatArray &answer, Load *, TimeStep *, ValueModeType mode);
    /**
     * Computes edge load vector contribution of receiver for given load (should has BoundaryLoad Base).
     * Each edge should have unique number assigned to identify it.
     * The default implementation does integration of load vetor in local edge space
     * (i.e. one dimensional integration is performed on line). This general implementation requires
     * that element must provide following services:
     * <UL>
     * <LI>
     * ComputeEgdeNMatrixAt - returns interpolation matrix of local edge DOFs in the local edge space.</LI>
     * <LI>
     * computeEdgeVolumeAround - returns volume corresponding to integration point at local edge.</LI>
     * <LI>
     * GiveEdgeDofMapping - returns integer array specifying local dof edge mapping to "global" element dofs.</LI>
     * </UL>
     * Integration rule is set up automatically, based on element interpolation order and load approximation.
     * Integration points are set-up using standart integration rule services (setUpIntegrationPoints method).
     * Gauss integration rule is used.
     * If derived class overrides this default implementation somehow, the above services
     * must not be implemented.
     *
     * @param answer computed load vector
     * @param load edge load object
     * @param tStep time step
     * @param mode determines response mode
     */
    virtual void   computeEdgeLoadVectorAt(FloatArray &answer, Load *, int, TimeStep *, ValueModeType mode);
    /**
     * Computes surface load vector contribution of receiver for given load (should has BoundaryLoad Base).
     * Each surface should have unique number assigned to identify it.
     * The default implementation does integration of load vetor in local surface space
     * (i.e. two dimensional integration is performed on triangle or square).
     * This general implementation requires
     * that element must provide following services:
     * <UL>
     * <LI>
     * GetSurfaceIntegrationRule - returns integration rule for surface for given polynomial order.</LI>
     * <LI>
     * ComputeSurfaceNMatrixAt - returns interpolation matrix of local surfare DOFs in the local edge space.</LI>
     * <LI>
     * computeSurfaceVolumeAround - returns volume corresponding to integration point of local surface.</LI>
     * <LI>
     * GiveSurfaceDofMapping - returns integer array specifying local dof surface mapping to "global" element dofs.</LI>
     * </UL>
     * Integration rule is set up automatically, based on element interpolation order and load approximation.
     * Integration points are set-up using standart integration rule services (setUpIntegrationPoints method).
     * Gauss integration rule is used.
     * If derived class overrides this default implementation somehow, the above services
     * must not be implemented.
     *
     * @param answer computed load vector
     * @param load surface load object
     * @param tStep time step
     * @param mode determines response mode
     */
    virtual void   computeSurfaceLoadVectorAt(FloatArray &answer, Load *, int, TimeStep *, ValueModeType mode);
    // interpolation matrices for nonzero unknowns
    // mapping to global element dofs is done using maps obtained
    // from GiveEdgeDofMapping or GiveSurfaceDofMapping functions with
    // particular edge or surface number as a parameter.

    /**
     * Computes Edge interpolation matrix. Interpolation matrix provide way, how to compute
     * local edge unknowns (nonzero element unknowns on edge) at any integration point of edge, based on
     * local edge unknowns in edge nodes.
     * The edge numbering and local edge coordinate system is element dependent.
     * The integration point is specified using one-dimensional iso coordinates.
     * @param answer interpolation matrix of edge
     * @param gp integration point
     */
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *)
    { answer.resize(0, 0); }
    /**
     * Computes surface interpolation matrix. Interpolation matrix provide way, how to compute
     * local surface unknowns (nonzero element unknowns on surface) at any integration point of surface, based on
     * local unknowns in surface nodes.
     * Local coordinate system of surfaceedge and element surface numbering is element dependent.
     * The integration point is specified using two-dimensional iso coordinates, or using area coordinates
     * for triangular surface.
     * @param answer interpolation matrix of edge
     * @param gp integration point
     */
    virtual void computeSurfaceNMatrixAt(FloatMatrix &answer, GaussPoint *)
    { answer.resize(0, 0); }
    /**
     * Assembles edge dof mapping mask, which provides mapping between edge local DOFs and "global" element
     * DOFs. Mask can be imagined as local edge code numbers used to localize local edge DOFs to
     * element DOFs.
     * @param answer edge DOF mask
     * @pram i edge number
     */
    virtual void  giveEdgeDofMapping(IntArray &answer, int) const
    { answer.resize(0); }
    /**
     * Assembles surface dof mapping mask, which provides mapping between surface local DOFs and "global" element
     * DOFs. Mask can be imagined as local surface code numbers used to localize local DOFs to
     * element DOFs.
     * @param answer surface DOF mask
     * @pram i surface number
     */
    virtual void  giveSurfaceDofMapping(IntArray &answer, int) const
    { answer.resize(0); }
    //virtual int    hasEdgeLoadSupport () {return 0;}
    //virtual int    hasSurfaceLoadSupport () {return 0;}
    /**
     * Returns integration rule for integration over element surface.
     * @param i order of integrated polynomial
     * @returns best integration rule to integrate polynomial of order i over element surface.
     */
    virtual IntegrationRule *GetSurfaceIntegrationRule(int) { return NULL; }
    /**
     * Computes volume related to integration point on local edge.
     * @param gp edge integration point
     * @param i edge number
     */
    virtual double        computeEdgeVolumeAround(GaussPoint *, int) { return 0.; }
    /**
     * Computes volume related to integration point on local surface.
     * @param gp surface integration point
     * @param i surface number
     */
    virtual double        computeSurfaceVolumeAround(GaussPoint *, int) { return 0.; }
    /**
     * Computes global coordinates of integration point on local edge.
     * @param answer global coordinates
     * @param gp edge integration point
     * @param i edge number
     */
    virtual void         computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *, int) { answer.resize(0); }
    /**
     * Computes global coordinates of integration point on  local surface.
     * @param answer global coordinates
     * @param gp surface integration point
     * @param i surface number
     */
    virtual void        computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *, int) { answer.resize(0); }





    // some nessesary transformation
    // Global to local element c.s transformation for load vector dofs
    /**
     * Returns transformation matrix from global coordinate system to local
     * element coordinate system for element load vector components.
     * If no transformation is necessary, answer is empty matrix (default);
     * @return nonzero if transformation matrix is not empty matrix, zero otherwise
     */
    virtual int  computeLoadGToLRotationMtrx(FloatMatrix &answer) { answer.beEmptyMtrx();
                                                                    return 0; }
    // Local edge (LE-local Edge c.s) or surface (LS-local surface c.s) c.s
    // to element local c.s for load vector dofs
    /**
     * Returns transformation matrix from local edge c.s  to element local coordinate system
     * of load vector components. Necessary, because integration must be done in local coordinate
     * system of entity (edge or surface).
     * If no transformation is necessary, answer is empty matrix (default);
     * @param i edge number
     * @param gp integration point (point, where transformation is computed, usefull for curved edges)
     * @return nonzero if transformation matrix is not empty matrix, zero otherwise
     */
    virtual int  computeLoadLEToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *) { answer.beEmptyMtrx();
                                                                                          return 0; }
    /**
     * Returns transformation matrix from local surface c.s  to element local coordinate system
     * of load vector components. Necessary, because integration must be done in local coordinate
     * system of entity (edge or surface).
     * If no transformation is necessary, answer is empty matrix (default);
     * @param i surface number
     * @param gp integration point (point, where transformation is computed, usefull for curved surfaces)
     * @return nonzero if transformation matrix is not empty matrix, zero otherwise
     */
    virtual int  computeLoadLSToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *) { answer.beEmptyMtrx();
                                                                                          return 0; }
    // END edge and surface load support
    //@}

    /**
     * Computes load vector due to prescribed strains. The load vector is obtained using numerical integration
     * (using default intagration rule of element) over
     * element volume \f$f=\int_V B^T D \varepsilon dV\f$, where \f$\varepsilon\f$ is stress independent strain vector
     * in particular untegration point, obtained using computeStressIndependentStrainVector service.
     * The load mode (Incremental or Total Load form) is  passed as parameter.
     * @param answer computed load vector contribution
     * @param tStep time step
     * @param mode load vector mode
     */
    void computePrescribedStrainLoadVectorAt(FloatArray & answer, TimeStep *, ValueModeType);
    virtual void computePrescribedStrainLocalLoadVectorAt(FloatArray & answer, TimeStep *, ValueModeType);

    // strains and stresses
    /**
     * Computes the stress vector of receiver at given integration point, at time step stepN.
     * The nature of these stresses depends on the element's type.
     * @param answer stress vector
     * @param gp integration point
     * @param stepN time step
     */
    virtual void   computeStressVector(FloatArray &answer, GaussPoint *, TimeStep *);

    /**
     * Computes the geometrical matrix of receiver in given integration point.
     * The product of this matrix (assembled at given integration point) and element displacement
     * vector is element strain vector. If lowerIndx and upperIndx parameters are specified,
     * answer is formed only for strains within this interval. This will affects the size of answer.
     *
     * @parm gp integration point for which answer is computed
     * @param answer geometric matrix of receiver
     * @param lowerIndx if specified, answer is formed only for strain with index equal and greater than  lowerIndx.
     * This parameter has default value 1 (answer is formed from first strain).
     * @param upperIndx if specified, answer is formed only for strain with index less and equal than  upperIndx.
     * This parameter has default value ALL_STRAINS (answer is formed for all strains).
     */
    virtual void  computeBmatrixAt(GaussPoint *, FloatMatrix &answer,
                                   int lowerIndx = 1, int upperIndx = ALL_STRAINS) = 0;
    /* computes interpolation matrix of unknowns which are taken into account
     * when integration of consiistent mass matrix is performed. */
    /**
     * Computes interpolation matrix for element unknowns.
     * The order and meaning of unknowns is element dependent.
     * @param gp integration point for which answer is assembled
     * @param answer interpolation matrix evaluated at gp
     */
    virtual void  computeNmatrixAt(GaussPoint *, FloatMatrix &)  = 0;

    /**
     * Updates rotation matrix r(l)=T r(g*) between  local and global coordinate system
     * taking into account also possible local - coordinate system in some elements
     * nodes.
     * Default implementation uses \ref computeGtoLRotationMatrix and
     \ref computeGNDofRotationMatrix  services to compute result.
     * Default implementation uses cached rotation matrix in
     * rotationMatrix attribute, so rotation matrix is computed only once.
     * @return nonzero if transformation is necessary.
     */
    virtual int   updateRotationMatrix();      //
    // give Transformation matrix from global coord. sysyt. to element-local c.s
    // i.e. r(l)=T r(h), if no trasformation necessary set anser to empty mtrx
    /**
     * Returns  transformation matrix from global coord. system to local element
     * coordinate system ( i.e. r(l)=T r(g)). if no trasformation is necessary
     * then answer is empty mtrx and zero value is returned.
     * @return nonzero if transformation is necessary, zero otherwise.
     */
    virtual int  computeGtoLRotationMatrix(FloatMatrix &answer) { answer.beEmptyMtrx();
                                                                  return 0; }
    // give Transformation matrix from global coord. syst. to local coordinate system in nodes.
    // i.e. r(n)=T r(g), if no trasformation necessary sets answer to empty mtrx.
    /**
     * Returns transformation matrix for DOFs from global coordinate system
     * to local coordinate system in nodes (i.e. r(n)=T r(g)) if mode == _toNodalCS.
     * If mode == _toGlobalCS, the transformation from local nodal cs to
     * global cs in node is returned. If no trasformation is
     * necessary sets answer to empty mtrx and returns zero value.
     * @return nonzero if transformation is necessary, zero otherwise.
     */
    virtual int  computeGNDofRotationMatrix(FloatMatrix &answer, DofManTrasfType mode);
    /**
     * Returns transformation matrix for loading from global coordinate system
     * to local coordinate system in nodes (i.e. r(n)=T r(g)) if mode == _toNodalCS.
     * If mode == _toGlobalCS, the transformation from local nodal cs to
     * global cs in node is returned. If no trasformation is
     * necessary sets answer to empty mtrx and returns zero value.
     * @return nonzero if transformation is necessary, zero otherwise.
     */
    virtual int  computeGNLoadRotationMatrix(FloatMatrix &answer, DofManTrasfType mode);

    // returns maximal approximation order of the receiver
    /**
     * Returns maximum approximation order used by receiver.
     * Must be implemented by derived classes
     */
    virtual int           giveApproxOrder() { return 0; }
    /**
     * Returns integration domain for receiver, used to initialize
     * integration point over receiver volume. Must be specialized.
     */
    virtual integrationDomain giveIntegrationDomain() { return _Unknown_integrationDomain; }
    /**
     * Return desired number of integration points for consistent mass matrix
     * computation, if required.
     */
    virtual int  giveNumberOfIPForMassMtrxIntegration() { return 0; }

    /**
     * General service for condenstaion of stiffness and optionally load vector and mass or initial stress matrices
     * of receiver.
     * @param stiff stifness matrix to be condensed. Must be specified.
     * @param mass mass or initial stress maatrix. If parameter is NULL, only stifness and/or load is condensed.
     * @param load load vector of receiver. If specified then it is condensed. If no load vector condensation is necessary
     * set parameter to NULL pointer.
     * @param what integar array. If at i-th position is nonzero, then i-th comoponet is condensed.
     */
    void condense(FloatMatrix *stiff, FloatMatrix *mass, FloatArray *load, IntArray *what);


    friend  class IDNLMaterial;
};

#endif // structuralelement_h








