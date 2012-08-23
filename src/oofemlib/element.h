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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef element_h
#define element_h

#include "femcmpnn.h"
#include "flotmtrx.h"
#include "flotarry.h"

#include "alist.h"
#include "intarray.h"
#include "error.h"
#include "integrationrule.h"
#include "classtype.h"
#include "chartype.h"
#include "elementgeometrytype.h"
#include "equationid.h"
#include "valuemodetype.h"
#include "internalstatetype.h"
#include "internalstatevaluetype.h"
#include "elementextension.h"
#include "entityrenumberingscheme.h"
#include "unknowntype.h"
#include "unknownnumberingscheme.h"

#ifdef __OOFEG
 #include "node.h"
 #include "engngm.h"
#endif

#include <cstdio>

namespace oofem {
class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class IntArray;
class CrossSection;
class ElementSide;
class FEInterpolation;

#ifdef __PARALLEL_MODE
class CommunicationBuffer;

/**
 * In parallel mode, this type indicates the mode of element.
 * In the case of element cut mode, the cut element is local on all partitions sharing it.
 * Some of such element nodes are local and some are remote. The local nodes are completely
 * surrounded by local element on particular partition.
 */
enum elementParallelMode {
    Element_local, ///< Element is local, there are no contributions from other domains to this element.
    // Element_shared, ///< Element is shared by neighboring partitions - not implemented.
    Element_remote, ///< Element in active domain is only mirror of some remote element.
};

#endif


/**
 * Abstract base class for all finite elements. Derived classes should be  base
 * classes for specific analysis type (for example base class for structural analysis,
 * thermal analysis or magnetostatics one). These derived classes then declare
 * analysis-specific part of interface and they provide default implementation
 * for these methods.
 * This abstract class declares (and possibly implements) general data and methods
 * common to all element types. General methods for obtaining characteristic vectors,
 * matrices and values are introduced and should be used instead of calling directly
 * specific member functions (these must be overloaded by derived analysis-specific
 * classes in order to invoke proper method according to type of component requested).
 * In general, element class member functions are called in following cases:
 * - When Engineering model assembles governing equation(s) from element's contributions.
 *   Typically when assembles some characteristic matrix of analyzed structure (method assemble),
 *   asks each element for its code numbers and for corresponding characteristic matrix, vector or
 *   value. This will typically uses above mentioned general method for obtaining characteristic
 *   matrices, vectors and values from elements.
 * - When Engineering model has some characteristic matrix stored in sparse form, then is necessary
 *   to build internal profile of such sparse matrix. Engineering model then typically calls
 *   buildInternalStructure member function of sparse matrix. This method then requests element code
 *   numbers from elements and builds its internal profile. This could look strange, because
 *   class sparse matrix should know "something about elements", but only sparse matrix knows, how
 *   to build its internal structure (this may require one or more loops over elements code number).
 * - When element computes its contribution, then it communicates with its cross-section (and cross
 *   section with corresponding material). For some cross-section or material models some further
 *   communication between these classes and element may be necessary (for example in cases of
 *   layered cross section model strains in each layer has to be evaluated from "integrated"
 *   strains (integrated means element-like strains including curvatures in case of beams and plates)
 *   by element, in cases of non local material and many other examples can be found).
 *
 * There are some general rules, that programmer must take into account.
 * - Element stores the numbers of its dofmanagers in dofManArray.
 *   These include nodes, element sides and internal DOFs that are not condensed at element level.
 *   Their order and meaning are determined by element definition.
 *   Local ordering of dofs for particular element is determined by local numbering of
 *   dofmanagers and their corresponding dofs. DOFS necessary for particular node/side is specified using node/side dof mask.
 *   Local DOF ordering must be taken into account when assembling various local characteristic
 *   vectors and matrices.
 */
class Element : public FEMComponent
{
protected:
    /// Number of dofmanagers
    int numberOfDofMans;
    /// Array containing dofmanager numbers.
    IntArray dofManArray;
    /// Number of associated material.
    int material;
    /// Number of associated cross section.
    int crossSection;
    /**
     * Array containing indexes of loads (body loads and boundary loads are kept separately),
     * that apply on receiver.
     */
    IntArray bodyLoadArray, boundaryLoadArray;
    /// Number of integration rules used by receiver.
    int numberOfIntegrationRules;
    /**
     * List of integration rules of receiver (each integration rule contains associated
     * integration points also). This list should contain only such integration rules,
     * that are used to integrate results depending on load time history. For all integration
     * points in these rules, history variables are stored and updated.
     * For integrations, where history stored in gauss points is not necessary
     * (mass matrix integration) and different integration rule is needed, one should preferably
     * use temporarily created integration rule.
     */
    IntegrationRule **integrationRulesArray;
    /// Array of code numbers of element.
    IntArray *locationArray;

    /// Transformation material matrix, used in orthotropic and anisotropic materials, global->local transformation
    FloatMatrix elemLocalCS;

#if defined ( __PARALLEL_MODE ) || defined ( __ENABLE_COMPONENT_LABELS )
    /**
     * In parallel mode, globalNumber contains globally unique DoFManager number.
     * The component number, inherited from FEMComponent class contains
     * local domain number.
     */
    int globalNumber;
#endif
#ifdef __PARALLEL_MODE
    elementParallelMode parallel_mode;
    /**
     * List of partition sharing the shared element or
     * remote partition containing remote element counterpart.
     */
    IntArray partitions;
#endif

public:
    /**
     * Constructor. Creates an element with number n belonging to domain aDomain.
     * @param n Element's number
     * @param aDomain Pointer to the domain to which element belongs.
     */
    Element(int n, Domain *aDomain);
    /// Virtual destructor.
    virtual ~Element();

    /**@name Methods referring to code numbers */
    //@{
    /**
     * Returns the location array (array of code numbers) of receiver for given numbering scheme.
     * Results are cached at receiver for default scheme in locationArray attribute.
     * The UnknownType parameter allows to distinguish between several possible governing equations, that
     * can be numbered separately. The default implementation assumes that location array will be assembled only for
     * one UnknownType value, and this array is cached on element level.
     */
    virtual void giveLocationArray(IntArray & locationArray, EquationID, const UnknownNumberingScheme & s) const;
    /**
     * Returns the location array for the boundary of the element.
     * The boundary is the corner nodes for 1D elements, the edges for a 2D element and the surfaces for a 3D element.
     */
    void giveBoundaryLocationArray(IntArray &locationArray, int boundary, EquationID eid, const UnknownNumberingScheme &s);
    /**
     * Invalidates location array in receiver. Each element stores its copy of location array(s), in order
     * to avoid time consuming assembly of code numbers every time when requested. Some engineering models
     * may support dynamic changes of static system (generally of boundary conditions during analysis),
     * then these models use this function to invalidate location array after finishing time step,
     * to enforce elements to update they code numbers, which may change. Changes of static system will lead
     * to different number of equations, which requires special attention (for example internal structure
     * of sparse matrices  must be reinitialized and characteristic vectors and matrices has to be
     * newly assembled).
     */
    void invalidateLocationArray();
    /**
     * @return Number of DOFs in element.
     */
    virtual int giveNumberOfDofs() { return 0; }
    /**
     * @return Number of internal DOF managers of element.
     */
    virtual int giveNumberOfInternalDofManagers() { return 0; }
    /**
     * Returns i-th internal element dof manager of the receiver
     * @param i Internal number of DOF.
     * @return DOF number i.
     */
    virtual DofManager *giveInternalDofManager(int i) const {
        _error2("No such DOF available on Element %d", number);
        return NULL;
    }
    //@}
    /**
     * @name General methods for obtaining element contributions
     * Note: These member functions have to  be overloaded by derived analysis-specific
     * classes in order to invoke proper method according to type of component requested.
     * Derived classes from these analysis-specific classes should not modify these functions.
     */
    //@{
    /**
     * Computes characteristic matrix of receiver of requested type in given time step.
     * @param answer Requested characteristic matrix.
     * If element has no capability to compute requested type of characteristic matrix
     * error function is invoked.
     * @param type   Id of characteristic component requested.
     * @param tStep  Time step when answer is computed.
     */
    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep);
    /**
     * Computes characteristic vector of receiver of requested type in given time step.
     * If element has no capability to compute requested type of characteristic vector
     * error function is invoked.
     * @param answer Requested characteristic vector.
     * @param type   Id  of characteristic component requested.
     * @param mode   Determines mode of answer.
     * @param tStep  Time step when answer is computed.
     */
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep);
    /**
     * Computes characteristic value of receiver of requested type in given time step.
     * If element has no capability to compute requested type of characteristic value
     * error function is invoked.
     * @param type  Id of characteristic component requested.
     * @param tStep Time step when answer is computed.
     * @return Requested value.
     */
    virtual double giveCharacteristicValue(CharType type, TimeStep *tStep);
    //@}

    //virtual MaterialMode giveMaterialMode() {return _Unknown;}

    /**@name General element functions */
    //@{
    /**
     * Returns local vector of unknowns. Local vector of unknowns is extracted from
     * engineering model global unknowns vector (if specific dof has assigned
     * corresponding equation number) and from boundary conditions (if dof has active boundary
     * and possibly initial condition). Because unknowns are obtained from engineering
     * model, this must support queries for given unknown and unknown mode. Also
     * engineering model must not be able to complete request for any TimeStep, because
     * keeping history of all unknowns may became impossible. But generally all engineering models
     * should be able return supported unknowns at current and previous time step. Consult
     * reference manual for particular engineering model.
     *
     * @param type   Identifies unknown type (eg. displacement or temperature vector).
     * @param u      Identifies mode of unknown (eg. total value or velocity of unknown).
     * @param stepN  Time step, when vector of unknowns is requested.
     * @param answer Local vector of unknowns.
     */
    virtual void computeVectorOf(EquationID type, ValueModeType u, TimeStep *stepN, FloatArray &answer);
    /**
     * Returns local vector of unknowns. Local vector of unknowns is extracted from
     * given field and from boundary conditions (if dof has active boundary
     * and possibly initial condition). Because unknowns are obtained from given field
     * model, this must support queries for given unknown.
     *
     * @param field  Source field (eg. displacement or temperature vector).
     * @param u      Value mode of unknown (incremental, total, ...).
     * @param stepN  Time step, when vector of unknowns is requested.
     * @param answer Local vector of unknowns.
     */
    virtual void computeVectorOf(PrimaryField &field, ValueModeType u, TimeStep *stepN, FloatArray &answer);
    /**
     * Returns local vector of prescribed unknowns. Local vector of prescribed unknowns is
     * extracted from nodal (and side - if they hold unknowns) boundary conditions.
     *
     * @param ut     Identifies mode of unknown (eg. total values or velocity of unknown).
     * @param type   Value mode of unknown (incremental, total, ...).
     * @param stepN  Time step, when vector of prescribed unknowns is requested.
     * @param answer Local vector of prescribed unknowns. If unknown is not prescribed,
     * zero value is placed on position of free dof.
     */
    void computeVectorOfPrescribed(EquationID ut, ValueModeType type, TimeStep *stepN, FloatArray &answer);

    /**
     * Computes or simply returns total number of element's local DOFs.
     * Must be defined by particular element.
     * @param eid Id of equation that DOFs belong to.
     * @return Number of DOFs belonging to ut.
     */
    virtual int computeNumberOfDofs(EquationID eid) { return 0; }
    /**
     * Computes the total number of element's global dofs.
     * The transitions from global c.s. to nodal c.s. should NOT be included.
     * @param eid Id of equation that DOFs belong to.
     * @return Total number of DOFs belonging to ut.
     */
    virtual int computeNumberOfGlobalDofs(EquationID eid);
    /**
     * Computes the total number of element's primary master DOFs.
     * @param eid ID of equation that DOFs belong to.
     * @return Total number of DOFs belonging to eid.
     */
    int computeNumberOfPrimaryMasterDofs(EquationID eid);
    /**
     * Returns transformation matrix from global c.s. to local element
     * c.s., i.e. @f$ r_l =T r_g @f$.
     * If no transformation is necessary then answer is empty matrix and zero value is returned.
     * @param answer Computed rotation matrix.
     * @return Nonzero if transformation is necessary, zero otherwise.
     */
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer) {
        answer.beEmptyMtrx();
        return false;
    }
    /**
     * Transformation matrices updates rotation matrix between element-local and primary DOFs,
     * taking into account nodal c.s. and master DOF weights.
     * @param eid Equation ID to compute rotation matrix for.
     * @return True if there is a rotation required, false otherwise.
     */
    virtual bool giveRotationMatrix(FloatMatrix &answer, EquationID eid);
    /**
     * Returns transformation matrix for DOFs from global coordinate system
     * to local coordinate system in nodes.
     * Also includes transformations to slave DOFs.
     * If no transformation is necessary sets answer to empty matrix and returns false.
     * Local stiffness matrix of element should be rotated with answer before assembly.
     * @note Function does most likely NOT need to be overridden.
     * @param answer Computed rotation matrix.
     * @param eid Equation ID.
     * @return True if transformation is necessary, false otherwise.
     */
    virtual bool computeDofTransformationMatrix(FloatMatrix &answer, EquationID eid);
    /**
     * Returns dofmanager dof mask for node. This mask defines the dofs which are used by element
     * in node. Mask influences the code number ordering for particular node. Code numbers are
     * ordered according to node order and dofs belonging to particular node are ordered
     * according to this mask. If element requests dofs using node mask which are not in node
     * then error is generated. This masking allows node to be shared by different elements with
     * different dofs in same node. Elements local code numbers are extracted from node using
     * this mask. Must be defined by particular element.
     *
     * @param inode  Mask is computed for local dofmanager with inode number.
     * @param ut     Equation DOFs belong to.
     * @param answer Mask for node.
     */
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const { answer.resize(0); }
    /**
     * Returns internal  dofmanager dof mask for node. This mask defines the dofs which are used by element
     * in node. Mask influences the code number ordering for particular node. Code numbers are
     * ordered according to node order and dofs belonging to particular node are ordered
     * according to this mask. If element requests dofs using node mask which are not in node
     * then error is generated. This masking allows node to be shared by different elements with
     * different dofs in same node. Elements local code numbers are extracted from node using
     * this mask. Must be defined by particular element.
     *
     * @param inode Mask is computed for local dofmanager with inode number.
     * @param ut Unknown type (support for several independent numberings within problem)
     * @param answer mask for node.
     */
    virtual void giveInternalDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
    { answer.resize(0); }
    /**
     * Returns element dof mask for node. This mask defines the dof ordering of the element interpolation.
     * Must be defined by particular element.
     *
     * @param ut Equation DOFs belong to.
     * @param answer DOF mask for receiver.
     */
    virtual void giveElementDofIDMask(EquationID ut, IntArray &answer) const { answer.resize(0); }
    /**
     * Returns volume related to given integration point. Used typically in subroutines,
     * that perform integration over element volume. Should be implemented by particular
     * elements.
     * @see GaussPoint
     * @param gp Integration point for which volume is computed.
     * @return Volume for integration point.
     */
    virtual double computeVolumeAround(GaussPoint *gp) { return 0.; }
    /// Computes the volume, area or length of the element depending on its spatial dimension.
    double computeVolumeAreaOrLength();
    /**
     * Computes the size of the element defined as its length.
     * @return Length, square root of area or cube root of volume (depending on spatial dimension).
     */
    double computeMeanSize();
    /**
     * Computes the volume.
     * @return Volume of element.
     */
    virtual double computeVolume();
    /**
     * Computes the area (zero for all but 2d geometries).
     * @return Element area.
     */
    virtual double computeArea();
    /**
     * Computes the length (zero for all but 1D geometries)
     * @return Element length.
     */
    virtual double computeLength();
    // If the need arises;
    /*
     * Computes the length of an edge.
     * @param iedge Edge number.
     * @return Edge length.
     */
    //virtual double computeEdgeLength(int iedge) { return 0.0; }
    /*
     * Computes the area of a surface.
     * @param isurf Surface number.
     * @param Surface area.
     */
    //virtual double computeSurfaceArea(int isurf) { return 0.0; }

    // data management
    /**
     * Translates local to global indices for dof managers.
     * @param i Local index of dof manager.
     * @return Global number of i-th dofmanager of element
     */
    int giveDofManagerNumber(int i) const { return dofManArray.at(i); }
    /// @return Receivers list of dof managers.
    IntArray &giveDofManArray() { return dofManArray; }
    /**
     * @param i Local index of the dof manager in element.
     * @return The i-th dofmanager of element.
     */
    DofManager *giveDofManager(int i) const;
    /**
     * Returns reference to the i-th node of element.
     * Default implementation returns i-th dofmanager of element converted to
     * Node class (check is made).
     * @param i Local index of node in element.
     * @return Requested node.
     */
    virtual Node *giveNode(int i) const;
    /**
     * Returns reference to the i-th side  of element.
     * Default implementation returns i-th dofmanager of element converted to
     * ElementSide class (check is made).
     * @param i Side number.
     * @return Requested element side.
     */
    virtual ElementSide *giveSide(int i) const;
    /// @return Interpolation of the element geometry, or NULL if none exist.
    virtual FEInterpolation *giveInterpolation() { return NULL; }
    /**
     * Returns the interpolation for the specific dof id.
     * Special elements which uses a mixed interpolation should reimplement this method.
     * @param id ID of the dof for the for the requested interpolation.
     * @return Appropriate interpolation, or NULL if none exists.
     */
    virtual FEInterpolation *giveInterpolation(DofIDItem id) { return giveInterpolation(); }
    /// @return Reference to the associated material of element.
    Material *giveMaterial();
    /// @return Reference to the associated crossSection of element.
    CrossSection *giveCrossSection();
    /**
     * Sets the material of receiver.
     * @param matIndx Index of new material.
     */
    void setMaterial(int matIndx) { this->material = matIndx; }
    /**
     * Sets the cross section model of receiver.
     * @param csIndx Index of new cross section.
     */
    void setCrossSection(int csIndx) { this->crossSection = csIndx; }

    /// @return Number of dofmanagers of receiver.
    int giveNumberOfDofManagers() const { return numberOfDofMans; }
    /**
     * Returns number of nodes of receiver.
     * Default implementation returns number of dofmanagers of element
     * @return Number of nodes of element.
     */
    virtual int giveNumberOfNodes() const { return numberOfDofMans; }
    /**
     * Sets receiver dofManagers.
     * @param dmans Array with dof manager indices.
     */
    void setDofManagers(const IntArray &dmans);

    /**
     * Sets integration rules.
     * @param irlist List of integration rules.
     */
    void setIntegrationRules(AList< IntegrationRule > *irlist);
    /**
     * Returns integration domain for receiver, used to initialize
     * integration point over receiver volume. Must be specialized.
     * @see IntegrationRule
     * @return Integration domain of element.
     */
    virtual integrationDomain giveIntegrationDomain() {
        IntegrationRule *ir = giveDefaultIntegrationRulePtr();
        return ir ? ir->giveIntegrationDomain() : _Unknown_integrationDomain;
    }
    /**
     * Returns material mode for receiver integration points. Should be specialized.
     * @return Material mode of element.
     */
    virtual MaterialMode giveMaterialMode() { return _Unknown; }
    /**
     * Assembles the code numbers of given integration element (sub-patch)
     * This is done by obtaining list of nonzero shape functions and
     * by collecting the code numbers of nodes corresponding to these
     * shape functions.
     * @return Nonzero if integration rule code numbers differ from element code numbers.
     */
    virtual int giveIntegrationRuleLocalCodeNumbers(IntArray &answer, IntegrationRule *ie, EquationID ut)
    { return 0; }

    // Returns number of sides (which have unknown dofs) of receiver
    //int giveNumberOfSides () {return numberOfSides;}

    /// @return Corresponding element region. Currently corresponds to cross section model number.
    int giveRegionNumber();

    /// Performs post initialization steps.
    void postInitialize();

    /**
     * Updates element state after equilibrium in time step has been reached.
     * Default implementation updates all integration rules defined by
     * integrationRulesArray member variable. Doing this, all integration points
     * and their material statuses are updated also. All temporary history variables,
     * which now describe equilibrium state are copied into equilibrium ones.
     * The existing internal state is used for update.
     * @param tStep Time step for newly reached state.
     * @see Material::updateYourself
     * @see IntegrationRule::updateYourself
     * @see GaussPoint::updateYourself
     * @see Element::updateInternalState
     */
    virtual void updateInternalState(TimeStep *tStep) { }
    /**
     * Updates element state after equilibrium in time step has been reached.
     * Default implementation updates all integration rules defined by
     * integrationRulesArray member variable. Doing this, all integration points
     * and their material statuses are updated also. All temporary history variables,
     * which now describe equilibrium state are copied into equilibrium ones.
     * The existing internal state is used for update.
     * @param tStep Time step for newly reached state.
     * @see Material::updateYourself
     * @see IntegrationRule::updateYourself
     * @see GaussPoint::updateYourself
     * @see Element::updateInternalState
     */
    virtual void updateYourself(TimeStep *tStep);
    // initialization to state given by initial conditions
    /** Initialization according to state given by initial conditions.
     * Some type of problems may require initialization of state variables
     * stored in integration points (in statuses related to material models)
     * according to initial conditions. Default implementation is not provided.
     * Typically, loop over all integration points, and call to some initialization
     * method of material (with necessary arguments) have to be made.
     * @param timeStepWhenICApply Time step when IC applies.
     */
    virtual void initializeYourself(TimeStep *timeStepWhenICApply) { }

    /**
     * Performs consistency check. This method is called at startup for all elements
     * in particular domain. This method is intended to check data compatibility.
     * Particular element types should test if compatible material
     * and crossSection both with required capabilities are specified.
     * Derived classes should provide their own analysis specific tests.
     * Some printed input if incompatibility is found should be provided
     * (error or warning member functions).
     * Method can be also used to initialize some variables, since
     * this is invoked after all domain components are instanciated.
     * @return Zero value if check fail, otherwise nonzero.
     */
    virtual int checkConsistency() { return 1; }

    // time step initialization (required for some non-linear solvers)
    /**
     * Initializes receivers state to new time step. It can be used also if
     * current time step must be re-started. Default implementation
     * invokes initForNewStep member function for all defined integrationRules.
     * Thus all state variables in all defined integration points are re initialized.
     * @see IntegrationRule::initForNewStep.
     */
    virtual void initForNewStep();
    /**
     * Returns the element geometry type.
     * This information is assumed to be of general interest, but
     * it is required only for some specialized tasks.
     * @return Geometry type of element.
     */
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_unknown; }
    /**
     * Returns the element spatial dimension (1, 2, or 3).
     * This is completely based on the geometrical shape, so a plane in space counts as 2 dimensions.
     * @return Number of spatial dimensions of element.
     */
    virtual int giveSpatialDimension() const;
    /**
     * @return Number of boundaries of element.
     */
    virtual int giveNumberOfBoundarySides() const;
    /**
     * Returns id of default integration rule. Various element types can use
     * different integration rules for implementation of selective or reduced
     * integration of selected components. One particular integration rule from
     * defined integration rules is default. There may be some operations (defined
     * by parent analysis type class) which use default integration rule.
     * @return Id of default integration rule. (index into integrationRulesArray).
     */
    virtual int giveDefaultIntegrationRule() const { return 0; }
    /**
     * Access method for default integration rule.
     * @return Pointer to default integration rule.
     * @see giveDefaultIntegrationRule
     */
    IntegrationRule *giveDefaultIntegrationRulePtr() {
        if ( this->giveNumberOfIntegrationRules() == 0 ) {
            return NULL;
        } else {
            return this->integrationRulesArray [ giveDefaultIntegrationRule() ];
        }
    }
    /// @return Number of integration rules for element.
    int giveNumberOfIntegrationRules() { return this->numberOfIntegrationRules; }
    /**
     * @param i Index of integration rule.
     * @return Requested integration rule.
     */
    IntegrationRule *giveIntegrationRule(int i) { return integrationRulesArray [ i ]; }
    /**
     * Tests if the element implements required extension. ElementExtension type defines
     * the list of all available element extensions.
     * @param ext Extension to be tested.
     * @return Nonzero if extension supported.
     * @see ElementExtension
     */
    virtual int testElementExtension(ElementExtension ext) { return 0; }
    //@}

    /**@name Methods required by some specialized models */
    //@{
    /**
     * Returns the integration point corresponding value in reduced form.
     * @param answer Contain corresponding integration point value, zero sized if not available.
     * @param gp Integration point to check.
     * @param type Determines the type of internal variable.
     * @param tStep Time step.
     * @return Nonzero if o.k, zero otherwise.
     */
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    /**
     * Returns the corresponding integration point value size in reduced form.
     * @param type determines the type of internal variable.
     * @param gp Integration point to check,
     * @return Nonzero if o.k, zero otherwise.
     */
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *gp);
    /**
     * Returns the type of internal variable (scalar, vector, tensor,...).
     * @param type Determines the type of internal variable.
     * @return Type of internal variable.
     */
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    /**
     * Returns the mask of reduced indexes of Internal Variable component .
     * @param answer Mask of full vector size, with components being the indices to reduced form vectors.
     * @param type Determines the internal variable requested.
     * @return Nonzero if o.k, zero otherwise.
     */
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type);

    // characteristic length in gp (for some material models)
    /**
     * Default implementation returns length of element projection into specified direction.
     * @return Element length in given direction.
     */
    virtual double giveLenghtInDir(const FloatArray &normalToCrackPlane);
    /**
     * Returns characteristic length of element in given integration point and in
     * given direction. Required by material models relying on crack-band approach to achieve
     * objectivity with respect to mesh size
     * @param gp Integration point.
     * @param normalToCrackPlane Normal to crack plane.
     * @return Characteristic length of element in given integration point and direction.
     */
    virtual double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane) { return 0.; }
    /**
     * Returns characteristic element size for a given integration point and
     * given direction. Required by material models relying on crack-band approach to achieve
     * objectivity with respect to mesh size.
     * Various techniques can be selected by changing the last parameter.
     * @param gp Integration point.
     * @param normalToCrackPlane Normal to assumed crack plane (some methods use it, some methods recompute it and return the new value).
     * @param method Selection of the specific method to be used.
     * @return Characteristic length of element in given integration point and direction.
     */
    virtual double giveCharacteristicSize(GaussPoint *gp, FloatArray &normalToCrackPlane, ElementCharSizeMethod method) { return giveCharacteristicLenght(gp, normalToCrackPlane); }
    /**
     * Updates internal element state (in all integration points of receiver)
     * before nonlocal averaging takes place. Used by so nonlocal materials,
     * because their response in particular point depends not only on state in this point, but
     * depends also on state in point's neighborhood. Nonlocal quantity is computed as nonlocal
     * average of local quantities. Therefore, before updating integration point state depending on
     * nonlocal quantity (or quantities), local quantities in all integration points must be updated
     * in advance. This function should be  overloaded by derived analysis-specific
     * class, which updates state in all receiver's integration points (using updateBeforeNonlocalAverage
     * member function declared at corresponding analysis specific material base class)
     * depending on driving variable (for example - strain vector in case of structural-analysis elements).
     * @param tStep Time step.
     */
    virtual void updateBeforeNonlocalAverage(TimeStep *tStep) { }
    /**
     * Computes the global coordinates from given element's local coordinates.
     * @param answer Requested global coordinates.
     * @param lcoords Local coordinates.
     * @return Nonzero if successful, zero otherwise.
     */
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    /**
     * Computes the element local coordinates from given global coordinates.
     * Should compute local coordinates even if point is outside element (for mapping purposes in adaptivity)
     * @return Nonzero if point is inside element; zero otherwise.
     * @
     */
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);
    /**
     * Returns local coordinate system of receiver
     * Required by material models with ortho- and anisotrophy.
     * Returns a unit vectors of local coordinate system at element stored row-wise.
     * If local system is equal to global one, set answer to empty matrix and return zero value.
     * @return nonzero if answer computed, zero value if answer is empty, i.e. no transformation is necessary.
     */
    virtual int giveLocalCoordinateSystem(FloatMatrix &answer);

    /**
     * Computes mid-plane normal of receiver at integration point.
     * Only for plane elements in space (3d)  (shells, plates, ....).
     * @param answer The mid plane normal.
     * @param gp The integration point at which to calculate the normal.
     */
    virtual void computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp);
    /**
     * Initializes the internal state variables stored in all IPs according to state in given domain.
     * Used in adaptive procedures.
     * @param oldd Old mesh reference.
     * @param tStep Time step.
     * @return Nonzero if o.k, otherwise zero.
     */
    virtual int adaptiveMap(Domain *oldd, TimeStep *tStep);
    /**
     * Updates the internal state variables stored in all IPs according to
     * already mapped state.
     * @param tStep Time step.
     * @return Nonzero if o.k, otherwise zero.
     */
    virtual int adaptiveUpdate(TimeStep *tStep) { return 1; }
    /**
     * Finishes the mapping for given time step.
     * @param tStep Time step.
     * @return Nonzero if o.k, otherwise zero.
     */
    virtual int adaptiveFinish(TimeStep *tStep);

    /**
     * Local renumbering support. For some tasks (parallel load balancing, for example) it is necessary to
     * renumber the entities. The various FEM components (such as nodes or elements) typically contain
     * links to other entities in terms of their local numbers, etc. This service allows to update
     * these relations to reflect updated numbering. The renumbering function is passed, which is supposed
     * to return an updated number of specified entity type based on old number.
     * @param f Decides the renumbering.
     */
    virtual void updateLocalNumbering(EntityRenumberingFunctor &f);

    /// Integration point evaluator, loops over receiver IP's and calls given function (passed as f parameter) on them. The IP is parameter to function f.
    template< class T >void ipEvaluator( T * src, void ( T :: *f )( GaussPoint * gp ) );
    /// Integration point evaluator, loops over receiver IP's and calls given function (passed as f parameter) on them. The IP is parameter to function f as well as additional array.
    template< class T, class S >void ipEvaluator(T * src, void ( T :: *f )( GaussPoint *, S & ), S & _val);

    //@}


#ifdef __OOFEG
    //
    // Graphics output
    //
    void         drawYourself(oofegGraphicContext &context);
    virtual void drawAnnotation(oofegGraphicContext &mode);
    virtual void drawRawGeometry(oofegGraphicContext &mode) { }
    virtual void drawDeformedGeometry(oofegGraphicContext &mode, UnknownType) { }
    virtual void drawScalar(oofegGraphicContext &context) { }
    virtual void drawSpecial(oofegGraphicContext &context) { }
    // added in order to hide IP element details from oofeg
    // to determine the max and min local values, when recovery does not takes place
    virtual void giveLocalIntVarMaxMin(oofegGraphicContext &context, TimeStep *, double &emin, double &emax) { emin = emax = 0.0; }

    //virtual void  drawInternalState (oofegGraphicContext& context) {}
    /**
     * Returns internal state variable (like stress,strain) at node of element in Reduced form,
     * the way how is obtained is dependent on InternalValueType.
     * The value may be local, or smoothed using some recovery technique.
     * Returns zero if element is unable to respond to request.
     * @param answer Contains result, zero sized if not supported.
     * @param type Determines the internal variable requested (physical meaning).
     * @param mode Determines the mode of variable (recovered, local, ...).
     * @param node Node number, for which variable is required.
     * @param atTime Time step.
     * @return Nonzero if o.k, zero otherwise.
     */
    virtual int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                        int node, TimeStep *atTime);
    /**
     * Returns internal state variable (like stress,strain) at side of element in Reduced form
     * If side is possessing DOFs, otherwise recover techniques will not work
     * due to absence of side-shape functions.
     * @param answer Contains result, zero sized if not supported.
     * @param type Determines the internal variable requested (physical meaning).
     * @param mode Determines the mode of variable (recovered, local, ...).
     * @param side Side number, for which variable is required.
     * @param atTime Time step.
     * @return Nonzero if o.k, zero otherwise.
     */
    virtual int giveInternalStateAtSide(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                        int side, TimeStep *atTime)
    {
        answer.resize(0);
        return 0;
    }

    /// Shows sparse structure
    virtual void showSparseMtrxStructure(CharType mtrx, oofegGraphicContext &gc, TimeStep *atTime) { }
    /// Shows extended sparse structure (for example, due to nonlocal interactions for tangent stiffness)
    virtual void showExtendedSparseMtrxStructure(CharType mtrx, oofegGraphicContext &gc, TimeStep *atTime) { }

#endif

#if defined ( __PARALLEL_MODE ) || defined ( __ENABLE_COMPONENT_LABELS )
    /**
     * @return Receivers globally unique number (label).
     */
    int giveLabel() const { return globalNumber; }
    /**
     * @return Receivers globally unique number.
     */
    int giveGlobalNumber() const { return globalNumber; }
    /**
     * Sets receiver globally unique number.
     * @param num New unique number.
     */
    void setGlobalNumber(int num) { globalNumber = num; }
#endif

#ifdef __PARALLEL_MODE
    /**
     * Return elementParallelMode of receiver. Defined for __Parallel_Mode only.
     */
    elementParallelMode giveParallelMode() const { return parallel_mode; }
    /// Sets parallel mode of element
    void setParallelMode(elementParallelMode _mode) { parallel_mode = _mode; }
 #ifdef __PARALLEL_MODE
    /**
     * Returns the parallel mode for particular knot span of the receiver.
     * The knot span identifies the sub-region of the finite element.
     */
    virtual elementParallelMode giveKnotSpanParallelMode(int) const { return parallel_mode; }
 #endif

    /**
     * Pack all necessary data of element (according to its parallel_mode) integration points
     * into given communication buffer. The corresponding cross section service is invoked, which in
     * turn should invoke material model service for particular integration point. The
     * nature of packed data is material model dependent.
     * Typically, for material of "local" response (response depends only on integration point local state)
     * no data are exchanged. For "nonlocal" constitutive models the send/receive of local values which
     * undergo averaging is performed between local and corresponding remote elements.
     * @param buff communication buffer
     * @param stepN solution step.
     */
    int packUnknowns(CommunicationBuffer &buff, TimeStep *stepN);
    /**
     * Unpack and updates all necessary data of element (according to its parallel_mode) integration points
     * into given communication buffer.
     * @see packUnknowns service.
     * @param buff communication buffer
     * @param stepN solution step.
     */
    int unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *stepN);
    /**
     * Estimates the necessary pack size to hold all packed data of receiver.
     * The corresponding cross section service is invoked, which in
     * turn should invoke material model service for particular integration point. The
     * nature of packed data is material model dependent.  */
    int estimatePackSize(CommunicationBuffer &buff);
    /**
     * Returns partition list of receiver.
     * @return partition array.
     */
    const IntArray *givePartitionList() const { return & partitions; }
    /**
     * Sets partition list of receiver
     */
    void setPartitionList(IntArray &pl) { partitions = pl; }
    /**
     * Returns the weight representing relative computational cost of receiver
     * The reference element is triangular plane stress element with
     * linear approximation, single integration point and linear isotropic material.
     * Its weight is equal to 1.0
     * Default implementation computes average computational cost of cross section model (this include material as well)
     * and multiplies it by element type weight (obtained by giveRelativeSelfComputationalCost())
     * The other elements should compare to this reference element.
     */
    virtual double predictRelativeComputationalCost();
    /**
     * Returns the weight representing relative computational cost of receiver
     * The reference element is triangular plane stress element.
     * Its weight is equal to 1.0
     * The other elements should compare to this reference element.
     */
    virtual double giveRelativeSelfComputationalCost() { return 1.0; }
    /**
     * Returns the relative redistribution cost of the receiver
     */
    virtual double predictRelativeRedistributionCost() { return 1.0; }
#endif

public:
    /// Returns array containing load numbers of loads acting on element
    IntArray *giveBodyLoadArray();
    /// Returns array containing load numbers of boundary loads acting on element.
    IntArray *giveBoundaryLoadArray();

    // Overloaded methods:
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual const char *giveClassName() const { return "Element"; }
    virtual classType giveClassID() const { return ElementClass; }

protected:
    /**
     * Initializes the array of integration rules and numberOfIntegrationRules member variable.
     * Element can have multiple integration rules for different tasks.
     * For example structural element family class uses this feature to implement
     * transparent support for reduced and selective integration of some strain components.
     * Must be defined by terminator classes.
     * @see IntegrationRule
     */
    virtual void computeGaussPoints() { }
};

template< class T >void
Element :: ipEvaluator( T *src, void ( T :: *f )( GaussPoint * gp ) )
{
    int ir, ip, nip;
    GaussPoint *gp;

    for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
        nip = integrationRulesArray [ ir ]->getNumberOfIntegrationPoints();
        for ( ip = 0; ip < nip; ip++ ) {
            gp = integrationRulesArray [ ir ]->getIntegrationPoint(ip);
            ( src->*f )(gp);
        }
    }
}

template< class T, class S >void
Element :: ipEvaluator(T *src, void ( T :: *f )( GaussPoint *, S & ), S &_val)
{
    int ir, ip, nip;
    GaussPoint *gp;

    for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
        nip = integrationRulesArray [ ir ]->getNumberOfIntegrationPoints();
        for ( ip = 0; ip < nip; ip++ ) {
            gp = integrationRulesArray [ ir ]->getIntegrationPoint(ip);
            ( src->*f )(gp, _val);
        }
    }
}
} // end namespace oofem
#endif //element_h
