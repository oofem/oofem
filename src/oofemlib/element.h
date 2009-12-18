/* $Header: /home/cvs/bp/oofem/oofemlib/src/element.h,v 1.27 2003/04/06 14:08:24 bp Exp $ */
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


//   *********************
//   *** CLASS ELEMENT ***
//   *********************


#ifndef element_h
#define element_h

#include "femcmpnn.h"
#include "domain.h"
#include "flotmtrx.h"

#include "primaryfield.h"
#include "integrationrule.h"
#include "elementheaders.h"
#include "equationid.h"
#include "valuemodetype.h"
#include "internalstatevaluetype.h"
#include "elementextension.h"
#include "entityrenumberingscheme.h"
#include "unknowntype.h"
#include "geometry.h"
#include "unknownnumberingscheme.h"

namespace oofem {

class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;
class CrossSection;
class ElementSide;
class CommunicationBuffer;

#ifdef __PARALLEL_MODE
class CommunicationBuffer;

/**
 * In parallel mode, this type indicates the mode of element.
 * <UL>
 * <LI>
 * Element_local mode - Element is local, there are no contribution from other domains to this Element.</LI>
 * <LI>
 * Element_shared mode - Element is shared by neighbouring partitions - not implemented.
 * In the case of eleemnt cut mode, the cutted element is local on all partitions sharing it.
 * Some of such element  nodes are local and some are remote. The local nodes are completely
 * surronded by local element on particular partition.</LI>
 * <LI>
 * Element_remote - Element in active domain is only mirror of some remote
 * Element. </LI>
 * </UL>
 */
enum elementParallelMode {
    Element_local,
    // Element_shared,
    Element_remote
};

#endif


/**
 * Abstract base class for all finite elements. Derived clases should be  base
 * clases for specific analysis type (for example base class for structural analysis,
 * thermal analysis or magnetostatics one). These derived classes then declare
 * analysis-specific part of interface and they provide default implementation
 * for these methods.
 * This abstract class declares (and possibly implements) general data and methods
 * common to all element types. General methods for obtaining characteristic vectors,
 * matrices and values are introduced and should be used instead of calling directly
 * specific member functions (these must be overloaded by derived analysis-specific
 * clases in order to invoke proper method acording to type of component requested).
 * In general, element class memeber functions are called in following cases:
 * <UL>
 * <LI>
 * When Engineering model assembles governing equation(s) from element's contributions.
 * Typically when assembles some characteristic matrix of analyzed structure (method assemble),
 * asks each element for its code numbers and for corresponding characteristic matrix, vector or
 * value. This will typically uses above mentioned general method for obtaining characteristic
 * matrices, vectors and values from elements.</LI>
 * <LI>
 * When Engineering model has some characteristic matrix stored in sparse form, then is necessary
 * to build internal profile of such sparse matrix. Engineering model then typically calls
 * buildInternalStructure member function of sparse matrix. This method then requests element code
 * numbers from elements and builds its internal profile. This could look strange, because
 * class sparse matrix should know "something about eleemnts", but only sparse matrix knows, how
 * to build its internal structure (this may require one or more loops over elements code number).</LI>
 * <LI>
 * When element computes its contribution, then it comunicates with its cross-section (and cross
 * section with corresponding material). For some cross-section or matrial models some further
 * communication between these classes and element may be necessary (for example in cases of
 * layered cross section model strains in each layer has to be evaluated from "integrated"
 * strains (integrated means element-like strains including curvatures in case of beams and plates)
 * by element, in cases of non local material and many other examples can be found).</LI>
 * </UL>
 * There are some general rules, that programmer must take into account.
 * <UL>
 * <LI>
 * Element stores the numbers of its dofnamagers in dofManArray.
 * These include nodes, element sides and internal DOFs that are not condensed at element level.
 * Their order and meaning are determined by element definition.
 * Local ordering of dofs for particular element is determined by local numbering of
 * dofnamagers and their corresponding dofs. DOFS necessary for particular node/side is specified using node/side dof mask.
 * Local DOF ordering must be taken into account when assembling various local characteristic
 * vectors and matrices.</LI>
 *
 *
 * </UL>
 *
 */
class Element : public FEMComponent, public ElementGeometry
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
     * Array containing indexes of loads (bodyloads and boundary loads are kept separately),
     * that apply on receiver.
     */
    IntArray bodyLoadArray, boundaryLoadArray;   /* load index into
                                                  * LoadList  */
    /// Number of integration rules used by receiver.
    int numberOfIntegrationRules;
    /** List of integration rules of receiver (each integration rule contains associated
     * integration points also). This list should contain only such integration rules,
     * that are used to integrate results depending on load time history. For all integration
     * points in these rules, history variables are stored and updated.
     * For integrations, where history stored in gauss points is not necessary
     * (mass matrix integration) and different integration rule is needed, one should preferably
     * use temporarly created integration rule.
     */
    IntegrationRule **integrationRulesArray;
    /// Array of code numbers of element.
    IntArray *locationArray;

    /**
     * Transformation material matrix, used in orthotropic and anisotropic materials, global->local transformation
     */
    FloatMatrix matLocalCS;


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
     * List of partition sharing the shared eleemnt or
     * remote partion containing remote element counterpart.
     */
    IntArray partitions;
#endif

public:
    /**
     * Constructor. Creates an element with number n belonging to domain aDomain.
     * @param n Element's number
     * @param aDomain Pointer to the domain to which element belongs.
     */
    Element(int n, Domain * aDomain);   // constructors
    /// Virtual destructor.
    virtual ~Element();                // destructor

    // assembly
    /**@name Methods refering to code numbers */
    //@{
    /**
     * Returns the location array (array of code numbers) of receiver for given numbering scheme.
     * Results are cached at receiver for default scheme in locationArray attribute.
     * The UnknownType parameter allows to distinguish between several possible governing equations, that
     * can be numbered separately. The default implementation assumes that location array will be assembled only for
     * one UnknownType value, and this array is cached on element level.
     */
    virtual void giveLocationArray(IntArray &locationArray, EquationID, const UnknownNumberingScheme &s) const;
    /**
     * Invalidates location array in receiver. Each element stores its copy of location array(s), in order
     * to avoid time consuming assembly of code numbers every time when requested. Some enginnering models
     * may support dynamic changes of static system (generally of boundary conditions during analysis),
     * then these models use this function to invalidate location array after finishing time step,
     * to enforce elements to update they code numbers, which may change. Changes of static system will lead
     * to different number of equations, which requires special attention (for example internal structure
     * of sparse matrices  must be reinitialized and characteristic vectors and matrices has to be
     * newly assembled).
     */
    void                  invalidateLocationArray();
    /**
     *  Returns the number of internal element dofs
     **/
    virtual int giveNumberOfDofs() { return 0; }
    /**
     * Returns i-th DOF of the receiver
     **/
    virtual Dof *giveDof(int i) {
        _error2("No such DOF available on Element %d", number);
        return NULL;
    }
    //@}
    /**@name General methods for obtaining element contributions
     * Note: These member functions have to  be overloaded by derived analysis-specific
     * clases in order to invoke proper method acording to type of component requested.
     * Derived classes from these analysis-specific clases should not modify these functions.
     */
    //@{
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
    virtual void     giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep);
    /**
     * Computes characteristic vector of receiver of requested type in given Timestep.
     * @param answer requested characteristic vector
     * @param type   id  of characteristic component requested.
     * @param mode   determines mode of answer.
     * @param tStep  time step when answer is computed.
     * @see CharType and ValueModeType types.
     * @return If element has no capability to compute requsted type of caharacteristic vector
     * error function is invoked.
     */
    virtual void     giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *);
    /**
     * Computes characteristic value of receiver of requested type in given Timestep.
     * @param mtrx   id  of characteristic component requested.
     * @param tStep  time step when answer is computed.
     * @see CharType type.
     * @return requested value. If element has no capability to compute requsted type of caharacteristic value
     * error function is invoked.
     */
    virtual double giveCharacteristicValue(CharType, TimeStep *);
    //@}
    /*
     *  // I have tryed name of next fuction  to be again giveCharacteristicVector
     *  // but there are strong disadvantages of it. (It may be compiler bug). This happen,
     *  // when you have the functions with same names and different arguments and both are
     *  // virtual.  When you overload one of then at parent level you must redefine also
     *  // the second one again. Without it, it does not work properly.
     *  virtual FloatArray*           GiveCharacteristicGeomVector (CharType, GaussPoint *);
     */
    //virtual MaterialMode          giveMaterialMode() {return _Unknown;}
    // vector of nodal unknowns

    /**
     * Gives transformation matrix Global=T*Local for material orientation on element. Only if defined by mlcs on element
     * @param answer transformation matrix 3x3
     */
    virtual void giveMatLocalCS(FloatMatrix &answer);


    /**@name General element functions */
    //@{
    /**
     * Returns local vector of unknowns. Local vector of unknowns is extracted from
     * engineering model global unknowns vector (if specific dof has assigned
     * corresponding equation number) and from boundary conditions (if dof has active boundary
     * and possibly initial condition). Because unknowns are obtained from engineering
     * model, this must support queries for given unknown and unknown mode. Also
     * engineering model must not be able to complete request for any TimeStep, because
     * keeping history of all unknowns may became unpossible. But generally all enginnering models
     * should be able return supported unknowns at current and previous time step. Consult
     * reference manual for particular engineering model.
     *
     * @param type Identifies unknown type (eg. displacement or temperature vector).
     * @param u    Identifies mode of unknown (eg. total value or velocity of unknown).
     * @param stepN Time step, when vector of unknowns is requested.
     * @param answer Local vector of unknowns.
     */
    virtual void          computeVectorOf(EquationID type, ValueModeType u,
                                          TimeStep *stepN, FloatArray &answer);
    /**
     * Returns local vector of unknowns. Local vector of unknowns is extracted from
     * given field and from boundary conditions (if dof has active boundary
     * and possibly initial condition). Because unknowns are obtained from given field
     * model, this must support queries for given unknown.
     * @param field source field (eg. displacement or temperature vector).
     * @param stepN Time step, when vector of unknowns is requested.
     * @param answer Local vector of unknowns.
     */
    virtual void          computeVectorOf(PrimaryField &field, ValueModeType u, TimeStep *stepN, FloatArray &answer);
    /**
     * Returns local vector of prescribed unknowns. Local vector of prescribed unknowns is
     * extracted from nodal (and side - if they hold unknowns) boundary conditions.
     * @param u    Identifies mode of unknown (eg. total values or velocity of unknown).
     * @param stepN Time step, when vector of prescribed unknowns is requested.
     * @param answer Local vector of prescribed unknowns. If unknown is not prescibed,
     * zero value is placed on position of free dof.
     */
    void                  computeVectorOfPrescribed(EquationID ut, ValueModeType type, TimeStep *stepN, FloatArray &answer);
    //void                  computeVectorOfPrescribed (UnknownType type, ValueModeType u,
    //                         TimeStep* stepN, FloatArray& answer) ;
    /**
     * Computes or simply returns total number of element's local dofs.
     * Must be defined by particular element.
     */
    virtual int            computeNumberOfDofs(EquationID ut) { return 0; }
    /**
     * Computes the total number of element's global dofs (the size of global element contribution).
     */
    virtual int            computeGlobalNumberOfDofs(EquationID ut);
    /**
     * Returns dofmanager dof mask for node. This mask defines the dofs which are used by element
     * in node. Mask influences the code number ordering for particular node. Code numbers are
     * ordered acording to node order and dofs belonging to particular node are ordered
     * according to this mask. If element requests dofs using node mask which are not in node
     * then error is generated. This masking allows node to be shared by different elements with
     * different dofs in same node. Elements local code numbers are extracted from node using
     * this mask. Must be defined by particular element.
     *
     * @param inode Mask is computed for local dofmanager with inode number.
     * @param ut unknown type (support for several independent numberings within problem)
     * @param answer mask for node.
     */
    virtual void           giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
    { answer.resize(0); }
    /**
     * Returns element dof mask for node. This mask defines the dof ordering of the element interpolation.
     * Must be defined by particular element.
     *
     * @param answer dof mask for receiver.
     */
    virtual void           giveElementDofIDMask(EquationID ut, IntArray &answer) const
    { answer.resize(0); }
    /**
     * Returns volume related to given integration point. Used typically in subroutines,
     * that perform integration over element volume. Should be implemented by particular
     * elements.
     * @param gp related volume for given integration point will be computed
     * @return related volume
     * @see GaussPoint class
     */
    virtual double        computeVolumeAround(GaussPoint *gp) { return 0.; }

    // data management
    ///Returns (global) number of i-th dofmanager of element
    int                    giveDofManagerNumber(int i) const { return dofManArray.at(i); }
    ///Rerurns reference to the i-th dofmanager of element.
    DofManager *giveDofManager(int i) const;
    /**Returns reference to the i-th node of element.
     * Default implementation returns i-th dofmanager of element converted to
     * Node class (check is made).
     */
    virtual Node *giveNode(int i) const;
    /**Returns reference to the i-th side  of element.
     * Default implementation returns i-th dofmanager of element converted to
     * ElementSide class (check is made).
     */
    virtual ElementSide *giveSide(int i) const;
    /// Returns reference to the associated geometry of element
    Geometry *giveGeometry() { return NULL; } // rch
    /// Returns interpolation of Element
    virtual FEInterpolation *giveInterpolation() { return NULL; } // rch
    ///Returns reference to the associated material of element.
    Material *giveMaterial();
    ///Returns reference to the associated crossSection of element.
    CrossSection *giveCrossSection();
    ///Sets the material of receiver
    void setMaterial(int matIndx) { this->material = matIndx; }
    ///Sets the cross section model of receiver
    void setCrossSection(int csIndx) { this->crossSection = csIndx; }

    ///Returns number of dofmanagers of receiver
    int                   giveNumberOfDofManagers() { return numberOfDofMans; }
    /**Returns number of nodes of receiver.
     * Default implementation returns number of dofmanagers of element
     */
    virtual int                   giveNumberOfNodes() { return numberOfDofMans; }
    /** Sets receiver dofManagers */
    void setDofManagers(const IntArray &_dmans);

    /** Sets integration rules */
    void setIntegrationRules(AList< IntegrationRule > *irlist); // rch
    /**
     * Returns integration domain for receiver, used to initialize
     * integration point over receiver volume. Must be specialized.
     */
    virtual integrationDomain giveIntegrationDomain() { return _Unknown_integrationDomain; }
    /**
     * Assembles the code numbers of given integration element (sub-patch)
     * This is done by obtaining list of nonzero shape functions and
     * by collecting the code numbers of nodes corresponding to these
     * shape functions
     * @returns returns nonzero if integration rule code numbers differ from element code numbers
     */
    virtual int giveIntegrationRuleLocalCodeNumbers(IntArray &answer, IntegrationRule *ie, EquationID ut)
    { return 0; }

    ///Returns number of sides (which have unknown dofs) of receiver
    //int                   giveNumberOfSides () {return numberOfSides;}
    /// Returns the corresponding element region. Currently corresponds to cross section model number.
    int                   giveRegionNumber();
    ///Initializes receiver acording to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);
    ///Performs a post initialization steps
    void         postInitialize();
    /**
     * Stores receiver state to output stream.
     * @exception throws an ContextIOERR exception if error encountered
     */
    contextIOResultType   saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores the receiver state previously written in stream.
     * @exception throws an ContextIOERR exception if error encountered
     */
    contextIOResultType   restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    // time step termination
    /**
     * Prints output of receiver to stream, for given time step.
     * Corresponding function for element gauss points is invoked
     * (gaussPoint::printOutputAt).
     */
    void                  printOutputAt(FILE *, TimeStep *);
    /**
     * Updates element state corresponding to newly reached solution.
     * Default is empty, derived classes should force the update of internal integration point values
     * acording to newly reached state.
     * @see Element::updateYourself
     */
    virtual void          updateInternalState(TimeStep *tStep) { }
    /**
     * Updates element state after equlibrium in time step has been reached.
     * Default implementation updates all integration rules defined by
     * integrationRulesArray member variable. Doing this, all integration points
     * and their material statuses are updated also. All temporary history variables,
     * which now describe equlibrium state are copied into equilibrium ones.
     * The existing internal state is used for update.
     * @see Material::updateYourself
     * @see IntegrationRule::updateYourself
     * @see gaussPoint::updateYourself
     * @see Element::updateInternalState
     */
    virtual void          updateYourself(TimeStep *tStep);
    // initialization to state given by initial conditions
    /** Initialization according to state given by initial conditions.
     * Some type of problems may require initialization of state variables
     * stored in integration points (in statuses related to material models)
     * according to initial conditions. Defaul implementation is not provided.
     * Typically, loop over all integration points, and call to some initialization
     * method of material (with necessary agruments) have to be made.
     */
    virtual void          initializeYourself(TimeStep *timeStepWhenICApply) { }

    // consistency check
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
     * @return zero value if check fail, othervise nonzero is returned.
     */
    virtual int    checkConsistency() { return 1; }

    // time step initialization (required for some non-linear solvers)
    /**
     * Initializes receivers state to new time step. It can be used also if
     * current time step must be re-started. Default implementation
     * invokes initForNewStep member function for all defined integrationRules.
     * Thus all state variables in all defined integration points are re initialized.
     * @see IntegrationRule::initForNewStep.
     */
    virtual void initForNewStep();

    // definition
    //       Element*              typed () ;
    /**
     * Returns a newly allocated element, with type depending on parameter.
     * Calls directly CreateUsrDefElementOfType global function to allocate
     * new instance of element of given type.
     * Calls global function CreateUsrDefElementOfType for creating appropriate
     * element instance. This function must be implemented by user.
     * @param aClass string with element name
     * @return newly allocated element with required type.
     * @see CreateUsrDefElementOfType function.
     */
    Element *ofType(char *aClass);
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Element"; }
    /** Returns classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    classType                giveClassID() const { return ElementClass; }
    /** Returns the element geometry type.
     * This information is assumed to be of general interest, but
     * it is required only for some specialized tasks.
     */
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_unknown; }
    /** Returns the element spacial dimension irrespectably to its position in space
     */
    int giveSpatialDimension(void);
    /** Returns the number of boundaries of dimension equal to element spatial dimension - 1
     */
    int giveNumberOfBoundarySides(void);
    /**
     * Returns id of default integration rule. Various element types can use
     * different integration rules for implementation of selective or reduced
     * integration of selected components. One particular integration rule from
     * defined integration rules is default. There may be some operations (defined
     * by parent analysis type class) which use default integration rule.
     * @return id of default integration rule. (index into integrationRulesArray).
     */
    virtual int giveDefaultIntegrationRule() { return 0; }
    IntegrationRule *giveDefaultIntegrationRulePtr() { return integrationRulesArray [ giveDefaultIntegrationRule() ]; }
    /**
     * Tests if the element implements required extension. ElementExtension type defines
     * the list of all available element extensions.
     * @param ext tested extension id
     * @return nonzero if extension supported.
     * @see ElementExtension
     */
    virtual int testElementExtension(ElementExtension ext) { return 0; }
    //@}
    /**@name Methods required by some specialized models */
    //@{
    //
    // functions required for special functionality
    //
    /**
     * Returns the integration point corresponding value in REDUCED form.
     * @param answer contain corresponding ip value, zero sized if not available.
     * @param aGaussPoint integration point
     * @param type determines the type of internal variable
     * @returns nonzero if ok, zero if var not supported
     */
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    /**
     * Returns the corresponding integration point value size in Reduced form.
     * @param type determines the type of internal variable
     * @param mat corresponding material
     * @return nonzero if o.k, zero otherwise
     */
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *);
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
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type);

    // characteristic length in gp (for some material models)
    /** Returns element length in given direction. Default implementation returns length
     * of element  projection into specified direction
     */
    virtual double        giveLenghtInDir(const FloatArray &normalToCrackPlane);
    /**
     * Returns characteristic length of element in given integration point and in
     * given direction. Required by material models relying on crack-band approach to achieve
     * objectivity with respect to mesh size
     */
    virtual double        giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane)
    { return 0.; }
    /**
     * Updates internal element state (in all integration points of receiver)
     * before nonlocal averaging takes place. Used by so nonlocal materials,
     * because their response in particular point depends not only on state in this point, but
     * depends also on state in point's neighbourhood. Nonlocal quatity is computed as nonlocal
     * average of local quantities. Therefore, before updating integration point state depending on
     * nonlocal quantity (or quantities), local quantities in all integration points must be updated
     * in advance. This function should be  overloaded by derived analysis-specific
     * class, which updates state in all receiver's integration points (using updateBeforeNonlocalAverage
     * member function deklared at corresponding analysis specific material base class)
     * depending on driving variable (for example - strain vector in case of structural-analysis elements).
     */
    virtual void updateBeforeNonlocalAverage(TimeStep *atTime) { }
    /**
     * Computes the global coordinates from given element's local coordinates.
     * @returns nonzero if successful; zero otherwise
     */
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) {
        answer.resize(0);
        return 0;
    }
    /**
     * Computes the element local coordinates from given global coordinates.
     * Should compute local coordinates even if point is outside element (for mapping purposes in adaptivity)
     * @returns nonzero if point is inside element; zero otherwise
     */
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords) {
        answer.resize(0);
        return 0;
    }
    // returns a unit vectors of local coordinate system at element
    // stored rowwise (mainly used by some materials with ortho and anisotrophy)
    // if local c.s == global c.s returns NULL
    /**
     * Returns local coordinate system of receiver. Required by material models with
     * ortho and anisotrophy. If local system is equal to global one, can set answer to empty mtrx
     * and return zero value.
     * @return nonzero if answer computed, zero value if answer is empty, i.e. no transformation is necessary.
     */
    virtual int  giveLocalCoordinateSystem(FloatMatrix &answer) {
        answer.beEmptyMtrx();
        return 0;
    }
    // mid-plane normal at gaussPoint - for materials with orthotrophy
    // valid only for plane elements in space (3d)  (shells, plates, ....)
    /**
     * Computes mid-plane normal of receiver at integration point.
     * Only for plane elements in space (3d)  (shells, plates, ....).
     */
    virtual FloatArray *ComputeMidPlaneNormal(GaussPoint *);
    /**
     * Initializes the internal state variables stored in all IPs according to
     * stete in given domain. Used in adaptive procedures.
     * @param oldd old mesh reference
     * @param tStep time step
     * @return nonzero if o.k.
     */
    virtual int adaptiveMap(Domain *oldd, TimeStep *tStep);
    /**
     * Updates the internal state variables stored in all IPs according to
     * already mapped state.
     * @param oldd old mesh reference
     * @param tStep time step
     * @return nonzero if o.k.
     */
    virtual int adaptiveUpdate(TimeStep *tStep) { return 1; }
    /**
     * Finishes the mapping for given time step.
     */
    virtual int adaptiveFinish(TimeStep *tStep);

    /**
     * Local renumbering support. For some tasks (parallel load balancing, for example) it is necessary to
     * renumber the entities. The various fem components (such as nodes or elements) typically contain
     * links to other entities in terms of their local numbers, etc. This service allows to update
     * these relations to reflext updated numbering. The renumbering funciton is passed, which is supposed
     * to return an updated number of specified entyty type based on old number.
     */
    virtual void updateLocalNumbering(EntityRenumberingFunctor &f);

    /// Integration point evaluator, loops over receiver IP's and calls given function (passed as f parameter) on them. The IP is parameter to function f.
    template< class T > void ipEvaluator( T *src, void ( T :: *f )( GaussPoint *gp ) );
    /// Integration point evaluator, loops over receiver IP's and calls given function (passed as f parameter) on them. The IP is parameter to function f as well as additional array.
    template< class T, class S > void ipEvaluator(T *src, void ( T :: *f )( GaussPoint *, S & ), S &_val);



    //@}


#ifdef __OOFEG
    //
    // Graphics output
    //
    void          drawYourself(oofegGraphicContext &context);
    virtual void  drawAnnotation(oofegGraphicContext &mode);
    virtual void  drawRawGeometry(oofegGraphicContext &mode) { }
    virtual void  drawDeformedGeometry(oofegGraphicContext &mode, UnknownType) { }
    virtual void  drawScalar(oofegGraphicContext &context) { }
    virtual void  drawSpecial(oofegGraphicContext &context) { }
    // added in order to hide IP element details from oofeg
    // to determine the max and min local values, when recovery does not takes place
    virtual void  giveLocalIntVarMaxMin(oofegGraphicContext &context, TimeStep *, double &emin, double &emax) { emin = emax = 0.0; }

    //virtual void  drawInternalState (oofegGraphicContext& context) {}
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
    /**
     * Returns internal state variable (like stress,strain) at side of element in Reduced form
     * (if side is possing DOFs, otherwise recover techniques will not work
     * due to absence of side-shape functions)
     * @param answer contains result, zero sized if not supported
     * @param type determines the internal variable requested (physical meaning)
     * @param mode determines the mode of variable (recovered, local, ...)
     * @param node node number, for which variable is required
     * @param atTime time step
     * @return nonzero if o.k, zero otherwise
     */
    virtual int   giveInternalStateAtSide(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                          int side, TimeStep *atTime)
    {
        answer.resize(0);
        return 0;
    }

    /// Shows sparse structure
    virtual void showSparseMtrxStructure(CharType mtrx, oofegGraphicContext &gc, TimeStep *atTime) { }
    /// Shows extended sparse structure (for example, due to nonlocal interactios for tangent stiffness)
    virtual void showExtendedSparseMtrxStructure(CharType mtrx, oofegGraphicContext &gc, TimeStep *atTime) { }

#endif

#if defined ( __PARALLEL_MODE ) || defined ( __ENABLE_COMPONENT_LABELS )
    /**
     * Returns receiver globally unique number (label).
     */
    int giveLabel() const { return globalNumber; }
    /**
     * Returns receiver globally unique number.
     */
    int giveGlobalNumber() const { return globalNumber; }
    /**
     * Sets receiver globally unique number.
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
    /**
     * Pack all necessary data of element (according to its parallel_mode) integration points
     * into given communication buffer. The corresponding cross section service is invoked, which in
     * turn should invoke material model service for particular integration point. The
     * nature of packed data is material model dependent.
     * Typically, for material of "local" response (response depeneds only on integration point local state)
     * no data are exchanged. For "nonlocal" constitutive models the send/receive of local values which
     * undergo averaging is performed between local and corressponding remote elements.
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
    //
    // public methods
    //

    ///Returns array containing load numbers of loads acting on element
    IntArray *giveBodyLoadArray();
    ///Returns array containing load numbers of boundary loads acting on element.
    IntArray *giveBoundaryLoadArray();

protected:
    //
    // protected methods
    //
    // interpolation, numerical integration
    /**
     * Initializes the array of integration rules and numberOfIntegrationRules member variable.
     * Element can have multiple integration rules for different tasks.
     * For example structural element family class uses this feature to implement
     * transparent support for reduced and selective integration of some strain components.
     * Must be defined by terminator clases.
     * @see IntegrationRule class.
     */
    virtual void          computeGaussPoints() { }
};

template< class T > void
Element :: ipEvaluator( T *src, void ( T :: *f )( GaussPoint *gp ) )
{
    int ir, ip, nir, nip;
    GaussPoint *gp;

    for ( ir = 0; ir < numberOfIntegrationRules; ir++ ) {
        nip = integrationRulesArray [ ir ]->getNumberOfIntegrationPoints();
        for ( ip = 0; ip < nip; ip++ ) {
            gp = integrationRulesArray [ ir ]->getIntegrationPoint(ip);
            ( src->*f )(gp);
        }
    }
}

template< class T, class S > void
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








