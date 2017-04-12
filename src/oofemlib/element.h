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

#ifndef element_h
#define element_h

#include "femcmpnn.h"
#include "error.h"
#include "chartype.h"
#include "domain.h"
#include "floatmatrix.h"
#include "integrationdomain.h"
#include "materialmode.h"
#include "elementgeometrytype.h"
#include "valuemodetype.h"
#include "internalstatemode.h"
#include "internalstatetype.h"
#include "elementextension.h"
#include "entityrenumberingscheme.h"
#include "matresponsemode.h"
#include "unknowntype.h"
#include "integrationrule.h"
#include "dofiditem.h"
#include "floatarray.h"

#include <cstdio>
#include <vector>
#include <memory>

///@name Input fields for general element.
//@{
#define _IFT_Element_mat "mat"
#define _IFT_Element_crosssect "crosssect"
#define _IFT_Element_nodes "nodes"
#define _IFT_Element_bodyload "bodyloads"
#define _IFT_Element_boundaryload "boundaryloads"
#define _IFT_Element_lcs "lcs"
#define _IFT_Element_partitions "partitions"
#define _IFT_Element_remote "remote"
#define _IFT_Element_activityTimeFunction "activityltf"
#define _IFT_Element_nip "nip"
//@}

namespace oofem {
class TimeStep;
class Node;
class Material;
class IntegrationRule;
class GaussPoint;
class FloatMatrix;
class IntArray;
class CrossSection;
class ElementSide;
class FEInterpolation;
class Load;
class BoundaryLoad;
class BodyLoad;
class SurfaceLoad;
class EdgeLoad;
class PrimaryField;
class UnknownNumberingScheme;

/**
 * In parallel mode, this type indicates the mode of element.
 * In the case of element cut mode, the cut element is local on all partitions sharing it.
 * Some of such element nodes are local and some are remote. The local nodes are completely
 * surrounded by local element on particular partition.
 */
enum elementParallelMode {
    Element_local, ///< Element is local, there are no contributions from other domains to this element.
    Element_remote, ///< Element in active domain is only mirror of some remote element.
};


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
class OOFEM_EXPORT Element : public FEMComponent
{
protected:
    /// Number of dofmanagers
    int numberOfDofMans; ///@todo We should remove this parameter. It's redundant (and therefore possibly wrong) since the dofManArray stores it's length internally. 
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
    /**
     * List of integration rules of receiver (each integration rule contains associated
     * integration points also). This list should contain only such integration rules,
     * that are used to integrate results depending on load time history. For all integration
     * points in these rules, history variables are stored and updated.
     * For integrations, where history stored in gauss points is not necessary
     * (mass matrix integration) and different integration rule is needed, one should preferably
     * use temporarily created integration rule.
     */
    std::vector< std :: unique_ptr< IntegrationRule > > integrationRulesArray;

    /// Transformation material matrix, used in orthotropic and anisotropic materials, global->local transformation
    FloatMatrix elemLocalCS;

    /// Element activity time function. If defined, nonzero value indicates active receiver, zero value inactive element.
    int activityTimeFunction;

    /**
     * In parallel mode, globalNumber contains globally unique DoFManager number.
     * The component number, inherited from FEMComponent class contains
     * local domain number.
     */
    int globalNumber;

    /**
     * Number of integration points as specified by nip.
     */
    int numberOfGaussPoints;

    /// Determines the parallel mode of the element
    elementParallelMode parallel_mode;

    /**
     * List of partition sharing the shared element or
     * remote partition containing remote element counterpart.
     */
    IntArray partitions;

public:
    /**
     * Constructor. Creates an element with number n belonging to domain aDomain.
     * @param n Element's number
     * @param aDomain Pointer to the domain to which element belongs.
     */
    Element(int n, Domain * aDomain);
    Element(const Element& src) = delete;
    Element &operator = (const Element &src) = delete;
    /// Virtual destructor.
    virtual ~Element();

    /**@name Methods referring to code numbers */
    //@{
    /**
     * Returns the location array (array of code numbers) of receiver for given numbering scheme.
     * Results are cached at receiver for default scheme in locationArray attribute.
     */
    void giveLocationArray(IntArray &locationArray, const UnknownNumberingScheme &s, IntArray *dofIds = NULL) const;
    void giveLocationArray(IntArray &locationArray, const IntArray &dofIDMask, const UnknownNumberingScheme &s, IntArray *dofIds = NULL) const;
    /**
     * Returns the location array for the boundary of the element.
     * Only takes into account nodes in the bNodes vector.
     */
    virtual void giveBoundaryLocationArray(IntArray &locationArray, const IntArray &bNodes, const UnknownNumberingScheme &s, IntArray *dofIds = NULL);
    virtual void giveBoundaryLocationArray(IntArray &locationArray, const IntArray &bNodes, const IntArray &dofIDMask, const UnknownNumberingScheme &s, IntArray *dofIds = NULL);
    /**
     * @return Number of DOFs in element.
     */
    virtual int giveNumberOfDofs() { return 0; }
    /**
     * @return Number of internal DOF managers of element.
     */
    virtual int giveNumberOfInternalDofManagers() const { return 0; }
    /**
     * Returns i-th internal element dof manager of the receiver
     * @param i Internal number of DOF.
     * @return DOF number i.
     */
    virtual DofManager *giveInternalDofManager(int i) const {
        OOFEM_ERROR("No such DOF available on Element %d", number);
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
     * @param answer Requested characteristic matrix (stiffness, tangent, ...).
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

    /**
     * @name General methods for computing the contribution from loads
     */
    //@{
    /**
     * Computes the contribution of the given body load (volumetric).
     * @param answer Requested contribution of load.
     * @param load   Load to compute contribution from.
     * @param type   Type of the contribution.
     * @param mode   Determines mode of answer.
     * @param tStep  Time step when answer is computed.
     */
    virtual void computeLoadVector(FloatArray &answer, BodyLoad *load, CharType type, ValueModeType mode, TimeStep *tStep);
    /**
     * Computes the contribution of the given load at the given boundary surface in global 
     * coordinate system. 
     * In general, the answer should include only relevant DOFs at the edge.
     * The related is giveBoundaryLocationArray method, which should return 
     * corresponding code numbers.
     * @note Elements which do not have an contribution should resize the vector to be empty.
     * @param answer Requested contribution of load.
     * @param load Load to compute contribution from.
     * @param boundary Boundary number.
     * @param type Type of the contribution.
     * @param mode Determines mode of answer.
     * @param tStep Time step when answer is computed.
     * @param global if true (default) then contribution is in global c.s., when false then contribution is in element local c.s.
     */
    virtual void computeBoundarySurfaceLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, bool global=true);
    /**
     * Computes the tangent contribution of the given load at the given boundary.
     * @note Elements which do not have an contribution should resize the vector to be empty.
     * @param answer Requested contribution of load.
     * @param load Load to compute contribution from.
     * @param boundary Surface number.
     * @param rmode Mode of the contribution.
     * @param tStep Time step when answer is computed.
     */
    virtual void computeTangentFromSurfaceLoad(FloatMatrix &answer, SurfaceLoad *load, int boundary, MatResponseMode rmode, TimeStep *tStep);
    /**
     * Computes the tangent contribution of the given load at the given boundary.
     * @note Elements which do not have an contribution should resize the vector to be empty.
     * @param answer Requested contribution of load.
     * @param load Load to compute contribution from.
     * @param boundary Surface number.
     * @param rmode Mode of the contribution.
     * @param tStep Time step when answer is computed.
     */
    virtual void computeTangentFromEdgeLoad(FloatMatrix &answer, EdgeLoad *load, int boundary, MatResponseMode rmode, TimeStep *tStep);
    /**
     * Computes the contribution of the given load at the given boundary edge. 
     * In general, the answer should include only relevant DOFs at the edge.
     * The related is giveBoundaryLocationArray method, which should return 
     * corresponding code numbers..
     * @note Elements which do not have an contribution should resize the vector to be empty.
     * @param answer Requested contribution of load (in Global c.s.).
     * @param load Load to compute contribution from.
     * @param edge Edge number.
     * @param type Type of the contribution.
     * @param mode Determines mode of answer.
     * @param tStep Time step when answer is computed.
     * @param global if true (default) then contribution is in global c.s., when false then contribution is in element local c.s.
     */
    virtual void computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode, TimeStep *tStep, bool global=true);

    /**
     * Returns receiver list of bodyloads
     *
     */
    const IntArray& giveBodyLoadList() const {return this->bodyLoadArray;}
    /**
     * Returns receiver list of boundary loads
     */
    const IntArray& giveBoundaryLoadList() const {return this->boundaryLoadArray;}

    //@}

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
     * @param dofIDMask Dof IDs for unknowns.
     * @param u      Identifies mode of unknown (eg. total value or velocity of unknown).
     * @param tStep  Time step, when vector of unknowns is requested.
     * @param answer Local vector of unknowns.
     */
    void computeVectorOf(ValueModeType u, TimeStep *tStep, FloatArray &answer);
    void computeVectorOf(const IntArray &dofIDMask, ValueModeType u, TimeStep *tStep, FloatArray &answer, bool padding = false);
    /**
     * Boundary version of computeVectorOf.
     * @param bNodes Boundary nodes.
     * @param dofIDMask Dof IDs for unknowns.
     * @param u Identifies mode of unknown (eg. total value or velocity of unknown).
     * @param tStep Time step, when vector of unknowns is requested.
     * @param answer Local vector of unknowns.
     */
    void computeBoundaryVectorOf(const IntArray &bNodes, const IntArray &dofIDMask, ValueModeType u, TimeStep *tStep, FloatArray &answer, bool padding = false);
    /**
     * Returns local vector of unknowns. Local vector of unknowns is extracted from
     * given field and from boundary conditions (if dof has active boundary
     * and possibly initial condition). Because unknowns are obtained from given field
     * model, this must support queries for given unknown.
     *
     * @param field  Source field (eg. displacement or temperature vector).
     * @param dofIDMask Dof IDs for unknowns.
     * @param u Value mode of unknown (incremental, total, ...).
     * @param tStep Time step, when vector of unknowns is requested.
     * @param answer Local vector of unknowns.
     */
    void computeVectorOf(PrimaryField &field, const IntArray &dofIDMask, ValueModeType u, TimeStep *tStep, FloatArray &answer, bool padding = false);
    /**
     * Returns local vector of prescribed unknowns. Local vector of prescribed unknowns is
     * extracted from nodal (and side - if they hold unknowns) boundary conditions.
     * @param u      Identifies mode of unknown (eg. total value or velocity of unknown).
     * @param tStep  Time step, when vector of unknowns is requested.
     * @param answer Local vector of unknowns.
     */
    void computeVectorOfPrescribed(ValueModeType u, TimeStep *tStep, FloatArray &answer);    
    /**
     * Returns local vector of prescribed unknowns. Local vector of prescribed unknowns is
     * extracted from nodal (and side - if they hold unknowns) boundary conditions.
     *
     * @param dofIDMask Dof IDs for unknowns.
     * @param ut Identifies mode of unknown (eg. total values or velocity of unknown).
     * @param tStep Time step, when vector of prescribed unknowns is requested.
     * @param answer Local vector of prescribed unknowns. If unknown is not prescribed,
     * zero value is placed on position of free dof.
     */
    void computeVectorOfPrescribed(const IntArray &dofIDMask, ValueModeType type, TimeStep *tStep, FloatArray &answer);

    /**
     * Computes or simply returns total number of element's local DOFs.
     * Must be defined by particular element.
     * @return Number of local DOFs of element.
     */
    virtual int computeNumberOfDofs() { return 0; }
    /**
     * Computes the total number of element's global dofs.
     * The transitions from global c.s. to nodal c.s. should NOT be included.
     * @return Total number of global DOFs of element.
     */
    virtual int computeNumberOfGlobalDofs();
    /**
     * Computes the total number of element's primary master DOFs.
     * @return Total number of DOFs belonging to eid.
     */
    int computeNumberOfPrimaryMasterDofs();
    /**
     * Returns transformation matrix from global c.s. to local element
     * c.s., i.e. @f$ r_l =T r_g @f$.
     * If no transformation is necessary then answer is empty matrix and zero value is returned.
     * @param answer Computed rotation matrix.
     * @return Nonzero if transformation is necessary, zero otherwise.
     */
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
    /**
     * Transformation matrices updates rotation matrix between element-local and primary DOFs,
     * taking into account nodal c.s. and master DOF weights.
     * @param answer Contains the rotation matrix on exit.
     * @return True if there is a rotation required, false otherwise.
     */
    virtual bool giveRotationMatrix(FloatMatrix &answer);
    /**
     * Returns transformation matrix for DOFs from global coordinate system
     * to local coordinate system in nodes.
     * Also includes transformations to slave DOFs.
     * If no transformation is necessary sets answer to empty matrix and returns false.
     * Local stiffness matrix of element should be rotated with answer before assembly.
     * @note Function does most likely NOT need to be overridden.
     * @param answer Computed rotation matrix.
     * @param nodes Nodes to include in element local ordering.
     * @param includeInternal Determines whether or not to include internal dof managers.
     * @return True if transformation is necessary, false otherwise.
     */
    virtual bool computeDofTransformationMatrix(FloatMatrix &answer, const IntArray &nodes, bool includeInternal);
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
     * @param answer Mask for node.
     */
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const { answer.clear(); }
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
     * @param answer mask for node.
     */
    virtual void giveInternalDofManDofIDMask(int inode, IntArray &answer) const
    { answer.clear(); }
    /**
     * Returns element dof mask for node. This mask defines the dof ordering of the element interpolation.
     * Default implementation for most elements, with noteable exceptions such as XFEM and some types of shell elements.
     *
     * @param ut Equation DOFs belong to.
     * @param answer DOF mask for receiver.
     */
    virtual void giveElementDofIDMask(IntArray &answer) const { this->giveDofManDofIDMask(1, answer); }
    /**
     * Computes the unknown vector interpolated at the specified local coordinates.
     * Used for exporting data and mapping fields.
     * @see giveElementDofIDMask The unknown vector should match the element field as specified by the element dof IDs.
     * @param mode Mode (total, increment, etc) of the output
     * @param tStep Time step to evaluate at
     * @param lcoords Local coordinates to evaluate at
     * @param answer Results
     */
    virtual void computeField(ValueModeType mode, TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
    { OOFEM_ERROR("Missing support for computing unknown vector at local element coordinates\n"); }
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
    virtual double computeVolumeAreaOrLength();
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
    /**
     * Returns list of receiver boundary nodes for given edge
     * @param bNodes list of boundary edge nodes 
     * @param boundary edge id
     */
    virtual void giveBoundaryEdgeNodes (IntArray& bNodes, int boundary);
    /**
     * Returns list of receiver boundary nodes for given surface
     * @param bNodes list of boundary surface nodes 
     * @param boundary surface id
     */
    virtual void giveBoundarySurfaceNodes (IntArray& bNodes, int boundary);
    /**
     * Returns boundary edge integration rule
     * @param order approximation order to integrate 
     * @param boundary boundary edge id
     * @note some elements may increase the order (like axusymmetric elements)
     */
    virtual IntegrationRule* giveBoundaryEdgeIntegrationRule (int order, int boundary);
    /**
     * Returns boundary surface integration rule
     * @param order approximation order to integrate 
     * @param boundary boundary surface id
     * @note some elements may increase the order (like axusymmetric elements)
      */
    virtual IntegrationRule* giveBoundarySurfaceIntegrationRule (int order, int boundary);



    
    // data management
    /**
     * Translates local to global indices for dof managers.
     * @param i Local index of dof manager.
     * @return Global number of i-th dofmanager of element
     */
    int giveDofManagerNumber(int i) const { return dofManArray.at(i); }
    /// @return Receivers list of dof managers.
    const IntArray &giveDofManArray() const { return dofManArray; }
    /**
     * @param i Local index of the dof manager in element.
     * @return The i-th dofmanager of element.
     */
    void addDofManager(DofManager *dMan);
    /**
     * @param dMan Pointer t a dof manager to add to the element.
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
    inline Node *giveNode(int i) const
    {
#ifdef DEBUG
        if ( ( i <= 0 ) || ( i > dofManArray.giveSize() ) ) {
            OOFEM_ERROR("Node is not defined");
        }
#endif
        return domain->giveNode( dofManArray.at(i) );
    }

    /**
     * Returns reference to the i-th side  of element.
     * Default implementation returns i-th dofmanager of element converted to
     * ElementSide class (check is made).
     * @param i Side number.
     * @return Requested element side.
     */
    virtual ElementSide *giveSide(int i) const;
    /// @return Interpolation of the element geometry, or NULL if none exist.
    virtual FEInterpolation *giveInterpolation() const { return NULL; }
    /**
     * Returns the interpolation for the specific dof id.
     * Special elements which uses a mixed interpolation should reimplement this method.
     * @param id ID of the dof for the for the requested interpolation.
     * @return Appropriate interpolation, or NULL if none exists.
     */
    virtual FEInterpolation *giveInterpolation(DofIDItem id) const { return giveInterpolation(); }
    /// @return Reference to the associated material of element.
    virtual Material *giveMaterial();
    /// @return Material number.
    int giveMaterialNumber() const {return material;}
    /// @return Reference to the associated crossSection of element.
    CrossSection *giveCrossSection();
    /**
     * Sets the cross section model of receiver.
     * @param csIndx Index of new cross section.
     */
    virtual void setCrossSection(int csIndx) { this->crossSection = csIndx; }

    /// @return Number of dofmanagers of receiver.
    virtual int giveNumberOfDofManagers() const { return numberOfDofMans; }
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
     * Sets receiver bodyLoadArray.
     * @param bodyLoads Array with body loads indices.
     */
    void setBodyLoads(const IntArray &bodyLoads);

    /**
     * Sets integration rules.
     * @param irlist List of integration rules.
     */
    void setIntegrationRules(std :: vector< std :: unique_ptr< IntegrationRule > > irlist);
    /**
     * Returns integration domain for receiver, used to initialize
     * integration point over receiver volume.
     * Default behavior is taken from the default interpolation.
     * @return Integration domain of element.
     */
    virtual integrationDomain giveIntegrationDomain() const;
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
     * @todo This is currently not used. It is intended for IGA elements? This seems redundant.
     */
    virtual int giveIntegrationRuleLocalCodeNumbers(IntArray &answer, IntegrationRule &ie)
    { return 0; }

    // Returns number of sides (which have unknown dofs) of receiver
    //int giveNumberOfSides () {return numberOfSides;}

    /// @return Corresponding element region. Currently corresponds to cross section model number.
    int giveRegionNumber();

    /// Performs post initialization steps.
    virtual void postInitialize();

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

    /**
     * @return True, if receiver is activated for given solution step, otherwise false.
     */
    virtual bool isActivated(TimeStep *tStep);

    /**
     * @return True, if the current time is higher than the casting time of the material, otherwise false.
     * Used from e.g. vtkxml export module to display only active elements
     * @note: The element can be activated (isActivated method) before its 
     * material is actually casted. This case has to be supported by 
     * the material and can be used to simulate the casting on deformed 
     * configuration, for example. 
     * In this case, the material has to define a small stifness 
     * for solution steps before is actually casted.
     */
    virtual bool isCast(TimeStep *tStep);

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
    virtual Element_Geometry_Type giveGeometryType() const;
    /**
     * Returns the element spatial dimension (1, 2, or 3).
     * This is completely based on the geometrical shape, so a plane in space counts as 2 dimensions.
     * @return Number of spatial dimensions of element.
     */
    virtual int giveSpatialDimension();
    /**
     * @return Number of boundaries of element.
     */
    virtual int giveNumberOfBoundarySides();
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
    virtual IntegrationRule *giveDefaultIntegrationRulePtr() {
        if ( this->giveNumberOfIntegrationRules() == 0 ) {
            return NULL;
        } else {
            return this->integrationRulesArray [ giveDefaultIntegrationRule() ].get();
        }
    }
    /// @return Number of integration rules for element.
    int giveNumberOfIntegrationRules() { return (int)this->integrationRulesArray.size(); }
    /**
     * @param i Index of integration rule.
     * @return Requested integration rule.
     */
    virtual IntegrationRule *giveIntegrationRule(int i) { return integrationRulesArray [ i ].get(); }
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
     * Returns the integration point corresponding value in full form.
     * @param answer Contain corresponding integration point value, zero sized if not available.
     * @param gp Integration point to check.
     * @param type Determines the type of internal variable.
     * @param tStep Time step.
     * @return Nonzero if o.k, zero otherwise.
     */
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    int giveGlobalIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    // characteristic length in gp (for some material models)
    /**
     * Default implementation returns length of element projection into specified direction.
     * @return Element length in given direction.
     */
    virtual double giveLengthInDir(const FloatArray &normalToCrackPlane) ;
    /**
     * Returns the size of element in the given direction, in some cases adjusted (e.g. if the direction is perpendicular to a planar element). 
     * Required by material models relying on the crack-band approach to achieve objectivity with respect to the mesh size. 
      * @param normalToCrackPlane Normal to the expected crack band.
     * @return Element size corresponding to the given direction (expected width of the crack band).
     */
    virtual double giveCharacteristicLength(const FloatArray &normalToCrackPlane) { OOFEM_ERROR("Function not overloaded, which probably means that the crack band approach should not be used for this element"); return 0.; }
    /**
     * Returns the size of element in the given direction if the direction is in the XY plane,
     * otherwise gives the mean size defined as the square root of the element area.
     * Required by material models relying on the crack-band approach to achieve objectivity with respect to the mesh size. 
      * @param normalToCrackPlane Normal to the expected crack band.
     * @return Element size corresponding to the given direction (expected width of the crack band).
     */
    double giveCharacteristicLengthForPlaneElements(const FloatArray &normalToCrackPlane) ;
    /**
     * Returns the size of an axisymmetric element in the given direction if the direction is in the XY plane,
     * otherwise gives the mean distance vrom the symmetry axis multiplied by pi.
     * Required by material models relying on the crack-band approach to achieve objectivity with respect to the mesh size. 
      * @param normalToCrackPlane Normal to the expected crack band.
     * @return Element size corresponding to the given direction (expected width of the crack band).
     */
    double giveCharacteristicLengthForAxisymmElements(const FloatArray &normalToCrackPlane) ;
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
    virtual double giveCharacteristicSize(GaussPoint *gp, FloatArray &normalToCrackPlane, ElementCharSizeMethod method) { return giveCharacteristicLength(normalToCrackPlane); }
    /**
     * Returns the size (length, area or volume depending on element type) of the parent
     * element. E.g. 4.0 for a quadrilateral.
     */
    virtual double giveParentElSize() const { return 0.0; }
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
     * @param answer Local coordinates.
     * @param gcoords Global coordinates.
     * @return Nonzero if point is inside element; zero otherwise.
     */
    virtual bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);
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
     * Maps the internal state variables stored in all IPs from the old domain to the new domain.
     * @param iOldDom Old domain.
     * @param iTStep Time step.
     * @return Nonzero if o.k, otherwise zero.
     */
    virtual int mapStateVariables(Domain &iOldDom, const TimeStep &iTStep);
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
    template< class T > void ipEvaluator( T *src, void ( T :: *f )( GaussPoint *gp ) );
    /// Integration point evaluator, loops over receiver IP's and calls given function (passed as f parameter) on them. The IP is parameter to function f as well as additional array.
    template< class T, class S > void ipEvaluator(T *src, void ( T :: *f )( GaussPoint *, S & ), S &_val);

    //@}


#ifdef __OOFEG
    //
    // Graphics output
    //
    virtual void drawYourself(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawAnnotation(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) { }
    virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) { }
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) { }
    virtual void drawSpecial(oofegGraphicContext &gc, TimeStep *tStep) { }
    // added in order to hide IP element details from oofeg
    // to determine the max and min local values, when recovery does not takes place
    virtual void giveLocalIntVarMaxMin(oofegGraphicContext &gc, TimeStep *tStep, double &emin, double &emax) { emin = emax = 0.0; }

    /**
     * Returns internal state variable (like stress,strain) at node of element in Reduced form,
     * the way how is obtained is dependent on InternalValueType.
     * The value may be local, or smoothed using some recovery technique.
     * Returns zero if element is unable to respond to request.
     * @param answer Contains result, zero sized if not supported.
     * @param type Determines the internal variable requested (physical meaning).
     * @param mode Determines the mode of variable (recovered, local, ...).
     * @param node Node number, for which variable is required.
     * @param tStep Time step.
     * @return Nonzero if o.k, zero otherwise.
     */
    virtual int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                        int node, TimeStep *tStep);
    /**
     * Returns internal state variable (like stress,strain) at side of element in Reduced form
     * If side is possessing DOFs, otherwise recover techniques will not work
     * due to absence of side-shape functions.
     * @param answer Contains result, zero sized if not supported.
     * @param type Determines the internal variable requested (physical meaning).
     * @param mode Determines the mode of variable (recovered, local, ...).
     * @param side Side number, for which variable is required.
     * @param tStep Time step.
     * @return Nonzero if o.k, zero otherwise.
     */
    virtual int giveInternalStateAtSide(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                        int side, TimeStep *tStep)
    {
        answer.clear();
        return 0;
    }

    /// Shows sparse structure
    virtual void showSparseMtrxStructure(CharType mtrx, oofegGraphicContext &gc, TimeStep *tStep) { }
    /// Shows extended sparse structure (for example, due to nonlocal interactions for tangent stiffness)
    virtual void showExtendedSparseMtrxStructure(CharType mtrx, oofegGraphicContext &gc, TimeStep *tStep) { }

#endif

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

    /**
     * Return elementParallelMode of receiver.
     */
    elementParallelMode giveParallelMode() const { return parallel_mode; }
    /// Sets parallel mode of element
    void setParallelMode(elementParallelMode _mode) { parallel_mode = _mode; }
    /**
     * Returns the parallel mode for particular knot span of the receiver.
     * The knot span identifies the sub-region of the finite element.
     */
    virtual elementParallelMode giveKnotSpanParallelMode(int) const { return parallel_mode; }
    /**
     * Pack all necessary data of element (according to its parallel_mode) integration points
     * into given communication buffer. The corresponding cross section service is invoked, which in
     * turn should invoke material model service for particular integration point. The
     * nature of packed data is material model dependent.
     * Typically, for material of "local" response (response depends only on integration point local state)
     * no data are exchanged. For "nonlocal" constitutive models the send/receive of local values which
     * undergo averaging is performed between local and corresponding remote elements.
     * @param buff communication buffer
     * @param tStep solution step.
     */
    int packUnknowns(DataStream &buff, TimeStep *tStep);
    /**
     * Unpack and updates all necessary data of element (according to its parallel_mode) integration points
     * into given communication buffer.
     * @see packUnknowns service.
     * @param buff communication buffer
     * @param tStep solution step.
     */
    int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep);
    /**
     * Estimates the necessary pack size to hold all packed data of receiver.
     * The corresponding cross section service is invoked, which in
     * turn should invoke material model service for particular integration point. The
     * nature of packed data is material model dependent.
     */
    int estimatePackSize(DataStream &buff);
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

public:
    /// Returns array containing load numbers of loads acting on element
    IntArray *giveBodyLoadArray();
    /// Returns array containing load numbers of boundary loads acting on element.
    IntArray *giveBoundaryLoadArray();

    // Overloaded methods:
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual const char *giveClassName() const { return "Element"; }

protected:
    /**
     * Initializes the array of integration rules member variable.
     * Element can have multiple integration rules for different tasks.
     * For example structural element family class uses this feature to implement
     * transparent support for reduced and selective integration of some strain components.
     * Must be defined by terminator classes.
     * @see IntegrationRule
     */
    virtual void computeGaussPoints() { }
};

template< class T > void
Element :: ipEvaluator( T *src, void ( T :: *f )( GaussPoint *gp ) )
{
    for ( auto &ir: integrationRulesArray ) {
        for ( GaussPoint *gp: *ir ) {
            ( src->*f )(gp);
        }
    }
}

template< class T, class S > void
Element :: ipEvaluator(T *src, void ( T :: *f )( GaussPoint *, S & ), S &_val)
{
    for ( auto &ir: integrationRulesArray ) {
        for ( GaussPoint *gp: *ir ) {
            ( src->*f )(gp, _val);
        }
    }
}

} // end namespace oofem
#endif //element_h
