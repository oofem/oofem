/* $Header: /home/cvs/bp/oofem/oofemlib/src/dofmanager.h,v 1.19.4.2 2004/05/14 13:45:27 bp Exp $ */
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


//
// class DofManager
//


#ifndef dofmanager_h
#define dofmanager_h

#include "alist.h"
#include "femcmpnn.h"
#include "domain.h"
#include "dof.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif
#include "classtype.h"
#include "equationid.h"
#include "valuemodetype.h"
#include "doftype.h"
#include "dofmantransftype.h"
#include "dofiditem.h"
#include "contextioresulttype.h"

#ifdef __PARALLEL_MODE
class CommunicationBuffer;
/**
 * In parallel mode, this type indicates the mode of DofManager.
 * <UL>
 * <LI>
 * DofManager_local mode - DofManager is local, there are no contribution from other domains to this DofManager.</LI>
 * <LI>
 * DofManager_shared mode - DofManager is shared by neighbouring partitions, it is necessary to
 * summ contributions from all contributing domains. Typical for node cut algorithm.</LI>
 * <LI>
 * DofManager_remote - DofManager in active domain is only mirror of some remote
 * DofManager. It is necessary to copy remote values into local ones. Typical for element cut.</LI>
 * <LI>
 * DofManager_null - DofManager in active domain is shared only by remote elements (these are only
 * introduced for nonlocal constitutive model to allow effective local averaging, so only local material
 * value to be averaged are transferred for these remote elements). Null nodes are therefore used only
 * for computing real integration point coordinates of remote elements and there is no reason to maintain
 * their unknowns (they have no equatiuon number assigned). </LI>
 *
 * </UL>
 */
enum dofManagerParallelMode { DofManager_local, DofManager_shared, DofManager_remote, DofManager_null };

#endif

class NodalLoad;
class TimeStep;
class FloatArray;
class IntArray;

/**
 * Base class for dof managers. Dof manager is an abstraction for object possessing degrees of
 * freedom. Dof managers (respectively derived clases like nodes or sides) are usually atributes of
 * elements and are maintained by domain. Degrees of freedom belonging to dof manager are
 * stored in 'dofArray'. Dof manager also
 * maintain a list of applied loads it is subjected to.
 * Number and physical meaning of dofs can be specified by user in input file
 * (see input file description). If it is not specified, default values are obtained from domain,
 * based on domain type of problem.
 */
class DofManager : public FEMComponent
{
    /*
     * This class implements absstract clas - a  dof manager in a finite element mesh.
     * It is a base class for FEM nodes and element sides, having dofs.
     * A dof manager  and its children (nodes, element sides) are
     * attributes of a domain.
     * They are  usually also attributes of a few elements.
     * DESCRIPTION
     * The DofManager possesses 'numberOfDofs' degrees of freedom, stored in 'dofArray'.
     * In 'loadArray' it stores the number of every  load it is subjected to
     * (typically, concentrated forces and moments).
     * In 'locationArray' the manager stores the equation number of each of its dofs.
     * This location array is used by the manager for assembling its load vector to
     * the right-hand side of the linear system ; it is also used by elements for
     * calculating their own location arrays.
     * TASKS
     * - managing its degrees of freedom (method 'giveDof') ;
     * - calculating its nodal  vector;
     * - printing and updating at end of step ;
     * - managing its swapping to and from disk.
     */

protected:
    /// total number of DOFs
    int numberOfDofs;
    /// array of DOFs
    Dof **dofArray;
    /// list of applied loads.
    IntArray loadArray;
    // IntArray*    locationArray ;
    /** Indicates if dofManager is boundary (true boundary or on boundary between regions) or interior.
     * This information is required by some recovery technigues. */
    bool isBoundaryFlag;
    /// flag indicating whether receiver has slave dofs
    bool hasSlaveDofs;
#if defined(__PARALLEL_MODE) || defined(__ENABLE_COMPONENT_LABELS)
    /**
     * In parallel mode, globalNumber contains globally unique DoFManager number.
     * The component number, inherited from FEMComponent class contains
     * local domain number.
     */
    int globalNumber;
#endif
#ifdef __PARALLEL_MODE
    dofManagerParallelMode parallel_mode;
    /**
     * List of partition sharing the shared dof manager or
     * remote partion containing remote dofmanager counterpart.
     */
    IntArray partitions;
#endif

public:
    /** Constructor. Creates DofManager with given number belonging to domain aDomain.
     * @param n DofManager's number in domain
     * @param aDomain reference to DofManager's domain
     */
    DofManager(int n, Domain *aDomain);                      // constructor
    /// Destructor.
    ~DofManager();                                 // destructor

    /**@name Dof management methods */
    //@{
    // dof management
    /** Returns reference (pointer) to i-th dof of receiver.
     * Index of Dof with required physical meaning can be obtained by invoking
     * method findDofWithDofId.
     * @see DofManager::findDofWithDofId
     */
    Dof *giveDof(int i) const;
    /// Returns DOF with given dofID; issues error if not present
    Dof *giveDofWithID(int dofID) const; // termitovo
    /** Returns total number of dofs managed by receiver */
    int          giveNumberOfDofs() const;

    /** Returns the number of primary dofs on which receiver dofs (given in dofArray) depend on.
     *  If receiver has only prinary dofs, the answer is the size of dofArray.
     */
    int giveNumberOfPrimaryMasterDofs(IntArray &dofArray) const;
    //virtual int          giveNumberOfDofs () const;
    /** Returns location array (array containing for each requested dof related equation number)
     * of receiver. Equation number of dof with active boundary condition is zero.
     * @param dofIDArry array containing dof mask. This mask containing DofIDItem values
     * (they describe physical meaning of dofs, see cltypes.h) is used to extract only
     * required values. If dof with requested physical meaning dos not exist in receiver,
     * an error is generated and execution exits.
     * @param locationArray - return parameter containing required equation numbers.
     * @see Element::giveDofManDofIDMask function.
     */
    virtual void         giveLocationArray(const IntArray &dofIDArry, IntArray &locationArray) const;
    /** Returns prescribed equations location  array
     * (array containing for each requested dof related equation number)
     * of receiver. Equation number of dof with inactive or no boundary condition is zero.
     * @param dofIDArry array containing dof mask. This mask containing DofIDItem values
     * (they describe physical meaning of dofs, see cltypes.h) is used to extract only
     * required values. If dof with requested physical meaning dos not exist in receiver,
     * an error is generated and execution exits.
     * @param locationArray - return parameter containing required equation numbers.
     * @see Element::giveDofManDofIDMask function.
     */
    virtual void         givePrescribedLocationArray(const IntArray &dofIDArry, IntArray &locationArray) const;
    /** Returns full location array of receiver containing equation numbers of all dofs
     * of receiver. Their order is specific to every DofManager. Mainly used at EngngModel level
     * to assemble DofManager contribution (typically load vector).
     * @param locationArray complete location array of receiver.
     */
    virtual void         giveCompleteLocationArray(IntArray &locationArray) const;
    /** Returns full location array of prescribed equations of receiver containing equation
     * numbers of all dofs of receiver. Their order is specific to every DofManager.
     * Mainly used at EngngModel level to assemble DofManager contribution (typically load vector).
     * @param locationArray complete location array of receiver.
     */
    virtual void         giveCompletePrescribedLocationArray(IntArray &locationArray) const;
    /** Returns DOFs numbers of receiver with required physical meaning.
     * @param dofIDArry array containing DofIDItem-type values (this is enumeration
     * identifying physical meaning of particular DOF, see cltypes.h).
     * @param answer array with DOF numbers. They are ordered according to dofIDArry.
     * @see cltypes.h
     */
    void         giveDofArray(const IntArray &dofIDArry, IntArray &answer) const;
    /**
     * Finds index of DOF with required physical meaning of receiver. This index can be different
     * for different DOFManagers (user can alter dof order and type in input file).
     * @param dofID physical meaning of DOF.
     * @return index of requested DOF. If such DOF with dofID does not exists, returns zero.
     */
    int          findDofWithDofId(DofID dofID) const;
    /**
     * Assembles the vector of unknowns in nodal c.s for given dofs of receiver.
     * This vector may have size different from number of dofs requested,
     * because some dofs may depend on other dofs. Default implementation uses
     * Dof::giveUnknown service.
     * @param answer result (in nodal cs.)
     * @param dofMask  dofIDArry array containing dof mask. This mask containing DofIDItem values
     * (they describe physical meaning of dofs, see cltypes.h) is used to extract only
     * required values. If dof with requested physical meaning dos not exist in receiver,
     * an error is generated and execution exits.
     * @param type physical meaning of  unknown.
     * @param mode mode of unknown (e.g, total value, velocity or acceleration of unknown).
     * @stepN time step when unknown requested. See documentation of particular EngngModel
     * class for valid StepN values (most implementaion can return only values for current
     * and possibly for previous time step).
     * @see Dof::giveUnknown service.
     */
    virtual void giveUnknownVector(FloatArray &answer, const IntArray &dofMask,
                                   EquationID type, ValueModeType mode, TimeStep *stepN);
    /**
     * Assembles the vector of unknowns of given filed in nodal c.s for given dofs of receiver.
     * This vector may have size different from number of dofs requested,
     * because some dofs may depend on other dofs. Default implementation uses
     * Dof::giveUnknown service.
     * @param answer result (in nodal cs.)
     * @param dofMask  dofIDArry array containing dof mask. This mask containing DofIDItem values
     * (they describe physical meaning of dofs, see cltypes.h) is used to extract only
     * required values. If dof with requested physical meaning dos not exist in receiver,
     * an error is generated and execution exits.
     * @param field primary filed
     * @stepN time step when unknown requested. See documentation of particular EngngModel
     * class for valid StepN values (most implementaion can return only values for current
     * and possibly for previous time step).
     * @see Dof::giveUnknown service.
     */
    virtual void giveUnknownVector(FloatArray &answer, const IntArray &dofMask,
                                   PrimaryField &field, ValueModeType mode, TimeStep *stepN);
    /**
     * Assembles the vector of prescribed unknowns in nodal c.s for given dofs of receiver.
     * This vector may have size different from number of dofs requested,
     * because some dofs may depend on other dofs. Default implementation uses
     * Dof::giveBcValue and Dof::hasBc service.
     * @param answer result (in nodal cs.)
     * @param dofMask  dofIDArry array containing dof mask. This mask containing DofIDItem values
     * (they describe physical meaning of dofs, see cltypes.h) is used to extract only
     * required values. If dof with requested physical meaning dos not exist in receiver,
     * an error is generated and execution exits.
     * @param mode mode of unknown (e.g, total value, velocity or acceleration of unknown).
     * @stepN time step when unknown requested. See documentation of particular EngngModel
     * class for valid StepN values (most implementaion can return only values for current
     * and possibly for previous time step).
     * @see Dof::giveBcValue, Dof::hasBc service.
     */
    virtual void givePrescribedUnknownVector(FloatArray &answer, const IntArray &dofMask,
                                             ValueModeType mode, TimeStep *stepN);
    //virtual void givePrescribedUnknownVector (FloatArray& answer, IntArray& dofMask,
    //                       UnknownType type, ValueModeType mode, TimeStep* stepN);



    //@}
    /**@name Transformation functions
     * The governing equations can be assembled not only in global coordinate system, but
     * also in user-defined local coordinate system of each dof manager. Methods in this section
     * introduce necessary transformation methods, allowing receiver dofs to be expressed in
     * their own local c.s. or to be dependent on other dofs on other dofManager (to implement
     * slave or rigid arm nodes etc.). The method for computing global c.s to receiver c.s
     * transformation matrix is provided.
     */
    //@{
    /** Computes receiver transformation matrix from global cs. to dofManager specific
     * coordinate system (in which governing equations are assembled, for example the
     * local coordinate system in node).
     * @param answer computed transformation matrix. It has generally dofIDArry.size rows and
     * if loc is obtained using giveLocationArray(dofIDArry, loc) call, loc.giveSize() columns.
     * This is because this transformation should generally include not only transformation to
     * dof manager local coordinate system, but receiver dofs can be expressed using
     * dofs of another dofManager (In this case, squre answer is produced anly if all
     * dof transformation is required).
     * @param dofIDArry array containing DofIDItem-type values (this is enumeration
     * identifying physical meaning of particular DOF, see cltypes.h) for which transfromation mtrx is
     * assembled. if dofIDArry is NULL, then all receiver dofs are assumed.
     */
    virtual void computeDofTransformation(FloatMatrix &answer, const IntArray *dofIDArry, DofManTransfType mode);
    virtual void computeLoadTransformation(FloatMatrix &answer, const IntArray *dofIDArry, DofManTransfType mode);
    /**
     * Indicates, whether dofManager requires the transformation from global c.s. to
     * dof manager specific coordinate system.
     * @return nonzero if transformation is necessary, even for single dof.
     */
    virtual int requiresTransformation() { return 0; }
    //@}
    /**@name Load related functions */
    //@{
    // nodal load vector
    /** Computes the load vector of receiver in given time.
     * @param answer load vector.
     * @param stepN time step when answer is computed.
     * @param mode determines response mode.
     */
    virtual void         computeLoadVectorAt(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
    /// Returns the array containing applied loadings of the receiver
    IntArray *giveLoadArray();
    /// Sets the array of applied loadings of the receiver
    void setLoadArray(IntArray &);
    //@}

    /**@name Position querry functions */
    //@{
    virtual bool        hasCoordinates() { return false; }
    /// Returns i-th coordinate of node.
    virtual double       giveCoordinate(int i) { return 0.0; }
    /// Returns pointer to node coordinate array.
    virtual FloatArray *giveCoordinates() { return NULL; }
    //@}



    // time step termination
    /** Prints output of receiver to stream, for given time step */
    void         printOutputAt(FILE *, TimeStep *);
    /// updates receiver after equlibrium in time step has been reached.
    void         updateYourself(TimeStep *);

    // miscellaneous
    bool         isBoundary() { return isBoundaryFlag; }
    void         setBoundaryFlag(bool _b) { this->isBoundaryFlag = _b; }
    /// Returns true if receiver contains slave dofs
    virtual int  hasAnySlaveDofs();
    /**
     * Returns true if the receiver is linked (its slave DOFs depend on master values) to some other dof managers.
     * In this case, the masters array should contain the list of masters. 
     * In both serial and parallel modes, local numbers are be provided. 
     * If the receiver contains only primary DOFs, false is returned.
     */
    virtual bool giveMasterDofMans (IntArray& masters);

    /**
     * Returns a newly allocated DofManager, with type depending on parameter.
     * Creates new object for following classes Node, ElementSide, RigidArmNode otherwise
     * calls global function CreateUsrDefDofManagerOfType for creating appropriate
     * instance. This function must be implemented by user.
     * @param aClass string with DofManager name
     * @return newly allocated DofManager with required type.
     */
    DofManager *ofType(char *);
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "DofManager"; }
    /** Returns classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    classType    giveClassID() const { return DofManagerClass; }
    ///Initializes receiver acording to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();
    /**
     * Stores receiver state to output stream.
     * @exception throws an ContextIOERR exception if error encountered
     */
    virtual contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores the receiver state previously written in stream.
     * @exception throws an ContextIOERR exception if error encountered
     */
    virtual contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    /// Returns true if dof of given type is allowed to be associated to receiver
    virtual bool isDofTypeCompatible(dofType type) const { return false; }
    /**
     * Checks internal data consistency in node.
     * Current implementation checks (when receiver has slave dofs) if receiver has the same
     * coordinate system as master dofManager of slave dof.
     * @return nonzero if receiver check is o.k.
     */
    virtual int    checkConsistency();
    /**
     * Local renumbering support. For some tasks (parallel load balancing, for example) it is necessary to
     * renumber the entities. The various fem components (such as nodes or elements) typically contain
     * links to other entities in terms of their local numbers, etc. This service allows to update
     * these relations to reflext updated numbering. The renumbering funciton is passed, which is supposed
     * to return an updated number of specified entyty type based on old number.
     */
    virtual void updateLocalNumbering( EntityRenumberingFunctor &f ) ;

    /**@name Advanced functions */
    //@{
    /** Sets number of dofs of the receiver; Dealocates existing DOFs;
     *  Resizes the dofArray accordingly */
    void setNumberOfDofs(int _ndofs);
    /** Sets i-th DOF of receiver to given DOF */
    void setDof(int i, Dof *dof);
    //@}
    /** Adds a Dof to i-th position in dofArray */
    void addDof(int i, Dof *dof);   // rch

#ifdef __OOFEG
    virtual void   drawYourself(oofegGraphicContext &context) { }
#endif

#if defined(__PARALLEL_MODE) || defined(__ENABLE_COMPONENT_LABELS)
    /**
     * Returns receiver globally unique number.
     */
    int giveGlobalNumber() const { return globalNumber; }
    int giveLabel() const {return globalNumber;}
    /** sets receiver global number */
    void setGlobalNumber(int _number) { globalNumber = _number; }
#endif
#ifdef __PARALLEL_MODE
    /**
     * Return dofManagerParallelMode of receiver. Defined for __Parallel_Mode only.
     */
    dofManagerParallelMode giveParallelMode() const { return parallel_mode; }
    /** Sets parallel mode of receiver */
    void setParallelMode(dofManagerParallelMode _mode) { parallel_mode = _mode; }
    /**
     * Packs specific  DOF Manager's dofs unknowns into communication buffer.
     * @param buff communication buffer to pack data.
     * @param type physical meaning of  unknown.
     * @param mode mode of unknown (e.g, total value, velocity or acceleration of unknown).
     * @stepN time step when unknown requested. See documentation of particular EngngModel
     * class for valid StepN values (most implementaion can return only values for current
     * and possibly for previous time step).
     * @return nonzero if succesfull
     */
    int packDOFsUnknowns(CommunicationBuffer &buff, EquationID type, ValueModeType mode, TimeStep *stepN);
    /**
     * Returns partition list of receiver.
     * @return partition array.
     */
    const IntArray *givePartitionList()  { return & partitions; }
    /** Sets receiver's partition list */
    void setPartitionList(const IntArray *_p) { partitions = * _p; }
    /// Removes given partition from receiver list
    void removePartitionFromList(int _part) {
        int _pos = partitions.findFirstIndexOf(_part);
        if ( _pos ) { partitions.erase(_pos); } }
    /// Merges receiver partition list with given lists
    void mergePartitionList(IntArray &_p);
    /**
     *  Returns number of partitions sharing given receiver (=number of shared partitions + local one)
     */
    const int givePartitionsConnectivitySize();
    /// Returns true if receiver is locally maintained
    bool isLocal();
    /// Returns true if receiver is shared
    bool isShared() { if ( parallel_mode == DofManager_shared ) { return true; } else { return false; } }

#endif
protected:
    virtual IRResultType resolveDofIDArray(InputRecord *ir, IntArray &dofIDArry);
    void computeSlaveLoadTransformation(FloatMatrix &answer, const IntArray *dofMask, DofManTransfType mode);
    void computeSlaveDofTransformation(FloatMatrix &answer, const IntArray *dofMask, DofManTransfType mode);
    IntArray *giveCompleteGlobalDofIDArray(void) const;
};


#endif // dofmanager_h






