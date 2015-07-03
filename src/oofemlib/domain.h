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

#ifndef domain_h
#define domain_h

#include "oofemcfg.h"
#include "domaintype.h"
#include "statecountertype.h"
#include "intarray.h"
#include "error.h"
#ifdef __PARALLEL_MODE
 #include "entityrenumberingscheme.h"
#endif

#include <memory>
#include <unordered_map>
#include <map>
#include <string>
#include <list>

///@name Input fields for domains
//@{
#define _IFT_Domain_type "domain" ///< This is trouble, will not work with dynamic input record
#define _IFT_Domain_ndofman "ndofman"
#define _IFT_Domain_nelem "nelem"
#define _IFT_Domain_nmat "nmat"
#define _IFT_Domain_ncrosssect "ncrosssect"
#define _IFT_Domain_nbc "nbc"
#define _IFT_Domain_nic "nic"
#define _IFT_Domain_nfunct "nltf"
#define _IFT_Domain_nset "nset"
#define _IFT_Domain_nbarrier "nbarrier"
#define _IFT_Domain_topology "topology"
#define _IFT_Domain_nxfemman "nxfemman" /// [in,optional] Specifies if there is an xfem-manager.
#define _IFT_Domain_ncontactman "ncontactman" /// [in,optional] Specifies if there is a contact manager.
#define _IFT_Domain_numberOfSpatialDimensions "nsd" ///< [in,optional] Specifies how many spatial dimensions the domain has.
#define _IFT_Domain_nfracman "nfracman" /// [in,optional] Specifies if there is a fracture manager.
#define _IFT_Domain_axisymmetric "axisymm" /// [optional] Specifies if the problem is axisymmetric.
//@}

namespace oofem {
class Element;
class Node;
class Material;
class GeneralBoundaryCondition;
class InitialCondition;
class Load;
class Function;
class CrossSection;
class ElementSide;
class DofManager;
class OutputManager;
class EngngModel;
class ConnectivityTable;
class ErrorEstimator;
class SpatialLocalizer;
class NodalRecoveryModel;
class NonlocalBarrier;
class DomainTransactionManager;
class XfemManager;
class TopologyDescription;
class DataReader;
class Set;
class FractureManager;
class oofegGraphicContext;
class ProcessCommunicator;
class ContactManager;
/**
 * Class and object Domain. Domain contains mesh description, or if program runs in parallel then it contains
 * description of domain associated to particular processor or thread of execution. Generally, it contain and
 * manages lists of Dof managers, elements, boundary conditions, cross sections and materials - these describe
 * the geometry of problem, its constitutive properties and applied boundary conditions. Services for accessing
 * these objects are provided. Domain is attribute of engineering model - which represent type of analysis
 * which should be performed.
 *
 * Domain also provides services for reading its description from
 * input stream and instantiating corresponding components accordingly. The basic Domain task are following
 * - Reading its description from input and creating corresponding objects.
 * - Provides services for accessing its particular components.
 * - Checking yourself.
 */
class OOFEM_EXPORT Domain
{
private:
    /// Element list.
    std :: vector< std :: unique_ptr< Element > > elementList;
    /// Dof manager list.
    std :: vector< std :: unique_ptr< DofManager > > dofManagerList;
    /// Material list.
    std :: vector< std :: unique_ptr< Material > > materialList;
    /// Cross section list.
    std :: vector< std :: unique_ptr< CrossSection > > crossSectionList;
    /// Boundary condition list.
    std :: vector< std :: unique_ptr< GeneralBoundaryCondition > > bcList;
    /// Initial condition list.
    std :: vector< std :: unique_ptr< InitialCondition > > icList;
    /// Load time function list.
    std :: vector< std :: unique_ptr< Function > > functionList;
    /// Set list.
    std :: vector< std :: unique_ptr< Set > > setList;
    /// Nonlocal barrier list.
    std :: vector< std :: unique_ptr< NonlocalBarrier > > nonlocalBarrierList;

    /// Default dofs for a node (depends on the domain type).
    IntArray defaultNodeDofIDArry;

    /**
     * Domain type. Determined by input data. It determines the problem type (like plane stress or plane strain mode).
     * According to this mode the default number of Dofs per node (or side) and their physical meaning are determined.
     * These default settings can be redefined by particular node or side. See related documentation for details.
     * @see Node
     * @see ElementSide
     */
    domainType dType;

    /**
     * Associated Engineering model. An abstraction for type of analysis which will be prformed.
     */
    EngngModel *engineeringModel;

    /**
     * Domain connectivity table. Table is build upon request.
     * Provides connectivity information of current domain.
     */
    std :: unique_ptr< ConnectivityTable > connectivityTable;
    /**
     * Spatial Localizer. It is build upon request.
     * Provides the spatial localization services.
     */
    std :: unique_ptr< SpatialLocalizer > spatialLocalizer;
    /// Output manager, allowing to filter the produced output.
    std :: unique_ptr< OutputManager > outputManager;
    /// Domain number.
    int number;
    /// Domain serial (version) number. Used for domain version identification during Adaptive computations.
    int serialNumber;
    /// Number of spatial dimensions
    int nsd;
    bool axisymm;
    /// nodal recovery object associated to receiver.
    std :: unique_ptr< NodalRecoveryModel > smoother; ///@todo I don't see why this has to be stored, and there is only one? /Mikael

    std :: string mDomainType;
    /**
     * For nonlocal models of integral type
     * it is necessary, mainly due to resulting efficiency, to compute variable(s)
     * which are nonlocally averaged in advance, before average process begins.
     * The loop over all  integration points is typically made to compute these variables.
     * To prevent doing this multiple times at the same solution state,
     * the modification time mark is kept.
     * This state counter could not be kept in static global variable,
     * because in case of multiple domains stateCounter should be kept independently for each domain.
     */
    StateCounterType nonlocalUpdateStateCounter;
    /// XFEM Manager
    std :: unique_ptr< XfemManager > xfemManager;

    /// Fracture Manager
    std :: unique_ptr< FractureManager > fracManager;

    /// Contact Manager
    std :: unique_ptr< ContactManager > contactManager;
    
    
    /**
     * Map from an element's global number to its place
     * in the element array. Added by ES 140326.
     */
    std::unordered_map< int, int > mElementPlaceInArray;

    /**
     * Map from a dofmans's global number to its place
     * in the dofman array.
     */
    std::unordered_map< int, int > mDofManPlaceInArray;

    /**
     * Map from material number to elements that have the
     * given material number. Added by ES 140718.
     */
    std::unordered_map< int, IntArray> mMapMaterialNum2El;

    /// Topology description
    std :: unique_ptr< TopologyDescription > topology;

    /// Keeps track of next free dof ID (for special Lagrange multipliers, XFEM and such)
    int freeDofID;

#ifdef __PARALLEL_MODE
    /**
     * Transaction manager. The purpose of this class is to
     * make the domain modification (in terms of adding and deleting components) versatile.
     */
    std :: unique_ptr< DomainTransactionManager > transactionManager;
    /// Global dof manager map (index is global of man number).
    std :: map< int, DofManager * >dmanMap;
    /// dmanMap init flag.
    bool dmanMapInitialized;
    /// Global element map (index is global of man number).
    std :: map< int, Element * >elementMap;
    /// dmanMap init flag.
    bool elementMapInitialized;


    /**@name Load Balancing data structures */
    //@{
    /// List of received elements.
    std :: list< Element * >recvElemList; ///@todo This seems like dead, old unused code. Remove it? / Mikael
    //@}
#endif

public:
    /**
     * Constructor. Creates empty n-th domain belonging to given engineering model.
     * @param n Domain number.
     * @param serNum Serial number
     * @param e Engineering model domain will belong to.
     * @see giveSerialNumber
     */
    Domain(int n, int serNum, EngngModel * e);

    Domain(const Domain& src) = delete;
    Domain &operator = (const Domain &src) = delete;

    /// Create a copy of the domain using the dynamic data reader.
    Domain *Clone();

    /// Destructor.
    ~Domain();

    /// Returns domain number.
    int giveNumber() { return this->number; }
    /// Returns domain number.
    void setNumber(int nn) { this->number = nn; }
    /// Returns domain serial (version) number.
    int giveSerialNumber() { return this->serialNumber; }

    // management of the mesh components
    /**
     * Service for accessing particular domain fe element.
     * Generates error if no such element is defined.
     * @param n Pointer to n-th element is returned.
     */
    Element *giveElement(int n);
    std :: vector< std :: unique_ptr< Element > > &giveElements() { return this->elementList; }
    /**
     * Service for accessing particular domain fe element.
     * Generates error if no such element is defined.
     * @param n Pointer to the element with id n
     */
    Element *giveGlobalElement(int n);
    /**
     * Returns the array index of the element with global
     * number iGlobalElNum, so that it can be fetched by
     * calling giveElement. Returns -1 if not found.
     */
    int giveElementPlaceInArray(int iGlobalElNum) const;
    /**
     * Returns the array index of the dofman with global
     * number iGlobalDofManNum, so that it can be fetched by
     * calling giveDofManager. Returns -1 if not found.
     */
    int giveDofManPlaceInArray(int iGlobalDofManNum) const;
    /**
     * Returns array with indices of elements that have a
     * given material number.
     */
    const IntArray &giveElementsWithMaterialNum(int iMaterialNum) const;
    /**
     * Returns engineering model to which receiver is associated.
     */
    EngngModel *giveEngngModel();
    /**
     * Service for accessing particular domain load.
     * Generates error if no such load is defined.
     * @param n Pointer to n-th load is returned.
     */
    Load *giveLoad(int n);
    /**
     * Service for accessing particular domain bc.
     * Generates error if no such bc is defined.
     * @param n Pointer to n-th bc is returned.
     */
    GeneralBoundaryCondition *giveBc(int n);
    std :: vector< std :: unique_ptr< GeneralBoundaryCondition > > &giveBcs() { return this->bcList; }
    /**
     * Service for accessing particular domain ic.
     * Generates error if no such ic is defined.
     * @param n Pointer to n-th ic is returned.
     */
    InitialCondition *giveIc(int n);
    std :: vector< std :: unique_ptr< InitialCondition > > &giveIcs() { return this->icList; }

    /**
     * Service for accessing particular domain load time function.
     * Generates error if no such load time function is defined.
     * @param n Pointer to n-th load time function is returned.
     */
    Function *giveFunction(int n);
    std :: vector< std :: unique_ptr< Function > > &giveFunctions() { return this->functionList; }
    /**
     * Service for accessing particular domain material model.
     * Generates error if no such material model is defined.
     * @param n Pointer to n-th material model is returned.
     */
    Material *giveMaterial(int n);
    std :: vector< std :: unique_ptr< Material > > &giveMaterials() { return this->materialList; }
    /**
     * Service for accessing particular domain cross section model.
     * Generates error if no such cross section  model is defined.
     * @param n Pointer to n-th cross section is returned.
     */
    CrossSection *giveCrossSection(int n);
    std :: vector< std :: unique_ptr< CrossSection > > &giveCrossSections() { return this->crossSectionList; }
    /**
     * Service for accessing particular domain nonlocal barrier representation.
     * Generates error if no such barrier  model is defined.
     * @param n Pointer to n-th barrier is returned.
     */
    NonlocalBarrier *giveNonlocalBarrier(int n);
    /**
     * Service for accessing particular domain set.
     * Generates error if no such set is defined.
     * @param n Pointer to n-th object is returned.
     */
    Set *giveSet(int n);
    std :: vector< std :: unique_ptr< Set > > &giveSets() { return this->setList; }
    /**
     * Service for accessing particular domain node.
     * Generates error if no such node is defined.
     * Note: nodes and element sides share common numbering (they are numbered as DofManagers).
     * @param n Pointer to n-th node is returned.
     */
    inline Node *giveNode(int n)
    {
    #ifdef DEBUG
        if ( n < 1 || n > (int)dofManagerList.size() ) {
            OOFEM_ERROR("undefined dofManager (%d)", n);
        }

        Node *node = reinterpret_cast< Node * >( dofManagerList[n-1].get() );

        return node;

    #else
        return reinterpret_cast< Node * >( dofManagerList[n-1].get() );
    #endif
    }
    /**
     * Service for accessing particular domain element side.
     * Generates error if no such element side is defined.
     * Note: nodes and element sides share common numbering (they are numbered as DofManagers).
     * @param n Pointer to n-th element side is returned.
     */
    ElementSide *giveSide(int n);
    /**
     * Service for accessing particular domain dof manager.
     * Generates error if no such dof manager is  defined.
     * Note: nodes and element sides share common numbering (they are numbered as DofManagers).
     * @param n Pointer to n-th dof manager is returned.
     */
    DofManager *giveDofManager(int n);
    std :: vector< std :: unique_ptr< DofManager > > &giveDofManagers() { return this->dofManagerList; }
    /**
     * Service for accessing particular domain dof manager.
     * Generates error if no such element is defined.
     * @param n Pointer to the element with id n
     */
    DofManager *giveGlobalDofManager(int n);
    /**
     * Reads receiver description from input stream and creates corresponding components accordingly.
     * It scans input file, each line is assumed to be single record describing type and parameters for
     * specific entity in domain. The record line is converted to lower case letters.
     * Corresponding component is created using classFactory.create* function of
     * corresponding base class, sending component name (extracted from corresponding record)
     * as parameter. After new object is created, its initializeFrom member function is
     * called with its record as parameter.
     * @param dr Input stream with domain description.
     * @return Nonzero if o.k.
     * @see FemComponent::initializeFrom
     */
    int instanciateYourself(DataReader *dr);
    /**
     * Performs post-initialization for all the domain contents (which is called after initializeFrom).
     * Currently, it only calls Element::postInitialize.
     */
    void postInitialize();
    /**
     * Automatically detects necessary nodal dofs and creates them accordingly.
     * Scans every element after its requested dof's and picks the union of all those dof types.
     * Intenal DOF managers are not affected, as those are created by the corresponding element/bc.
     */
    void createDofs();
    /// Returns number of dof managers in domain.
    int giveNumberOfDofManagers() const { return (int)dofManagerList.size(); }
    /// Returns number of elements in domain.
    int giveNumberOfElements() const { return (int)elementList.size(); }
    /// Returns number of material models in domain.
    int giveNumberOfMaterialModels() const { return (int)materialList.size(); }
    /// Returns number of cross section models in domain.
    int giveNumberOfCrossSectionModels() const { return (int)crossSectionList.size(); }
    /// Returns number of boundary conditions in domain.
    int giveNumberOfBoundaryConditions() const { return (int)bcList.size(); }
    /// Returns number of initial conditions in domain.
    int giveNumberOfInitialConditions() const { return (int)icList.size(); }
    /// Returns number of load time functions in domain.
    int giveNumberOfFunctions() const { return (int)functionList.size(); }
    /// Returns number of regions. Currently regions corresponds to cross section models.
    int giveNumberOfRegions() const { return this->giveNumberOfCrossSectionModels(); }
    /// Returns number of nonlocal integration barriers
    int giveNumberOfNonlocalBarriers() const { return (int)nonlocalBarrierList.size(); }
    /// Returns number of sets
    int giveNumberOfSets() const { return (int)setList.size(); }

    /// Returns number of spatial dimensions.
    int giveNumberOfSpatialDimensions();
    /// Returns true of axisymmetry is in effect.
    bool isAxisymmetric();
    /**
     * @name Advanced domain manipulation methods.
     */
    //@{
    /// Resizes the internal data structure to accommodate space for _newSize dofManagers.
    void resizeDofManagers(int _newSize);
    /// Resizes the internal data structure to accommodate space for _newSize elements.
    void resizeElements(int _newSize);
    /// Resizes the internal data structure to accommodate space for _newSize cross section models.
    void resizeCrossSectionModels(int _newSize);
    /// Resizes the internal data structure to accommodate space for _newSize materials.
    void resizeMaterials(int _newSize);
    /// Resizes the internal data structure to accommodate space for _newSize nonlocal barriers.
    void resizeNonlocalBarriers(int _newSize);
    /// Resizes the internal data structure to accommodate space for _newSize boundary conditions.
    void resizeBoundaryConditions(int _newSize);
    /// Resizes the internal data structure to accommodate space for _newSize initial conditions.
    void resizeInitialConditions(int _newSize);
    /// Resizes the internal data structure to accommodate space for _newSize load time functions.
    void resizeFunctions(int _newSize);
    /// Resizes the internal data structure to accommodate space for _newSize sets.
    void resizeSets(int _newSize);

    /// Sets i-th component. The component will be further managed and maintained by domain object.
    void setDofManager(int i, DofManager *obj);
    /// Sets i-th component. The component will be further managed and maintained by domain object.
    void setElement(int i, Element *obj);
    /// Sets i-th component. The component will be further managed and maintained by domain object.
    void setCrossSection(int i, CrossSection *obj);
    /// Sets i-th component. The component will be further managed and maintained by domain object.
    void setMaterial(int i, Material *obj);
    /// Sets i-th component. The component will be further managed and maintained by domain object.
    void setNonlocalBarrier(int i, NonlocalBarrier *obj);
    /// Sets i-th component. The component will be further managed and maintained by domain object.
    void setBoundaryCondition(int i, GeneralBoundaryCondition *obj);
    /// Sets i-th component. The component will be further managed and maintained by domain object.
    void setInitialCondition(int i, InitialCondition *obj);
    /// Sets i-th component. The component will be further managed and maintained by domain object.
    void setFunction(int i, Function *obj);
    /// Sets i-th component. The component will be further managed and maintained by domain object.
    void setSet(int i, Set *obj);

    /// Temporary function, sets xfemManager.
    void setXfemManager(XfemManager *ipXfemManager);

    XfemManager *giveXfemManager();
    bool hasXfemManager();

    ContactManager *giveContactManager();
    bool hasContactManager();
    
    FractureManager *giveFractureManager();
    bool hasFractureManager();

    /**
     * Sets receiver's associated topology description.
     * @param topo New topology description for receiver.
     * @param destroyOld Determines if any preexisting topology description should be deleted.
     */
    void setTopology(TopologyDescription *topo, bool destroyOld = true);
    /// Clear all boundary conditions.
    void clearBoundaryConditions();
    /// Clear all elements
    void clearElements();
    /// Clear receiver.
    void clear();
    //@}
    /**
     * Stores the domain state to output stream. Stores recursively the state of all
     * managed objects, like DofManagers and Elements.
     * Stored context is associated with current time step. One time step can have only
     * one associated context. Multiple call to saveContext within same time step
     * override previously saved context for this step.
     * By default the stream parameter is used to store data and is not closed.
     * If stream is NULL, new file descriptor is created and this must be also closed at the end.
     * @param stream Context stream. If NULL then new file descriptor will be opened and closed
     * at the end else the stream given as parameter will be used and not closed at the end.
     * @param mode Determines amount of info in stream.
     * @param obj Void pointer to an int array containing two values:time step number and
     * version of a context file to be restored.
     * @return contextIOResultType.
     * @exception ContextIOERR If error encountered.
     */
    contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores the domain state from output stream. Restores recursively the state of all
     * managed objects, like DofManagers and Elements.
     * Each context is associated with unique time step. Only one context per time step is
     * allowed. Restore context function will restore such context, which is related
     * (through its step number) to time step number and version given in obj parameter.
     * Restoring context will change current time step in order to correspond to newly restored
     * context.
     * @param stream Context file.
     * @param mode Determines amount of info in stream.
     * @param obj Void pointer to an int array containing two values:time step number and
     * version of a context file to be restored.
     * @return contextIOResultType.
     * @exception ContextIOERR exception if error encountered.
     */
    contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    /**
     * Returns default DofID array which defines physical meaning of particular DOFs.
     * of nodal dofs. Default values are determined using current domain type.
     */
    const IntArray &giveDefaultNodeDofIDArry();
    /// Returns domain type.
    domainType giveDomainType() { return dType; }
    /// Sets domain type
    void setDomainType(domainType _dType) { this->dType = _dType; }
    /**
     * Checks internal consistency of domain and all domain components.
     * The checkConsistency of all domain components is invoked.
     * @return nonzero if test is o.k.
     * @see FEMComponent::checkConsistency
     */
    int checkConsistency();

    /**
     * Gives the sum of the area of all elements.
     * @return Total area.
     */
    double giveArea();
    /**
     * Gives the sum of the volume of all elements.
     * @return Total volume.
     */
    double giveVolume();
    /**
     * Gives the sum of the volume or area of all elements.
     * @return Total volume.
     */
    double giveSize();

    /**
     * Gives the next free dof ID.
     * Useful for XFEM and other boundary conditions that introduce other unique Lagrange multipliers.
     * @return The next free dof ID.
     */
    int giveNextFreeDofID(int increment = 1);
    /**
     * Resets the free dof IDs.
     */
    void resetFreeDofID();
    /**
     * Gives the current maximum dof ID used.
     */
    int giveMaxDofID() { return this->freeDofID - 1; }

    /**
     * Returns receiver's associated connectivity table.
     */
    ConnectivityTable *giveConnectivityTable();
    /**
     * Returns receiver's associated spatial localizer.
     */
    SpatialLocalizer *giveSpatialLocalizer();
    /**
     * Returns domain output manager.
     */
    OutputManager *giveOutputManager();
    /**
     * Returns Error Estimator associated to receiver.
     * Calls corresponding EngngModel Service.
     */
    ErrorEstimator *giveErrorEstimator();
    /**
     * Returns the actual Smoother associated to receiver.
     * Creates the default, if no one associated.
     */
    NodalRecoveryModel *giveSmoother();
    /**
     * Returns receiver's associated topology description.
     */
    TopologyDescription *giveTopology();
    /**
     * Sets the given smoother as an actual smoother for receiver.
     * @param smoother New smoother for receiver.
     * @param destroyOld Determines if any preexisting smoother should be deleted.
     */
    void setSmoother(NodalRecoveryModel *newSmoother, bool destroyOld = true);

#ifdef __PARALLEL_MODE
    /**@name Domain transaction support methods.
     * The purpose of these methods is to provide a unified approach
     * for changing domain at runtime (meaning mainly adding and deleting dofmanagers and elements).
     * The changes are recorded in transaction manager and until the are committed,
     * no change is reflected in domain itself.
     */
    //@{
    /**
     * Returns domain transaction manager.
     */
    DomainTransactionManager *giveTransactionManager();
    /**
     * Commits transactions recorded in transaction manager. The purpose of transaction manager is to
     * make the domain modification (adding and deleting components) possible and versatile.
     *
     * The changes are recorded in transaction manager and until the are committed,
     * no change is reflected in domain itself. After transactions are committed, the local numbering can change.
     * A message to the system is sent to update the numbering.
     * @param tm Manager to commit transactions to.
     */
    int commitTransactions(DomainTransactionManager *tm);

    /**
     * Initializes global dof man map according to domain dofman list.
     */
    void initGlobalDofManMap(bool forceinit = false);
    void initGlobalElementMap(bool forceinit = false);
    /**
     * Assigns new local number (stored as dofmanager number, so it can be requested)
     * to all dofManagers available in domanMap.
     */
    void renumberDofManagers();
    void renumberDofManData(DomainTransactionManager *tm);
    void renumberElements();
    void renumberElementData(DomainTransactionManager *tm);
    /** Return updated local entity number after load balancing */
    int LB_giveUpdatedLocalNumber(int oldnum, EntityRenumberingScheme scheme);
    /** Return updated local entity number after load balancing */
    int LB_giveUpdatedGlobalNumber(int oldnum, EntityRenumberingScheme scheme);
    //@}


    /**@name Load Balancing support methods */
    //@{
    int dofmanGlobal2Local(int _globnum);
    int elementGlobal2Local(int _globnum);
    //@}
#endif

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Domain"; }

    /// Returns the value of nonlocalUpdateStateCounter
    StateCounterType giveNonlocalUpdateStateCounter() { return this->nonlocalUpdateStateCounter; }
    /// sets the value of nonlocalUpdateStateCounter
    void setNonlocalUpdateStateCounter(StateCounterType val) { this->nonlocalUpdateStateCounter = val; }

private:
    void resolveDomainDofsDefaults(const char *);

    /// Returns string for prepending output (used by error reporting macros).
    std :: string errorInfo(const char *func) const;

private:
    /**
     * Construct map from an element's global number to
     * its place the element array.
     */
    void BuildElementPlaceInArrayMap();

    /**
     * Construct map from a dofman's global number to
     * its place the dofman array.
     */
    void BuildDofManPlaceInArrayMap();

    /**
     * Construct map from a material number to
     * elements with the given material number.
     */
    void BuildMaterialToElementMap();
};
} // end namespace oofem
#endif // domain_h
