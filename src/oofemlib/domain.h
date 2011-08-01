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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef domain_h
#define domain_h

#include "alist.h"
#include "datareader.h"
#include "domaintype.h"
#include "statecountertype.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <time.h>
 #include <map>
 #ifdef __PARALLEL_MODE
  #include <list>
  #include "entityrenumberingscheme.h"
 #endif
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
class Element;
class Node;
class Material;
class TimeStep;
class GeneralBoundaryCondition;
class InitialCondition;
class Load;
class LoadTimeFunction;
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
class RandomFieldGenerator;
class XfemManager;
class TopologyDescription;

#ifdef __PARALLEL_MODE
class ProcessCommunicator;
class LoadBalancer;
#endif
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
class Domain
{
private:
    /// Element list.
    AList< Element > *elementList;
    /// Dof manager list.
    AList< DofManager > *dofManagerList;
    /// Material list.
    AList< Material > *materialList;
    /// Cross section list.
    AList< CrossSection > *crossSectionList;
    /// Boundary condition list.
    AList< GeneralBoundaryCondition > *bcList;
    /// Initial condition list.
    AList< InitialCondition > *icList;
    /// Load time function list.
    AList< LoadTimeFunction > *loadTimeFunctionList;
    /// Nonlocal barrier list.
    AList< NonlocalBarrier > *nonlocalBarierList;
    /// List of Random generators.
    AList< RandomFieldGenerator > *randomFieldGeneratorList;

    // numberOfDefaultDofsPerNode specifies default number of dofs per node
    // for current domain type. The defaultDofMask describes the physical meaning of these
    // dofs.
    int numberOfDefaultDofsPerNode;
    int numberOfDefaultDofsPerSide;
    IntArray defaultNodeDofIDArry;
    IntArray defaultSideDofIDArry;

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
    ConnectivityTable *connectivityTable;
    /**
     * Spatial Localizer. It is build upon request.
     * Provides the spatial localization services.
     */
    SpatialLocalizer *spatialLocalizer;
    /// Output manager, allowing to filter the produced output.
    OutputManager *outputManager;
    /// Domain number.
    int number;
    /// Domain serial (version) number. Used for domain version identification during Adaptive computations.
    int serialNumber;
    /// nodal recovery object associated to receiver.
    NodalRecoveryModel *smoother;
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
    /// XFEM Manager.
    XfemManager *xfemManager;
    /// Topology description
    TopologyDescription *topology;

#ifdef __PARALLEL_MODE
    /**
     * Transaction manager. The purpose of this class is to
     * make the domain modification (in terms of adding and deleting components) versatile.
     */
    DomainTransactionManager *transactionManager;
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
    std :: list< Element * >recvElemList;
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
    Domain(int n, int serNum, EngngModel *e);
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
    /**
     * Service for accessing particular domain ic.
     * Generates error if no such ic is defined.
     * @param n Pointer to n-th ic is returned.
     */
    InitialCondition *giveIc(int n);

    /**
     * Service for accessing particular domain load time function.
     * Generates error if no such load time function is defined.
     * @param n Pointer to n-th load time function is returned.
     */
    LoadTimeFunction *giveLoadTimeFunction(int n);
    /**
     * Service for accessing particular domain material model.
     * Generates error if no such material model is defined.
     * @param n Pointer to n-th material model is returned.
     */
    Material *giveMaterial(int n);
    /**
     * Service for accessing particular domain cross section model.
     * Generates error if no such cross section  model is defined.
     * @param n Pointer to n-th cross section is returned.
     */
    CrossSection *giveCrossSection(int n);
    /**
     * Service for accessing particular domain nonlocal barrier representation.
     * Generates error if no such barrier  model is defined.
     * @param n Pointer to n-th barrier is returned.
     */
    NonlocalBarrier *giveNonlocalBarrier(int n);
    /**
     * Service for accessing particular domain random field generator.
     * Generates error if no such generator is defined.
     * @param n Pointer to n-th object is returned.
     */
    RandomFieldGenerator *giveRandomFieldGenerator(int n);

    /**
     * Service for accessing particular domain node.
     * Generates error if no such node is defined.
     * Note: nodes and element sides share common numbering (they are numbered as DofManagers).
     * @param n Pointer to n-th node is returned.
     */
    Node *giveNode(int n);
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
    /**
     * Reads receiver description from input stream and creates corresponding components accordingly.
     * It scans input file, each line is assumed to be single record describing type and parameters for
     * specific entity in domain. The record line is converted to lower case letters.
     * Corresponding component is created using ofType member function of
     * corresponding base class, sending component name (extracted from corresponding record)
     * as parameter. After new object is created, its initializeFrom member function is
     * called with its record as parameter. The ofType member functions are using global
     * user modifiable functions to allow simple extension of library. See base classes documentation for details.
     * @param dr Input stream with domain description.
     * @return Nonzero if o.k.
     * @see FemComponent::initializeFrom
     */
    int instanciateYourself(DataReader *dr);
    //int giveNumberOfNodes () {return nodeList->giveSize();}
    //int giveNumberOfSides () {return elementSideList->giveSize();}
    /// Returns number of dof managers in domain.
    int giveNumberOfDofManagers() { return dofManagerList->giveSize(); }
    /// Returns number of elements in domain.
    int giveNumberOfElements() { return elementList->giveSize(); }
    /// Returns number of material models in domain.
    int giveNumberOfMaterialModels() { return materialList->giveSize(); }
    /// Returns number of cross section models in domain.
    int giveNumberOfCrossSectionModels() { return crossSectionList->giveSize(); }
    /// Returns number of boundary conditions in domain.
    int giveNumberOfBoundaryConditions() { return bcList->giveSize(); }
    /// Returns number of initial conditions in domain.
    int giveNumberOfInitialConditions() { return icList->giveSize(); }
    /// Returns number of load time functions in domain.
    int giveNumberOfLoadTimeFunctions() { return loadTimeFunctionList->giveSize(); }
    /// Returns number of regions. Currently regions corresponds to cross section models.
    int giveNumberOfRegions() { return this->giveNumberOfCrossSectionModels(); }
    /// Returns number of nonlocal integration barriers
    int giveNumberOfNonlocalBarriers() { return nonlocalBarierList->giveSize(); }
    int giveCorrespondingCoordinateIndex(int);
    /// Returns number of spatial dimensions.
    int giveNumberOfSpatialDimensions();
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
    void resizeLoadTimeFunctions(int _newSize);

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
    void setLoadTimeFunction(int i, LoadTimeFunction *obj);
    /// Temporary function, sets xfemManager.
    void setXfemManager(XfemManager *xfemManager);
    /**
     * Sets receiver's associated topology description.
     * @param topo New topology description for receiver.
     * @param destroyOld Determines if any preexisting topology description should be deleted.
     */
    void setTopology(TopologyDescription *topo, bool destroyOld = true);
    /// Clear all boundary conditions.
    void clearBoundaryConditions();
    //@}
    /**
     * Returns default DofID array which defines physical meaning of particular DOFs.
     * of nodal dofs. Default values are determined using current domain type.
     */
    const IntArray &giveDefaultNodeDofIDArry();
    /**
     * Returns default DofID array which defines physical meaning of particular DOFs.
     * of element side dofs. Default values are determined using current domain type.
     */
    const IntArray &giveDefaultSideDofIDArry();
    /**
     * Returns default number of dofs per node. Determined using current domain type.
     */
    int giveNumberOfDefaultNodeDofs();
    /**
     * Returns default number of dofs per side. Determined using current domain type.
     */
    int giveNumberOfDefaultSideDofs();

    /// Returns domain type.
    domainType giveDomainType() { return dType; }
    /// Sets domain type
    void setDomainType(domainType _dType) { this->dType = _dType; }
    // consistency check
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
    OutputManager *giveOutputManager() { return outputManager; }
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
    TopologyDescription *giveTopology() { return topology; }
    /**
     * Sets the given smoother as an actual smoother for receiver.
     * @param smoother New smoother for receiver.
     * @param destroyOld Determines if any preexisting smoother should be deleted.
     */
    void setSmoother(NodalRecoveryModel *smoother, bool destroyOld = true);

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

    void initGlobalDofManMap(bool forceinit = false);
    void initGlobalElementMap(bool forceinit = false);
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

#ifdef __OOFEG
    void drawYourself(oofegGraphicContext &context);
    void drawElements(oofegGraphicContext &context);
    void drawNodes(oofegGraphicContext &context);
#endif
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Domain"; }

    /// Returns the value of nonlocalUpdateStateCounter
    StateCounterType giveNonlocalUpdateStateCounter() { return this->nonlocalUpdateStateCounter; }
    /// sets the value of nonlocalUpdateStateCounter
    void setNonlocalUpdateStateCounter(StateCounterType val) { this->nonlocalUpdateStateCounter = val; }

private:
    void resolveDomainDofsDefaults(char *);

    void error(const char *file, int line, const char *format, ...);
    void warning(const char *file, int line, const char *format, ...);
};
} // end namespace oofem
#endif // domain_h

