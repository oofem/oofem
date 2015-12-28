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

#ifndef subdivision_h
#define subdivision_h

#include "mesherinterface.h"
#include "floatarray.h"
#include "intarray.h"
#include "element.h"
#include "dofmanager.h"
#include "connectivitytable.h"

#include <memory>
#include <vector>
#include <queue>
#include <list>

namespace oofem {

#define SHARED_IRREGULAR_DATA_TAG 7654
#define SUBDIVISION_SHARED_IRREGULAR_REC_TAG 7655
#define SUBDIVISION_END_DATA 7656
#define SUBDIVISION_MIGRATE_REMOTE_ELEMENTS_TAG 7657
#define SHARED_EDGE_DATA_TAG 7658
#define SUBDIVISION_SHARED_EDGE_REC_TAG 7659

class TimeStep;
class ProcessCommunicator;

/**
 * This class represents the Rivara Subdivision algorithm for triangular meshes.
 * based on M.C. Rivara. Local modification of meshes for adaptive and/or multigrid
 * finite-element methods. J. Comput. Appl. Math., 36:79–89, 1991.
 */
class OOFEM_EXPORT Subdivision : public MesherInterface
{
protected:
    class RS_Mesh;
    class RS_Node
    {
protected:
        FloatArray coords;
        double requiredDensity;
        int number;
        int parent; // number of parent node or zero for new nodes
        bool boundary;
        // mesh
        RS_Mesh *mesh;
        // Sorted list of local elements sharing the node
        // if it is shared (in such a case the list is built a priori)
        // if it is on boundary of 3D region and incident to edge subjected to bisection
        // (in such a case the list is built on the fly)
        IntArray connectedElements;
        int globalNumber;
#ifdef __PARALLEL_MODE
        dofManagerParallelMode parallel_mode;
        /**
         * List of partition sharing the shared dof manager or
         * remote partion containing remote dofmanager counterpart.
         */
        IntArray partitions;
#endif
public:
        RS_Node(int n, Subdivision :: RS_Mesh * m, int parent, FloatArray & c, double rd, bool boundary) {
            this->number = n;
            this->mesh = m;
            this->coords = c;
            this->requiredDensity = rd;
            this->boundary = boundary;
            this->parent = parent;
            this->globalNumber  = 0;
#ifdef __PARALLEL_MODE
            this->parallel_mode = DofManager_local;
#endif
        }
        virtual ~RS_Node() { }
        double giveRequiredDensity() { return requiredDensity; }
        FloatArray *giveCoordinates() { return & coords; }
        double giveCoordinate(int i) { return coords.at(i); }
        int giveParent() { return this->parent; }
        bool isBoundary() { return this->boundary; }
        void setBoundary(bool b) { this->boundary = b; }
        int giveNumber() { return this->number; }
        void setNumber(int _n) { this->number = _n; }
        int buildTopLevelNodeConnectivity(ConnectivityTable *ct);
        const IntArray *giveConnectedElements()  { return & connectedElements; }
        void setConnectedElements(IntArray _conn) { connectedElements = std :: move(_conn); }
        void insertConnectedElement(int num) { connectedElements.insertSorted(num, 10); }
        void eraseConnectedElement(int num) { connectedElements.eraseSorted(num); }
        void preallocateConnectedElements(int size) { connectedElements.preallocate(size); }
        int giveGlobalNumber() { return globalNumber; }
        void setGlobalNumber(int gn) { this->globalNumber = gn; }
        virtual bool isIrregular() { return false; }

#ifdef __PARALLEL_MODE
        void numberSharedEdges();
        dofManagerParallelMode giveParallelMode() const { return parallel_mode; }
        void setParallelMode(dofManagerParallelMode _mode) { parallel_mode = _mode; }
        const IntArray *givePartitions()  { return & partitions; }
        void setPartitions(IntArray _p) { partitions = std :: move(_p); }
        int importConnectivity(ConnectivityTable *ct);
#endif
#ifdef __OOFEG
        void  drawGeometry();
#endif
    };

    class RS_IrregularNode : public RS_Node
    {
protected:
        int iNode, jNode; // parent edge nodes
public:
        RS_IrregularNode(int n, Subdivision :: RS_Mesh * mesh, int parent, FloatArray & c, double rd, bool boundary) : RS_Node(n, mesh, parent, c, rd, boundary) { }
        void setEdgeNodes(int i, int j) {
            iNode = i;
            jNode = j;
        }
        void giveEdgeNodes(int &i, int &j) {
            i = iNode;
            j = jNode;
        }
        virtual bool isIrregular() { return true; }
    };

    class RS_Element
    {
protected:
        int number;
        // element regular nodes
        IntArray nodes;
        // element neighbours (on the same level of refinement)
        IntArray neghbours_base_elements;
        // irregular nodes associated to corresponding element edge
        IntArray irregular_nodes;
        // children
        IntArray children;
        // parent
        int parent;
        // longest edge index
        int leIndex;
        // mesh
        RS_Mesh *mesh;
        // flag whether element is in bisection queue
        bool queue_flag;
        int globalNumber;
#ifdef __PARALLEL_MODE
        elementParallelMode parallel_mode;
        // numbers of shared edges
        IntArray shared_edges;
#endif
public:
        RS_Element(int number, Subdivision :: RS_Mesh * m, int parent, IntArray & nodes) {
            this->number = number;
            this->nodes = nodes;
            this->leIndex = 0;
            this->mesh = m;
            this->parent = parent;
            this->queue_flag = false;
            this->globalNumber = -1;
#ifdef __PARALLEL_MODE
            this->parallel_mode = Element_local;
#endif
        }
        virtual ~RS_Element() { }

        /// Returns true if element has some irregular nodes
        bool hasIrregulars() { return !irregular_nodes.containsOnlyZeroes(); }
        /// Returns true if receiver is terminal (not further subdivided)
        bool isTerminal() { return children.isEmpty(); }

        int giveIrregular(int iedge) { return irregular_nodes.at(iedge); }
        void setIrregular(int iedge, int ir) { this->irregular_nodes.at(iedge) = ir; }

        virtual int evaluateLongestEdge() { return 0; }
        virtual void bisect(std :: queue< int > &subdivqueue, std :: list< int > &sharedIrregularsQueue) { }
        virtual void generate(std :: list< int > &sharedEdgesQueue) { }
        virtual void update_neighbours() { }
        virtual double giveDensity() { return 0.0; }
        virtual double giveRequiredDensity();
        const IntArray *giveChildren() { return & this->children; }
        virtual bool isNeighborOf(Subdivision :: RS_Element *elem) = 0;
        const IntArray *giveNeighbors() { return & this->neghbours_base_elements; }
        int giveNeighbor(int iside) { return neghbours_base_elements.at(iside); }
        void setNeighbor(int iside, int nb) { this->neghbours_base_elements.at(iside) = nb; }
        bool containsNode(int _node) { return nodes.findFirstIndexOf(_node) > 0; }
        virtual void giveSideNodes(int iside, IntArray &snodes) = 0;
        int giveParent() { return this->parent; }
        int giveTopParent();
        bool giveQueueFlag() { return this->queue_flag; }
        void setQueueFlag(bool _qf) { this->queue_flag = _qf; }
        void buildTopLevelNodeConnectivity(Subdivision :: RS_Node *node);
        virtual void importConnectivity(ConnectivityTable *ct) = 0;
        const IntArray *giveNodes() { return & this->nodes; }
        int giveNode(int i) { return this->nodes.at(i); }
        int giveNumber() { return this->number; }
        void setNumber(int newNum) { this->number = newNum; }
        virtual int giveEdgeIndex(int iNode, int jNode) = 0;
        /// Returns the longest edge index of the receiver
        int giveLeIndex() { return this->leIndex; }
        /// Sets the longest edge index
        void setLeIndex(int _n) { this->leIndex = _n; }

#ifdef __OOFEG
        virtual void  drawGeometry() { }
#endif
        int giveGlobalNumber() { return globalNumber; }
        void setGlobalNumber(int gn) { this->globalNumber = gn; }

#ifdef __PARALLEL_MODE
        virtual void numberSharedEdges(int iNode, IntArray &connNodes) = 0;
        const IntArray *giveSharedEdges()  { return & shared_edges; }
        virtual void makeSharedEdges() = 0;
        int giveSharedEdge(int iedge) { return shared_edges.at(iedge); }
        void setSharedEdge(int iedge, int num) { this->shared_edges.at(iedge) = num; }
        elementParallelMode giveParallelMode() const { return parallel_mode; }
        // Sets parallel mode of element
        void setParallelMode(elementParallelMode _mode) { parallel_mode = _mode; }
#endif
    };

    class RS_Triangle : public Subdivision :: RS_Element
    {
public:
        RS_Triangle(int number, Subdivision :: RS_Mesh * mesh, int parent, IntArray & nodes);
        int evaluateLongestEdge();
        void bisect(std :: queue< int > &subdivqueue, std :: list< int > &sharedIrregularsQueue);
        void generate(std :: list< int > &sharedEdgesQueue);
        void update_neighbours();
        double giveDensity();
        bool isNeighborOf(Subdivision :: RS_Element *elem);
        void giveSideNodes(int iside, IntArray &snodes);
        int giveEdgeIndex(int iNode, int jNode);
        virtual void importConnectivity(ConnectivityTable *ct);
#ifdef __OOFEG
        void drawGeometry();
#endif
#ifdef __PARALLEL_MODE
        void numberSharedEdges(int iNode, IntArray &connNodes);
        void makeSharedEdges() {
            shared_edges.resize(3);
            shared_edges.zero();
        }
#endif
    };

    class RS_Tetra : public Subdivision :: RS_Element
    {
protected:
        IntArray side_leIndex;
public:
        RS_Tetra(int number, Subdivision :: RS_Mesh * mesh, int parent, IntArray & nodes);
        int evaluateLongestEdge();
        void bisect(std :: queue< int > &subdivqueue, std :: list< int > &sharedIrregularsQueue);
        void generate(std :: list< int > &sharedEdgesQueue);
        void update_neighbours();
        double giveDensity();
        bool isNeighborOf(Subdivision :: RS_Element *elem);
        void giveSideNodes(int iside, IntArray &snodes);
        int giveEdgeIndex(int iNode, int jNode);
        virtual void importConnectivity(ConnectivityTable *ct);
#ifdef __OOFEG
        void drawGeometry();
#endif
#ifdef __PARALLEL_MODE
        void numberSharedEdges(int iNode, IntArray &connNodes);
        void makeSharedEdges() {
            shared_edges.resize(6);
            shared_edges.zero();
        }
#endif
    };

#ifdef __PARALLEL_MODE
    class RS_SharedEdge
    {
protected:
        // parent edge nodes
        int iNode, jNode;
        // shared partitions except the local one
        IntArray partitions;
        // mesh
        RS_Mesh *mesh;
public:
        RS_SharedEdge(Subdivision :: RS_Mesh * m) {
            this->mesh = m;
        }
        void setEdgeNodes(int i, int j) {
            iNode = i;
            jNode = j;
        }
        void giveEdgeNodes(int &i, int &j) {
            i = iNode;
            j = jNode;
        }
        const IntArray *givePartitions()  { return & partitions; }
        void setPartitions(IntArray _p) { partitions = std :: move(_p); }
        void addPartition(int _p, int allocChunk) { partitions.followedBy(_p, allocChunk); }
        void removePartitions() { partitions.clear(); }
        int giveSharedPartitions(IntArray &partitions);
    };
#endif

    class RS_Mesh
    {
        // HUHU protected? private?
        std :: vector< std :: unique_ptr< Subdivision :: RS_Node > >nodes;
        std :: vector< std :: unique_ptr< Subdivision :: RS_Element > >elements;
#ifdef __PARALLEL_MODE
        std :: vector< std :: unique_ptr< Subdivision :: RS_SharedEdge > >edges;
#endif
        Subdivision *subdivision;

#ifdef __PARALLEL_MODE
        /// Global shared node map (index is global shared node number)
        std :: map< int, Subdivision :: RS_Node * >sharedNodeMap;
        /// sharedNodeMap init flag
        bool sharedNodeMapInitialized;
#endif

public:
#ifdef __PARALLEL_MODE
        RS_Mesh(Subdivision * s) : nodes(), elements(), edges() {
            this->subdivision = s;
            sharedNodeMapInitialized = false;
        }
#else
        RS_Mesh(Subdivision * s) : nodes(), elements() {
            this->subdivision = s;
        }
#endif
        ~RS_Mesh() {}

        Subdivision :: RS_Node *giveNode(int i) { return nodes[i-1].get(); }
        Subdivision :: RS_Element *giveElement(int i) { return elements[i-1].get(); }
        int giveNumberOfNodes() { return (int)nodes.size(); }
        int giveNumberOfElements() { return (int)elements.size(); }
        void addNode(Subdivision :: RS_Node *obj) { nodes.emplace_back(obj); }
        void addElement(Subdivision :: RS_Element *obj) { elements.emplace_back(obj); }
#ifdef __PARALLEL_MODE
        Subdivision :: RS_SharedEdge *giveEdge(int i) { return edges[i-1].get(); }
        int giveNumberOfEdges() { return (int)edges.size(); }
        void addEdge(Subdivision :: RS_SharedEdge *obj) { edges.emplace_back(obj); }
        void initGlobalSharedNodeMap() { sharedNodeMap.clear(); }
        void insertGlobalSharedNodeMap(Subdivision :: RS_Node *node);
        int sharedNodeGlobal2Local(int _globnum);
#endif

        Subdivision *giveSubdivision() { return this->subdivision; }
    };

    class RS_CompareNodePositions
    {
        RS_Mesh *m;
public:
        RS_CompareNodePositions(RS_Mesh * _m) {
            m = _m;
        }
        int operator() (int i, int j);
    };

    struct RS_packRemoteElemsStruct {
        Domain *d;                      // HUHU proc, subdivision ma pristup na domain
        IntArray *parentElemMap;
    };

    //RS_Mesh *oldMesh, *newMesh;
    RS_Mesh *mesh;
    std :: queue< int >subdivqueue;
    std :: list< int >sharedIrregularsQueue;
    std :: list< int >sharedEdgesQueue;
    // smoothing flag
    bool smoothingFlag;

public:
    /// Constructor
    Subdivision(Domain * d) : MesherInterface(d) {
        mesh = 0;
        smoothingFlag = false;
    }
    virtual ~Subdivision() {
        if ( mesh ) {
            delete mesh;
        }
    }

    /// Runs the mesh generation, mesh will be written to corresponding domain din file
    virtual returnCode createMesh(TimeStep *tStep, int domainNumber, int domainSerNum, Domain **dNew);
    const char *giveClassName() { return "Subdivision"; }
    Domain *giveDomain() { return domain; }

protected:
    Subdivision :: RS_Mesh *giveMesh() { return mesh; }
    void bisectMesh();
    void smoothMesh();

    bool isNodeLocalIrregular(Subdivision :: RS_Node *node, int myrank);
    void assignGlobalNumbersToElements(Domain *d);

#ifdef __PARALLEL_MODE
    /**
     * Exchanges the shared irregulars between partitions. Returns true if any shared irregular has been
     * exchanged between any partitions.
     */
    bool exchangeSharedIrregulars();
    int packSharedIrregulars(Subdivision *s, ProcessCommunicator &pc);
    int unpackSharedIrregulars(Subdivision *s, ProcessCommunicator &pc);
    void assignGlobalNumbersToSharedIrregulars();
    int packIrregularSharedGlobnums(Subdivision *s, ProcessCommunicator &pc);
    int unpackIrregularSharedGlobnums(Subdivision *s, ProcessCommunicator &pc);
    void exchangeSharedEdges();
    int packSharedEdges(Subdivision *s, ProcessCommunicator &pc);
    int unpackSharedEdges(Subdivision *s, ProcessCommunicator &pc);

    /// Returns true if receiver is irregular, shared node locally maintatined
    bool isNodeLocalSharedIrregular(Subdivision :: RS_Node *node, int myrank);

    int giveRank();
    int giveNumberOfProcesses();

    void exchangeRemoteElements(Domain *d, IntArray &);
    int packRemoteElements(RS_packRemoteElemsStruct *s, ProcessCommunicator &pc);
    int unpackRemoteElements(Domain *d, ProcessCommunicator &pc);

#endif
};
} // end namespace oofem
#endif // subdivision_h
