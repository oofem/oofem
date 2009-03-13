/* $Header$ */
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

//   ******************************************
//   *** CLASS Rivara Subdivision algorithm ***
//   ******************************************

#ifndef subdivision_h
#define subdivision_h

#include "mesherinterface.h"
#include "flotarry.h"
#include "intarray.h"
#include "alist.h"
#include "conTable.h"
#ifndef __MAKEDEPEND
#include <queue>
#include <list>
#endif



#define SHARED_IRREGULAR_DATA_TAG 7654
#define SUBIVISION_IRREGULAR_REC_TAG 7655
#define SUBDIVISION_END_DATA 7656
#define SUBDIVISION_MIGRATE_REMOTE_ELEMENTS_TAG 7657

class TimeStep;

/**
 * This class represents the Rivara Subdivision algorithm for trinagular meshes.
 * based on M.C. Rivara. Local modification of meshes for adaptive and/or multigrid
 * finite-element methods. J. Comput. Appl. Math., 36:79–89, 1991.
 */
class Subdivision : public MesherInterface
{
 protected:
  class RS_Mesh;
  class RS_Node {
  protected:
    FloatArray coords;
    double requiredDensity;
    int number;
    int parent; // number of parent node or zero for new nodes
		bool boundary;
#ifdef __PARALLEL_MODE
    int globalNumber;
    dofManagerParallelMode parallel_mode;
    /**
     * List of partition sharing the shared dof manager or
     * remote partion containing remote dofmanager counterpart.
     */
    IntArray partitions;
#endif
  public:
    RS_Node (int n, int parent, FloatArray& c, double rd, bool boundary) {
      this->number=n; this->coords=c; this->requiredDensity = rd; this->boundary=boundary;
      this->parent=parent;
#ifdef __PARALLEL_MODE
      this->parallel_mode = DofManager_local;
      this->globalNumber  = 0;
#endif
    }
    virtual ~RS_Node() {}
    double giveRequiredDensity() {return requiredDensity;}
    FloatArray* giveCoordinates() {return &coords;}
    double giveCoordinate(int i) {return coords.at(i);}
    int giveParent () {return this->parent;}
		bool isBoundary () {return this->boundary;}
		void setBoundary (bool b) {this -> boundary=b;}
    int giveNumber() {return this->number;}
    void setNumber(int _n) {this->number = _n;}
#ifdef __PARALLEL_MODE
    dofManagerParallelMode giveParallelMode() const { return parallel_mode; }
    void setParallelMode(dofManagerParallelMode _mode) { parallel_mode = _mode; }
    void setPartition (int p) {this->partitions.resize(1); this->partitions.at(1)=p;}
    const IntArray *givePartitions ()  { return & partitions; }
    /** Sets receiver's partition list */
    void setPartitions (const IntArray &_p) { partitions = _p; }
    int giveGlobalNumber() {return globalNumber;}
    void setGlobalNumber(int gn) {this->globalNumber=gn;}
    virtual bool isIrregular() {return false;}
#endif
#ifdef __OOFEG
    void  drawGeometry ();
#endif
  };

  class RS_IrregularNode : public RS_Node {
  protected:
    int iNode, jNode; // paren edge nodes
  public:
    RS_IrregularNode (int n, int parent, FloatArray& c, double rd, bool boundary) : RS_Node (n, parent, c, rd, boundary) {}
    void setEdgeNodes (int i, int j) {iNode=i; jNode=j;}
    void giveEdgeNodes (int& i, int& j) {i=iNode; j=jNode;}
    virtual bool isIrregular() {return true;}
  };

  class RS_Element {
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
    RS_Mesh* mesh;
		// flag whether element is in bisection queue
		bool queue_flag;
#ifdef __PARALLEL_MODE
    elementParallelMode parallel_mode;
    int globalNumber;
#endif
  public:
    RS_Element(int number, Subdivision::RS_Mesh* m, int parent, IntArray &nodes) {
      this->number = number; this->nodes = nodes; this->leIndex=0; 
			this->mesh = m; this->parent = parent; this->queue_flag = false; 
#ifdef __PARALLEL_MODE
      this->globalNumber = -1;
      this->parallel_mode = Element_local;
			children.resize(0);
#endif
    }

    /// Returns true if element has some irregular nodes
    bool hasIrregulars() {return !irregular_nodes.containsOnlyZeroes();}
    /// Returns true if receiver is terminal (not further subdivided)
    bool isTerminal() {return children.isEmpty();}
    /**
       Add Irregular node to receiver, irregular identified by boundary entity identification
       (typically by an edge identified by its nodes 
    */
    void addIrregularOn (int iNum, const IntArray& bEnttity);
    /**
       Returns Irregular node of receiver identified by boundary entity identification
       (typically by an edge identified by its nodes 
       Returns zero if not found.
    */
    int giveIrregularOn (const IntArray& bEnttity);
		int giveIrregular(int iedge) {return irregular_nodes.at(iedge);}
		void setIrregular (int iedge, int ir) {this->irregular_nodes.at(iedge) = ir;}

		virtual int evaluateLongestEdge() {}
    virtual void bisect(std::queue<int> &subdivqueue, std::list<int> &sharedIrregularsQueue) {}
    virtual void generate () {}
    virtual void update_neighbours() {}
    virtual double giveDensity() {return 0.0;}
    virtual double giveRequiredDensity();
    void giveChildren (IntArray& c) {c=this->children;}
    virtual bool isNeighborOf (Subdivision::RS_Element* elem) = 0;
		const IntArray *giveNeighbors () {return &this->neghbours_base_elements;}
    int giveNeighbor (int iside) {return neghbours_base_elements.at(iside);}
    void setNeighbor (int iside, int nb) {this->neghbours_base_elements.at(iside) = nb;}
    bool containsNode (int _node) {return nodes.findFirstIndexOf (_node);}
		virtual void giveSideNodes (int iside, IntArray& snodes) = 0;
    int giveParent() {return this->parent;}
		int giveTopParent();
		bool giveQueueFlag() {return this->queue_flag;}
		void setQueueFlag(bool _qf) {this->queue_flag=_qf;}
    virtual void importConnectivities(ConnectivityTable* ct) = 0;
    const IntArray* giveNodes () {return &this->nodes;}
    int giveNode(int i) {return this->nodes.at(i);}
    int giveNumber() {return this->number;}
    void setNumber (int newNum) {this->number=newNum;}
		virtual int giveEdgeIndex (int iNode, int jNode) = 0;
    /// Returns the longest edge index of the receiver
    int giveLeIndex() {return this->leIndex;}
    /// Sets the longest edge index 
    void setLeIndex (int _n) {this->leIndex = _n;}

#ifdef __OOFEG
    virtual void  drawGeometry () {}
#endif

#ifdef __PARALLEL_MODE
    // returns number of shared remote partitions (returned in partitions param)
    virtual int giveIrregularSharedPartitions(int eIndex, int iNode, int jNode, IntArray& partitions) = 0;
    elementParallelMode giveParallelMode() const { return parallel_mode; }
    /// Sets parallel mode of element
    void setParallelMode(elementParallelMode _mode) { parallel_mode = _mode; }
    // return local number of given global counterpart, if contained within receiver, zero otherwise
    int  containsGlobalNode (int gnode);
    int giveGlobalNumber() {return globalNumber;}
    void setGlobalNumber(int gn) {this->globalNumber=gn;}
#endif

  };
  
  
  class RS_Triangle : public Subdivision::RS_Element {
  public:
    RS_Triangle (int number, Subdivision::RS_Mesh* mesh, int parent, IntArray& nodes) ;
		int evaluateLongestEdge();
    void bisect(std::queue<int> &subdivqueue, std::list<int> &sharedIrregularsQueue);
    void generate ();
    void update_neighbours();
    double giveDensity();
    bool isNeighborOf (Subdivision::RS_Element* elem);
		void giveSideNodes (int iside, IntArray& snodes);
		int giveEdgeIndex (int iNode, int jNode);
		virtual void importConnectivities(ConnectivityTable* ct);
    /** 
        Add Irregular node and associate it to corresponding element edge, that is shared by given neighborElement
    */
    void addIrregular (int neighbour_base_element, int node);
#ifdef __OOFEG
    void drawGeometry() ;
#endif
#ifdef __PARALLEL_MODE
    // returns number of shared remote partitions (returned in partitions param)
    int giveIrregularSharedPartitions(int eIndex, int iNode, int jNode, IntArray& partitions);
#endif
  };

	class RS_Tetra : public Subdivision::RS_Element {
  protected:
		IntArray side_leIndex;
  public:
    RS_Tetra (int number, Subdivision::RS_Mesh* mesh, int parent, IntArray& nodes) ;
		int evaluateLongestEdge();
    void bisect(std::queue<int> &subdivqueue, std::list<int> &sharedIrregularsQueue);
    void generate ();
    void update_neighbours();
    double giveDensity();
    bool isNeighborOf (Subdivision::RS_Element* elem);
		void giveSideNodes (int iside, IntArray& snodes);
		int giveEdgeIndex (int iNode, int jNode);
    virtual void importConnectivities(ConnectivityTable* ct);
#ifdef __OOFEG
    void drawGeometry() ;
#endif
#ifdef __PARALLEL_MODE
    // returns number of shared remote partitions (returned in partitions param)
    int giveIrregularSharedPartitions(int eIndex, int iNode, int jNode, IntArray& partitions);
#endif
  };

  class RS_Mesh {
    AList<Subdivision::RS_Node> nodes;
    AList<Subdivision::RS_Element> elements;
    Subdivision *subdivision;
    
  public:
		// increase the number of at once allocated fields if the fine (multilevel) subdivision 
		// with many new nodes and elements is too slow
    RS_Mesh (Subdivision* s) : nodes(0,20), elements(0,20) {this->subdivision=s;} // HUHU
    ~RS_Mesh () {nodes.clear(); elements.clear();}


    Subdivision::RS_Node* giveNode(int i) {return nodes.at(i);}
    Subdivision::RS_Element* giveElement(int i) {return elements.at(i);}
    int giveNumberOfNodes() {return nodes.giveSize();}
    int giveNumberOfElements() {return elements.giveSize();}
    void addNode (int num, Subdivision::RS_Node* obj) {nodes.put(num, obj);}
    void addElement (int num, Subdivision::RS_Element* obj) {elements.put(num, obj);}
    Subdivision *giveSubdivision() {return this->subdivision;}
  };
  
  class RS_CompareNodePositions {
    RS_Mesh *m;
  public:
    RS_CompareNodePositions(RS_Mesh *_m) {m=_m;}
    int operator()(int i, int j);
  };

  struct RS_packRemoteElemsStruct {
    Domain* d;
    IntArray* parentElemMap;
  };

  //RS_Mesh *oldMesh, *newMesh;
  RS_Mesh *mesh;
  std::queue<int> subdivqueue;
  std::list<int> sharedIrregularsQueue;
  // smoothing flag
  bool smoothingFlag;

 public:
  /// Constructor
  Subdivision(Domain* d) : MesherInterface(d) {mesh=0; smoothingFlag=true;}
  ~Subdivision() {if (mesh) delete mesh;}
  
  /// Runs the mesh generation, mesh will be written to corresponding domain din file
  virtual returnCode createMesh(TimeStep *tStep, int domainNumber, int domainSerNum, Domain** dNew);
  const char* giveClassName () {return "Subdivision";}
  Domain* giveDomain() {return domain;}

 protected:
  Subdivision::RS_Mesh* giveMesh() {return mesh;}
  void bisectMesh (); 
  void smoothMesh (); 
#ifdef __PARALLEL_MODE
  /** Exchanges the shared irregulars between partitions. Returns true if any shared irregular has been
      exchanged between any partitions.
  */
  bool exchangeSharedIrregulars ();
  int  packSharedIrregulars (Subdivision *s, ProcessCommunicator &pc);
  int  unpackSharedIrregulars (Subdivision *s, ProcessCommunicator &pc);
  void assignGlobalNumbersToSharedIrregulars ();
  int packIrregularSharedGlobnums (Subdivision *s, ProcessCommunicator &pc);
  int unpackIrregularSharedGlobnums (Subdivision *s, ProcessCommunicator &pc);
  

  bool isNodeLocalIrregular (Subdivision::RS_Node* node, int myrank);
  /// Returns true if receiver is irregular, shared node locally maintatined 
  bool isNodeLocalSharedIrregular (Subdivision::RS_Node* node, int myrank);

  int giveRank() ;
  int giveNumberOfProcesses();

  void exchangeRemoteElements (Domain* d, IntArray&);
  int packRemoteElements (RS_packRemoteElemsStruct *s, ProcessCommunicator &pc);
  int unpackRemoteElements (Domain *d, ProcessCommunicator &pc);

  void assignGlobalNumbersToElements(Domain* d);
#endif
};

#endif // subdivision_h
