/* $Header:$ */
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

#include "subdivision.h"
#include "errorestimator.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "engngm.h"
#include "mathfem.h"
#include "masterdof.h"
#include "nonlocalbarrier.h"
#include "initial.h"
#include "usrdefsub.h"
#include "oofemtxtinputrecord.h"
#include "outputmanager.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#include <queue>
#include <set>
#endif

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#endif


#ifdef __PARALLEL_MODE
#include "parallel.h"
#include "communicator.h"
#include "datastream.h"
#include "domaintransactionmanager.h"
#endif

#ifdef __OOFEG
void Subdivision::RS_Node::drawGeometry()
//
// draws graphics representation of receiver
//
{
    GraphicObj *go;

    WCRec p [ 1 ]; /* point */
    p [ 0 ].x = ( FPNum ) this->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) this->giveCoordinate(3);
    
    EASValsSetMType(FILLED_CIRCLE_MARKER);
    //EASValsSetColor( gc.getNodeColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    EASValsSetMSize(8);
    go = CreateMarker3D(p);
    EGWithMaskChangeAttributes(LAYER_MASK | MTYPE_MASK | MSIZE_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);


    char num [ 6 ];
    //EASValsSetColor( gc.getNodeColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) this->giveCoordinate(3);
    sprintf( num, "%d", this->number );
    go = CreateAnnText3D(p, num);
    EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}

#endif


void
Subdivision::RS_Element::addIrregular (int neighborElement, int node)
{
  int ind;
  if ((ind = neghbours_base_elements.findFirstIndexOf(neighborElement))) {
    irregular_nodes.at(ind) = node;
  } else {
    OOFEM_ERROR3("Subdivision::RS_Element::addIrregular: element %d is not a neighbor of element %d", neighborElement, this->number);
  }
}

double
Subdivision::RS_Element::giveRequiredDensity()
{
  int i, nnodes = nodes.giveSize();
  double answer=0.0;
  for (i=1; i<=nnodes; i++) {
    answer+=mesh->giveNode(nodes.at(i))->giveRequiredDensity();
  }
  answer = answer/nnodes;
  return answer;
}


Subdivision::RS_Triangle::RS_Triangle (int number, RS_Mesh* mesh, int parent, IntArray& nodes) : 
  Subdivision::RS_Element (number, mesh, parent, nodes) {

  irregular_nodes.resize(3); irregular_nodes.zero();
  neghbours_base_elements.resize(3); 
  neghbours_base_elements.zero();

  double elength1, elength2, elength3;
  elength1 = mesh->giveNode(this->nodes.at(1))->giveCoordinates()->distance(*(mesh->giveNode(this->nodes.at(2))->giveCoordinates()));
  elength2 = mesh->giveNode(this->nodes.at(2))->giveCoordinates()->distance(*(mesh->giveNode(this->nodes.at(3))->giveCoordinates()));
  elength3 = mesh->giveNode(this->nodes.at(3))->giveCoordinates()->distance(*(mesh->giveNode(this->nodes.at(1))->giveCoordinates()));
  leIndex = 1;
  if (elength2>elength1) leIndex = 2;
  if (elength3>elength1) leIndex = 3;

  
}


void
Subdivision::RS_Triangle::bisect(std::queue<int> &subdivqueue, std::queue<int> &sharedIrregularsQueue) {
  /* this is symbolic bisection - no new elements are added, only irregular nodes are introduced */
  int inode, jnode;
  double density;
  FloatArray coords;

#ifdef __PARALLEL_MODE
  if (this->giveParallelMode() != Element_local) return;
#endif

  // test if element has other irregular nodes than the one from longest edge bissection
  if (irregular_nodes.at(leIndex)) {
    // irregular on the longest edge already exist
    return;
  } else {
    int iNum;
    RS_IrregularNode* irregular;
    // introduce new irregular node on the longest edge
    inode = leIndex;
    jnode = (leIndex<3)?inode+1:1;
    // compute coordinates of new irregular
    coords = *(mesh->giveNode(nodes.at(inode))->giveCoordinates());
    coords.add (mesh->giveNode(nodes.at(jnode))->giveCoordinates());
    coords.times(0.5);
    // compute required density of a new node
    density = 0.5*(mesh->giveNode(nodes.at(inode))->giveRequiredDensity()+
		   mesh->giveNode(nodes.at(jnode))->giveRequiredDensity());
    // add irregular to receiver and its neighbour
    // 
    iNum = mesh->giveNumberOfNodes() + 1;
    irregular = new Subdivision::RS_IrregularNode (iNum, 0, coords, density);
    mesh->addNode (iNum, irregular);
#ifdef __OOFEG
    irregular->drawGeometry();
#endif

    this->irregular_nodes.at(leIndex) = iNum;
    if (this->neghbours_base_elements.at(leIndex)) {
      mesh->giveElement(this->neghbours_base_elements.at(leIndex))->addIrregular(this->number, iNum);
      // add neighbour to list of elements for subdivision
      subdivqueue.push(this->neghbours_base_elements.at(leIndex));
    }

#ifdef __PARALLEL_MODE
    int partition;
    // test if new node is on shared interpartition boundary
    if (this->isIrregularShared (leIndex, inode, jnode, partition)) {
      // mark irregular as shared (can be shared by only one partition (in 2D))
      irregular->setParallelMode (DofManager_shared);
      irregular->setPartition (partition);
      // and put its number into queue of shared irregulars that is later used to inform remote partition about this fact
      irregular->setEdgeNodes(nodes.at(inode), nodes.at(jnode));
      sharedIrregularsQueue.push(iNum);
#ifdef __VERBOSE_PARALLEL
      OOFEM_LOG_INFO("RS_Triangle::bisect: Shared Irregular detected, number %d [%d[%d] %d[%d]], elem %d, partition %d\n", iNum, nodes.at(inode), mesh->giveNode(nodes.at(inode))->giveGlobalNumber(), nodes.at(jnode), mesh->giveNode(nodes.at(jnode))->giveGlobalNumber(), this->number, partition);
#endif
    }
#endif

  }
}

void
Subdivision::RS_Triangle::addIrregularOn (int iNum, const IntArray& bEntity)
{
  if (bEntity.giveSize() != 2) OOFEM_ERROR ("Subdivision::RS_Triangle::addIrregularOn: bEntity size mismatch");
  int anode = bEntity.at(1), bnode = bEntity.at(2), inode, jnode, ie;
  
  for (ie=1; ie<=3; ie++) { // loop over edges
    inode = ie;
    jnode = (ie<3)?ie+1: 1;
    
    if (((this->nodes.at(inode) == anode) && (this->nodes.at(jnode) == bnode)) ||
	((this->nodes.at(inode) == bnode) && (this->nodes.at(jnode) == anode))) {
      this->irregular_nodes.at(ie) = iNum;
      return;
    }
  }
  OOFEM_ERROR ("Subdivision::RS_Triangle::addIrregularOn: no such bEntity found");
}

int
Subdivision::RS_Triangle::giveIrregular (const IntArray& bEntity)
{
  if (bEntity.giveSize() != 2) OOFEM_ERROR ("Subdivision::RS_Triangle::giveIrregular: bEntity size mismatch");
  int anode = bEntity.at(1), bnode = bEntity.at(2), inode, jnode, ie;
  
  for (ie=1; ie<=3; ie++) { // loop over edges
    inode = ie;
    jnode = (ie<3)?ie+1: 1;
    
    if (((this->nodes.at(inode) == anode) && (this->nodes.at(jnode) == bnode)) ||
	((this->nodes.at(inode) == bnode) && (this->nodes.at(jnode) == anode))) {
      return this->irregular_nodes.at(ie);
    }
  }
  return 0;
}


#ifdef __PARALLEL_MODE
bool
Subdivision::RS_Element::isIrregularShared(int leIndex, int inode, int jnode, int& partition)
{
  bool found = false;
  if ((mesh->giveNode(nodes.at(inode))->giveParallelMode() == DofManager_shared) &&
      (mesh->giveNode(nodes.at(jnode))->giveParallelMode() == DofManager_shared)) {
    // test if neigbor element exist (if yes and is local then the edge is not shared
    if (this->neghbours_base_elements.at(leIndex)) {
      if (mesh->giveElement(this->neghbours_base_elements.at(leIndex))->giveParallelMode() != Element_local) {
      
	// test if parent vertex nodes are shared by the same partition
	const IntArray *ipartitions = mesh->giveNode(nodes.at(inode))->givePartitions();
	const IntArray *jpartitions = mesh->giveNode(nodes.at(jnode))->givePartitions();
	// find common partitions 
	int i, ipsize = ipartitions->giveSize();
	for (i=1; i<=ipsize; i++) {
	  if (jpartitions->contains (ipartitions->at(i)) &&
	      (ipartitions->at(i) != this->mesh->giveSubdivision()->giveRank())) { // 
	    partition = ipartitions->at(i);
	    if (!found) found = true; else OOFEM_ERROR ("Subdivision::RS_Triangle::isIrregularShared: topology error");
	    
	  }
	}
      }
    }
  }
  return found;
}

int
Subdivision::RS_Element::containsGlobalNode (int gnum)
{
  int in, nnodes = nodes.giveSize();
  for (in=1; in<=nnodes; in++) {
    if (this->mesh->giveNode(nodes.at(in))->giveGlobalNumber() == gnum) return in;
  }
  return 0;
}
#endif



void
Subdivision::RS_Triangle::generate()
{
  /* generates the children elements of already bisected element and adds them into mesh.
     Also children array of receiver is updated to contain children numbers */
#ifdef __PARALLEL_MODE
  if ((irregular_nodes.containsOnlyZeroes()) || (this->giveParallelMode() == Element_remote)) {
#else
  if (irregular_nodes.containsOnlyZeroes()) {
#endif
    // no subdivision of receiver required 
    children.resize(0);
  } else {
    int i, nIrregulars = 0;
    int childNum, iedge, jedge, kedge, inode, jnode, knode;
    IntArray _nodes(3);
    Subdivision::RS_Triangle *child;
    
    for (i=1; i<=irregular_nodes.giveSize(); i++) if (irregular_nodes.at(i)) nIrregulars++;
    this->children.resize(nIrregulars+1);
    
    // leIndex determines primary edge to be subdivided
    /*

                               inode
                                 o
                                

                   kedge x                x jedge


                o                x                o
              jnode            iedge               knode
                              =leIndex

    */
    iedge = leIndex;
    jedge = (iedge <3)?iedge+1:1;
    kedge = (jedge <3)?jedge+1:1;

    inode = (iedge>1)?iedge-1:3;
    jnode = (inode<3)?inode+1:1;
    knode = (jnode<3)?jnode+1:1;
    
    if (irregular_nodes.at(iedge) && irregular_nodes.at(jedge) && irregular_nodes.at(kedge)) {
      _nodes.at(1)=nodes.at(inode); _nodes.at(2)=irregular_nodes.at(kedge); _nodes.at(3)=irregular_nodes.at(iedge);
      childNum = mesh->giveNumberOfElements() + 1;
      child = new Subdivision::RS_Triangle (childNum, mesh, this->number, _nodes); // number, parent, coords
      mesh->addElement (childNum, child);
      children.at(1) = childNum;
      // set neigbour info
      child->setNeighbor(1,-this->giveNeighbor(kedge));
      child->setNeighbor(2, childNum+2);
      child->setNeighbor(3, childNum+1);

      _nodes.at(1)=nodes.at(inode); _nodes.at(2)=irregular_nodes.at(iedge); _nodes.at(3)=irregular_nodes.at(jedge);
      childNum = mesh->giveNumberOfElements() + 1;
      child = new Subdivision::RS_Triangle (childNum, mesh, this->number, _nodes); // number, parent, coords
      mesh->addElement (childNum, child);
      children.at(2) = childNum;
      // set neigbour info
      child->setNeighbor(1, childNum-1);
      child->setNeighbor(2, childNum+2);
      child->setNeighbor(3,-this->giveNeighbor(jedge));

      _nodes.at(1)=nodes.at(jnode); _nodes.at(2)=irregular_nodes.at(iedge); _nodes.at(3)=irregular_nodes.at(kedge);
      childNum = mesh->giveNumberOfElements() + 1;
      child = new Subdivision::RS_Triangle (childNum, mesh, this->number, _nodes); // number, parent, coords
      mesh->addElement (childNum, child);
      children.at(3) = childNum;
      // set neigbour info
      child->setNeighbor(1,-this->giveNeighbor(iedge));
      child->setNeighbor(2, childNum-2);
      child->setNeighbor(3,-this->giveNeighbor(kedge));

      _nodes.at(1)=nodes.at(knode); _nodes.at(2)=irregular_nodes.at(jedge); _nodes.at(3)=irregular_nodes.at(iedge);
      childNum = mesh->giveNumberOfElements() + 1;
      child = new Subdivision::RS_Triangle (childNum, mesh, this->number, _nodes); // number, parent, coords
      mesh->addElement (childNum, child);
      children.at(4) = childNum;
      // set neigbour info
      child->setNeighbor(1,-this->giveNeighbor(jedge)); 
      child->setNeighbor(2, childNum-2);
      child->setNeighbor(3,-this->giveNeighbor(iedge));      

    } else if (irregular_nodes.at(iedge) && irregular_nodes.at(jedge)) {
      _nodes.at(1)=nodes.at(inode); _nodes.at(2)=nodes.at(jnode); _nodes.at(3)=irregular_nodes.at(iedge);
      childNum = mesh->giveNumberOfElements() + 1;
      child = new Subdivision::RS_Triangle (childNum, mesh, this->number, _nodes); // number, parent, coords
      mesh->addElement (childNum, child);
      children.at(1) = childNum;
      // set neigbour info
      child->setNeighbor(1,-this->giveNeighbor(kedge)); 
      child->setNeighbor(2,-this->giveNeighbor(iedge));
      child->setNeighbor(3, childNum+1);

      _nodes.at(1)=nodes.at(inode); _nodes.at(2)=irregular_nodes.at(iedge); _nodes.at(3)=irregular_nodes.at(jedge);
      childNum = mesh->giveNumberOfElements() + 1;
      child = new Subdivision::RS_Triangle (childNum, mesh, this->number, _nodes); // number, parent, coords
      mesh->addElement (childNum, child);
      children.at(2) = childNum;
      // set neigbour info
      child->setNeighbor(1, childNum-1);
      child->setNeighbor(2, childNum+1);
      child->setNeighbor(3,-this->giveNeighbor(jedge));      

      _nodes.at(1)=nodes.at(knode); _nodes.at(2)=irregular_nodes.at(jedge); _nodes.at(3)=irregular_nodes.at(iedge);
      childNum = mesh->giveNumberOfElements() + 1;
      child = new Subdivision::RS_Triangle (childNum, mesh, this->number, _nodes); // number, parent, coords
      mesh->addElement (childNum, child);
      children.at(3) = childNum;
      // set neigbour info
      child->setNeighbor(1,-this->giveNeighbor(jedge)); 
      child->setNeighbor(2, childNum-1);
      child->setNeighbor(3,-this->giveNeighbor(iedge));      

    } else if (irregular_nodes.at(iedge) && irregular_nodes.at(kedge)) {

      _nodes.at(1)=nodes.at(inode); _nodes.at(2)=irregular_nodes.at(kedge); _nodes.at(3)=irregular_nodes.at(iedge);
      childNum = mesh->giveNumberOfElements() + 1;
      child = new Subdivision::RS_Triangle (childNum, mesh, this->number, _nodes); // number, parent, coords
      mesh->addElement (childNum, child);
      children.at(1) = childNum;
      // set neigbour info
      child->setNeighbor(1,-this->giveNeighbor(kedge)); 
      child->setNeighbor(2, childNum+1);
      child->setNeighbor(3, childNum+2);

      _nodes.at(1)=nodes.at(jnode); _nodes.at(2)=irregular_nodes.at(iedge); _nodes.at(3)=irregular_nodes.at(kedge);
      childNum = mesh->giveNumberOfElements() + 1;
      child = new Subdivision::RS_Triangle (childNum, mesh, this->number, _nodes); // number, parent, coords
      mesh->addElement (childNum, child);
      children.at(2) = childNum;
      // set neigbour info
      child->setNeighbor(1,-this->giveNeighbor(iedge)); 
      child->setNeighbor(2, childNum-1);
      child->setNeighbor(3,-this->giveNeighbor(kedge));      

      _nodes.at(1)=nodes.at(inode); _nodes.at(2)=irregular_nodes.at(iedge); _nodes.at(3)=nodes.at(knode);
      childNum = mesh->giveNumberOfElements() + 1;
      child = new Subdivision::RS_Triangle (childNum, mesh, this->number, _nodes); // number, parent, coords
      mesh->addElement (childNum, child);
      children.at(3) = childNum;
      // set neigbour info
      child->setNeighbor(1, childNum-2);
      child->setNeighbor(2,-this->giveNeighbor(iedge));
      child->setNeighbor(3,-this->giveNeighbor(jedge));      

    } else if (irregular_nodes.at(iedge)) {

      _nodes.at(1)=nodes.at(inode); _nodes.at(2)=nodes.at(jnode); _nodes.at(3)=irregular_nodes.at(iedge);
      childNum = mesh->giveNumberOfElements() + 1;
      child = new Subdivision::RS_Triangle (childNum, mesh, this->number, _nodes); // number, parent, coords
      mesh->addElement (childNum, child);
      children.at(1) = childNum;
      // set neigbour info
      child->setNeighbor(1,-this->giveNeighbor(kedge));
      child->setNeighbor(2,-this->giveNeighbor(iedge));
      child->setNeighbor(3, childNum+1);
      

      _nodes.at(1)=nodes.at(inode); _nodes.at(2)=irregular_nodes.at(iedge); _nodes.at(3)=nodes.at(knode);
      childNum = mesh->giveNumberOfElements() + 1;
      child = new Subdivision::RS_Triangle (childNum, mesh, this->number, _nodes); // number, parent, coords
      mesh->addElement (childNum, child);
      children.at(2) = childNum;
      // set neigbour info
      child->setNeighbor(1, childNum-1);
      child->setNeighbor(2,-this->giveNeighbor(iedge));
      child->setNeighbor(3,-this->giveNeighbor(jedge));
     
    } else {
      OOFEM_ERROR2("Subdivision::RS_Triangle::generate - element %d internal data inconsistency", this->number);
    }
    
  }
}





void
Subdivision::RS_Triangle::update_neighbours()
{
  bool found;
  int i, iside,inode,jnode,parentNeighbor;
  IntArray ca;

  /*
    neghbours_base_elements pre-set in generate using this rule:
    value > 0 .. real neighbour known, value equals to real neighbor on new mesh
    value = 0 .. boundary
    value < 0 .. neighbour set to parent element, actual neighbour has to be determined from parent children 
                 or if no child exist then to parent new counterpart.

   */


  //this->neghbours_base_elements.resize(3);
  for (iside=1; iside<=3; iside++) {
    if (this->neghbours_base_elements.at(iside) < 0) {
      inode = iside;
      jnode = (iside<3)?iside+1: 1;
      parentNeighbor = -this->neghbours_base_elements.at(iside);
      // test if parentNeighbor has been subdivided
#ifdef __PARALLEL_MODE
      if ((mesh->giveElement(parentNeighbor)->giveParallelMode() == Element_local) &&
	  (mesh->giveElement(parentNeighbor)->hasIrregulars())) {
#else
      if (mesh->giveElement(parentNeighbor)->hasIrregulars()) {
#endif
	// old neighbour bisected -> we loop over its children to find an appropriate new neighbor
	mesh->giveElement(parentNeighbor)->giveChildren(ca);
	found = false;
	for (i=1; i<=ca.giveSize(); i++) {
	  if (this->isNeighborOf(mesh->giveElement(ca.at(i)))) {
	    this->neghbours_base_elements.at(iside) = ca.at(i);
	    found = true;
	    break;
	  }
	}
	if (!found) {
	  OOFEM_ERROR2 ("Subdivision::RS_Triangle::update_neighbours failed for element %d", this->number);
	}
      } else {
	  // parent neighbour remains actual neighbour
	  this->neghbours_base_elements.at(iside) = parentNeighbor;
      }

    } else if (this->neghbours_base_elements.at(iside) > 0) {
      // neighbor element already set
    }
  } // end loop over element sides
}


bool
Subdivision::RS_Triangle::isNeighborOf (Subdivision::RS_Element* elem)
{
  // Simplified implementation, considering only one element type - triangles
  int i, _c = 0;
  for (i=1; i<=3; i++) {
    _c +=  elem->containsNode(this->nodes.at(i));
  }
  return (_c==2);
}

double
Subdivision::RS_Triangle::giveDensity() 
{
  double x1,x2,x3,y1,y2,y3;
  x1 = mesh->giveNode(nodes.at(1))->giveCoordinate(1);
  x2 = mesh->giveNode(nodes.at(2))->giveCoordinate(1);
  x3 = mesh->giveNode(nodes.at(3))->giveCoordinate(1);
  
  y1 = mesh->giveNode(nodes.at(1))->giveCoordinate(2);
  y2 = mesh->giveNode(nodes.at(2))->giveCoordinate(2);
  y3 = mesh->giveNode(nodes.at(3))->giveCoordinate(2);
  
  return sqrt (fabs( 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 )));
}


void
Subdivision::RS_Triangle::importConnectivities(ConnectivityTable* ct)
{
  IntArray me(1), conn;
  int iside, inode, jnode, el, i;
  bool _found;
  me.at(1) = this->number;

  neghbours_base_elements.resize(3);
  neghbours_base_elements.zero();
  ct->giveElementNeighbourList(conn, me);

  for (iside=1; iside<4; iside++) {
    inode=nodes.at(iside);
    jnode=nodes.at((iside==3)?1:iside+1);
    _found = false;

    // select right neighbour
    for (i=1; i<=conn.giveSize(); i++) {
      el = conn.at(i);
      if (el != this->number) {
	if (mesh->giveElement(el)->containsNode(inode) &&
	    mesh->giveElement(el)->containsNode(jnode)) {
	  _found = true;
	  neghbours_base_elements.at(iside) = el;
	}
      }
    }
    if (!_found) neghbours_base_elements.at(iside) = 0;// no neighbour (boundary edge)
  }
}

#ifdef __OOFEG
void
Subdivision::RS_Triangle::drawGeometry()
{
    WCRec p [ 3 ];
    GraphicObj *go;

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    //EASValsSetColor( gc.getElementColor() );
    //EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(TRUE);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    //    EASValsSetShrink(0.8);
    p [ 0 ].x = ( FPNum ) mesh->giveNode(nodes.at(1))->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) mesh->giveNode(nodes.at(1))->giveCoordinate(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) mesh->giveNode(nodes.at(2))->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) mesh->giveNode(nodes.at(2))->giveCoordinate(2);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) mesh->giveNode(nodes.at(3))->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) mesh->giveNode(nodes.at(3))->giveCoordinate(2);
    p [ 2 ].z = 0.;

    go =  CreateTriangle3D(p);
    //EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | SHRINK_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGWithMaskChangeAttributes(WIDTH_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}
#endif


void
Subdivision::bisectMesh () {
  
  int ie, nelems = mesh->giveNumberOfElements();
  int in, nnodes=mesh->giveNumberOfNodes();
  double iedensity, rdensity;
  //std::queue<int> subdivqueue;
#ifdef __PARALLEL_MODE
  int myrank=this->giveRank();
  int problem_size=this->giveNumberOfProcesses();
#endif
  // first select all candidates for local bisection based on required mesh density

  for (ie=1; ie<=nelems; ie++) {
#ifdef __PARALLEL_MODE
    if (mesh->giveElement(ie)->giveParallelMode() == Element_remote) continue;
#endif
    if (!mesh->giveElement(ie)->isTerminal()) continue;
    iedensity = mesh->giveElement(ie)->giveDensity();
    rdensity = mesh->giveElement(ie)->giveRequiredDensity();
    
    if (rdensity < iedensity) {
      subdivqueue.push(ie);
    }
  }
#ifdef __PARALLEL_MODE
  do {
#endif
    // loop over subdivision queue to bisect all local elements there
    while( !subdivqueue.empty() ) {
      mesh->giveElement(subdivqueue.front())->bisect(subdivqueue, sharedIrregularsQueue);
      subdivqueue.pop();
    }   
    
#ifdef __PARALLEL_MODE
   // in parallel communicate with neighbours the irregular nodes on shared bondary
  } while (!exchangeSharedIrregulars ());

  // assign global numbers to newly introduced irregulars while
  // keeping global numbering of existing (master) nodes
  
  // idea: first determine the max globnum already assigned
  // and start global numbering of new nodes from this value up
  // 1. determine max global number of local nodes
  int maxlocalglobal = 0, maxglobalnumber;
  int localIrregulars = 0;
  nnodes=mesh->giveNumberOfNodes();
  for (in=1; in<=nnodes; in++) {
    maxlocalglobal = max (maxlocalglobal, mesh->giveNode(in)->giveGlobalNumber());
  }
  // determine max global number on all partitions
  MPI_Allreduce (&maxlocalglobal, &maxglobalnumber, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  // determine number of local irregulars
  // this is needed to determine offsets on each partition to start global numbering
  // the numbering of irregulars starts on rank 0 partition by numbering subsequntly all its irregulars
  // then we proceed with rank 1, etc. 
  // The shared irregulars receive their number from partition with the lowest rank.
  
  // ok. count local irregulars that receive their global number from this partition
  for (in=1; in<=nnodes; in++) {
    if (this->isNodeLocalIrregular(mesh->giveNode(in), myrank))
      localIrregulars ++;
  }

#ifdef __VERBOSE_PARALLEL
  OOFEM_LOG_INFO ("[%d] Subdivision::bisectMesh: number of new local irregulars is %d\n", myrank, localIrregulars);
#endif
  int irank, localOffset, partitionsIrregulars[problem_size];
  // gather number of local irregulars from all partitions
  MPI_Allgather (&localIrregulars, 1, MPI_INT, partitionsIrregulars, 1, MPI_INT, MPI_COMM_WORLD);
  // compute local offset
  for (localOffset=0, irank=0; irank <myrank; irank++) localOffset += partitionsIrregulars[irank];
  // start to assign global numbers to local irregulars
  int availGlobNum = maxglobalnumber+localOffset;
  for (in=1; in<=nnodes; in++) {
    if (this->isNodeLocalIrregular(mesh->giveNode(in), myrank))
      mesh->giveNode(in)->setGlobalNumber(++availGlobNum);
  }
  // finally, communicate global numbers assigned to shared irregulars
  this->assignGlobalNumbersToSharedIrregulars();
#endif

   // symbolic bisection is finished. Now we need to create new mesh. Also there is a need to update 
   // element connectivities
   for (ie=1; ie<=nelems; ie++) { mesh->giveElement(ie)->generate ();  }

   int nelems_new = mesh->giveNumberOfElements();
   for (ie=1; ie<=nelems_new; ie++) { mesh->giveElement(ie)->update_neighbours();}

}


MesherInterface::returnCode
Subdivision :: createMesh(TimeStep *stepN, int domainNumber, int domainSerNum, Domain** dNew)
{
  // import data from old domain
  int i, j, parent, nnodes = domain->giveNumberOfDofManagers(), nelems=domain->giveNumberOfElements();
  int inode, idof, ielem, ndofs, num;
  IntArray enodes;
  Subdivision::RS_Node *_node;
  Subdivision::RS_Triangle *_element;
  IRResultType result;                            // Required by IR_GIVE_FIELD macro

  oofem_timeval st, dt;
  ::getUtime(st);

  if (this->mesh) delete mesh;
  mesh = new Subdivision::RS_Mesh(this);

  // import nodes
  for (i=1; i<=nnodes; i++) {
    _node = new Subdivision::RS_Node(i, i, *(domain->giveNode(i)->giveCoordinates()), 
				     domain->giveErrorEstimator()->giveRemeshingCrit()->giveRequiredDofManDensity(i, stepN));
#ifdef __PARALLEL_MODE
    _node->setGlobalNumber(domain->giveNode(i)->giveGlobalNumber());
    _node->setParallelMode(domain->giveNode(i)->giveParallelMode());
    _node->setPartitions (*domain->giveNode(i)->givePartitionList());
#endif
    this->mesh->addNode (i, _node);
  }
  
  // import elements
  for (i=1; i<=nelems; i++) {
    if (domain->giveElement(i)->giveGeometryType() == EGT_triangle_1) {
      enodes.resize(3);
      for (j=1; j<=3; j++) enodes.at(j)=domain->giveElement(i)->giveDofManagerNumber(j);
      _element = new Subdivision::RS_Triangle (i,mesh,i,enodes);
      this->mesh->addElement (i, _element);

    } else {
      OOFEM_ERROR2 ("Subdivision::createMesh: Unsupported element geometry (element %d)", i);
    }
#ifdef __PARALLEL_MODE
    _element->setGlobalNumber(domain->giveElement(i)->giveGlobalNumber());
    _element->setParallelMode(domain->giveElement(i)->giveParallelMode());
#endif
  }
  // import connectivities
  for (i=1; i<=nelems; i++) {
    this->mesh->giveElement(i)->importConnectivities(domain->giveConnectivityTable());
  }
  
#ifdef __OOFEG
  //nelems=oldMesh->giveNumberOfElements();
  //for (i=1; i<=nelems; i++) {
  //  oldMesh->giveElement(i)->drawGeometry();
  // }
  //ESIEventLoop (YES, "Subdivision Bisection Finished; Press Ctrl-p to continue"); 
#endif

  // bisect mesh
  this->bisectMesh();
#ifdef __OOFEG
  /*
  nelems=newMesh->giveNumberOfElements();
  for (i=1; i<=nelems; i++) {
    newMesh->giveElement(i)->drawGeometry();
    //ESIEventLoop (YES, "Subdivision Bisection Finished; Press Ctrl-p to continue");
  }
  ESIEventLoop (YES, "Subdivision Bisection Finished; Press Ctrl-p to continue");
  */
#endif


  Dof *idofPtr, *dof;
  DofManager *parentNodePtr, *node;
  Element *parentElementPtr, *elem;
  CrossSection* crossSection;
  Material* mat;
  NonlocalBarrier* barrier;
  GeneralBoundaryCondition *bc;
  InitialCondition *ic;
  LoadTimeFunction *ltf;
  char name [ MAX_NAME_LENGTH ];
  const char *__proc = "createMesh"; // Required by IR_GIVE_FIELD macro
  
  // create new mesh (missing param for new mesh!)
  nnodes = mesh->giveNumberOfNodes();
  (*dNew) = new Domain (2, domain->giveSerialNumber()+1, domain->giveEngngModel());
  (*dNew)->setDomainType (domain->giveDomainType());
  
  // copy dof managers
  (*dNew) -> resizeDofManagers (nnodes);
  const IntArray dofIDArrayPtr = domain->giveDefaultNodeDofIDArry();
  for (inode=1; inode<=nnodes; inode++) {
    parent = mesh->giveNode(inode)->giveParent();
    if (parent) {
      parentNodePtr = domain->giveNode(parent);
      // inherit all data from parent (bc, ic, load, etc.)
      node = ::CreateUsrDefDofManagerOfType (parentNodePtr->giveClassID(), inode, *dNew);
      ndofs = parentNodePtr->giveNumberOfDofs();
      node->setNumberOfDofs (ndofs);
      node->setLoadArray (*parentNodePtr->giveLoadArray());
      if (parentNodePtr->giveClassID() == NodeClass) {
	((Node*) node)->setCoordinates(*parentNodePtr->giveCoordinates());
      }
      // create individual DOFs
      for (idof=1; idof<=ndofs; idof++) {
	idofPtr = parentNodePtr->giveDof(idof);
	if (idofPtr->giveClassID() == MasterDofClass) {
	  dof = new MasterDof (idof, node, idofPtr->giveBcId(), idofPtr->giveIcId(), idofPtr->giveDofID());
	} else OOFEM_ERROR ("Subdivision :: createMesh: unsupported DOF type");
	node->setDof (idof, dof);
      }
#ifdef __PARALLEL_MODE
      node->setGlobalNumber(parentNodePtr->giveGlobalNumber());
      node->setParallelMode(parentNodePtr->giveParallelMode());
      node->setPartitionList (parentNodePtr->givePartitionList());
#endif      
    } else {
      // newly created "bisected" node
      node = ::CreateUsrDefDofManagerOfType (NodeClass, inode, *dNew);
      //create new node with default DOFs 

      ndofs = dofIDArrayPtr.giveSize();
      node->setNumberOfDofs (ndofs);
      ((Node*) node)->setCoordinates(*mesh->giveNode(inode)->giveCoordinates());
      // create individual DOFs
      for (idof=1; idof<=ndofs; idof++) {
	dof = new MasterDof (idof, node, 0, 0, dofIDArrayPtr.at(idof));
	node->setDof (idof, dof);
      }
#ifdef __PARALLEL_MODE
      // ?????   node->setGlobalNumber(parentPtr->giveGlobalNumber());
      node->setParallelMode(mesh->giveNode(inode)->giveParallelMode());
      node->setGlobalNumber(mesh->giveNode(inode)->giveGlobalNumber());
      node->setPartitionList (mesh->giveNode(inode)->givePartitions());
#endif      
      
    }
    (*dNew) -> setDofManager (inode, node);
  } // end creating dof managers

  // create elements
  // count number of local terminal elements first
  int count, nterminals = 0;
  nelems=mesh->giveNumberOfElements();
  for (ielem=1; ielem<=nelems; ielem++) {
#ifdef __PARALLEL_MODE
    if (mesh->giveElement(ielem)->giveParallelMode()==Element_remote) continue;
#endif
    if (mesh->giveElement(ielem)->isTerminal()) nterminals++;
  }

  (*dNew) -> resizeElements (nterminals);
  int eNum = 0;
  for (ielem=1; ielem<=nelems; ielem++) {
#ifdef __PARALLEL_MODE
    if (mesh->giveElement(ielem)->giveParallelMode()==Element_remote) continue;
#endif
    if (!mesh->giveElement(ielem)->isTerminal()) continue;
    eNum ++;
    parent = mesh->giveElement(ielem)->giveParent();
    if (parent) {
      parentElementPtr = domain->giveElement(parent);
      elem = ::CreateUsrDefElementOfType (parentElementPtr->giveClassID(), eNum, *dNew);
      (*dNew)->setElement (eNum, elem);
      elem->setDofManagers(*mesh->giveElement(ielem)->giveNodes());
      elem->setMaterial (parentElementPtr->giveMaterial()->giveNumber());
      elem->setCrossSection(parentElementPtr->giveCrossSection()->giveNumber());
#ifdef __PARALLEL_MODE
      elem->setParallelMode(parentElementPtr->giveParallelMode());
      elem->setGlobalNumber(mesh->giveElement(ielem)->giveGlobalNumber());
#endif
      elem->postInitialize();
    } else OOFEM_ERROR ("Subdivision :: createMesh: parent element missing");
  } // end loop over elements
  OOFEMTXTInputRecord _ir, *irPtr = &_ir;
  std::string irString;
  // create the rest of the model description (BCs, CrossSections, Materials, etc)
  // Cross sections
  int ncrosssect = domain->giveNumberOfCrossSectionModels();
  (*dNew) -> resizeCrossSectionModels(ncrosssect);
  for (i=1; i<=ncrosssect; i++) {
    domain->giveCrossSection(i)->giveInputRecordString (irString);
    irPtr->setRecordString (irString.c_str());
    IR_GIVE_RECORD_KEYWORD_FIELD(irPtr, name, num, MAX_NAME_LENGTH);
    
    ( crossSection  = ( CrossSection * )
      ( CrossSection(i, *dNew).ofType(name) ) )->initializeFrom(irPtr);
    (*dNew) -> setCrossSection (i, crossSection);
  }
  // materials
  int nmat = domain->giveNumberOfMaterialModels();
  (*dNew) -> resizeMaterials(nmat);
  for (i=1; i<=nmat; i++) {
    domain->giveMaterial(i)->giveInputRecordString (irString);
    irPtr->setRecordString (irString.c_str());
    IR_GIVE_RECORD_KEYWORD_FIELD(irPtr, name, num, MAX_NAME_LENGTH);
    
    ( mat  = ( Material * )
      ( Material(i, *dNew).ofType(name) ) )->initializeFrom(irPtr);
    (*dNew) -> setMaterial (i, mat);
  }
  // barriers
  int nbarriers = domain->giveNumberOfNonlocalBarriers();
  (*dNew) -> resizeNonlocalBarriers(nbarriers);
  for (i=1; i<=nbarriers; i++) {
    domain->giveNonlocalBarrier(i)->giveInputRecordString (irString);
    irPtr->setRecordString (irString.c_str());
    IR_GIVE_RECORD_KEYWORD_FIELD(irPtr, name, num, MAX_NAME_LENGTH);
    
    barrier = CreateUsrDefNonlocalBarrierOfType(name, i, *dNew);
    barrier->initializeFrom(irPtr);
    (*dNew) -> setNonlocalBarrier (i, barrier);
  }
  // boundary conditions
  int nbc = domain->giveNumberOfBoundaryConditions();
  (*dNew) -> resizeBoundaryConditions(nbc);
  for (i=1; i<=nbc; i++) {
    domain->giveBc(i)->giveInputRecordString (irString);
    irPtr->setRecordString (irString.c_str());
    IR_GIVE_RECORD_KEYWORD_FIELD(irPtr, name, num, MAX_NAME_LENGTH);
    
    ( bc = ( GeneralBoundaryCondition * )
      ( GeneralBoundaryCondition(i, *dNew).ofType(name) ) )->initializeFrom(irPtr);
    (*dNew) -> setBoundaryCondition (i, bc);
  }
  
  // initial conditions
  int nic = domain->giveNumberOfInitialConditions();
  (*dNew) -> resizeInitialConditions(nic);
  for (i=1; i<=nic; i++) {
    domain->giveIc(i)->giveInputRecordString (irString);
    irPtr->setRecordString (irString.c_str());
    IR_GIVE_RECORD_KEYWORD_FIELD(irPtr, name, num, MAX_NAME_LENGTH);
    
    ic = new InitialCondition(i, *dNew);
    ic->initializeFrom(irPtr);
    (*dNew) -> setInitialCondition (i, ic);
  }
  // load time functions
  int nltf = domain->giveNumberOfLoadTimeFunctions();
  (*dNew) -> resizeLoadTimeFunctions(nltf);
  for (i=1; i<=nltf; i++) {
    domain->giveLoadTimeFunction(i)->giveInputRecordString (irString);
    irPtr->setRecordString (irString.c_str());
    IR_GIVE_RECORD_KEYWORD_FIELD(irPtr, name, num, MAX_NAME_LENGTH);
    
    ( ltf  = ( LoadTimeFunction * )
      ( LoadTimeFunction(i, *dNew).ofType(name) ) )->initializeFrom(irPtr);
    (*dNew) -> setLoadTimeFunction (i, ltf);
  }
  
  // copy output manager settings
  (*dNew)->giveOutputManager()->beCopyOf (domain->giveOutputManager());
  :: getRelativeUtime(dt, st);
  
  OOFEM_LOG_INFO( "Subdivision: created new mesh (%d nodes and %d elements) in %.2fs\n",
		  nnodes, eNum, (double)(dt.tv_sec+dt.tv_usec/(double)OOFEM_USEC_LIM));
  
#ifdef __PARALLEL_MODE
#ifdef __VERBOSE_PARALLEL
  nnodes = (*dNew)->giveNumberOfDofManagers();
  for (inode=1; inode<=nnodes; inode++) {
    if ((*dNew)->giveDofManager(inode)->giveParallelMode()==DofManager_shared) {
      OOFEM_LOG_INFO ("[%d] Shared Node %d[%d]\n", this->giveRank(), inode, (*dNew)->giveDofManager(inode)->giveGlobalNumber());
    }
  }
#endif
  
  


    // we need to assign global numbers to newly generated elements 
    this->assignGlobalNumbersToElements(*dNew);

    int im;
    bool nonloc = false;
    nmat = (*dNew)->giveNumberOfMaterialModels();
    for (im=1; im<=nmat; im++)
      if ((*dNew)->giveMaterial(im)->giveInterface(NonlocalMaterialExtensionInterfaceType)) nonloc=true;
    if (nonloc) exchangeRemoteElements (*dNew);
    
    
    (*dNew)->commitTransactions( (*dNew)->giveTransactionManager() );
    
    // print some statistics
    nelems=(*dNew)->giveNumberOfElements();
    int localVals[2], globalVals[2];
    for (localVals[0] = 0, ielem=1; ielem<=nelems; ielem++) {
      if ((*dNew)->giveElement(ielem)->giveParallelMode() == Element_local) localVals[0]++ ;
    }
    nnodes = (*dNew)->giveNumberOfDofManagers();
    for (localVals[1]=0, inode=1; inode<=nnodes; inode++) {
      if ((*dNew)->giveDofManager(inode)->isLocal()) localVals[1]++;
    }
    MPI_Reduce(localVals, globalVals, 2, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (this->giveRank()==0) {
      OOFEM_LOG_INFO("Subdivision: new mesh info: %d nodes, %d elements in total\n", globalVals[1], globalVals[0]);
    }
    
    
#endif
    
    return MI_OK;
}


#ifdef __PARALLEL_MODE
bool
Subdivision :: exchangeSharedIrregulars () {
  // loop over local sharedIrregularQueue 
  int globalIrregularQueueEmpty, localSharedIrregularQueueEmpty = this->sharedIrregularsQueue.empty();
#ifdef __VERBOSE_PARALLEL
  OOFEM_LOG_INFO ("[%d]Subdivision :: exchangeSharedIrregulars: localSharedIrregularQueueEmpty %d",
		  this->giveRank(), localSharedIrregularQueueEmpty);
#endif
  MPI_Allreduce (&localSharedIrregularQueueEmpty, &globalIrregularQueueEmpty, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD); 
#ifdef __VERBOSE_PARALLEL
  OOFEM_LOG_INFO ("[%d]Subdivision :: exchangeSharedIrregulars:  globalIrregularQueueEmpty %d",
		  this->giveRank(), globalIrregularQueueEmpty);
#endif
  if (globalIrregularQueueEmpty) {
    return true; 
  } else {
#ifdef __VERBOSE_PARALLEL
    OOFEM_LOG_INFO ("[%d]Subdivision :: exchangeSharedIrregulars: started\n", this->giveRank());
#endif
    
    // there are some shared irregulars -> data exchange
    CommunicatorBuff cb(this->giveNumberOfProcesses(), CBT_dynamic);
    Communicator com(domain->giveEngngModel(), &cb, this->giveRank(), this->giveNumberOfProcesses(), CommMode_Dynamic);
    
    com.packAllData(this, this, &Subdivision::packSharedIrregulars);
    com.initExchange(SHARED_IRREGULAR_DATA_TAG);
    com.unpackAllData(this, this, &Subdivision::unpackSharedIrregulars);
    
    return false;
  }
}


int
Subdivision::packSharedIrregulars (Subdivision *s, ProcessCommunicator &pc)
{
  int pi, inode, jnode, rproc = pc.giveRank();
  int myrank = this->giveRank();
  const IntArray *sharedPartitions;
  IntArray edgeInfo (2);
  if ( rproc == myrank ) return 1; // skip local partition

  // query process communicator to use
  ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();
  while( !sharedIrregularsQueue.empty() ) {
    pi = sharedIrregularsQueue.front(); sharedIrregularsQueue.pop();
    sharedPartitions = this->mesh->giveNode(pi)->givePartitions();
    if (sharedPartitions->contains(rproc)) {
      // the info about new local shared irregular node needs to be sent to remote partition
      // the new ireegular on remote partition is identified using two nodes defining 
      // an edge on which irregular node is introduced
      ((RS_IrregularNode*)this->mesh->giveNode(pi))->giveEdgeNodes(inode, jnode);
      edgeInfo.at(1) = this->mesh->giveNode(inode)->giveGlobalNumber();
      edgeInfo.at(2) = this->mesh->giveNode(jnode)->giveGlobalNumber();
      pcbuff->packInt (SUBIVISION_IRREGULAR_REC_TAG);
      pcbuff->packIntArray (edgeInfo);
#ifdef __VERBOSE_PARALLEL
      OOFEM_LOG_INFO("[%d] Subdivision::packSharedIrregulars: packing shared node %d [%d %d]\n",
		     myrank,pi, edgeInfo.at(1), edgeInfo.at(2));
#endif
    }
  }
  pcbuff->packInt (SUBDIVISION_END_DATA);
}

int 
Subdivision::unpackSharedIrregulars (Subdivision *s, ProcessCommunicator &pc)
{
  int myrank = this->giveRank();
  int ie, _type, inode, jnode, iNum, iproc = pc.giveRank();
  int nelems=mesh->giveNumberOfElements();
  double density;
  IntArray edgeInfo(2), bEnt(2);
  FloatArray coords;
  RS_Element *elem;
  RS_IrregularNode *irregular;
  bool found = false;

  if ( iproc == myrank ) {
    return 1;                // skip local partition
  }
  
  // query process communicator to use
  ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();
    
  pcbuff->unpackInt(_type);
  // unpack dofman data
  while ( _type != SUBDIVISION_END_DATA ) {
    if (_type == SUBIVISION_IRREGULAR_REC_TAG) {
      pcbuff->unpackIntArray (edgeInfo);
#ifdef __VERBOSE_PARALLEL
      OOFEM_LOG_INFO("[%d] Subdivision::unpackSharedIrregulars: received shared node record [%d %d]...",
		     myrank, edgeInfo.at(1), edgeInfo.at(2));
#endif
      // look for local element containing an edge identical to the one just received
      for (found=false, ie=1; ie<=nelems; ie++) { 
	elem = mesh->giveElement(ie);
	if (elem->giveParallelMode() == Element_remote) continue;
	if (!elem->isTerminal()) continue;
	if ((inode=elem->containsGlobalNode(edgeInfo.at(1))) &&
	    (jnode=elem->containsGlobalNode(edgeInfo.at(2)))) {

	  /* first check if Irregular is already there
	     this can happen, when both partitions introduce locally shared node on the same edge
	     and this infor is broadcasted mutually between them
	  */
	  bEnt.at(1) = elem->giveNode(inode); bEnt.at(2) = elem->giveNode(jnode);
	  if (elem->giveIrregular (bEnt)) {
	    // entry already exists 
#ifdef __VERBOSE_PARALLEL
	    OOFEM_LOG_INFO ("already exists as %d on %d elem\n", elem->giveIrregular (bEnt), ie);
#endif
	    found = true;
	    break;
	  } else { 
	    // irregular does not exist:
	    // element found => introduce new Irregular node there
	    // compute coordinates of new irregular
	    coords = *(mesh->giveNode(elem->giveNode(inode))->giveCoordinates());
	    coords.add (mesh->giveNode(elem->giveNode(jnode))->giveCoordinates());
	    coords.times(0.5);
	    // compute required density of a new node
	    density = 0.5*(mesh->giveNode(elem->giveNode(inode))->giveRequiredDensity()+
			   mesh->giveNode(elem->giveNode(jnode))->giveRequiredDensity());
	    // add irregular to receiver 
	    iNum = mesh->giveNumberOfNodes() + 1;
	    irregular = new Subdivision::RS_IrregularNode (iNum, 0, coords, density);
	    irregular->setPartition (iproc);
	    irregular->setParallelMode(DofManager_shared);
	    mesh->addNode (iNum, irregular);
	    elem->addIrregularOn (iNum, bEnt);
	    irregular->setEdgeNodes(bEnt.at(1), bEnt.at(2));
	    // schedule elem for bisection
	    subdivqueue.push(ie);
	    found = true;
#ifdef __VERBOSE_PARALLEL
	    OOFEM_LOG_INFO ("added as %d on %d elem\n", iNum, ie);
#endif
	    break;
	  } 
	}
      }
      if (!found) OOFEM_ERROR ("Subdivision::unpackSharedIrregulars: no element found\n");
    } else {
      OOFEM_ERROR ("Subdivision::unpackSharedIrregulars: unknown tag received");
    }
    // get type of the next record
    pcbuff->unpackInt(_type);
  }; // while (_type != LOADBALANCER_END_DATA);
  
}

void
Subdivision :: assignGlobalNumbersToSharedIrregulars () {
  // there are some shared irregulars -> data exchange
  CommunicatorBuff cb(this->giveNumberOfProcesses(), CBT_dynamic);
  Communicator com(domain->giveEngngModel(), &cb, this->giveRank(), 
		   this->giveNumberOfProcesses(), CommMode_Dynamic);
  
  com.packAllData(this, this, &Subdivision::packIrregularSharedGlobnums);
  com.initExchange(SHARED_IRREGULAR_DATA_TAG);
  com.unpackAllData(this, this, &Subdivision::unpackIrregularSharedGlobnums);
  
}  

int
Subdivision::packIrregularSharedGlobnums (Subdivision *s, ProcessCommunicator &pc)
{
  int rproc = pc.giveRank();
  int i, inode, jnode, nnodes=mesh->giveNumberOfNodes();
  int myrank = this->giveRank();
  IntArray edgeInfo(3);
  const IntArray *sharedPartitions;
  if ( rproc == myrank ) return 1; // skip local partition
  
  // query process communicator to use
  ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();
  for (i=1; i<=nnodes; i++) {
    if (this->isNodeLocalSharedIrregular(mesh->giveNode(i), myrank)) {
      sharedPartitions = this->mesh->giveNode(i)->givePartitions();
      if (sharedPartitions->contains(rproc)) {
	
	// the info about new local shared irregular node needs to be sent to remote partition
	// the new ireegular on remote partition is identified using two nodes defining 
	// an edge on which irregular node is introduced
	((RS_IrregularNode*)this->mesh->giveNode(i))->giveEdgeNodes(inode,jnode);
	edgeInfo.at(1) = this->mesh->giveNode(i)->giveGlobalNumber();
	edgeInfo.at(2) = this->mesh->giveNode(inode)->giveGlobalNumber();
	edgeInfo.at(3) = this->mesh->giveNode(jnode)->giveGlobalNumber();
#ifdef __VERBOSE_PARALLEL
	OOFEM_LOG_INFO("[%d]packIrregularSharedGlobnums: sending [%d][%d %d]\n",myrank, edgeInfo.at(1), edgeInfo.at(2), edgeInfo.at(3) );
#endif
	pcbuff->packInt (SUBIVISION_IRREGULAR_REC_TAG);
	pcbuff->packIntArray (edgeInfo);
      }
    }
  }
  pcbuff->packInt (SUBDIVISION_END_DATA);
}


int 
Subdivision::unpackIrregularSharedGlobnums (Subdivision *s, ProcessCommunicator &pc)
{
  bool found;
  int myrank = this->giveRank();
  int ie, _type, iproc = pc.giveRank();
  int inode, jnode, irregNum, nelems=mesh->giveNumberOfElements();
  IntArray edgeInfo(3), be(2);
  RS_Element *elem;
  RS_Node *node;

  if ( iproc == myrank ) {
    return 1;                // skip local partition
  }
  
  // query process communicator to use
  ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();
  
  pcbuff->unpackInt(_type);
  // unpack dofman data
  while ( _type != SUBDIVISION_END_DATA ) {
    if (_type == SUBIVISION_IRREGULAR_REC_TAG) {
      pcbuff->unpackIntArray (edgeInfo);
      
      // look for local element containing an edge identical to the one just received
      for (found = false, ie=1; ie<=nelems; ie++) { 
	elem = mesh->giveElement(ie);
	if (elem->giveParallelMode() == Element_remote) continue;
	if (!elem->isTerminal()) continue;
	if ((inode=elem->containsGlobalNode(edgeInfo.at(2))) &&
	    (jnode=elem->containsGlobalNode(edgeInfo.at(3)))) {
	  // element found => set globnum to corresponding irregular
	  be.at(1) = elem->giveNode(inode); be.at(2) = elem->giveNode(jnode);
	  irregNum = elem->giveIrregular (be);
	  if (irregNum) {
	    node = mesh->giveNode(irregNum);
	    node->setGlobalNumber(edgeInfo.at(1));
#ifdef __VERBOSE_PARALLEL
	    OOFEM_LOG_INFO("[%d]unpackIrregularSharedGlobnums: received [%d][%d %d] as %d\n",myrank, edgeInfo.at(1), edgeInfo.at(2), edgeInfo.at(3), irregNum);
#endif
	    found = true;
	    break;
	  } else {
	    OOFEM_ERROR ("unpackIrregularSharedGlobnums: internal inconsistency, element irregular not found\n");
	  }
	}
      }
      if (!found) {
	char buff[1024];
	sprintf (buff, "[%d] unpackIrregularSharedGlobnums: nonmatching record received %d [%d %d]", myrank, edgeInfo.at(1), edgeInfo.at(2), edgeInfo.at(3));
	OOFEM_ERROR (buff);
      }

    } else {
      OOFEM_ERROR ("Subdivision::unpackIrregularSharedGlobnums: unknown tag received");
    }
    pcbuff->unpackInt(_type);
  }
}
 
bool
Subdivision::isNodeLocalSharedIrregular (Subdivision::RS_Node* node, int myrank) {
  if (node->isIrregular()) {
    if (node->giveParallelMode() == DofManager_shared) {
      int i, minpart, npart;
      const IntArray *partitions = node->givePartitions();
      npart = partitions->giveSize();
      minpart = myrank;
      for (i=1; i<=npart; i++) minpart = min (minpart, partitions->at(i));
      if (minpart == myrank) return true; else return false;
    } else return false;
  } else {
    return false;
  }
}
 

bool
Subdivision::isNodeLocalIrregular (Subdivision::RS_Node* node, int myrank) {
  if (node->isIrregular()) {
    if (node->giveParallelMode() == DofManager_local) {
      return true;
    } else if (node->giveParallelMode() == DofManager_shared) {
      int i, minpart, npart;
      const IntArray *partitions = node->givePartitions();
      npart = partitions->giveSize();
      minpart = myrank;
      for (i=1; i<=npart; i++) minpart = min (minpart, partitions->at(i));
      if (minpart == myrank) return true; else return false;
    } else return false;
  } else {
    return false;
  }
}


int
Subdivision::giveRank() 
{
  return this->domain->giveEngngModel()->giveRank();
}

int 
Subdivision::giveNumberOfProcesses()
{
  return this->domain->giveEngngModel()->giveNumberOfProcesses();
}


/* 
   Idea is to rebuild the remote elements from scratch after subdivision
   Reason is that sbdivision in mot cases followed by smoothing, 
   so that it becomes difficult to reconstruct the remot elements from old remote elements
   
   the algorith should be summarized as follows:
   1. For each remote partition:
      For each node (INODE) shared with remote partition:
           Add all elements sharing INODE into the queue:
           For each element in the queue:
	       Mark element as remote 
	       For each node (JNODE) of this element:
	          Compute its distance from INODE
	          If less than nonlocal radius:
		      Add all elements sharing JNODE to queue

    This is an aproximate algrithm, that works correctly for reasonably shaped grids,
    for grids with badly shaped elements, some interactions can be neglected, due to the fact
    that real distances between intagration points are not considered directly.
*/
void
Subdivision::exchangeRemoteElements (Domain* d)
{
  int nproc = this->giveNumberOfProcesses(), myrank = this->giveRank();
  CommunicatorBuff cb(nproc, CBT_dynamic);
  Communicator com(d->giveEngngModel(), &cb, myrank, nproc, CommMode_Dynamic);
  
  // move existing dofmans and elements, that will be local on current partition,
  // into local map
  com.packAllData(this, d, & Subdivision :: packRemoteElements);
  com.initExchange(SUBDIVISION_MIGRATE_REMOTE_ELEMENTS_TAG);

  // remove existing remote elements and null nodes
  int i, nelem = d->giveNumberOfElements();
  int nnodes = d->giveNumberOfDofManagers();
  DomainTransactionManager *dtm = d->giveTransactionManager();
  for (i=1; i<=nnodes; i++) {
    if (d->giveDofManager(i)->giveParallelMode() == DofManager_null) {
      dtm->addTransaction(DomainTransactionManager :: DTT_Remove, DomainTransactionManager :: DCT_DofManager, d->giveDofManager(i)->giveGlobalNumber(), NULL);
    }
  }

  for (i=1; i<=nelem; i++) {
    if (d->giveElement(i)->giveParallelMode() == Element_remote) {
      dtm->addTransaction(DomainTransactionManager :: DTT_Remove, DomainTransactionManager :: DCT_Element, d->giveElement(i)->giveGlobalNumber(), NULL);
    }
  }


  // receive remote data
  com.unpackAllData(this, d, & Subdivision :: unpackRemoteElements);
}

int
Subdivision::packRemoteElements (Domain *d, ProcessCommunicator &pc)
{
  int rproc = pc.giveRank();
  int i, inode, ec, elem, nnodes=d->giveNumberOfDofManagers();
  int j, elemNodes, nn, in, nmat;
  int myrank = this->giveRank();
  bool nonlocal = false;
  double radius, _radius, dist;
  classType dtype;
  ConnectivityTable *ct = d->giveConnectivityTable();
  DofManager *jnodePtr, *inodePtr;
  Element *elemPtr, *relemPtr;
  NonlocalMaterialExtensionInterface *itrfc;
  const IntArray *inodeConnectivity, *jnodeConnectivity;
  std::queue<int> elemCandidates;
  std::set<int> remoteElements, processedElements, nodesToSend;
  std::set<int>::const_iterator si;

  if ( rproc == myrank ) return 1; // skip local partition

  // determine the maximum interaction radius
  radius = 0.0;
  nmat = d->giveNumberOfMaterialModels();
  for (i=1; i<=nmat; i++) {
    itrfc = (NonlocalMaterialExtensionInterface*) 
      d->giveMaterial(i)->giveInterface(NonlocalMaterialExtensionInterfaceType);
    if (itrfc) {
      itrfc->giveSupportRadius(_radius);
      radius = max (radius, _radius);
      nonlocal=true;
    }
  }
  
  // assemble the list of elements to send 
  // step 1
  for (i=1; i<=nnodes; i++) {
    inodePtr = d->giveDofManager(i);
    if (d->giveDofManager(i)->isShared() && d->giveDofManager(i)->givePartitionList()->contains(rproc)) {
      // found INODE
      while (!elemCandidates.empty()) elemCandidates.pop();  // clear candiate queue
      processedElements.clear();
      inodeConnectivity = ct->giveDofManConnectivityArray(i);
      for (ec = 1; ec<= inodeConnectivity->giveSize(); ec++) { // loop over element candidates
	  elemCandidates.push(inodeConnectivity->at(ec));
      }
      
      // step 2:
      while (!elemCandidates.empty()) {
	elem = elemCandidates.front(); elemCandidates.pop();  
	elemPtr = d->giveElement(elem);
	if (elemPtr->giveParallelMode() != Element_local) continue;

	remoteElements.insert(elem);
	elemNodes = d->giveElement(elem)->giveNumberOfDofManagers();
	for (j=1; j<= elemNodes; j++) {
	  jnodePtr = elemPtr->giveDofManager(j);
	  if (jnodePtr->hasCoordinates()) {
	    dist = jnodePtr->giveCoordinates()->distance(inodePtr->giveCoordinates());
	    if (dist < radius) {
	      jnodeConnectivity = ct->giveDofManConnectivityArray(jnodePtr->giveNumber());
	      for (ec = 1; ec<= jnodeConnectivity->giveSize(); ec++) { // loop over element candidates
		  if (processedElements.find(jnodeConnectivity->at(ec)) == processedElements.end()) {
		    elemCandidates.push(jnodeConnectivity->at(ec));
		    processedElements.insert(jnodeConnectivity->at(ec));
		  }
	      }
	    }
	  }
	}
      } // end while (!elemCandidates.empty())
    }
  }  // end for (i=1; i<=nnodes; i++)

  // now the list of elements to became remote on given remote partition is in remoteElements set

  // need to pack elements definition and corresponding nodes!
  /* 
     mark shecheduled nodes:
     loop over elements in remoteElements set and add all their nodes (except those that are shared)
  */
  for (si = remoteElements.begin(); si != remoteElements.end(); si++) {
    relemPtr = d->giveElement(*si);
    nn = relemPtr->giveNumberOfNodes();
    for (in = 1; in <= nn; in ++) {
      inode=relemPtr->giveDofManagerNumber(in);
      if (d->giveDofManager(inode)->giveParallelMode() == DofManager_local) {
	nodesToSend.insert(inode);
      }
    }
  }


  //-----------------end here-------------------

  // query process communicator to use
  ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();
  ProcessCommDataStream pcDataStream(pcbuff);
  // send nodes that define remote_elements gometry
  for (si=nodesToSend.begin(); si!=nodesToSend.end(); si++) {
    inodePtr = d->giveDofManager(*si);
    dtype = inodePtr->giveClassID();
 
    pcbuff->packInt(dtype);
    pcbuff->packInt(inodePtr->giveGlobalNumber() );
    inodePtr->saveContext(& pcDataStream, CM_Definition);
  }

  // pack end-of-element-record
  pcbuff->packInt(SUBDIVISION_END_DATA);


  // send elements
  for (si=remoteElements.begin(); si != remoteElements.end(); si++) {
    elemPtr = d->giveElement(*si);
    // pack local element (node numbers shuld be global ones!!!)
    // pack type
    pcbuff->packInt(elemPtr->giveClassID() );
    // nodal numbers shuld be packed as global !!
    elemPtr->saveContext(& pcDataStream, CM_Definition | CM_DefinitionGlobal);
    //OOFEM_LOG_INFO ("[%d] Sending Remote elem %d[%d] to rank %d\n", myrank,*si, elemPtr->giveGlobalNumber(), rproc );
  }
  
  // pack end-of-element-record
  pcbuff->packInt(SUBDIVISION_END_DATA);

  return 1;
}


int 
Subdivision::unpackRemoteElements (Domain *d, ProcessCommunicator &pc)
{
  int myrank = d->giveEngngModel()->giveRank();
  int iproc = pc.giveRank();
  int _globnum, _type;
  bool _newentry;
  classType _etype;
  DofManager *dofman;
  DomainTransactionManager *dtm = d->giveTransactionManager();

  if ( iproc == myrank ) {
    return 1;                // skip local partition
  }
  
  // query process communicator to use
  ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();
  ProcessCommDataStream pcDataStream(pcbuff);
  
  pcbuff->unpackInt(_type);
  // unpack dofman data
  while ( _type != SUBDIVISION_END_DATA ) {
    _etype = ( classType ) _type;
    // receiving new local dofManager
    pcbuff->unpackInt(_globnum);

    _newentry = false;
    if ( ( dofman = dtm->giveDofManager(_globnum) ) == NULL ) {
      // data not available -> create a new one
      _newentry = true;
      dofman = CreateUsrDefDofManagerOfType(_etype, 0, d);
    }

    dofman->setGlobalNumber(_globnum);
    // unpack dofman state (this is the local dofman, not available on remote)
    dofman->restoreContext(& pcDataStream, CM_Definition );
    dofman->setParallelMode(DofManager_null);
    // add transaction if new entry allocated; otherwise existing one has been modified via returned dofman
    if ( _newentry ) {
      dtm->addTransaction(DomainTransactionManager :: DTT_ADD, DomainTransactionManager :: DCT_DofManager, _globnum, dofman);
    }
    pcbuff->unpackInt(_type);
  }

  Element *elem;
  IntArray elemPartitions(1); elemPartitions.at(1) = iproc;
  int nrecv = 0;
  
  do {
    pcbuff->unpackInt(_type);
    if ( _type == SUBDIVISION_END_DATA ) {
      break;
    }

    _etype = ( classType ) _type;
    elem = CreateUsrDefElementOfType(_etype, 0, d);
    elem->restoreContext(& pcDataStream, CM_Definition | CM_DefinitionGlobal);
    elem->setParallelMode (Element_remote);
    elem->setPartitionList (elemPartitions);
    dtm->addTransaction(DomainTransactionManager :: DTT_ADD, DomainTransactionManager :: DCT_Element, elem->giveGlobalNumber(), elem);
    nrecv++;
    //OOFEM_LOG_INFO ("[%d] Received Remote elem [%d] to rank %d\n", myrank, elem->giveGlobalNumber(), iproc );
    //recvElemList.push_back(elem);
  } while ( 1 );

  return 1;
}

void
Subdivision::assignGlobalNumbersToElements (Domain *d)
{
  int problem_size=this->giveNumberOfProcesses();
  int myrank=this->giveRank();
  int i, nelems, numberOfLocalElementsToNumber = 0, partitionNumberOfElements[problem_size];
  int localMaxGlobnum = 0, globalMaxGlobnum;

  // idea: first determine the number of local elements waiting for new global id
  // and also determine max global number assigned up to now
  nelems=d->giveNumberOfElements();
  for (i=1; i<=nelems; i++) {
    localMaxGlobnum = max (localMaxGlobnum, d->giveElement(i)->giveGlobalNumber());
    if ((d->giveElement(i)->giveParallelMode() == Element_local) &&
	(d->giveElement(i)->giveGlobalNumber() <= 0)) {
      numberOfLocalElementsToNumber++;
    }
  }
  // determine number of elements across all partitions
  MPI_Allgather (&numberOfLocalElementsToNumber, 1, MPI_INT, 
		 partitionNumberOfElements, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allreduce (&localMaxGlobnum, &globalMaxGlobnum, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#ifdef __VERBOSE_PARALLEL
  OOFEM_LOG_INFO("[%d]Subdivision::assignGlobalNumbersToElements: max globnum %d, new elements %d\n",myrank,globalMaxGlobnum,numberOfLocalElementsToNumber);
#endif

  // compute local offset
  int startOffset = globalMaxGlobnum, availGlobNum;
  for (i=0; i<myrank; i++) startOffset+= partitionNumberOfElements[i];
  // ok. lets assign global numbers on each partition to local elements
  availGlobNum = startOffset;
  for (i=1; i<=nelems; i++) {
    if ((d->giveElement(i)->giveParallelMode() == Element_local) && 
	(d->giveElement(i)->giveGlobalNumber() <= 0)) {
      d->giveElement(i)->setGlobalNumber (++availGlobNum);
    }
  }

#ifdef __VERBOSE_PARALLEL
  char *elemmode[] = {"local", "remote"};
  for (i=1; i<=nelems; i++) {
    OOFEM_LOG_INFO ("[%d] Element %d[%d], %s\n", myrank, i,d->giveElement(i)->giveGlobalNumber(),  elemmode[(int) d->giveElement(i)->giveParallelMode()]);
  }
#endif
}

#endif
