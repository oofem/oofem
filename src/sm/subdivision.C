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

//#define __VERBOSE_PARALLEL

//#define DEBUG_INFO
#define DEBUG_CHECK

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
Subdivision::RS_Triangle::addIrregular (int neighborElement, int node)
{
  int ind;
  if ((ind = neghbours_base_elements.findFirstIndexOf(neighborElement))) {
    irregular_nodes.at(ind) = node;
  } else {
    OOFEM_ERROR4("Subdivision::RS_Triangle::addIrregular: element %d is not neighbor %d of element %d", 
								 neighborElement, neighborElement, this->number);
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


int
Subdivision::RS_Element::giveTopParent()
{
	RS_Element *elem=this;
	while(elem->giveNumber()!=elem->giveParent())elem=mesh->giveElement(elem->giveParent());
	return(elem->giveNumber());
}


Subdivision::RS_Triangle::RS_Triangle (int number, RS_Mesh* mesh, int parent, IntArray& nodes) : 
  Subdivision::RS_Element (number, mesh, parent, nodes) {

  irregular_nodes.resize(3); irregular_nodes.zero();
  neghbours_base_elements.resize(3); 
  neghbours_base_elements.zero();

  double elength1, elength2, elength3;
  elength1 = mesh->giveNode(this->nodes.at(1))->giveCoordinates()->distance_square(*(mesh->giveNode(this->nodes.at(2))->giveCoordinates()));
  elength2 = mesh->giveNode(this->nodes.at(2))->giveCoordinates()->distance_square(*(mesh->giveNode(this->nodes.at(3))->giveCoordinates()));
  elength3 = mesh->giveNode(this->nodes.at(3))->giveCoordinates()->distance_square(*(mesh->giveNode(this->nodes.at(1))->giveCoordinates()));
  leIndex = 1;
  if (elength2>elength1) {
    leIndex = 2;
    if (elength3>elength2) leIndex = 3;
  } else if (elength3>elength1) {
    leIndex = 3;
  }
}


Subdivision::RS_Tetra::RS_Tetra (int number, RS_Mesh* mesh, int parent, IntArray& nodes) : 
  Subdivision::RS_Element (number, mesh, parent, nodes) {

  irregular_nodes.resize(6); irregular_nodes.zero();
  neghbours_base_elements.resize(6); 
  neghbours_base_elements.zero();

#ifdef DEBUG_CHECK
	if(nodes.findFirstIndexOf(0)){
		OOFEM_ERROR2("Subdivision::RS_Tetra::RS_Tetra: 0 node of element %d", this->number);
  }
#endif

  double elength1, elength2, elength3, elength4, elength5, elength6;
	double maxlength, mlength;
  elength1 = mesh->giveNode(this->nodes.at(1))->giveCoordinates()->distance_square(*(mesh->giveNode(this->nodes.at(2))->giveCoordinates()));
  elength2 = mesh->giveNode(this->nodes.at(2))->giveCoordinates()->distance_square(*(mesh->giveNode(this->nodes.at(3))->giveCoordinates()));
  elength3 = mesh->giveNode(this->nodes.at(3))->giveCoordinates()->distance_square(*(mesh->giveNode(this->nodes.at(1))->giveCoordinates()));
  elength4 = mesh->giveNode(this->nodes.at(1))->giveCoordinates()->distance_square(*(mesh->giveNode(this->nodes.at(4))->giveCoordinates()));
  elength5 = mesh->giveNode(this->nodes.at(2))->giveCoordinates()->distance_square(*(mesh->giveNode(this->nodes.at(4))->giveCoordinates()));
  elength6 = mesh->giveNode(this->nodes.at(3))->giveCoordinates()->distance_square(*(mesh->giveNode(this->nodes.at(4))->giveCoordinates()));
	// get largest edge
  leIndex = 1; maxlength = elength1;
  if (elength2>maxlength) {leIndex = 2; maxlength=elength2;}
  if (elength3>maxlength) {leIndex = 3; maxlength=elength3;}
  if (elength4>maxlength) {leIndex = 4; maxlength=elength4;}
  if (elength5>maxlength) {leIndex = 5; maxlength=elength5;}
  if (elength6>maxlength) {leIndex = 6; maxlength=elength6;}
	// get second largest edge
	mlength = 0.0;
	if(leIndex != 1 && elength1 > mlength)mlength=elength1;
	if(leIndex != 2 && elength2 > mlength)mlength=elength2;
	if(leIndex != 3 && elength3 > mlength)mlength=elength3;
	if(leIndex != 4 && elength4 > mlength)mlength=elength4;
	if(leIndex != 5 && elength5 > mlength)mlength=elength5;
	if(leIndex != 6 && elength6 > mlength)mlength=elength6;

	// IMPORTANT: ambiguity must be handled !!!

	// in order to achieve that shared tetra faces are subdivided in compatible way, always the same longest edge must be identified;
	// this may be problem if there are two or more edges of the "same" length
	// then the longest is dependent on the order in which edges are processed;
	// also in order to get always the same distance of two nodes, the nodes should be processes in the same order;
	
	// however it is not clear whether one may rely that the same distance is always obtained for the same pair of nodes
	// (in proper order) especially if processed on different processors

	// if the seccond largest edge is not equal to the largest edge and if the largest and second largest edges
	// are opposite edges, no action is theoretically required;
  // but because it is not clear whether different result would be obtained if edge lengths are evaluated
	// using node ordering, this possibility (of no action) is not considered

	if(mlength > maxlength * 0.9999999){
		int i, l_nd[4];
		int nd[4] = {0, 1, 2, 3}, ed[4][4] = {{0, 1, 3, 4}, {1, 0, 2, 5}, {3, 2, 0, 6}, {4, 5, 6, 0}};
#ifdef __PARALLEL_MODE
		int g_nd[4];
#endif
		bool swap;

		// sort node ids
		for(i = 0; i < 4; i++){
			l_nd[i] = nodes.at(i+1);
#ifdef __PARALLEL_MODE
			g_nd[i] = mesh->giveNode(l_nd[i])->giveGlobalNumber();
#endif
		}
		swap = true;
		while(swap){
			swap = false;
			for(i = 0; i < 3; i++){
#ifdef __PARALLEL_MODE
				if(g_nd[i] > g_nd[i+1])
#else
				if(l_nd[i] > l_nd[i+1])
#endif
				{
					int tmp;
#ifdef __PARALLEL_MODE
					tmp = g_nd[i]; g_nd[i] = g_nd[i+1]; g_nd[i+1] = tmp;
#endif
					tmp = l_nd[i]; l_nd[i] = l_nd[i+1]; l_nd[i+1] = tmp;
					tmp = nd[i]; nd[i] = nd[i+1]; nd[i+1] = tmp;
					swap = true;
				}
			}
		}

		// evaluate again all edges because not only the largest two may be of the same length
		// process them from smallest node id
		elength1 = mesh->giveNode(l_nd[0])->giveCoordinates()->distance_square(*(mesh->giveNode(l_nd[1])->giveCoordinates()));
		elength2 = mesh->giveNode(l_nd[0])->giveCoordinates()->distance_square(*(mesh->giveNode(l_nd[2])->giveCoordinates()));
		elength3 = mesh->giveNode(l_nd[0])->giveCoordinates()->distance_square(*(mesh->giveNode(l_nd[3])->giveCoordinates()));
		elength4 = mesh->giveNode(l_nd[1])->giveCoordinates()->distance_square(*(mesh->giveNode(l_nd[2])->giveCoordinates()));
		elength5 = mesh->giveNode(l_nd[2])->giveCoordinates()->distance_square(*(mesh->giveNode(l_nd[3])->giveCoordinates()));
		elength6 = mesh->giveNode(l_nd[3])->giveCoordinates()->distance_square(*(mesh->giveNode(l_nd[3])->giveCoordinates()));

		// get largest edge
		maxlength = 0.0;
		if (elength1>maxlength) {maxlength=elength1; leIndex = ed[nd[0]][nd[1]];}
		if (elength2>maxlength) {maxlength=elength2; leIndex = ed[nd[0]][nd[2]];}
		if (elength3>maxlength) {maxlength=elength3; leIndex = ed[nd[0]][nd[3]];}
		if (elength4>maxlength) {maxlength=elength4; leIndex = ed[nd[1]][nd[2]];}
		if (elength5>maxlength) {maxlength=elength5; leIndex = ed[nd[1]][nd[3]];}
		if (elength6>maxlength) {maxlength=elength6; leIndex = ed[nd[2]][nd[3]];}
	}
}


int
Subdivision::RS_Triangle::giveEdgeIndex(int inode, int jnode) {
	/* returns zero, if triangle does not have inode or jnode */
	int in=0, jn=0;

	int k;
	for(k=1;k<=3;k++){
		if(nodes.at(k)==inode)in=k;
		if(nodes.at(k)==jnode)jn=k;
	}

	if(in==0 || jn==0)return 0;

	if(in < jn){
		if(jn == in+1)return(in);
		return(3);
	} else {
		if(in == jn+1)return(jn);
		return(3);
	}
}


int
Subdivision::RS_Tetra::giveEdgeIndex(int inode, int jnode) {
	/* returns zero, if tetra does not have inode or jnode */
	int in=0, jn=0;

	int k;
	for(k=1;k<=4;k++){
		if(nodes.at(k)==inode)in=k;
		if(nodes.at(k)==jnode)jn=k;
	}

	if(in==0 || jn==0)return 0;

	if(in < jn){
		if(jn==4)return(in+3);
		if(jn == in+1)return(in);
		return(3);
	} else {
		if(in==4)return(jn+3);
		if(in == jn+1)return(jn);
		return(3);
	}
}


void
Subdivision::RS_Triangle::bisect(std::queue<int> &subdivqueue, std::list<int> &sharedIrregularsQueue) {
  /* this is symbolic bisection - no new elements are added, only irregular nodes are introduced */
  int inode, jnode;
  double density;
	bool boundary=false;
  FloatArray coords;
  Subdivision::RS_Triangle* triangleElem;

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

		/*
		// check whether new node is boundary
		if(this->neghbours_base_elements.at(leIndex)){
#ifdef __PARALLEL_MODE
			if(mesh->giveElement(this->neghbours_base_elements.at(leIndex))->giveParallelMode() == Element_local){
#endif
        Domain *dorig = mesh->giveSubdivision()->giveDomain();
				if(mesh->giveNode(nodes.at(inode))->isBoundary() || mesh->giveNode(nodes.at(jnode))->isBoundary()){
					if(dorig->giveElement(this->giveTopParent())->giveRegionNumber() != dorig->giveElement(mesh->giveElement(this->neghbours_base_elements.at(leIndex))->giveTopParent())->giveRegionNumber()) boundary=true;
				} 
#ifdef __PARALLEL_MODE
			} else {
				boundary=true;
			}
#endif
		} else {
			boundary=true;
		}
		*/

		// quick hack to simulate shared boundary in sequential analysis
		// important: a single element can have only two original nodes marked as boundary !!!
		if(mesh->giveNode(nodes.at(inode))->isBoundary() && mesh->giveNode(nodes.at(jnode))->isBoundary())boundary=true; // HUHU

    // add irregular to receiver and its neighbour
    iNum = mesh->giveNumberOfNodes() + 1;
    irregular = new Subdivision::RS_IrregularNode (iNum, 0, coords, density, boundary);
    mesh->addNode (iNum, irregular);

#ifdef __PARALLEL_MODE
      OOFEM_LOG_INFO ("[%d] RS_Triangle::bisecting %d[%d] nodes %d %d %d, lindex %d, new irregular %d\n",mesh->giveSubdivision()->giveRank(),this->number, this->giveGlobalNumber(),nodes.at(1), nodes.at(2), nodes.at(3), leIndex, iNum);
#else
      OOFEM_LOG_INFO ("RS_Triangle::bisecting %d, new irregular %d\n",this->number, iNum);
#endif

#ifdef __OOFEG
    irregular->drawGeometry();
#endif

    this->irregular_nodes.at(leIndex) = iNum;
    if (this->neghbours_base_elements.at(leIndex)) {
#ifdef __PARALLEL_MODE
			if (mesh->giveElement(this->neghbours_base_elements.at(leIndex))->giveParallelMode() == Element_local){
#endif
        Subdivision::RS_Triangle* triangleElem;
        triangleElem = dynamic_cast< Subdivision::RS_Triangle*> (mesh->giveElement(this->neghbours_base_elements.at(leIndex)));
        if (triangleElem) {
          triangleElem ->addIrregular(this->number, iNum);
        } else {
          OOFEM_ERROR2 ("RS_Triangle::bisect: dynamic cast failed on %d element", this->neghbours_base_elements.at(leIndex));
        }
				// add neighbour to list of elements for subdivision
				subdivqueue.push(this->neghbours_base_elements.at(leIndex));
#ifdef __PARALLEL_MODE
			}
#endif
		}

#ifdef __PARALLEL_MODE
    IntArray partitions;
    // test if new node is on shared interpartition boundary
    if (this->giveIrregularSharedPartitions (leIndex, inode, jnode, partitions)) {
      // mark irregular as shared (can be shared by only one partition (in 2D))
      irregular->setParallelMode (DofManager_shared);
      irregular->setPartition (partitions.at(1));
      // and put its number into queue of shared irregulars that is later used to inform remote partition about this fact
      irregular->setEdgeNodes(nodes.at(inode), nodes.at(jnode));
      sharedIrregularsQueue.push_back(iNum);
#ifdef __VERBOSE_PARALLEL
      OOFEM_LOG_INFO("RS_Triangle::bisect: Shared Irregular detected, number %d [%d[%d] %d[%d]], elem %d, partition %d\n", iNum, nodes.at(inode), mesh->giveNode(nodes.at(inode))->giveGlobalNumber(), nodes.at(jnode), mesh->giveNode(nodes.at(jnode))->giveGlobalNumber(), this->number, partitions.at(1));
#endif
    }
#endif
		
  }
}


void
Subdivision::RS_Tetra::bisect(std::queue<int> &subdivqueue, std::list<int> &sharedIrregularsQueue) {
  /* this is symbolic bisection - no new elements are added, only irregular nodes are introduced */
  int inode, jnode, ngb1, ngb2, eIndex, reg;
  double density;
	bool boundary=false, iboundary, jboundary;
	Subdivision::RS_Element *elem1, *elem2;
  FloatArray coords;
  Domain *dorig;

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
		if(leIndex <= 3){
			inode = leIndex;
			jnode = (leIndex<3)?inode+1:1;
			ngb1 = 1;
			ngb2 = inode+1;
		}
		else{
			inode = leIndex-3;
			jnode = 4;
			ngb1 = inode+1;
			ngb2 = (inode>1)?inode:4;
		}
    // compute coordinates of new irregular
    coords = *(mesh->giveNode(nodes.at(inode))->giveCoordinates());
    coords.add (mesh->giveNode(nodes.at(jnode))->giveCoordinates());
    coords.times(0.5);
    // compute required density of a new node
    density = 0.5*(mesh->giveNode(nodes.at(inode))->giveRequiredDensity()+
									 mesh->giveNode(nodes.at(jnode))->giveRequiredDensity());
    // add irregular to receiver and its neighbour
    iNum = mesh->giveNumberOfNodes() + 1;
    irregular = new Subdivision::RS_IrregularNode (iNum, 0, coords, density, boundary);
    mesh->addNode (iNum, irregular);
#ifdef __OOFEG
    irregular->drawGeometry();
#endif

    this->irregular_nodes.at(leIndex) = iNum;

#ifdef DEBUG_INFO
		fprintf(stderr, "Irregular %d added on %d (edge %d, nodes %d %d, nds %d %d %d %d, ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
						iNum, this->number, leIndex, nodes.at(inode), nodes.at(jnode), 
						nodes.at(1), nodes.at(2), nodes.at(3), nodes.at(4), 
						neghbours_base_elements.at(1), neghbours_base_elements.at(2), 
						neghbours_base_elements.at(3), neghbours_base_elements.at(4),
						irregular_nodes.at(1), irregular_nodes.at(2), irregular_nodes.at(3), 
						irregular_nodes.at(4), irregular_nodes.at(5), irregular_nodes.at(6));
#endif

		IntArray bEnt(2);
		bEnt.at(1)=nodes.at(inode);
		bEnt.at(2)=nodes.at(jnode);

		iboundary=mesh->giveNode(nodes.at(inode))->isBoundary();
		jboundary=mesh->giveNode(nodes.at(jnode))->isBoundary();
		boundary=false;
		// traverse neighbours in one direction
		elem1=this;
		if(iboundary == true || jboundary == true) {
      dorig = mesh->giveSubdivision()->giveDomain();
      reg=dorig->giveElement(this->giveTopParent())->giveRegionNumber();
    }
		while(elem1->giveNeighbor(ngb1)){
			elem2=mesh->giveElement(elem1->giveNeighbor(ngb1));
			if(elem2==this){
				ngb2=0; // there is no need to traverse in the other direction (loop was completed)
				break;
			}
			eIndex=elem2->giveEdgeIndex(bEnt.at(1),bEnt.at(2));
			// do not stop traversal if the neighbour is remote
			// there still may be local elements further
#ifdef __PARALLEL_MODE
			if(elem2->giveParallelMode() == Element_local){
#endif
				elem2->addIrregularOn(iNum,bEnt);
				// add neighbour to list of elements for subdivision
				subdivqueue.push(elem2->giveNumber());

#ifdef DEBUG_INFO
				fprintf(stderr, "Irregular %d added on (ngb1) %d (edge %d, nodes %d %d, nds %d %d %d %d, ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
								iNum, elem2->giveNumber(), eIndex, bEnt.at(1), bEnt.at(2), 
								elem2->giveNode(1), elem2->giveNode(2), elem2->giveNode(3), elem2->giveNode(4), 
								elem2->giveNeighbor(1), elem2->giveNeighbor(2), elem2->giveNeighbor(3), elem2->giveNeighbor(4),
								elem2->giveIrregular(1), elem2->giveIrregular(2), elem2->giveIrregular(3), 
								elem2->giveIrregular(4), elem2->giveIrregular(5), elem2->giveIrregular(6));
#endif

				if(!boundary){
					if(iboundary == true || jboundary == true){
						if(dorig->giveElement(elem2->giveTopParent())->giveRegionNumber()!=reg)boundary=true;
					}
				}
#ifdef __PARALLEL_MODE
			} else {
				boundary=true;
			}
#endif

			if(eIndex <= 3){
				if(elem2->giveNeighbor(1)==elem1->giveNumber())
					ngb1=eIndex+1;
				else
					ngb1=1;
			}
			else{
				if(elem2->giveNeighbor(eIndex-2)==elem1->giveNumber())
					ngb1=(eIndex>4)?eIndex-3:4;
				else
					ngb1=eIndex-2;
			}
			elem1=elem2;
		}

		if(ngb2){
			boundary=true;
			// traverse neighbours in other direction
			elem1=this;
			while(elem1->giveNeighbor(ngb2)){
				elem2=mesh->giveElement(elem1->giveNeighbor(ngb2));
				eIndex=elem2->giveEdgeIndex(bEnt.at(1),bEnt.at(2));
				// do not stop traversal if the neighbour is remote
				// there still may be local elements further
#ifdef __PARALLEL_MODE
				if(elem2->giveParallelMode() == Element_local){
#endif
					elem2->addIrregularOn(iNum,bEnt);
					// add neighbour to list of elements for subdivision
					subdivqueue.push(elem2->giveNumber());

#ifdef DEBUG_INFO
					fprintf(stderr, "Irregular %d added on (ngb2) %d (edge %d, nodes %d %d, nds %d %d %d %d, ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
									iNum, elem2->giveNumber(), eIndex, bEnt.at(1), bEnt.at(2), 
									elem2->giveNode(1), elem2->giveNode(2), elem2->giveNode(3), elem2->giveNode(4), 
									elem2->giveNeighbor(1), elem2->giveNeighbor(2), elem2->giveNeighbor(3), elem2->giveNeighbor(4),
									elem2->giveIrregular(1), elem2->giveIrregular(2), elem2->giveIrregular(3), 
									elem2->giveIrregular(4), elem2->giveIrregular(5), elem2->giveIrregular(6));
#endif

#ifdef __PARALLEL_MODE
				}
#endif
				
				if(eIndex <= 3){
					if(elem2->giveNeighbor(1)==elem1->giveNumber())
						ngb2=eIndex+1;
					else
						ngb2=1;
				}
				else{
					if(elem2->giveNeighbor(eIndex-2)==elem1->giveNumber())
						ngb2=(eIndex>4)?eIndex-3:4;
					else
						ngb2=eIndex-2;
				}
				elem1=elem2;
			}
		}
		if(boundary)irregular->setBoundary(boundary);
			
#ifdef __PARALLEL_MODE
    IntArray partitions;
    // test if new node is on shared interpartition boundary
    if (this->giveIrregularSharedPartitions (leIndex, inode, jnode, partitions)) {
      // mark irregular as shared (can be shared by only more partitions (in 3D))
      irregular->setParallelMode (DofManager_shared);
      irregular->setPartitions (partitions);
      // and put its number into queue of shared irregulars that is later used to inform remote partition about this fact
      irregular->setEdgeNodes(nodes.at(inode), nodes.at(jnode));
      sharedIrregularsQueue.push_back(iNum);
#ifdef __VERBOSE_PARALLEL
			char buffer[1024];
			int i, len=0;
			for(i=1;i<=partitions.giveSize();i++){
				sprintf(buffer+len, " %d", partitions.at(i));
				len=strlen(buffer);
				if(len>=1024)OOFEM_ERROR ("Subdivision::RS_Tetra::bisect: too small buffer");
			}
      OOFEM_LOG_INFO("RS_Tetra::bisect: Shared Irregular detected, number %d [%d[%d] %d[%d]], elem %d, partitions%s\n", iNum, nodes.at(inode), mesh->giveNode(nodes.at(inode))->giveGlobalNumber(), nodes.at(jnode), mesh->giveNode(nodes.at(jnode))->giveGlobalNumber(), this->number, buffer);
#endif
    }
#endif

  }
}



void
Subdivision::RS_Element::addIrregularOn (int iNum, const IntArray& bEntity)
{
  if (bEntity.giveSize() != 2) OOFEM_ERROR ("Subdivision::RS_Element::addIrregularOn: bEntity size mismatch");
	int eIndex;

	eIndex=this->giveEdgeIndex(bEntity.at(1),bEntity.at(2));
	if(!eIndex)OOFEM_ERROR4 ("Subdivision::RS_Element::addIrregularOn: no bEntity connecting nodes %d and %d found on elem %d", 
													 bEntity.at(1), bEntity.at(2), this->number);
	this->irregular_nodes.at(eIndex) = iNum;
}


/*
void
Subdivision::RS_Triangle::addIrregularOn (int iNum, const IntArray& bEntity)  // HUHU
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


void
Subdivision::RS_Tetra::addIrregularOn (int iNum, const IntArray& bEntity)  // HUHU
{
  if (bEntity.giveSize() != 2) OOFEM_ERROR ("Subdivision::RS_Tetra::addIrregularOn: bEntity size mismatch");
  int anode = bEntity.at(1), bnode = bEntity.at(2), inode, jnode, ie;
  
  for (ie=1; ie<=3; ie++) { // loop over edges
    inode = ie;
    jnode = (ie<3)?ie+1: 1;
    
    if (((this->nodes.at(inode) == anode) && (this->nodes.at(jnode) == bnode)) ||
				((this->nodes.at(inode) == bnode) && (this->nodes.at(jnode) == anode))) {
      this->irregular_nodes.at(ie) = iNum;
      return;
    }
    if (((this->nodes.at(inode) == anode) && (this->nodes.at(4) == bnode)) ||
				((this->nodes.at(inode) == bnode) && (this->nodes.at(4) == anode))) {
      this->irregular_nodes.at(ie) = iNum;
      return;
    }
  }
  OOFEM_ERROR ("Subdivision::RS_Tetra::addIrregularOn: no such bEntity found");
}
*/



int
Subdivision::RS_Element::giveIrregularOn (const IntArray& bEntity)
{
  if (bEntity.giveSize() != 2) OOFEM_ERROR ("Subdivision::RS_Element::giveIrregularOn: bEntity size mismatch");
	int eIndex;

	eIndex=this->giveEdgeIndex(bEntity.at(1),bEntity.at(2));
	if(eIndex)return this->irregular_nodes.at(eIndex);
  return 0;
}


/*
int
Subdivision::RS_Triangle::giveIrregularOn (const IntArray& bEntity)  // HUHU
{
  if (bEntity.giveSize() != 2) OOFEM_ERROR ("Subdivision::RS_Triangle::giveIrregularOn: bEntity size mismatch");
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


int
Subdivision::RS_Tetra::giveIrregularOn (const IntArray& bEntity)  // HUHU
{
  if (bEntity.giveSize() != 2) OOFEM_ERROR ("Subdivision::RS_Tetra::giveIrregularOn: bEntity size mismatch");
  int anode = bEntity.at(1), bnode = bEntity.at(2), inode, jnode, ie;
  
  for (ie=1; ie<=3; ie++) { // loop over edges
    inode = ie;
    jnode = (ie<3)?ie+1: 1;
    
    if (((this->nodes.at(inode) == anode) && (this->nodes.at(jnode) == bnode)) ||
				((this->nodes.at(inode) == bnode) && (this->nodes.at(jnode) == anode))) {
      return this->irregular_nodes.at(ie);
    }
    if (((this->nodes.at(inode) == anode) && (this->nodes.at(4) == bnode)) ||
				((this->nodes.at(inode) == bnode) && (this->nodes.at(4) == anode))) {
      return this->irregular_nodes.at(ie);
    }
  }
  return 0;
}
*/

#ifdef __PARALLEL_MODE
int
Subdivision::RS_Element::giveIrregularSharedPartitions(int leIndex, int inode, int jnode, IntArray& partitions)
{
  int j=0;
  if ((mesh->giveNode(nodes.at(inode))->giveParallelMode() == DofManager_shared) &&
      (mesh->giveNode(nodes.at(jnode))->giveParallelMode() == DofManager_shared)) {
    // test if neigbor element exist (if yes and is local then the edge is not shared
    bool has_local_neighbor = false;
    int base_neighbor = this->neghbours_base_elements.at(leIndex);
    if (base_neighbor) {
      if (mesh->giveElement(base_neighbor)->giveParallelMode() == Element_local) has_local_neighbor = true;
    }
    if (!has_local_neighbor) { 
      // test if parent vertex nodes are shared by the same partition
      const IntArray *ipartitions = mesh->giveNode(nodes.at(inode))->givePartitions();
      const IntArray *jpartitions = mesh->giveNode(nodes.at(jnode))->givePartitions();
      // find common partitions 
      int i, ipsize = ipartitions->giveSize();
      for (i=1; i<=ipsize; i++) {
        if (jpartitions->contains (ipartitions->at(i)) &&
            (ipartitions->at(i) != this->mesh->giveSubdivision()->giveRank())) { //
          partitions.followedBy(ipartitions->at(i),2);
					j++;
				}
			}
			if(j==0)OOFEM_ERROR ("Subdivision::RS_Element::giveIrregularPartitions: topology error");
		}
	}
  return j;
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
     Also children array of receiver is update to contain children numbers */
#ifdef __PARALLEL_MODE
  if ((irregular_nodes.containsOnlyZeroes()) || (this->giveParallelMode() == Element_remote)) {
#else
  if (irregular_nodes.containsOnlyZeroes()) {
#endif
    // no subdivision of receiver required 
    children.resize(0);
  } else {
    int i, ind, nIrregulars = 0;
    int childNum, iedge, jedge, kedge, inode, jnode, knode;
    IntArray _nodes(3);
    Subdivision::RS_Triangle *child;
		Subdivision::RS_Element *ngb;
    
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

		// if there is neighbor (of "this") not designated for bisection
		// its neihgbor corresponding to "this" must be made negative to enforce update_neighbors
		for(i=1;i<=3;i++){
			if(this->neghbours_base_elements.at(i)){

#ifdef DEBUG_CHECK
				if(this->neghbours_base_elements.at(i) < 0){
					OOFEM_ERROR2("Subdivision::RS_Triangle::generate - negative neighbor of %d not expected", this->number);
				}
#endif

				ngb=mesh->giveElement(this->neghbours_base_elements.at(i));
#ifdef __PARALLEL_MODE
				if(ngb->giveParallelMode() == Element_local){
#endif
					if (!ngb->hasIrregulars()) {
						ind=ngb->giveNeighbors()->findFirstIndexOf(this->number);
						if(ind)ngb->setNeighbor(ind, -this->number);
					}
#ifdef __PARALLEL_MODE
				}
#endif
			}
		}
  }
}


void
Subdivision::RS_Tetra::generate()
{
  /* generates the children elements of already bisected element and adds them into mesh.
		 This is done recursively.
     Also children array of receiver is updated to contain children numbers */

	if(this->number < 0)return;    // skip tmp tetras in bisection loop 2 and more

#ifdef __PARALLEL_MODE
  if ((irregular_nodes.containsOnlyZeroes()) || (this->giveParallelMode() == Element_remote)) {
#else
  if (irregular_nodes.containsOnlyZeroes()) {
#endif
    // no subdivision of receiver required 
    children.resize(0);
  } else {
		int irregulars1 = 0, irregulars2 = 0;
    int childNum, iedge, jedge, kedge, iiedge, jjedge, kkedge, inode, jnode, knode, nnode, iface, jface, kface, nface;
		double maxlength, elength;
		int i, ind, leIndex1, leIndex2, ed1[3][3], ed2[3][3], l_nd[3], nd[3], nd1[3], nd2[3];
		int ed_1[3][3], ed_2[3][3];
    bool swap;
#ifdef __PARALLEL_MODE
		int g_nd[3];
#endif
    IntArray _nodes(4);
    Subdivision::RS_Tetra *child;
		Subdivision::RS_Element *ngb;
    
    // leIndex determines primary edge to be subdivided
		// it is identified either with iedge or iiedge
    /*

                               inode
                                 o

                              iiedge=leIndex(if>3)
                                 x

                   kedge x     nnode      x jedge
                                 o

                  jjedge x                x kkedge

                o                x                o
              jnode            iedge               knode
                              =leIndex(if<=3)

		 children keep all nodes in the same order except for that one which is replaced by irregular on leIndex edge !!!
    */
		if(leIndex <=3)
			iedge = leIndex;
		else
			iedge = (leIndex==6)?1:leIndex-2;

		jedge = (iedge <3)?iedge+1:1;
		kedge = (jedge <3)?jedge+1:1;

		inode = (iedge>1)?iedge-1:3;
		jnode = (inode<3)?inode+1:1;
		knode = (jnode<3)?jnode+1:1;
		nnode = 4;

		iiedge = inode+3;
		jjedge = jnode+3;
		kkedge = knode+3;

		iface = iedge+1;
		jface = jedge+1;
		kface = kedge+1;
		nface = 1;

#ifdef DEBUG_CHECK
		// check neighbours
		for(i=1;i<=4;i++){
			if(this->neghbours_base_elements.at(i) < 0){
				OOFEM_ERROR2("Subdivision::RS_Tetra::generate - negative neighbor of %d not expected", this->number);
			}
		}
#endif

		this->children.resize(2);

		if(leIndex <=3){
			// count number of irregulars on each child;
			// if there is one irregular only, leIndex is set
			if(irregular_nodes.at(iiedge)){
				irregulars1++;
				irregulars2++;
				leIndex1 = leIndex2 = 4;
			}
			if(irregular_nodes.at(kedge)){irregulars1++; leIndex1 = 1;}
			if(irregular_nodes.at(jjedge)){irregulars1++; leIndex1 = 5;}
			if(irregular_nodes.at(jedge)){irregulars2++; leIndex2 = 3;}
			if(irregular_nodes.at(kkedge)){irregulars2++; leIndex2 = 6;}

			if(irregulars1 > 1){
				nd1[0] = inode;
				nd1[1] = nnode;
				nd1[2] = jnode;

				// these are relevant edge indices of child1 with respect to nd1 array
				ed1[0][1] = ed1[1][0] = 4;
				ed1[0][2] = ed1[2][0] = 1;
				ed1[1][2] = ed1[2][1] = 5;

				// these are relevant edge indices of "this" with respect to nd1 array
				ed_1[0][1] = ed_1[1][0] = iiedge;
				ed_1[0][2] = ed_1[2][0] = kedge;
				ed_1[1][2] = ed_1[2][1] = jjedge;
			}
			if(irregulars2 > 1){
				nd2[0] = inode;
				nd2[1] = nnode;
				nd2[2] = knode;

				// these are relevant edge indices of child2 with respect to nd2 array
				ed2[0][1] = ed2[1][0] = 4;
				ed2[0][2] = ed2[2][0] = 3;
				ed2[1][2] = ed2[2][1] = 6;

				// these are relevant edge indices of "this" with respect to nd2 array
				ed_2[0][1] = ed_2[1][0] = iiedge;
				ed_2[0][2] = ed_2[2][0] = jedge;
				ed_2[1][2] = ed_2[2][1] = kkedge;
			}

			_nodes.at(1)=nodes.at(inode); _nodes.at(2)=nodes.at(jnode); 
			_nodes.at(3)=irregular_nodes.at(iedge); _nodes.at(4)=nodes.at(nnode); 
			childNum = mesh->giveNumberOfElements() + 1;
			child = new Subdivision::RS_Tetra (childNum, mesh, this->number, _nodes); // number, parent, coords
			mesh->addElement (childNum, child);
			children.at(1) = childNum;

			// set neigbour info
			if(irregulars1){
				child->setNeighbor(1, this->giveNeighbor(nface));
				child->setNeighbor(2, this->giveNeighbor(kface));
				child->setNeighbor(3, this->giveNeighbor(iface));

				child->setIrregular(1, irregular_nodes.at(kedge));
				child->setIrregular(5, irregular_nodes.at(jjedge));
				child->setIrregular(4, irregular_nodes.at(iiedge));
			}
			else{
				child->setNeighbor(1, -this->giveNeighbor(nface));
				child->setNeighbor(2, -this->giveNeighbor(kface));
				child->setNeighbor(3, -this->giveNeighbor(iface));
			}
			// neihgbor4 of child1 is changed to negative during subdivision (if any) of child2
			child->setNeighbor(4, childNum+1);

#ifdef DEBUG_INFO
			fprintf(stderr, "Child %d generated on parent %d (leIndex %d, nds %d %d %d %d, ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
							childNum, this->number, this->leIndex, _nodes.at(1), _nodes.at(2), _nodes.at(3), _nodes.at(4),
							child->giveNeighbor(1), child->giveNeighbor(2), child->giveNeighbor(3), child->giveNeighbor(4),
							child->giveIrregular(1), child->giveIrregular(2), child->giveIrregular(3), 
							child->giveIrregular(4), child->giveIrregular(5), child->giveIrregular(6));
#endif

			_nodes.at(1)=nodes.at(inode); _nodes.at(2)=irregular_nodes.at(iedge); 
			_nodes.at(3)=nodes.at(knode); _nodes.at(4)=nodes.at(nnode); 
			childNum = mesh->giveNumberOfElements() + 1;
			child = new Subdivision::RS_Tetra (childNum, mesh, this->number, _nodes); // number, parent, coords
			mesh->addElement (childNum, child);
			children.at(2) = childNum;
			// set neigbour info
			if(irregulars2){
				child->setNeighbor(1, this->giveNeighbor(nface));
				child->setNeighbor(3, this->giveNeighbor(iface));
				child->setNeighbor(4, this->giveNeighbor(jface));

				child->setIrregular(3, irregular_nodes.at(jedge));
				child->setIrregular(6, irregular_nodes.at(kkedge));
				child->setIrregular(4, irregular_nodes.at(iiedge));
			}
			else{
				child->setNeighbor(1, -this->giveNeighbor(nface));
				child->setNeighbor(3, -this->giveNeighbor(iface));
				child->setNeighbor(4, -this->giveNeighbor(jface));
			}
			// neihgbor2 of child2 is changed to negative during subdivision (if any) of child1
			child->setNeighbor(2, childNum-1);

#ifdef DEBUG_INFO
			fprintf(stderr, "Child %d generated on parent %d (leIndex %d, nds %d %d %d %d, ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
							childNum, this->number, this->leIndex, _nodes.at(1), _nodes.at(2), _nodes.at(3), _nodes.at(4),
							child->giveNeighbor(1), child->giveNeighbor(2), child->giveNeighbor(3), child->giveNeighbor(4),
							child->giveIrregular(1), child->giveIrregular(2), child->giveIrregular(3), 
							child->giveIrregular(4), child->giveIrregular(5), child->giveIrregular(6));
#endif

		} else {
			// count number of irregulars on each child;
			// if there is one irregular only, leIndex is set
			if(irregular_nodes.at(iedge)){
				irregulars1++;
				irregulars2++;
				leIndex1 = leIndex2 = 2;
			}
			if(irregular_nodes.at(jedge)){irregulars1++; leIndex1 = 3;}
			if(irregular_nodes.at(kedge)){irregulars1++; leIndex1 = 1;}
			if(irregular_nodes.at(jjedge)){irregulars2++; leIndex2 = 5;}
			if(irregular_nodes.at(kkedge)){irregulars2++; leIndex2 = 6;}

			if(irregulars1 > 1){
				nd1[0] = jnode;
				nd1[1] = knode;
				nd1[2] = inode;

				// these are relevant edge indices of child1 with respect to nd1 array
				ed1[0][1] = ed1[1][0] = 2;
				ed1[0][2] = ed1[2][0] = 1;
				ed1[1][2] = ed1[2][1] = 3;

				// these are relevant edge indices of "this" with respect to nd1 array
				ed_1[0][1] = ed_1[1][0] = iedge;
				ed_1[0][2] = ed_1[2][0] = kedge;
				ed_1[1][2] = ed_1[2][1] = jedge;
			}
			if(irregulars2 > 1){
				nd2[0] = jnode;
				nd2[1] = knode;
				nd2[2] = nnode;

				// these are relevant edge indices of child2 with respect to nd2 array
				ed2[0][1] = ed2[1][0] = 2;
				ed2[0][2] = ed2[2][0] = 5;
				ed2[1][2] = ed2[2][1] = 6;

				// these are relevant edge indices of "this" with respect to nd2 array
				ed_2[0][1] = ed_2[1][0] = iedge;
				ed_2[0][2] = ed_2[2][0] = jjedge;
				ed_2[1][2] = ed_2[2][1] = kkedge;
			}

			_nodes.at(1)=nodes.at(inode); _nodes.at(2)=nodes.at(jnode); 
			_nodes.at(3)=nodes.at(knode); _nodes.at(4)=irregular_nodes.at(iiedge); 
			childNum = mesh->giveNumberOfElements() + 1;
			child = new Subdivision::RS_Tetra (childNum, mesh, this->number, _nodes); // number, parent, coords
			mesh->addElement (childNum, child);
			children.at(1) = childNum;
			// set neigbour info
			if(irregulars1){
				child->setNeighbor(1, this->giveNeighbor(nface));
				child->setNeighbor(2, this->giveNeighbor(kface));
				child->setNeighbor(4, this->giveNeighbor(jface));

				child->setIrregular(1, irregular_nodes.at(kedge));
				child->setIrregular(3, irregular_nodes.at(jedge));
				child->setIrregular(2, irregular_nodes.at(iedge));
			}
			else{
				child->setNeighbor(1, -this->giveNeighbor(nface));
				child->setNeighbor(2, -this->giveNeighbor(kface));
				child->setNeighbor(4, -this->giveNeighbor(jface));
			}
			// neihgbor3 of child1 is changed to negative during subdivision (if any) of child2
			child->setNeighbor(3, childNum+1);

#ifdef DEBUG_INFO
			fprintf(stderr, "Child %d generated on parent %d (leIndex %d, nds %d %d %d %d, ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
							childNum, this->number, this->leIndex, _nodes.at(1), _nodes.at(2), _nodes.at(3), _nodes.at(4),
							child->giveNeighbor(1), child->giveNeighbor(2), child->giveNeighbor(3), child->giveNeighbor(4),
							child->giveIrregular(1), child->giveIrregular(2), child->giveIrregular(3), 
							child->giveIrregular(4), child->giveIrregular(5), child->giveIrregular(6));
#endif

			_nodes.at(1)=irregular_nodes.at(iiedge); _nodes.at(2)=nodes.at(jnode); 
			_nodes.at(3)=nodes.at(knode); _nodes.at(4)=nodes.at(nnode); 
			childNum = mesh->giveNumberOfElements() + 1;
			child = new Subdivision::RS_Tetra (childNum, mesh, this->number, _nodes); // number, parent, coords
			mesh->addElement (childNum, child);
			children.at(2) = childNum;
			// set neigbour info
			if(irregulars2){
				child->setNeighbor(2, this->giveNeighbor(kface));
				child->setNeighbor(3, this->giveNeighbor(iface));
				child->setNeighbor(4, this->giveNeighbor(jface));

				child->setIrregular(5, irregular_nodes.at(jjedge));
				child->setIrregular(6, irregular_nodes.at(kkedge));
				child->setIrregular(2, irregular_nodes.at(iedge));
			}
			else{
				child->setNeighbor(2, -this->giveNeighbor(kface));
				child->setNeighbor(3, -this->giveNeighbor(iface));
				child->setNeighbor(4, -this->giveNeighbor(jface));
			}
			// neihgbor1 of child2 is changed to negative during subdivision (if any) of child1
			child->setNeighbor(1, childNum-1);

#ifdef DEBUG_INFO
			fprintf(stderr, "Child %d generated on parent %d (leIndex %d, nds %d %d %d %d, ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
							childNum, this->number, this->leIndex, _nodes.at(1), _nodes.at(2), _nodes.at(3), _nodes.at(4),
							child->giveNeighbor(1), child->giveNeighbor(2), child->giveNeighbor(3), child->giveNeighbor(4),
							child->giveIrregular(1), child->giveIrregular(2), child->giveIrregular(3), 
							child->giveIrregular(4), child->giveIrregular(5), child->giveIrregular(6));
#endif

		}

		// IMPORTANT: ambiguity must be handled !!!

		// in order to achieve that shared tetra faces are subdivided in compatible way, always the same longest edge must be identified;
		// this may be problem if there are two or more edges of the "same" length
		// then the longest is dependent on the order in which edges are processed;
		// also in order to get always the same distance of two nodes, the nodes should be processes in the same order;
		
		// note this special handling is not needed if there is just one irregular in the child
		// or if lengths of edges with irregular nodes are not the same
		// (the later condition is not considered in this implementation)

		// however it is not clear whether one may rely that the same distance is always obtained for the same pair of nodes
		// (in proper order) especially if processed on different processors

		// recursively subdivide children
		if(irregulars1){

#ifdef DEBUG_CHECK
			ed1[0][0] = ed1[1][1] = ed1[2][2] = 0;
			ed_1[0][0] = ed_1[1][1] = ed_1[2][2] = 0;
#endif
		
			if(irregulars1 != 1){
				// sort node ids
				for(i = 0; i < 3; i++){
					nd[i] = i;
					l_nd[i] = nodes.at(nd1[i]);
#ifdef __PARALLEL_MODE
					g_nd[i] = mesh->giveNode(l_nd[i])->giveGlobalNumber();
#endif
				}
				swap = true;
				while(swap){
					swap = false;
					for(i = 0; i < 2; i++){
#ifdef __PARALLEL_MODE
						if(g_nd[i] > g_nd[i+1])
#else
						if(l_nd[i] > l_nd[i+1])
#endif
						{
							int tmp;
#ifdef __PARALLEL_MODE
							tmp = g_nd[i]; g_nd[i] = g_nd[i+1]; g_nd[i+1] = tmp;
#endif
							tmp = l_nd[i]; l_nd[i] = l_nd[i+1]; l_nd[i+1] = tmp;
							tmp = nd[i]; nd[i] = nd[i+1]; nd[i+1] = tmp;
							swap = true;
						}
					}
				}

				// get largest edge
				// process the edges from smallest node id
				maxlength = 0.0;
				if(irregular_nodes.at(ed_1[nd[0]][nd[1]])){
					elength = mesh->giveNode(l_nd[0])->giveCoordinates()->distance_square(*(mesh->giveNode(l_nd[1])->giveCoordinates()));
					if (elength>maxlength) {maxlength=elength; leIndex1 = ed1[nd[0]][nd[1]];}
				}
				if(irregular_nodes.at(ed_1[nd[0]][nd[2]])){
					elength = mesh->giveNode(l_nd[0])->giveCoordinates()->distance_square(*(mesh->giveNode(l_nd[2])->giveCoordinates()));
					if (elength>maxlength) {maxlength=elength; leIndex1 = ed1[nd[0]][nd[2]];}
				}
				if(irregular_nodes.at(ed_1[nd[1]][nd[2]])){
					elength = mesh->giveNode(l_nd[1])->giveCoordinates()->distance_square(*(mesh->giveNode(l_nd[2])->giveCoordinates()));
					if (elength>maxlength) {maxlength=elength; leIndex1 = ed1[nd[1]][nd[2]];}
				}
			}
				
			mesh->giveElement(children.at(1))->setLeIndex(leIndex1);
			mesh->giveElement(children.at(1))->generate();
			mesh->giveElement(children.at(1))->setNumber(-children.at(1));            // mark tmp tetra
		}
		if(irregulars2){

#ifdef DEBUG_CHECK
			ed2[0][0] = ed2[1][1] = ed2[2][2] = 0;
			ed_2[0][0] = ed_2[1][1] = ed_2[2][2] = 0;
#endif

			if(irregulars2 != 1){
				// sort node ids
				for(i = 0; i < 3; i++){
					nd[i] = i;
					l_nd[i] = nodes.at(nd2[i]);
#ifdef __PARALLEL_MODE
					g_nd[i] = mesh->giveNode(l_nd[i])->giveGlobalNumber();
#endif
				}
				swap = true;
				while(swap){
					swap = false;
					for(i = 0; i < 2; i++){
#ifdef __PARALLEL_MODE
						if(g_nd[i] > g_nd[i+1])
#else
						if(l_nd[i] > l_nd[i+1])
#endif
						{
							int tmp;
#ifdef __PARALLEL_MODE
							tmp = g_nd[i]; g_nd[i] = g_nd[i+1]; g_nd[i+1] = tmp;
#endif
							tmp = l_nd[i]; l_nd[i] = l_nd[i+1]; l_nd[i+1] = tmp;
							tmp = nd[i]; nd[i] = nd[i+1]; nd[i+1] = tmp;
							swap = true;
						}
					}
				}

				// get largest edge
				// process the edges from smallest node id
				maxlength = 0.0;
				if(irregular_nodes.at(ed_2[nd[0]][nd[1]])){
					elength = mesh->giveNode(l_nd[0])->giveCoordinates()->distance_square(*(mesh->giveNode(l_nd[1])->giveCoordinates()));
					if (elength>maxlength) {maxlength=elength; leIndex2 = ed2[nd[0]][nd[1]];}
				}
				if(irregular_nodes.at(ed_2[nd[0]][nd[2]])){
					elength = mesh->giveNode(l_nd[0])->giveCoordinates()->distance_square(*(mesh->giveNode(l_nd[2])->giveCoordinates()));
					if (elength>maxlength) {maxlength=elength; leIndex2 = ed2[nd[0]][nd[2]];}
				}
				if(irregular_nodes.at(ed_2[nd[1]][nd[2]])){
					elength = mesh->giveNode(l_nd[1])->giveCoordinates()->distance_square(*(mesh->giveNode(l_nd[2])->giveCoordinates()));
					if (elength>maxlength) {maxlength=elength; leIndex2 = ed2[nd[1]][nd[2]];}
				}
			}

			mesh->giveElement(children.at(2))->setLeIndex(leIndex2);
			mesh->giveElement(children.at(2))->generate();
			mesh->giveElement(children.at(2))->setNumber(-children.at(2));            // mark tmp tetra
		}

		// if there is neighbor (of "this") not designated for bisection
		// its neihgbor corresponding to "this" must be made negative to enforce update_neighbors
		for(i=1;i<=4;i++){
			if(this->neghbours_base_elements.at(i)){

#ifdef DEBUG_CHECK
				if(this->neghbours_base_elements.at(i) < 0){
					OOFEM_ERROR2("Subdivision::RS_Tetra::generate - negative neighbor of %d not expected", this->number);
				}
#endif

				ngb=mesh->giveElement(this->neghbours_base_elements.at(i));
#ifdef __PARALLEL_MODE
				if(ngb->giveParallelMode() == Element_local){
#endif
					if (!ngb->hasIrregulars()) {
						ind=ngb->giveNeighbors()->findFirstIndexOf(this->number);
						if(ind)ngb->setNeighbor(ind, -this->number);
					}
#ifdef __PARALLEL_MODE
				}
#endif
			}
		}
  }
}



void
Subdivision::RS_Triangle::update_neighbours()
{
  bool found;
  int i, iside,parentNeighbor;
  IntArray ca;

  /*
    neghbours_base_elements pre-set in generate using this rule:
    value > 0 .. real neighbour known, value equals to real neighbor on new mesh
    value = 0 .. boundary
    value < 0 .. neighbour set to parent element, actual neighbour has to be determined from parent children 
                 or if no child exist then to parent new counterpart.
   */
#ifdef __PARALLEL_MODE
  if (this->giveParallelMode() == Element_remote) return;
#endif

  //this->neghbours_base_elements.resize(3);
  for (iside=1; iside<=3; iside++) {
    if (this->neghbours_base_elements.at(iside) < 0) {
      parentNeighbor = -this->neghbours_base_elements.at(iside);
      // test if parentNeighbor has been subdivided
#ifdef __PARALLEL_MODE
      if ((mesh->giveElement(parentNeighbor)->giveParallelMode() == Element_local) &&
          (mesh->giveElement(parentNeighbor)->hasIrregulars())) {      // HUHU replace by !isTerminal
#else
			if (mesh->giveElement(parentNeighbor)->hasIrregulars()) {        // HUHU replace by !isTerminal
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


void
Subdivision::RS_Tetra::update_neighbours()
{
  bool found;
  int i, j, k, iside,parentNeighbor;
  IntArray ca1, ca2, ca3;

	if(this->number < 0)return;    // skip tmp tetras in bisection pass 2 and more

  /*
    neghbours_base_elements pre-set in generate using this rule:
    value > 0 .. real neighbour known, value equals to real neighbor on new mesh
    value = 0 .. boundary
    value < 0 .. neighbour set to parent element, actual neighbour has to be determined from terminal children of parent
                 or if no child exist then to parent new counterpart.
   */

  //this->neghbours_base_elements.resize(4);

#ifdef __PARALLEL_MODE
  if (this->giveParallelMode() == Element_remote) return;
#endif

  for (iside=1; iside<=4; iside++) {
    if (this->neghbours_base_elements.at(iside) < 0) {
      parentNeighbor = -this->neghbours_base_elements.at(iside);
      // test if parentNeighbor has been subdivided
#ifdef __PARALLEL_MODE
      if ((mesh->giveElement(parentNeighbor)->giveParallelMode() == Element_local) &&
					(mesh->giveElement(parentNeighbor)->hasIrregulars())) {      // HUHU replace by !isTerminal
#else
			if (mesh->giveElement(parentNeighbor)->hasIrregulars()) {      // HUHU replace by !isTerminal
#endif
				// old neighbour bisected -> we loop over its terminal children to find an appropriate new neighbor
				// there may be at maximum 3 levels of parent tetra subdivision
				mesh->giveElement(parentNeighbor)->giveChildren(ca1);
				found = false;
				for (i=1; i<=ca1.giveSize(); i++) {
					if(mesh->giveElement(ca1.at(i))->isTerminal()){
						if (this->isNeighborOf(mesh->giveElement(ca1.at(i)))) {
							this->neghbours_base_elements.at(iside) = ca1.at(i);
							found = true;
							break;
						}
					} else{
						mesh->giveElement(ca1.at(i))->giveChildren(ca2);
						for (j=1; j<=ca2.giveSize(); j++) {
							if(mesh->giveElement(ca2.at(j))->isTerminal()){
								if (this->isNeighborOf(mesh->giveElement(ca2.at(j)))) {
									this->neghbours_base_elements.at(iside) = ca2.at(j);
									found = true;
									break;
								}
							} else{
								mesh->giveElement(ca2.at(j))->giveChildren(ca3);
								for (k=1; k<=ca3.giveSize(); k++) {
									if(mesh->giveElement(ca3.at(k))->isTerminal()){
										if (this->isNeighborOf(mesh->giveElement(ca3.at(k)))) {
											this->neghbours_base_elements.at(iside) = ca3.at(k);
											found = true;
											break;
										}
#ifdef DEBUG_CHECK
									} else {
										OOFEM_ERROR2 ("Subdivision::RS_Tetra::update_neighbours: too many subdivision levels for element %d", this->number);
									}
#endif
								}
							}
							if (found) break;
						}
					}
					if (found) break;
				}
				if (!found) {
					OOFEM_ERROR4 ("Subdivision::RS_Tetra::update_neighbours failed for element %d (side %d, elem %d)", 
												this->number, iside, -this->neghbours_base_elements.at(iside));
				}
			} else {
				// parent neighbour remains actual neighbour
				this->neghbours_base_elements.at(iside) = parentNeighbor;
			}
		} else if (this->neghbours_base_elements.at(iside) > 0) {
      // neighbor element already set
    }
  } // end loop over element side faces

#ifdef DEBUG_CHECK
	// check updated neighbors
	IntArray snodes1, snodes2;
	RS_Element *ngb;
	for(i=1;i<=4;i++){
		if(this->neghbours_base_elements.at(i)){
			if(this->neghbours_base_elements.at(i) > 0){
				this->giveSideNodes(i, snodes1);
				ngb=mesh->giveElement(this->neghbours_base_elements.at(i));
				if(ngb->giveNumber() < 0){
					OOFEM_ERROR3 ("Subdivision::RS_Tetra::update_neighbours neighbour %d of %d is temporary", 
												this->neghbours_base_elements.at(i), this->number);
				}
				if(!ngb->isTerminal()){
					OOFEM_ERROR3 ("Subdivision::RS_Tetra::update_neighbours neighbour %d of %d not terminal", 
												this->neghbours_base_elements.at(i), this->number);
				}
				if(ngb->giveNumber() < this->number){
					// ngb has been already processed
					j=ngb->giveNeighbors()->findFirstIndexOf(this->number);
					if(j){
						ngb->giveSideNodes(j, snodes2);
						for(k=1;k<=3;k++){
							if(snodes2.findFirstIndexOf(snodes1.at(k)))continue;
							OOFEM_ERROR4 ("Subdivision::RS_Tetra::update_neighbours: element %d is not neighbor (%d) of element %d", 
														this->number, i, this->neghbours_base_elements.at(i));
						}
					} else {
						OOFEM_ERROR4 ("Subdivision::RS_Tetra::update_neighbours: element %d is not neighbor (%d) of element %d", 
													this->number, i, this->neghbours_base_elements.at(i));
					}
				}
				else{
					// ngb has not been yet processed ==> I cannot check particular side of ngb
					for(k=1;k<=3;k++){
						if(ngb->giveNodes()->findFirstIndexOf(snodes1.at(k)))continue;
						OOFEM_ERROR3 ("Subdivision::RS_Tetra::update_neighbours: element %d is not neighbor of element %d", 
													this->number, this->neghbours_base_elements.at(i));
					}
				}
			}
			else{
				OOFEM_ERROR3 ("Subdivision::RS_Tetra::update_neighbours negative neighbour %d of %d not expected", 
											this->neghbours_base_elements.at(i), this->number);
			}
		}
	}
#endif

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


bool
Subdivision::RS_Tetra::isNeighborOf (Subdivision::RS_Element* elem)
{
  // Simplified implementation, considering only one element type - tetras
  int i, _c = 0;
  for (i=1; i<=4; i++) {
    _c +=  elem->containsNode(this->nodes.at(i));
  }
  return (_c==3);
}


void
Subdivision::RS_Triangle::giveSideNodes (int iside, IntArray& snodes)
{
	int inode,jnode;

	inode=iside;
	jnode=(iside<3)?iside+1:1;

	snodes.resize(2);
	snodes.at(1)=nodes.at(inode);
	snodes.at(2)=nodes.at(jnode);
}
	

void
Subdivision::RS_Tetra::giveSideNodes (int iside, IntArray& snodes)
{
	int inode,jnode,knode;

	if(iside==1){
		inode=1;
		jnode=2;
		knode=3;
	} else {
		inode=iside-1;
		jnode=(iside<4)?iside:1;
		knode=4;
	}

	snodes.resize(3);
	snodes.at(1)=nodes.at(inode);
	snodes.at(2)=nodes.at(jnode);
	snodes.at(3)=nodes.at(knode);
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


double
Subdivision::RS_Tetra::giveDensity()
{
  double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
	double dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3;
	double vol, d;
  x1 = mesh->giveNode(nodes.at(1))->giveCoordinate(1);
  x2 = mesh->giveNode(nodes.at(2))->giveCoordinate(1);
  x3 = mesh->giveNode(nodes.at(3))->giveCoordinate(1);
  x4 = mesh->giveNode(nodes.at(4))->giveCoordinate(1);
  
  y1 = mesh->giveNode(nodes.at(1))->giveCoordinate(2);
  y2 = mesh->giveNode(nodes.at(2))->giveCoordinate(2);
  y3 = mesh->giveNode(nodes.at(3))->giveCoordinate(2);
  y4 = mesh->giveNode(nodes.at(4))->giveCoordinate(2);
  
  z1 = mesh->giveNode(nodes.at(1))->giveCoordinate(3);
  z2 = mesh->giveNode(nodes.at(2))->giveCoordinate(3);
  z3 = mesh->giveNode(nodes.at(3))->giveCoordinate(3);
  z4 = mesh->giveNode(nodes.at(4))->giveCoordinate(3);
  
	dx1=x2-x1; dy1=y2-y1; dz1=z2-z1;
	dx2=x3-x1; dy2=y3-y1; dz2=z3-z1;
	dx3=x4-x1; dy3=y4-y1; dz3=z4-z1;

	vol=(dx3*(dy1*dz2-dz1*dy2)+dy3*(dz1*dx2-dx1*dz2)+dz3*(dx1*dy2-dy1*dx2)) / 6.0;
	d=exp(log(vol)/3.0);

  return d;
}


void
Subdivision::RS_Triangle::importConnectivities(ConnectivityTable* ct)
{
  IntArray me(1), conn;
  int iside, inode, jnode, el, i;
  me.at(1) = this->number;

  neghbours_base_elements.resize(3);
  neghbours_base_elements.zero();             // initialized to have no neighbour
  ct->giveElementNeighbourList(conn, me);

	// HUHU is conn ordered? if yes, function can be optimized

  for (iside=1; iside<=3; iside++) {
    inode=nodes.at(iside);
    jnode=nodes.at((iside==3)?1:iside+1);

    // select right neighbour
    for (i=1; i<=conn.giveSize(); i++) {
      el = conn.at(i);
      if (el != this->number) {
        if (mesh->giveElement(el)->containsNode(inode) &&
            mesh->giveElement(el)->containsNode(jnode)) {
          neghbours_base_elements.at(iside) = el;
					break;
        }
      }
    }
  }
}


void
Subdivision::RS_Tetra::importConnectivities(ConnectivityTable* ct)
{
  IntArray me(1), conn;
  int iside, inode, jnode, knode, el, i;
  me.at(1) = this->number;

  neghbours_base_elements.resize(4);
  neghbours_base_elements.zero();             // initialized to have no neighbour
  ct->giveElementNeighbourList(conn, me);

	// HUHU is conn ordered? if yes, function can be optimized

	inode=nodes.at(1);
	jnode=nodes.at(2);
	knode=nodes.at(3);

    // select right base neighbour
	for (i=1; i<=conn.giveSize(); i++) {
		el = conn.at(i);
		if (el != this->number) {
			if (mesh->giveElement(el)->containsNode(inode) &&
					mesh->giveElement(el)->containsNode(jnode) &&
					mesh->giveElement(el)->containsNode(knode)) {
				neghbours_base_elements.at(1) = el;
				break;
			}
		}
	}

	knode=nodes.at(4);
  for (iside=1; iside<=3; iside++) {
    inode=nodes.at(iside);
    jnode=nodes.at((iside==3)?1:iside+1);

    // select right side neighbour
    for (i=1; i<=conn.giveSize(); i++) {
      el = conn.at(i);
      if (el != this->number) {
				if (mesh->giveElement(el)->containsNode(inode) &&
						mesh->giveElement(el)->containsNode(jnode) &&
						mesh->giveElement(el)->containsNode(knode)) {
					neghbours_base_elements.at(iside+1) = el;
					break;
				}
      }
    }
  }

#ifdef DEBUG_INFO
	fprintf(stderr, "Elem %d connectivity imported (nds %d %d %d %d ngbs %d %d %d %d)\n", this->number, 
					nodes.at(1), nodes.at(2), nodes.at(3), nodes.at(4),
					neghbours_base_elements.at(1), neghbours_base_elements.at(2), 
					neghbours_base_elements.at(3), neghbours_base_elements.at(4));
#endif
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


void
Subdivision::RS_Tetra::drawGeometry()
{
    WCRec p [ 4 ];
    GraphicObj *go;

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    //EASValsSetColor( gc.getElementColor() );
    //EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(TRUE);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    //    EASValsSetShrink(0.8);
    p [ 0 ].x = ( FPNum ) mesh->giveNode(nodes.at(1))->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) mesh->giveNode(nodes.at(1))->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) mesh->giveNode(nodes.at(1))->giveCoordinate(3);
    p [ 1 ].x = ( FPNum ) mesh->giveNode(nodes.at(2))->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) mesh->giveNode(nodes.at(2))->giveCoordinate(2);
    p [ 1 ].z = ( FPNum ) mesh->giveNode(nodes.at(2))->giveCoordinate(3);
    p [ 2 ].x = ( FPNum ) mesh->giveNode(nodes.at(3))->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) mesh->giveNode(nodes.at(3))->giveCoordinate(2);
    p [ 2 ].z = ( FPNum ) mesh->giveNode(nodes.at(3))->giveCoordinate(3);
    p [ 3 ].x = ( FPNum ) mesh->giveNode(nodes.at(4))->giveCoordinate(1);
    p [ 3 ].y = ( FPNum ) mesh->giveNode(nodes.at(4))->giveCoordinate(2);
    p [ 3 ].z = ( FPNum ) mesh->giveNode(nodes.at(4))->giveCoordinate(3);

    go =  CreateTetra(p);
    //EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | SHRINK_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGWithMaskChangeAttributes(WIDTH_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}
#endif


void
Subdivision::bisectMesh () {
  
  int ie, nelems = mesh->giveNumberOfElements(), nelems_old=0;
  int in, nnodes=mesh->giveNumberOfNodes();
  double iedensity, rdensity,rate,min_rate;
	int repeat=1, loop=1;
  //std::queue<int> subdivqueue;
#ifdef __PARALLEL_MODE
  int myrank=this->giveRank();
  int problem_size=this->giveNumberOfProcesses();
#endif
  // first select all candidates for local bisection based on required mesh density

	// repeat bisection until no new element is created
	while(repeat){
    OOFEM_LOG_INFO ("Subdivision::bisectMesh: entering bisection loop %d\n", loop++);
    repeat = 0;
		min_rate = 1.0;
		for (ie=nelems_old+1; ie<=nelems; ie++) {
#ifdef __PARALLEL_MODE
			if(nelems_old == 0){
				if (mesh->giveElement(ie)->giveParallelMode() == Element_remote) continue;    // skip remote elems in first bisection pass
			}
#endif
			if (!mesh->giveElement(ie)->isTerminal()) continue;

			iedensity = mesh->giveElement(ie)->giveDensity();
			rdensity = mesh->giveElement(ie)->giveRequiredDensity();
			
			if (rdensity < iedensity) {
				rate = rdensity / iedensity;
				if(rate < min_rate)min_rate = rate;

				subdivqueue.push(ie);
        repeat = 1;
				/*
					#ifdef __PARALLEL_MODE
					OOFEM_LOG_INFO("[%d] Subdivision: scheduling element %d[%d] for bisection, dens=%lf rdens=%lf\n", myrank, ie, mesh->giveElement(ie)->giveGlobalNumber(), iedensity, rdensity);
					#else
					OOFEM_LOG_INFO("Subdivision: scheduling element %d for bisection, dens=%lf rdens=%lf\n", ie, iedensity, rdensity);
					#endif
				*/
			}
		}
		fprintf(stderr, "min refinement %f\n", min_rate);
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
#endif
    
#ifdef __PARALLEL_MODE
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
      if (this->isNodeLocalIrregular(mesh->giveNode(in), myrank) && (mesh->giveNode(in)->giveGlobalNumber() == 0)) 
        localIrregulars ++;
    }
    
#ifdef __VERBOSE_PARALLEL
    OOFEM_LOG_INFO ("[%d] Subdivision::bisectMesh: number of new local irregulars is %d\n", myrank, localIrregulars);
#endif
    int irank, localOffset, gnum, partitionsIrregulars[problem_size];
    // gather number of local irregulars from all partitions
    MPI_Allgather (&localIrregulars, 1, MPI_INT, partitionsIrregulars, 1, MPI_INT, MPI_COMM_WORLD);
    // compute local offset
    for (localOffset=0, irank=0; irank <myrank; irank++) localOffset += partitionsIrregulars[irank];
    // start to assign global numbers to local irregulars
    int availGlobNum = maxglobalnumber+localOffset;
    for (in=1; in<=nnodes; in++) {
      if (this->isNodeLocalIrregular(mesh->giveNode(in), myrank) && (mesh->giveNode(in)->giveGlobalNumber() == 0))
        // set negative glognum to mark newly assigned nodes with glognum to participate in mutual globnum data exchange 
        mesh->giveNode(in)->setGlobalNumber(-(++availGlobNum));
    }
    // finally, communicate global numbers assigned to shared irregulars
    this->assignGlobalNumbersToSharedIrregulars();
    for (in=1; in<=nnodes; in++) {
      gnum = mesh->giveNode(in)->giveGlobalNumber();
      if (gnum < 0) {
        // turn all globnums to positive
        mesh->giveNode(in)->setGlobalNumber(abs(mesh->giveNode(in)->giveGlobalNumber()));
      } else if (gnum == 0) {
        OOFEM_ERROR3 ("[%d] Subdivision::bisectMesh: zero globnum identified on node %d", myrank, in);
      }
    }
    
#endif


		// symbolic bisection is finished. Now we need to create new mesh. Also there is a need to update 
		// element connectivities
		for (ie=1; ie<=nelems; ie++) { 
			if (!mesh->giveElement(ie)->isTerminal()) continue;
			mesh->giveElement(ie)->generate ();                      
		}
		
		nelems_old=nelems;
		nelems = mesh->giveNumberOfElements();
		for (ie=1; ie<=nelems; ie++) {
			if (!mesh->giveElement(ie)->isTerminal()) continue;
			mesh->giveElement(ie)->update_neighbours();
		}

#ifdef __PARALLEL_MODE
		int global_repeat=0;

		// determine whether any partition needs additional bisection pass
		MPI_Allreduce (&repeat, &global_repeat, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);  
		repeat=global_repeat;
#endif
    // debug only
    // repeat = 0;
	}


}


void
Subdivision::smoothMesh ()
{
	int nnodes, nelems, i, j, in, ie;
	int pos, number, reg, nd, cycles=6;
	IntArray snodes;
  FloatArray *coords;
	RS_Element *elem;
	//bool fixed;
	IntArray node_num_elems, node_con_elems;
	IntArray node_num_nodes, node_con_nodes;

#ifdef SMOOTH_DEBUG
	char buffer[1024];
	int len;
#endif
  OOFEM_LOG_INFO("Subdivision::smoothMesh\n");
	nnodes = mesh->giveNumberOfNodes();
	nelems = mesh->giveNumberOfElements();

	// count number of elements incident to nodes
	node_num_elems.resize(nnodes+1);
	node_num_elems.zero();

	for(ie=1;ie<=nelems;ie++){
		elem=mesh->giveElement(ie);
		if (!elem->isTerminal()) continue;
#ifdef __PARALLEL_MODE
    if(elem->giveParallelMode() == Element_remote)continue;
#endif

		for(i=1;i<=elem->giveNodes()->giveSize();i++)node_num_elems.at(elem->giveNode(i))++;
	}

	// recalculate the number of elements to current addresses
	pos=1;
	for(in=1;in<=nnodes;in++){
		number = node_num_elems.at(in);
		node_num_elems.at(in) = pos;
		pos += number;
	}
		
	node_num_elems.at(nnodes+1)=pos;
	node_con_elems.resize(pos);
		
	// store element numbers incident to nodes
	for(ie=1;ie<=nelems;ie++){
		elem=mesh->giveElement(ie);
		if (!elem->isTerminal()) continue;
#ifdef __PARALLEL_MODE
    if(elem->giveParallelMode() == Element_remote)continue;
#endif

		for(i=1;i<=elem->giveNodes()->giveSize();i++)node_con_elems.at(node_num_elems.at(elem->giveNode(i))++)=ie;
	}
		
	// recalculate the addresses to address of the first element
	pos=1;
	for(in=1;in<=nnodes;in++){
		number = node_num_elems.at(in) - pos;
		node_num_elems.at(in) = pos;
		pos += number;
	}

#ifdef SMOOTH_DEBUG
#ifdef __PARALLEL_MODE
	OOFEM_LOG_INFO("[%d] Subdivision::smoothMesh: connectivity node_id: elem_ids\n", this->giveRank());
#else
	OOFEM_LOG_INFO("Subdivision::smoothMesh: connectivity node_id: elem_ids\n");
#endif
	for(in=1;in<=nodes;in++){
		len=0;
		for(i=node_num_elems.at(in);i<node_num_elems.at(in+1);i++){
			sprintf(buffer+len, " %d", node_con_elems.at(i));
			len=strlen(buffer);
			if(len>=1024)OOFEM_ERROR ("Subdivision::smoothMesh: too small buffer");
		}
		OOFEM_LOG_INFO ("%d:%s\n", in, buffer);
	}
#endif

	// count number of nodes incident to nodes
	node_num_nodes.resize(nnodes+1);
	node_num_nodes.zero();

	for(in=1;in<=nnodes;in++){
		for(i=node_num_elems.at(in);i<node_num_elems.at(in+1);i++){
			elem=mesh->giveElement(node_con_elems.at(i));
			for(j=1;j<=elem->giveNodes()->giveSize();j++){
				nd=elem->giveNode(j);
				if(nd != in && mesh->giveNode(nd)->giveNumber() > 0){
					mesh->giveNode(nd)->setNumber(-nd);      // mark connected node
					node_num_nodes.at(in)++;
				}
			}
		}
			
		// unmark connected nodes
		for(i=node_num_elems.at(in);i<node_num_elems.at(in+1);i++){
			elem=mesh->giveElement(node_con_elems.at(i));
			for(j=1;j<=elem->giveNodes()->giveSize();j++){
				nd=elem->giveNode(j);
				mesh->giveNode(nd)->setNumber(nd);
			}
		}
	}
		
	// recalculate the number of nodes to current addresses
	pos=1;
	for(in=1;in<=nnodes;in++){
		number = node_num_nodes.at(in);
		node_num_nodes.at(in) = pos;
		pos += number;
	}
			
	node_num_nodes.at(nnodes+1) = pos;
	node_con_nodes.resize(pos);
			
	// store nodes incident to nodes
	for(in=1;in<=nnodes;in++){
		for(i=node_num_elems.at(in);i<node_num_elems.at(in+1);i++){
			elem=mesh->giveElement(node_con_elems.at(i));
			for(j=1;j<=elem->giveNodes()->giveSize();j++){
				nd=elem->giveNode(j);
				if(nd != in && mesh->giveNode(nd)->giveNumber() > 0){
					mesh->giveNode(nd)->setNumber(-nd);      // mark connected node
					node_con_nodes.at(node_num_nodes.at(in)++)=nd;
				}
			}
		}
			
		// unmark connected nodes
		for(i=node_num_elems.at(in);i<node_num_elems.at(in+1);i++){
			elem=mesh->giveElement(node_con_elems.at(i));
			for(j=1;j<=elem->giveNodes()->giveSize();j++){
				nd=elem->giveNode(j);
				mesh->giveNode(nd)->setNumber(nd);
			}
		}
	}
			
	// recalculate the addresses to address of the first node
	pos=1;
	for(in=1;in<=nnodes;in++){
		number = node_num_nodes.at(in) - pos;
		node_num_nodes.at(in) = pos;
		pos += number;
	}

#ifdef SMOOTH_DEBUG
#ifdef __PARALLEL_MODE
	OOFEM_LOG_INFO("[%d] Subdivision::smoothMesh: connectivity node_id: node_ids\n", this->giveRank());
#else
	OOFEM_LOG_INFO("Subdivision::smoothMesh: connectivity node_id: node_ids\n");
#endif
	for(in=1;in<=nodes;in++){
		len=0;
		for(i=node_num_nodes.at(in);i<node_num_nodes.at(in+1);i++){
			sprintf(buffer+len, " %d", node_con_nodes.at(i));
			len=strlen(buffer);
			if(len>=1024)OOFEM_ERROR ("Subdivision::smoothMesh: too small buffer");
		}
		OOFEM_LOG_INFO ("%d:%s\n", in, buffer);
	}
#endif

	/* this is not needed any more, because boundary flag is used !!! */  // HUHU

	// identify fixed nodes which should not be subjected to smoothing;
	// these are nodes on boundary (geometrical or material)
	bool fixed;
	for(in=1;in<=nnodes;in++){
		fixed=false;
		elem=mesh->giveElement(node_con_elems.at(node_num_elems.at(in)));
		// check for geometrical boundary
		for(j=1;j<=elem->giveNeighbors()->giveSize();j++){
			if(elem->giveNeighbor(j)==0){
				elem->giveSideNodes(j,snodes);
				if(snodes.findFirstIndexOf(in)!=0){
					mesh->giveNode(in)->setNumber(-in);     // mark fixed node
					fixed=true;
					break;
				}
			}
		}
		if(!fixed){
			reg=domain->giveElement(elem->giveTopParent())->giveRegionNumber();
			for(i=node_num_elems.at(in)+1;i<node_num_elems.at(in+1);i++){
				elem=mesh->giveElement(node_con_elems.at(i));
				// check for material boundary
				if(domain->giveElement(elem->giveTopParent())->giveRegionNumber() != reg){
					mesh->giveNode(in)->setNumber(-in);      //mark fixed node
					fixed=true;
					break;
				}
				// check for geometrical boundary
				for(j=1;j<=elem->giveNeighbors()->giveSize();j++){
					if(elem->giveNeighbor(j)==0){
						elem->giveSideNodes(j,snodes);
						if(snodes.findFirstIndexOf(in)!=0){
							mesh->giveNode(in)->setNumber(-in);      //mark fixed node
							fixed=true;
							break;
						}
					}
				}
				if(fixed)break;
			}
		}
	}

	// quick hack
	int jn;
	IntArray orderedNodes;
  Subdivision::RS_CompareNodePositions cmp (mesh);
	orderedNodes.resize(nnodes);
	for(jn=1;jn<=nnodes;jn++)orderedNodes.at(jn)=jn;
	sort (orderedNodes, cmp);

	while(cycles--){
		for(jn=1;jn<=nnodes;jn++){
			//for(in=1;in<=nnodes;in++){
			in=orderedNodes.at(jn);
#ifdef __PARALLEL_MODE
      if((mesh->giveNode(in)->giveParallelMode() == DofManager_shared) || 
				 (mesh->giveNode(in)->giveParallelMode() == DofManager_null))continue;    // skip shared and remote
#endif
			if(mesh->giveNode(in)->giveNumber() < 0)continue;         // skip fixed node
			if(mesh->giveNode(in)->isBoundary())continue;             // skip fixed node
			coords = mesh->giveNode(in)->giveCoordinates();

#ifdef DEBUG_CHECK
      if (coords) {
				int count=0;
        coords -> zero();
        for(i=node_num_nodes.at(in);i<node_num_nodes.at(in+1);i++) {
          if (mesh->giveNode(node_con_nodes.at(i))) {
            if (mesh->giveNode(node_con_nodes.at(i))->giveCoordinates()) {
              coords->add(mesh->giveNode(node_con_nodes.at(i))->giveCoordinates());
							count++;
            } else {
              OOFEM_ERROR2 ("Smooth: node %d without coordinates", in);
            }
          } else {
            OOFEM_ERROR2 ("Smooth: undefined node %d", in);
          }
        }
				if(!count){
					OOFEM_ERROR2 ("Smooth: node %d without connectivity", in);
				}
        coords->times(1.0/(node_num_nodes.at(in+1)-node_num_nodes.at(in)));
      } else {
        OOFEM_ERROR2 ("Smooth: node %d without coordinates", in);
      }
#else
			coords -> zero();
			for(i=node_num_nodes.at(in);i<node_num_nodes.at(in+1);i++) {
				coords->add(mesh->giveNode(node_con_nodes.at(i))->giveCoordinates());
			}
			coords->times(1.0/(node_num_nodes.at(in+1)-node_num_nodes.at(in)));
#endif

		}
	}
			
	// unmark fixed nodes
	for(in=1;in<=nnodes;in++){
		if(mesh->giveNode(in)->giveNumber() < 0)mesh->giveNode(in)->setNumber(in);
	}
}


MesherInterface::returnCode
Subdivision :: createMesh(TimeStep *stepN, int domainNumber, int domainSerNum, Domain** dNew)
{
  // import data from old domain
  int i, j, parent, nnodes = domain->giveNumberOfDofManagers(), nelems=domain->giveNumberOfElements();
  int inode, idof, ielem, ndofs, num;
  IntArray enodes;
  Subdivision::RS_Node *_node;
  Subdivision::RS_Element *_element;
  IRResultType result;                            // Required by IR_GIVE_FIELD macro

  oofem_timeval st, dt;
  ::getUtime(st);

  if (this->mesh) delete mesh;
  mesh = new Subdivision::RS_Mesh(this);

  // import nodes
  for (i=1; i<=nnodes; i++) {
    _node = new Subdivision::RS_Node(i, i, *(domain->giveNode(i)->giveCoordinates()), 
                                     domain->giveErrorEstimator()->giveRemeshingCrit()->giveRequiredDofManDensity(i, stepN),
                                     domain->giveNode(i)->isBoundary());
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
		} else if (domain->giveElement(i)->giveGeometryType() == EGT_tetra_1) {
      enodes.resize(4);
      for (j=1; j<=4; j++) enodes.at(j)=domain->giveElement(i)->giveDofManagerNumber(j);
      _element = new Subdivision::RS_Tetra (i,mesh,i,enodes);
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
	// smooth mesh
  if (smoothingFlag) this->smoothMesh();
#ifdef __OOFEG
  nelems=mesh->giveNumberOfElements();
  for (i=1; i<=nelems; i++) {
    if (mesh->giveElement(i)->isTerminal()) mesh->giveElement(i)->drawGeometry();
    //ESIEventLoop (YES, "Subdivision Bisection Finished; Press Ctrl-p to continue");
  }
  ESIEventLoop (YES, "Subdivision Bisection & smoothing Finished; Press Ctrl-p to continue");
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
    // set node coordinates
    ((Node*) node)->setCoordinates(*mesh->giveNode(inode)->giveCoordinates());
		node->setBoundaryFlag (mesh->giveNode(inode)->isBoundary());
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
    parent = mesh->giveElement(ielem)->giveTopParent();
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
  OOFEM_LOG_INFO ("[%d] Subdivision :: exchangeSharedIrregulars: localSharedIrregularQueueEmpty %d\n",
                  this->giveRank(), localSharedIrregularQueueEmpty);
#endif
  MPI_Allreduce (&localSharedIrregularQueueEmpty, &globalIrregularQueueEmpty, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD); 
#ifdef __VERBOSE_PARALLEL
  OOFEM_LOG_INFO ("[%d] Subdivision :: exchangeSharedIrregulars:  globalIrregularQueueEmpty %d\n",
                  this->giveRank(), globalIrregularQueueEmpty);
#endif
  if (globalIrregularQueueEmpty) {
    return true; 
  } else {
#ifdef __VERBOSE_PARALLEL
    OOFEM_LOG_INFO ("[%d] Subdivision :: exchangeSharedIrregulars: started\n", this->giveRank());
#endif
    
    // there are some shared irregulars -> data exchange
    CommunicatorBuff cb(this->giveNumberOfProcesses(), CBT_dynamic);
    Communicator com(domain->giveEngngModel(), &cb, this->giveRank(), this->giveNumberOfProcesses(), CommMode_Dynamic);
    
    com.packAllData(this, this, &Subdivision::packSharedIrregulars);
    com.initExchange(SHARED_IRREGULAR_DATA_TAG);
    com.unpackAllData(this, this, &Subdivision::unpackSharedIrregulars);
    
    this->sharedIrregularsQueue.clear();
    return false;
  }
}


int
Subdivision::packSharedIrregulars (Subdivision *s, ProcessCommunicator &pc)
{
  int pi, inode, jnode, rproc = pc.giveRank();
  int myrank = this->giveRank();
  const IntArray *sharedPartitions;
  std::list<int>::const_iterator sharedIrregQueueIter;
  IntArray edgeInfo (2);
  if ( rproc == myrank ) return 1; // skip local partition

  // query process communicator to use
  ProcessCommunicatorBuff *pcbuff = pc.giveProcessCommunicatorBuff();
  for (sharedIrregQueueIter = sharedIrregularsQueue.begin(); 
       sharedIrregQueueIter != sharedIrregularsQueue.end(); 
       sharedIrregQueueIter++) {
    pi = (*sharedIrregQueueIter);
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
  RS_Element *elem,  *elem1, *elem2;
  RS_IrregularNode *irregular;
  bool found = false;
	int eIndex, ngb1, ngb2;

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
			// HUHU do not search over all elements, use connectivity if available
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
					found=true;
          if (elem->giveIrregularOn (bEnt)) {
            // entry already exists 
#ifdef __VERBOSE_PARALLEL
            OOFEM_LOG_INFO ("already exists as %d on %d elem\n", elem->giveIrregularOn (bEnt), ie);
#endif
						// search over remaining local elements can be safely skipped
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
            irregular = new Subdivision::RS_IrregularNode (iNum, 0, coords, density, true);
            irregular->setPartition (iproc);
            irregular->setParallelMode(DofManager_shared);
            irregular->setEdgeNodes(bEnt.at(1), bEnt.at(2));
            mesh->addNode (iNum, irregular);
						// simplified implementation (either purely 2D or 3D mesh assumed)   // HUHU
						if (domain->giveElement(elem->giveTopParent())->giveGeometryType() == EGT_triangle_1) {
							// there may be only one local element to add irregular on
							elem->addIrregularOn (iNum, bEnt);
							// schedule elem for bisection
							subdivqueue.push(ie);
#ifdef __VERBOSE_PARALLEL
							OOFEM_LOG_INFO ("added as %d on %d elem\n", iNum, ie);
#endif
						} else if (domain->giveElement(elem->giveTopParent())->giveGeometryType() == EGT_tetra_1) {
							// add iregular to all local elements sharing bEnt edge
							elem->addIrregularOn (iNum, bEnt);
							// schedule elem for bisection
							subdivqueue.push(ie);

#ifdef __VERBOSE_PARALLEL
							OOFEM_LOG_INFO ("added as %d on elems: %d", iNum, ie);
#endif
							eIndex = elem->giveEdgeIndex(bEnt.at(1), bEnt.at(2));
							if(eIndex <= 3){
								ngb1 = 1;
								ngb2 = eIndex+1;
							}
							else{
								ngb1 = eIndex-2;
								ngb2 = (eIndex>4)?eIndex-3:4;
							}

#ifdef DEBUG_INFO
							fprintf(stderr, "[%d] Irregular %d added on %d (edge %d, nodes %d %d, nds %d %d %d %d, ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
											myrank, iNum, elem->giveNumber(), eIndex, bEnt.at(1), bEnt.at(2), 
											elem->giveNode(1), elem->giveNode(2), elem->giveNode(3), elem->giveNode(4), 
											elem->giveNeighbor(1), elem->giveNeighbor(2), elem->giveNeighbor(3), elem->giveNeighbor(4),
											elem->giveIrregular(1), elem->giveIrregular(2), elem->giveIrregular(3), 
											elem->giveIrregular(4), elem->giveIrregular(5), elem->giveIrregular(6));
#endif

							elem1=elem;
							while(elem1->giveNeighbors()->at(ngb1)){
								elem2=mesh->giveElement(elem1->giveNeighbors()->at(ngb1));
								if(elem2==elem){
									ngb2=0; // there is no need to traverse in the other direction (loop was completed)
									break;
								}
								// do not stop traversal if the neighbour is remote
								// there still may be local elements further
								if(elem2->giveParallelMode() == Element_local){
									elem2->addIrregularOn(iNum,bEnt);
									eIndex=elem2->giveEdgeIndex(bEnt.at(1),bEnt.at(2));
									// add neighbour to list of elements for subdivision
									subdivqueue.push(elem2->giveNumber());
#ifdef __VERBOSE_PARALLEL
									OOFEM_LOG_INFO (" %d", elem2->giveNumber());
#endif

#ifdef DEBUG_INFO
									fprintf(stderr, "[%d] Irregular %d added on (ngb1) %d (edge %d, nodes %d %d, nds %d %d %d %d, ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
													myrank, iNum, elem2->giveNumber(), eIndex, bEnt.at(1), bEnt.at(2), 
													elem2->giveNode(1), elem2->giveNode(2), elem2->giveNode(3), elem2->giveNode(4), 
													elem2->giveNeighbor(1), elem2->giveNeighbor(2), elem2->giveNeighbor(3), elem2->giveNeighbor(4),
													elem2->giveIrregular(1), elem2->giveIrregular(2), elem2->giveIrregular(3), 
													elem2->giveIrregular(4), elem2->giveIrregular(5), elem2->giveIrregular(6));
#endif

								}
								if(eIndex <= 3){
									if(elem2->giveNeighbor(1)==elem1->giveNumber())
										ngb1=eIndex+1;
									else
										ngb1=1;
								}
								else{
									if(elem2->giveNeighbor(eIndex-2)==elem1->giveNumber())
										ngb1=(eIndex>4)?eIndex-3:4;
									else
										ngb1=eIndex-2;
								}
								elem1=elem2;
							}

							if(ngb2){
								// traverse neighbours in other direction
								elem1=elem;
								while(elem1->giveNeighbor(ngb2)){
									elem2=mesh->giveElement(elem1->giveNeighbor(ngb2));
									// do not stop traversal if the neighbour is remote
									// there still may be local elements further
									if(elem2->giveParallelMode() == Element_local){
										elem2->addIrregularOn(iNum,bEnt);
										eIndex=elem2->giveEdgeIndex(bEnt.at(1),bEnt.at(2));
										// add neighbour to list of elements for subdivision
										subdivqueue.push(elem2->giveNumber());
#ifdef __VERBOSE_PARALLEL
										OOFEM_LOG_INFO (" %d", elem2->giveNumber());
#endif

#ifdef DEBUG_INFO
										fprintf(stderr, "[%d] Irregular %d added on (ngb2) %d (edge %d, nodes %d %d, nds %d %d %d %d, ngbs %d %d %d %d, irr %d %d %d %d %d %d)\n",
														myrank, iNum, elem2->giveNumber(), eIndex, bEnt.at(1), bEnt.at(2), 
														elem2->giveNode(1), elem2->giveNode(2), elem2->giveNode(3), elem2->giveNode(4), 
														elem2->giveNeighbor(1), elem2->giveNeighbor(2), elem2->giveNeighbor(3), elem2->giveNeighbor(4),
														elem2->giveIrregular(1), elem2->giveIrregular(2), elem2->giveIrregular(3), 
														elem2->giveIrregular(4), elem2->giveIrregular(5), elem2->giveIrregular(6));
#endif

									}
									if(eIndex <= 3){
										if(elem2->giveNeighbor(1)==elem1->giveNumber())
											ngb2=eIndex+1;
										else
											ngb2=1;
									}
									else{
										if(elem2->giveNeighbor(eIndex-2)==elem1->giveNumber())
											ngb2=(eIndex>4)?eIndex-3:4;
										else
											ngb2=eIndex-2;
									}
									elem1=elem2;
								}
							}
#ifdef __VERBOSE_PARALLEL
							OOFEM_LOG_INFO ("\n");
#endif
 						} else {
							OOFEM_ERROR2 ("Subdivision::unpackSharedIrregulars: Unsupported element geometry (element %d)", ie);
						}
						// search over remaining local elements can be safely skipped
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
    if (this->isNodeLocalSharedIrregular(mesh->giveNode(i), myrank) && (mesh->giveNode(i)->giveGlobalNumber() < 0) ) {
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
        OOFEM_LOG_INFO("[%d] packIrregularSharedGlobnums: sending [%d][%d %d]\n",myrank, edgeInfo.at(1), edgeInfo.at(2), edgeInfo.at(3) );
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
			// HUHU do not search over all elements, use connectivity if available
      for (found = false, ie=1; ie<=nelems; ie++) { 
        elem = mesh->giveElement(ie);
        if (elem->giveParallelMode() == Element_remote) continue;
        if (!elem->isTerminal()) continue;
        if ((inode=elem->containsGlobalNode(edgeInfo.at(2))) &&
            (jnode=elem->containsGlobalNode(edgeInfo.at(3)))) {
          // element found => set globnum to corresponding irregular
          be.at(1) = elem->giveNode(inode); be.at(2) = elem->giveNode(jnode);
          irregNum = elem->giveIrregularOn (be);
          if (irregNum) {
            node = mesh->giveNode(irregNum);
            node->setGlobalNumber(edgeInfo.at(1));
#ifdef __VERBOSE_PARALLEL
            OOFEM_LOG_INFO("[%d] unpackIrregularSharedGlobnums: received [%d][%d %d] as %d\n",myrank, edgeInfo.at(1), edgeInfo.at(2), edgeInfo.at(3), irregNum);
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
  OOFEM_LOG_INFO("[%d] Subdivision::assignGlobalNumbersToElements: max globnum %d, new elements %d\n",myrank,globalMaxGlobnum,numberOfLocalElementsToNumber);
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


int 
Subdivision::RS_CompareNodePositions::operator() (int i, int j)
{
	double icoord, jcoord;
	
	icoord=m->giveNode(i)->giveCoordinates()->at(1);
	jcoord=m->giveNode(j)->giveCoordinates()->at(1);

  if (icoord < jcoord) return -1;
	if (icoord > jcoord) return 1;
	
 	icoord=m->giveNode(i)->giveCoordinates()->at(2);
	jcoord=m->giveNode(j)->giveCoordinates()->at(2);

  if (icoord < jcoord) return -1;
	if (icoord > jcoord) return 1;
	
	icoord=m->giveNode(i)->giveCoordinates()->at(3);
	jcoord=m->giveNode(j)->giveCoordinates()->at(3);

  if (icoord < jcoord) return -1;
	if (icoord > jcoord) return 1;
	
	return 0;
}
