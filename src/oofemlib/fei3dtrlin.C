/* $Header: /home/cvs/bp/oofem/oofemlib/src/fei3dtrlin.C,v 1.1.4.1 2004/04/05 15:19:43 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/

#include "fei3dtrlin.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "node.h"
#include "mathfem.h"

void
FEI3dTrLin::evalN (FloatArray& answer, const FloatArray& lcoords, double time)
{

 answer.resize (4);
 
 answer.at(1) = lcoords.at(1);
 answer.at(2) = lcoords.at(2);
 answer.at(3) = lcoords.at(3);
 answer.at(4) = 1. - lcoords.at(1) - lcoords.at(2) - lcoords.at(3);
 
 return;
}

void
FEI3dTrLin :: evaldNdx (FloatMatrix&answer, const FloatArray** coords, const FloatArray& lcoords, double time)
{
 double x1,x2,x3,x4, y1,y2,y3,y4, z1,z2,z3,z4, vol;
 answer.resize(4,3);

 x1 = coords[0] -> at(1) ;
 x2 = coords[1] -> at(1) ;
 x3 = coords[2] -> at(1) ;
 x4 = coords[3] -> at(1) ;
 
 y1 = coords[0] -> at(2) ;
 y2 = coords[1] -> at(2) ;
 y3 = coords[2] -> at(2) ;
 y4 = coords[3] -> at(2) ;

 z1 = coords[0] -> at(3) ;
 z2 = coords[1] -> at(3) ;
 z3 = coords[2] -> at(3) ;
 z4 = coords[3] -> at(3) ;
 
 vol = ((x4-x1)*(y2-y1)*(z3-z1) - (x4-x1)*(y3-y1)*(z2-z1) +
     (x3-x1)*(y4-y1)*(z2-z1) - (x2-x1)*(y4-y1)*(z3-z1) +
     (x2-x1)*(y3-y1)*(z4-z1) - (x3-x1)*(y2-y1)*(z4-z1))/6.;
 
 if (vol <= 0.0) {
  OOFEM_ERROR ("FEI3dTrLin :: evaldNdx: negative volume");
 }

 answer.at(1,1) = -((y3-y2)*(z4-z2)-(y4-y2)*(z3-z2));
 answer.at(2,1) = (y4-y3)*(z1-z3)-(y1-y3)*(z4-z3);
 answer.at(3,1) = -((y1-y4)*(z2-z4)-(y2-y4)*(z1-z4));
 answer.at(4,1) = (y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
 
 answer.at(1,2) = -((x4-x2)*(z3-z2)-(x3-x2)*(z4-z2));
 answer.at(2,2) = (x1-x3)*(z4-z3)-(x4-x3)*(z1-z3);
 answer.at(3,2) = -((x2-x4)*(z1-z4)-(x1-x4)*(z2-z4));
 answer.at(4,2) = (x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
 
 answer.at(1,3) = -((x3-x2)*(y4-y2)-(x4-x2)*(y3-y2));
 answer.at(2,3) = (x4-x3)*(y1-y3)-(x1-x3)*(y4-y3);
 answer.at(3,3) = -((x1-x4)*(y2-y4)-(x2-x4)*(y1-y4));
 answer.at(4,3) = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
 
 answer.times(1./(6.*vol));
}

void
FEI3dTrLin :: local2global (FloatArray& answer, const FloatArray** coords, const FloatArray& lcoords, double time)
{
 int i;
 FloatArray l(4);
 answer.resize (3); 
 answer.zero();
 
 l.at(1) = lcoords.at(1);
 l.at(2) = lcoords.at(2);
 l.at(3) = lcoords.at(3);
 l.at(4) = 1.0 - l.at(1) - l.at(2) - l.at(3);
 
 for (i=1; i<=4; i++) { 
  answer.at(1) += l.at(i)*coords[i-1]->at(1);
  answer.at(2) += l.at(i)*coords[i-1]->at(2);
  answer.at(3) += l.at(i)*coords[i-1]->at(3);
 }

}

#define POINT_TOL 1.e-3

int
FEI3dTrLin :: global2local (FloatArray& answer, const FloatArray** nc, const FloatArray& coords, double time)
{
  double       x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,xp,yp,zp, volume;
  answer.resize(4);

  x1 = nc[0] -> at(1) ;
  x2 = nc[1] -> at(1) ;
  x3 = nc[2] -> at(1) ;
  x4 = nc[3] -> at(1) ;
  
  y1 = nc[0] -> at(2) ;
  y2 = nc[1] -> at(2) ;
  y3 = nc[2] -> at(2) ;
  y4 = nc[3] -> at(2) ;
  
  z1 = nc[0] -> at(3) ;
  z2 = nc[1] -> at(3) ;
  z3 = nc[2] -> at(3) ;
  z4 = nc[3] -> at(3) ;
  
  xp= coords.at(1);
  yp= coords.at(2);
  zp= coords.at(3);
  
  volume = ((x4-x1)*(y2-y1)*(z3-z1) - (x4-x1)*(y3-y1)*(z2-z1) +
     (x3-x1)*(y4-y1)*(z2-z1) - (x2-x1)*(y4-y1)*(z3-z1) +
     (x2-x1)*(y3-y1)*(z4-z1) - (x3-x1)*(y2-y1)*(z4-z1))/6.;

  answer.resize(4);
  
  answer.at(1) = ((x3-x2)*(yp-y2)*(z4-z2) - (xp-x2)*(y3-y2)*(z4-z2) +
          (x4-x2)*(y3-y2)*(zp-z2) - (x4-x2)*(yp-y2)*(z3-z2) +
          (xp-x2)*(y4-y2)*(z3-z2) - (x3-x2)*(y4-y2)*(zp-z2))/6./volume;
  
  answer.at(2) = ((x4-x1)*(yp-y1)*(z3-z1) - (xp-x1)*(y4-y1)*(z3-z1) +
          (x3-x1)*(y4-y1)*(zp-z1) - (x3-x1)*(yp-y1)*(z4-z1) +
          (xp-x1)*(y3-y1)*(z4-z1) - (x4-x1)*(y3-y1)*(zp-z1))/6./volume;
  
  answer.at(3) = ((x2-x1)*(yp-y1)*(z4-z1) - (xp-x1)*(y2-y1)*(z4-z1) +
         (x4-x1)*(y2-y1)*(zp-z1) - (x4-x1)*(yp-y1)*(z2-z1) +
          (xp-x1)*(y4-y1)*(z2-z1) - (x2-x1)*(y4-y1)*(zp-z1))/6./volume;
  
  answer.at(4) = 1.0-answer.at(1)-answer.at(2)-answer.at(3);
  
  // test if inside
  for (int i=1; i<=4; i++) {
   if (answer.at(i)<(0.-POINT_TOL)) return 0;
   if (answer.at(i)>(1.+POINT_TOL)) return 0;
  }
  return 1;
 }


double
FEI3dTrLin :: giveTransformationJacobian (const FloatArray** coords, const FloatArray& lcoords, double time)
{
 double       volume, x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;

 x1 = coords[0] -> at(1) ;
 x2 = coords[1] -> at(1) ;
 x3 = coords[2] -> at(1) ;
 x4 = coords[3] -> at(1) ;
 
 y1 = coords[0] -> at(2) ;
 y2 = coords[1] -> at(2) ;
 y3 = coords[2] -> at(2) ;
 y4 = coords[3] -> at(2) ;
 
 z1 = coords[0] -> at(3) ;
 z2 = coords[1] -> at(3) ;
 z3 = coords[2] -> at(3) ;
 z4 = coords[3] -> at(3) ;
 
 volume = ((x4-x1)*(y2-y1)*(z3-z1) - (x4-x1)*(y3-y1)*(z2-z1) +
      (x3-x1)*(y4-y1)*(z2-z1) - (x2-x1)*(y4-y1)*(z3-z1) +
      (x2-x1)*(y3-y1)*(z4-z1) - (x3-x1)*(y2-y1)*(z4-z1))/6.;
 
 if (volume <= 0.0) {
  OOFEM_ERROR ("FEI3dTrLin :: giveTransformationJacobian: negative volume encountered");
 }
 return volume;
}


void 
FEI3dTrLin :: edgeEvalN (FloatArray& answer, const FloatArray& lcoords, double time)
{
 double ksi = lcoords.at(1);
 answer.resize(2);
 
 answer.at(1) = (1. - ksi) * 0.5 ;
 answer.at(2) = (1. + ksi) * 0.5 ;
}

void 
FEI3dTrLin :: edgeEvaldNdx (FloatMatrix&answer, int iedge, 
			    const FloatArray** nc, const FloatArray& lcoords, double time)
{
 double coeff,l,x1,x2,y1,y2,z1,z2;
 IntArray edgeNodes;
 this->computeLocalEdgeMapping (edgeNodes, iedge);
 l = this->edgeComputeLength(edgeNodes, nc);
 coeff = 1.0/l/l;

 x1 = nc[edgeNodes.at(1)-1]->at(1);
 y1 = nc[edgeNodes.at(1)-1]->at(2);
 z1 = nc[edgeNodes.at(1)-1]->at(3);
 x2 = nc[edgeNodes.at(2)-1]->at(1);
 y2 = nc[edgeNodes.at(2)-1]->at(2);
 z2 = nc[edgeNodes.at(2)-1]->at(3);

 answer.resize (2,3);
 answer.at(1,1) = (x1-x2)*coeff;
 answer.at(1,2) = (y1-y2)*coeff;
 answer.at(1,3) = (z1-z2)*coeff;

 answer.at(2,1) = (x2-x1)*coeff;
 answer.at(2,2) = (y2-y1)*coeff;
 answer.at(2,3) = (z2-z1)*coeff;
}

void
FEI3dTrLin :: edgeLocal2global (FloatArray& answer, int iedge, 
				const FloatArray** nc, const FloatArray& lcoords, double time) 
{
 IntArray edgeNodes;
 FloatArray n;
 this->computeLocalEdgeMapping (edgeNodes, iedge);
 this->edgeEvalN (n, lcoords, time);
 
 answer.resize(3); 
 answer.at(1) = (n.at(1)* nc[edgeNodes.at(1)-1]->at(1) +
		 n.at(2)* nc[edgeNodes.at(2)-1]->at(1));
 answer.at(2) = (n.at(1)* nc[edgeNodes.at(1)-1]->at(2) +
		 n.at(2)* nc[edgeNodes.at(2)-1]->at(2));
 answer.at(3) = (n.at(1)* nc[edgeNodes.at(1)-1]->at(3) +
		 n.at(2)* nc[edgeNodes.at(2)-1]->at(3));
  
}


double
FEI3dTrLin :: edgeGiveTransformationJacobian (int iedge, const FloatArray** nc, const FloatArray& lcoords, double time) 
{
 IntArray edgeNodes;
 this->computeLocalEdgeMapping (edgeNodes, iedge);
 return 0.5 * this->edgeComputeLength(edgeNodes, nc);
}


void 
FEI3dTrLin :: computeLocalEdgeMapping (IntArray& edgeNodes, int iedge)
{
 int aNode=0, bNode=0;
 edgeNodes.resize(2);
 
  if (iedge == 1) { // edge between nodes 1 2
  aNode = 1;
  bNode = 2;
  } else if (iedge == 2) { // edge between nodes 2 3
  aNode = 2;
  bNode = 3;
  } else if (iedge == 3) { // edge between nodes 3 1
  aNode = 3;
  bNode = 1;
 } else if (iedge == 4) { // edge between nodes 1 4
  aNode = 1;
  bNode = 4;
 } else if (iedge == 5) { // edge between nodes 2 4
  aNode = 2;
  bNode = 4;
 } else if (iedge == 6) { // edge between nodes 3 4
  aNode = 3;
  bNode = 4;
  } else {
  OOFEM_ERROR2 ("FEI3dTrLin :: computeEdgeMapping: wrong egde number (%d)", iedge);
  }
 
 edgeNodes.at(1) = (aNode);
 edgeNodes.at(2) = (bNode);
}

double
FEI3dTrLin :: edgeComputeLength (IntArray& edgeNodes, const FloatArray** nc)
{
  double dx,dy,dz;
  int nodeA,nodeB ;

  nodeA   = (edgeNodes.at(1)-1) ;
  nodeB   = (edgeNodes.at(2)-1) ;
 
  dx      = nc[nodeB]->at(1) - nc[nodeA]->at(1) ;
  dy      = nc[nodeB]->at(2) - nc[nodeA]->at(2) ;
  dz      = nc[nodeB]->at(3) - nc[nodeA]->at(3) ;
  return (sqrt(dx*dx + dy*dy + dz*dz));
}

void
FEI3dTrLin :: surfaceEvalN (FloatArray& answer, const FloatArray& lcoords, double time)
{
 answer.resize (3);
 
 answer.at(1) = lcoords.at(1);
 answer.at(2) = lcoords.at(2);
 answer.at(3) = 1. - lcoords.at(1) - lcoords.at(2);
 
 return;
}

void 
FEI3dTrLin :: surfaceLocal2global (FloatArray& answer, int iedge,
				   const FloatArray** nc, const FloatArray& lcoords, double time)
{
  double l1,l2,l3;
  answer.resize (3);
  IntArray nodes(3);
 
  computeLocalSurfaceMapping (nodes, iedge);
  
  l1 = lcoords.at(1);
  l2 = lcoords.at(2);
  l3 = 1.0 - l1 - l2;
  
  answer.at(1) = (l1*nc[nodes.at(1)-1]->at(1)+
		  l2*nc[nodes.at(2)-1]->at(1)+
		  l3*nc[nodes.at(3)-1]->at(1));
  answer.at(2) = (l1*nc[nodes.at(1)-1]->at(2)+
		  l2*nc[nodes.at(2)-1]->at(2) +
		  l3*nc[nodes.at(3)-1]->at(2));
  answer.at(3) = (l1*nc[nodes.at(1)-1]->at(3)+
		  l2*nc[nodes.at(2)-1]->at(3) +
		  l3*nc[nodes.at(3)-1]->at(3));
}

double 
FEI3dTrLin :: surfaceGiveTransformationJacobian (int isurf, const FloatArray** nc, const FloatArray& lcoords,
						 double time)
{
 int n1,n2,n3;
 IntArray snodes (3);
 this->computeLocalSurfaceMapping (snodes, isurf);

 n1 = (snodes.at(1)-1);
 n2 = (snodes.at(2)-1);
 n3 = (snodes.at(3)-1);

 FloatArray a(3),b(3),c(3);
 for (int i=1; i<=3; i++) {
  a.at(i) = nc[n2]->at(i)-nc[n1]->at(i);
  b.at(i) = nc[n3]->at(i)-nc[n1]->at(i);
 }
 
 c.beVectorProductOf (a,b);
 return 0.5 * sqrt(dotProduct (c,c,3));
}

void 
FEI3dTrLin :: computeLocalSurfaceMapping (IntArray& surfNodes, int isurf)
{
 int aNode=0, bNode=0, cNode=0;
 surfNodes.resize(3);
 
  if (isurf == 1) { // surface 1 - nodes 1 3 2
  aNode = 1;
  bNode = 3;
  cNode = 2;
  } else if (isurf == 2) { // surface 2 - nodes 1 2 4
  aNode = 1;
  bNode = 2;
  cNode = 4;
  } else if (isurf == 3) { // surface 3  - nodes 2 3 4
  aNode = 2;
  bNode = 3;
  cNode = 4;
 } else if (isurf == 4) { // surface 4 - nodes 1 4 3
  aNode = 1;
  bNode = 4;
  cNode = 3;
 } else {
  OOFEM_ERROR2 ("FEI3dTrLin :: computeSurfaceMapping: wrong surface number (%d)", isurf);
 }
 
 surfNodes.at(1) = (aNode);
 surfNodes.at(2) = (bNode);
 surfNodes.at(3) = (cNode);
}



