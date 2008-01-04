/* $Header: /home/cvs/bp/oofem/tm/src/quad1_ht.C,v 1.2 2003/04/23 14:22:15 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2002   Borek Patzak                                       



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

#include "quad1_ht.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "engngm.h"
#include "structuralms.h"
#include "load.h"
#ifndef __MAKEDEPEND
#include <math.h>
#include <stdio.h>
#endif

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#include "oofegutils.h"
#include "conTable.h"
#endif


Quad1_ht :: Quad1_ht (int n, Domain* aDomain, ElementMode em)
: TransportElement (n,aDomain,em)
   // Constructor.
{
 numberOfDofMans  = 4 ;
}

Quad1_ht :: ~Quad1_ht ()
// Destructor
{}

void
Quad1_ht :: giveDerivativeKsi (FloatArray& answer, double eta)
{
 answer.resize(4);

   answer.at(1) =  0.25*(1.+eta);
   answer.at(2) = -0.25*(1.+eta);
   answer.at(3) = -0.25*(1.-eta);
   answer.at(4) =  0.25*(1.-eta);

}

void
Quad1_ht :: giveDerivativeEta (FloatArray& answer, double ksi)
{
 answer.resize(4);
 
 answer.at(1) =  0.25*(1.+ksi);
 answer.at(2) =  0.25*(1.-ksi);
 answer.at(3) = -0.25*(1.-ksi);
 answer.at(4) = -0.25*(1.+ksi);

}



void
Quad1_ht :: computeNSubMatrixAt (FloatMatrix& answer, FloatArray* coords)
// Returns the displacement interpolation matrix {N} of the receiver,
// evaluated at aGaussPoint.
{
 double       ksi,eta;

 ksi = coords->at(1) ;
 eta = coords->at(2) ;

 answer.resize(1,4);

 answer.at(1,1) = (1. + ksi) * (1. + eta) * 0.25 ;
 answer.at(1,2) = (1. - ksi) * (1. + eta) * 0.25 ;
 answer.at(1,3) = (1. - ksi) * (1. - eta) * 0.25 ;
 answer.at(1,4) = (1. + ksi) * (1. - eta) * 0.25 ;

 return ;
}

void  
Quad1_ht :: computeNmatrixAt (FloatMatrix& answer, FloatArray* coords)
{
 if (emode==HeatTransferEM) this->computeNSubMatrixAt (answer, coords);
 else {
  FloatMatrix n;
  int i,j;

  this->computeNSubMatrixAt (n,coords);
  answer.resize(2,8);
  for (i=1; i<=2; i++)
   for (j=1;j<=4;j++) answer.at(i,(j-1)*2+i) = n.at(1,j);
 }
 
}

void  
Quad1_ht :: computeGradientMatrixAt (FloatMatrix& answer, GaussPoint* aGaussPoint)
{
 int i;
 FloatArray nk,ne;
 FloatMatrix jacobianMtrx(2,2),inv;
 double ksi = aGaussPoint -> giveCoordinate(1);
 double eta = aGaussPoint -> giveCoordinate(2);
 
 this->giveDerivativeKsi (nk,eta);
 this->giveDerivativeEta (ne,ksi);
 
 this->computeJacobianMatrix(jacobianMtrx,aGaussPoint);
 inv.beInverseOf(jacobianMtrx);

 answer.resize(2,4);
 answer.zero();
 
 for (i=1;i<=4;i++){
  answer.at(1,i) = nk.at(i)*inv.at(1,1) + ne.at(i)*inv.at(1,2);
  answer.at(2,i) = nk.at(i)*inv.at(2,1) + ne.at(i)*inv.at(2,2);
 }
}


void
Quad1_ht :: computeJacobianMatrix (FloatMatrix& answer, GaussPoint* aGaussPoint)
// Returns the jacobian matrix  J (x,y)/(ksi,eta)  of the receiver.
// Computes it if it does not exist yet.
{
 int        i;
 double     ksi,eta,x,y;
 FloatArray nx, ny;

 answer.resize(2,2);
 answer.zero();

 ksi = aGaussPoint -> giveCoordinate(1);
 eta = aGaussPoint -> giveCoordinate(2);

 this->giveDerivativeKsi(nx,eta);
 this->giveDerivativeEta(ny,ksi);

 for (i=1;i<=4;i++) {
  x = this -> giveNode(i) -> giveCoordinate(1);
  y = this -> giveNode(i) -> giveCoordinate(2);

  answer.at(1,1) += nx.at(i)*x;
  answer.at(1,2) += nx.at(i)*y;
  answer.at(2,1) += ny.at(i)*x;
  answer.at(2,2) += ny.at(i)*y;
 }
}

void
Quad1_ht :: computeGaussPoints ()
// Sets up the array containing the four Gauss points of the receiver.
{
 MaterialMode mmode;
 
 if (emode==HeatTransferEM) mmode = _2dHeat;
 else mmode = _2dHeMo;

 numberOfIntegrationRules = 1 ;
 integrationRulesArray = new IntegrationRule*[1];
 integrationRulesArray[0] = new GaussIntegrationRule (1,this, 1, 3);
 integrationRulesArray[0]->setUpIntegrationPoints (_Square, numberOfGaussPoints, mmode);

}

void
Quad1_ht ::   giveDofManDofIDMask  (int inode, EquationID, IntArray& answer) const {
// returns DofId mask array for inode element node.
// DofId mask array determines the dof ordering requsted from node.
// DofId mask array contains the DofID constants (defined in cltypes.h)
// describing physical meaning of particular DOFs.
 if (emode == HeatTransferEM) {
  answer.resize (1);
  answer.at(1) = T_f;
 } else if (emode == HeatMass1TransferEM) {
  answer.resize (2);
  answer.at(1) = T_f;
  answer.at(2) = C_1;
 } else {
  _error ("Unknown ElementMode");
 }
}

IRResultType
Quad1_ht :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

 this->Element :: initializeFrom (ir);
 numberOfGaussPoints = 4;
 IR_GIVE_OPTIONAL_FIELD (ir, numberOfGaussPoints, IFT_Quad1_ht_nip, "nip"); // Macro
 
 if (!((numberOfGaussPoints == 4) ||
   (numberOfGaussPoints == 9) ||
   (numberOfGaussPoints == 16)))
   numberOfGaussPoints = 4;

 this -> computeGaussPoints();
 return IRRT_OK;
}


double
Quad1_ht :: computeVolumeAround (GaussPoint* aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
   double       determinant,weight,thickness,volume;
  FloatMatrix jacobianMtrx;
  this->computeJacobianMatrix(jacobianMtrx,aGaussPoint);
   determinant = fabs (jacobianMtrx.giveDeterminant());
   weight      = aGaussPoint -> giveWeight();
   thickness   = this -> giveCrossSection() -> give('t');
   volume      = determinant * weight * thickness ;
  
   return volume;
}

void 
Quad1_ht :: computeEgdeNMatrixAt(FloatMatrix& answer, GaussPoint* gp)
{
  /*
  
  computes interpolation matrix for element edge.
  we assemble locally this matrix for only nonzero 
  shape functions. 
  (for example only two nonzero shape functions for 2 dofs are
  necessary for linear plane stress tringle edge).
  These nonzero shape functions are then mapped to 
  global element functions.
  
  Using mapping technique will allow to assemble shape functions 
  without regarding particular side
  */
 
  double ksi,n1,n2;
 answer.resize (1,2);
  answer.zero();

  ksi = gp -> giveCoordinate(1) ;
  n1  = (1. - ksi) * 0.5 ;
  n2  = (1. + ksi) * 0.5 ;

  answer.at(1,1) = n1 ;
  answer.at(1,2) = n2 ;
  
  return  ;
}


double
Quad1_ht :: computeEdgeVolumeAround(GaussPoint* gp, int iEdge)
{ 
 double dx,dy, length, thick ;
  Node   *nodeA,*nodeB ;
  int aNode=0, bNode=0;
 
  if (iEdge == 1) { // edge between nodes 1 2
  aNode = 1;
  bNode = 2;
  } else if (iEdge == 2) { // edge between nodes 2 3
  aNode = 2;
  bNode = 3;
  } else if (iEdge == 3) { // edge between nodes 3 4
  aNode = 3;
  bNode = 4;
  } else if (iEdge == 4) { // edge between nodes 4 1
  aNode = 4;
  bNode = 1;
  } else {
  _error ("computeEdgeVolumeAround: wrong egde number");
  }
 
  nodeA   = this->giveNode(aNode) ;
  nodeB   = this->giveNode(bNode) ;
  
  dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1) ;
  dy      = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2) ;
  length = sqrt(dx*dx + dy*dy);
 thick = this -> giveCrossSection() -> give('t');
  return 0.5 * length * thick * gp -> giveWeight() ;
 
}


void
Quad1_ht :: giveEdgeDofMapping (IntArray& answer, int iEdge)
{
  /* 
   provides dof mapping of local edge dofs (only nonzero are taken into account)
   to global element dofs
  */

  answer.resize(2);
  if (iEdge == 1)  { // edge between nodes 1,2
  answer.at(1) = 1;
  answer.at(2) = 2;
  } else if (iEdge == 2) { // edge between nodes 2 3
  answer.at(1) = 2;
  answer.at(2) = 3;
  } else if (iEdge == 3) {// edge between nodes 3 4
  answer.at(1) = 3;
  answer.at(2) = 4;
  } else if (iEdge == 4) {// edge between nodes 4 1
  answer.at(1) = 4;
  answer.at(2) = 1;
  } else {
  _error ("giveEdgeDofMapping: wrong edge number");
  }

  return ;
}

void      
Quad1_ht :: computeEdgeIpGlobalCoords (FloatArray& answer, GaussPoint* gp, int iEdge) 
{

  double n1, n2, ksi = gp->giveCoordinate(1) ;
  Node   *nodeA,*nodeB ;
  int aNode=0, bNode=0;

  if (iEdge == 1) { // edge between nodes 1 2
  aNode = 1;
  bNode = 2;
  } else if (iEdge == 2) { // edge between nodes 2 3
  aNode = 2;
  bNode = 3;
  } else if (iEdge == 3) { // edge between nodes 3 4
  aNode = 3;
  bNode = 4;
  } else if (iEdge == 4) { // edge between nodes 4 1
  aNode = 4;
  bNode = 1;
  } else {
  _error ("computeEdgeIpGlobalCoords: wrong egde number");
  }

 n1  = (1. - ksi) * 0.5 ;
 n2  = (1. + ksi) * 0.5 ;

  nodeA   = this->giveNode(aNode) ;
  nodeB   = this->giveNode(bNode) ;
  
 answer.resize (2);
 answer.at(1) = n1*nodeA->giveCoordinate(1)+n2*nodeB->giveCoordinate(1);
 answer.at(2) = n1*nodeA->giveCoordinate(2)+n2*nodeB->giveCoordinate(2);
}

void 
Quad1_ht :: computeInternalSourceRhsVectorAt (FloatArray& answer, TimeStep* atTime, ValueModeType mode)
{
 if (emode == HeatTransferEM) {
  this->computeInternalSourceRhsSubVectorAt (answer, atTime, mode, 1) ;
 } else if (emode == HeatMass1TransferEM) {
  FloatArray subAnswer;
  int i;
  
  for (i=1; i<=2; i++) {
   this->computeInternalSourceRhsSubVectorAt (subAnswer, atTime, mode, i);
   if (subAnswer.isNotEmpty()) {
    if (answer.isEmpty()) {answer.resize(8); answer.zero();}
    this->assembleLocalContribution (answer, subAnswer, 2, i, 1.0);
   }
  }    
 } else {
  _error ("Unknown ElementMode");
 }
}

int 
Quad1_ht :: computeGlobalCoordinates (FloatArray& answer, const FloatArray& lcoords)
{
 double ksi, eta, n1, n2, n3, n4;
 
 ksi = lcoords.at(1) ;
 eta = lcoords.at(2) ;
 
 n1 = (1. + ksi) * (1. + eta) * 0.25 ;
 n2 = (1. - ksi) * (1. + eta) * 0.25 ;
 n3 = (1. - ksi) * (1. - eta) * 0.25 ;
 n4 = (1. + ksi) * (1. - eta) * 0.25 ;
 
 answer.resize(2);
 answer.at(1) = n1*this->giveNode(1)->giveCoordinate(1) + n2*this->giveNode(2)->giveCoordinate(1)+
  n3 * this->giveNode(3)->giveCoordinate(1) + n4 * this->giveNode(4)->giveCoordinate(1);
 answer.at(2) = n1*this->giveNode(1)->giveCoordinate(2) + n2*this->giveNode(2)->giveCoordinate(2)+
  n3 * this->giveNode(3)->giveCoordinate(2) + n4 * this->giveNode(4)->giveCoordinate(2);

 return 1;
}


Interface* 
Quad1_ht :: giveInterface (InterfaceType interface)
{
 if (interface == SpatialLocalizerInterfaceType) return (SpatialLocalizerInterface*) this;
 else if (interface == EIPrimaryFieldInterfaceType) return (EIPrimaryFieldInterface*) this;
 return NULL;
}


int 
Quad1_ht :: SpatialLocalizerI_containsPoint (const FloatArray& coords) {

 FloatArray lcoords;
 return this->computeLocalCoordinates (lcoords, coords);
}


double 
Quad1_ht::SpatialLocalizerI_giveDistanceFromParametricCenter (const FloatArray& coords)
{
 FloatArray lcoords(2), gcoords;
 double dist;
 int size, gsize;

 lcoords.at(1) = lcoords.at(2) = 0.0;
 this -> computeGlobalCoordinates (gcoords, lcoords);

 if((size = coords.giveSize()) < (gsize = gcoords.giveSize()))
  _error("SpatialLocalizerI_giveDistanceFromParametricCenter: coordinates size mismatch");

 if(size == gsize){
  dist = coords.distance(gcoords);
 }
 else{
  FloatArray helpCoords = coords;

  helpCoords.resize(gsize);
  dist = helpCoords.distance(gcoords);
 }

 return dist;
}


#define POINT_TOL 1.e-6

int
Quad1_ht::computeLocalCoordinates (FloatArray& answer, const FloatArray& coords)
{
  Node    *node1,*node2,*node3,*node4;
  double  x1,x2,x3,x4,y1,y2,y3,y4,a1,a2,a3,a4,b1,b2,b3,b4;
  double a,b,c, ksi1, ksi2, ksi3, eta1=0.0, eta2=0.0, denom;
  int nroot;
  
  answer.resize(2);

  node1 = this -> giveNode(1) ;
  node2 = this -> giveNode(2) ;
  node3 = this -> giveNode(3) ;
  node4 = this -> giveNode(4) ;
  

  x1 = node1 -> giveCoordinate(1) ;
  x2 = node2 -> giveCoordinate(1) ;
  x3 = node3 -> giveCoordinate(1) ;
  x4 = node4 -> giveCoordinate(1) ;
  
  y1 = node1 -> giveCoordinate(2) ;
  y2 = node2 -> giveCoordinate(2) ;
  y3 = node3 -> giveCoordinate(2) ;
  y4 = node4 -> giveCoordinate(2) ;
  
  a1 = x1+x2+x3+x4;
  a2 = x1-x2-x3+x4;
  a3 = x1+x2-x3-x4;
  a4 = x1-x2+x3-x4;
  
  b1 = y1+y2+y3+y4;
  b2 = y1-y2-y3+y4;
  b3 = y1+y2-y3-y4;
  b4 = y1-y2+y3-y4;

  a = a2*b4-b2*a4;
  b = a1*b4+a2*b3-a3*b2-b1*a4-b4*4.0*coords.at(1)+a4*4.0*coords.at(2);
  c = a1*b3-a3*b1-4.0*coords.at(1)*b3+4.0*coords.at(2)*a3;

  // solve quadratic equation for ksi
  cubic (0.0, a, b, c, &ksi1, &ksi2, &ksi3, &nroot);

  if (nroot == 0) return 0;
  if(nroot){
   denom = (b3+ksi1*b4);
   if (fabs(denom) <= 1.0e-10) 
    eta1 = (4.0*coords.at(1)-a1-ksi1*a2) / (a3+ksi1*a4);
   else eta1 = (4.0*coords.at(2)-b1-ksi1*b2) / denom;
  }

  if(nroot > 1){
   double diff_ksi1, diff_eta1, diff_ksi2, diff_eta2, diff1, diff2;

   denom = b3+ksi2*b4;
   if (fabs(denom) <= 1.0e-10) 
    eta2 = (4.0*coords.at(1)-a1-ksi2*a2) / (a3+ksi2*a4);
   else eta2 = (4.0*coords.at(2)-b1-ksi2*b2) / denom;
   
   // choose the one which seems to be closer to the parametric space (square <-1;1>x<-1;1>)
   diff_ksi1 = 0.0;
   if(ksi1 > 1.0)diff_ksi1 = ksi1 - 1.0;
   if(ksi1 < -1.0)diff_ksi1 = ksi1 + 1.0;

   diff_eta1 = 0.0;
   if(eta1 > 1.0)diff_eta1 = eta1 - 1.0;
   if(eta1 < -1.0)diff_eta1 = eta1 + 1.0;
   
   diff_ksi2 = 0.0;
   if(ksi2 > 1.0)diff_ksi2 = ksi2 - 1.0;
   if(ksi2 < -1.0)diff_ksi2 = ksi2 + 1.0;

   diff_eta2 = 0.0;
   if(eta2 > 1.0)diff_eta2 = eta2 - 1.0;
   if(eta2 < -1.0)diff_eta2 = eta2 + 1.0;
   
   diff1 = diff_ksi1 * diff_ksi1 + diff_eta1 * diff_eta1;
   diff2 = diff_ksi2 * diff_ksi2 + diff_eta2 * diff_eta2;

   // ksi2, eta2 seems to be closer
   if(diff1 > diff2){
    ksi1 = ksi2;
    eta1 = eta2;
   }
  }

  answer.at(1) = ksi1;
  answer.at(2) = eta1;

  // test if inside
  for (int i=1; i<=2; i++) {
   if (answer.at(i)<(-1.-POINT_TOL)) return 0;
   if (answer.at(i)>(1.+POINT_TOL)) return 0;
  }

  return 1;
}

/*
#define _SMALLNUM 1.e-6
int
Quad1_ht::computeLocalCoordinates (FloatArray& answer, const FloatArray& coords)
{
   Node    *node1,*node2,*node3,*node4;
   double  x1,x2,x3,x4,y1,y2,y3,y4,a1,a2,a3,a4,b1,b2,b3,b4;
  double a,b,c, r1, r2, r3, eta, denom1, denom2;
  int nroot;

  answer.resize(2);

   node1 = this -> giveNode(1) ;
   node2 = this -> giveNode(2) ;
   node3 = this -> giveNode(3) ;
   node4 = this -> giveNode(4) ;


   x1 = node1 -> giveCoordinate(1) ;
   x2 = node2 -> giveCoordinate(1) ;
   x3 = node3 -> giveCoordinate(1) ;
   x4 = node4 -> giveCoordinate(1) ;

   y1 = node1 -> giveCoordinate(2) ;
   y2 = node2 -> giveCoordinate(2) ;
   y3 = node3 -> giveCoordinate(2) ;
   y4 = node4 -> giveCoordinate(2) ;

  a1 = x1+x2+x3+x4;
  a2 = x1-x2-x3+x4;
  a3 = x1+x2-x3-x4;
  a4 = x1-x2+x3-x4;

  b1 = y1+y2+y3+y4;
  b2 = y1-y2-y3+y4;
  b3 = y1+y2-y3-y4;
  b4 = y1-y2+y3-y4;

  a = a2*b4-b2*a4;
  b = a1*b4+a2*b3-a3*b2-b1*a4-b4*4.0*coords.at(1)+a4*4.0*coords.at(2);
  c = a1*b3-a3*b1-4.0*coords.at(1)*b3+4.0*coords.at(2)*a3;

  // solve quadratic equation for ksi
  cubic (0.0, a, b, c, &r1, &r2, &r3, &nroot);
  if (nroot) {
   if ((r1>=(-1.0-_SMALLNUM))&&(r1<=(1.0+_SMALLNUM))) {
    denom1 = (b3+r1*b4);
    denom2 = (a3+r1*a4);
    if (fabs(denom2) > fabs(denom1)) 
     eta = (4.0*coords.at(1)-a1-r1*a2) / denom2;
    else eta = (4.0*coords.at(2)-b1-r1*b2) / denom1;
    
    if ((eta>=(-1.0-_SMALLNUM))&&(eta<=(1.0+_SMALLNUM))) {
     answer.at(1) = r1;
     answer.at(2) = eta;
     return 1;
    }
   } 
  }
  if (nroot > 1) {
   if ((r2>=(-1.0-_SMALLNUM))&&(r2<=(1.0+_SMALLNUM))) {
    eta = (4.0*coords.at(2)-b1-r2*b2) / (b3+r2*b4);
    fprintf (stderr, "%e\n", eta);
    if ((eta>=(-1.0-_SMALLNUM))&&(eta<=(1.0+_SMALLNUM))) {
     answer.at(1) = r2;
     answer.at(2) = eta;
     return 1;
    }
   }
  }
   
  // given point not inside receiver volume
  return 0;
}
*/



#ifdef __OOFEG
void Quad1_ht :: drawRawGeometry (oofegGraphicContext &gc)
{
  WCRec p[4];
  GraphicObj *go;

 if (!gc.testElementGraphicActivity(this)) return; 

  EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
  EASValsSetColor(gc.getElementColor());
  EASValsSetEdgeColor(gc.getElementEdgeColor());
  EASValsSetEdgeFlag(TRUE);
  EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
 EASValsSetFillStyle (FILL_HOLLOW);
  p[0].x = (FPNum) this->giveNode(1)->giveCoordinate(1);
  p[0].y = (FPNum) this->giveNode(1)->giveCoordinate(2);
  p[0].z = 0.;
  p[1].x = (FPNum) this->giveNode(2)->giveCoordinate(1);
  p[1].y = (FPNum) this->giveNode(2)->giveCoordinate(2);
  p[1].z = 0.;
  p[2].x = (FPNum) this->giveNode(3)->giveCoordinate(1);
  p[2].y = (FPNum) this->giveNode(3)->giveCoordinate(2);
  p[2].z = 0.;
  p[3].x = (FPNum) this->giveNode(4)->giveCoordinate(1);
  p[3].y = (FPNum) this->giveNode(4)->giveCoordinate(2);
  p[3].z = 0.;
   
  go =  CreateQuad3D(p);
  EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
  EGAttachObject(go, (EObjectP) this);
  EMAddGraphicsToModel(ESIModel(), go);
}

void Quad1_ht :: drawScalar   (oofegGraphicContext& context)
{
 int i, indx, result=0;
 WCRec p[4];
 GraphicObj *tr;
 TimeStep* tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
 double s[4];
 FloatArray v[4];
 IntArray map;
 InternalStateType itype = context.giveIntVarType();

 if (!context.testElementGraphicActivity(this)) return; 

 EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);

 if (itype == IST_HydrationDegree) { 
   for (i=1; i<=4; i++)
     result+= this->giveInternalStateAtNode (v[i-1], context.giveIntVarType(), context.giveIntVarMode(), i, tStep);
   if (result !=4) return;
   
   this->giveIntVarCompFullIndx (map, context.giveIntVarType());
   if ((indx = map.at(context.giveIntVarIndx())) == 0) return;
   for (i=1; i<=4; i++) s[i-1]=v[i-1].at(indx);
   
   if (context.getScalarAlgo() == SA_ISO_SURF) {
     
     for (i=0; i< 4; i++) {
       p[i].x = (FPNum) this->giveNode(i+1)->giveCoordinate(1);
       p[i].y = (FPNum) this->giveNode(i+1)->giveCoordinate(2);
       p[i].z = 0.;
     }

   context.updateFringeTableMinMax (s, 4);
   tr =  CreateQuadWD3D(p, s[0], s[1], s[2], s[3]);
   EGWithMaskChangeAttributes(LAYER_MASK, tr);
   EMAddGraphicsToModel(ESIModel(), tr);
   }
 } else  if ((((emode == HeatTransferEM)||(emode == HeatMass1TransferEM)) && (itype == IST_Temperature)) || 
             ((emode == HeatMass1TransferEM) && (itype == IST_MassConcentration_1))   ) {
   IntArray dofMask (1) ; 
   if (itype == IST_Temperature) dofMask.at(1) = T_f;
   else dofMask.at(1) = C_1;
   
   FloatArray r;
   for (i=0; i< 4; i++) {
     this->giveNode(i+1)->giveUnknownVector (r, dofMask, EID_ConservationEquation, VM_Total, tStep);
     s[i]=r.at(1);
     
     p[i].x = (FPNum) this->giveNode(i+1)->giveCoordinate(1);
     p[i].y = (FPNum) this->giveNode(i+1)->giveCoordinate(2);
     p[i].z = 0.;
   }
   
   context.updateFringeTableMinMax (s, 4);
   tr =  CreateQuadWD3D(p, s[0], s[1], s[2], s[3]);
   EGWithMaskChangeAttributes(LAYER_MASK, tr);
   EMAddGraphicsToModel(ESIModel(), tr);
 }
 return;
}



#endif
