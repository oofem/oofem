/* $Header: /home/cvs/bp/oofem/sm/src/ltrelemppde.C,v 1.4 2003/04/06 14:08:30 bp Exp $ */
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

//   file ltrelemppde.CC

 
#include "ltrelemppde.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
// #include "polynoxy.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "verbose.h"
#include "engngm.h"
#ifndef __MAKEDEPEND
#include <math.h>
#include <stdio.h>
#endif
#include "mathfem.h"

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#include "conTable.h"
#endif

LTrElementPPDE :: LTrElementPPDE (int n, Domain* aDomain)
      : PPdeElement (n,aDomain)
   // Constructor.
{
   numberOfDofMans  = 3 ;
   area = -1;
}


double
LTrElementPPDE :: giveArea ()
// returns the area occupied by the receiver
{
 if (area > 0) return area;  // check if previously computed

   Node    *node1,*node2,*node3;
   double  x1,x2,x3,y1,y2,y3;

   node1 = this -> giveNode(1) ;
   node2 = this -> giveNode(2) ;
   node3 = this -> giveNode(3) ;

   x1 = node1 -> giveCoordinate(1) ;
   x2 = node2 -> giveCoordinate(1) ;
   x3 = node3 -> giveCoordinate(1) ;

   y1 = node1 -> giveCoordinate(2) ;
   y2 = node2 -> giveCoordinate(2) ;
   y3 = node3 -> giveCoordinate(2) ;

   return (area = 0.5*(x2*y3+x1*y2+y1*x3-x2*y1-x3*y2-x1*y3)) ;
 
}


FloatArray* LTrElementPPDE :: GivebCoeff ()
//
// Returns coefficients of partial derivatives of shape functions
// with respect to x
//
{
   double  y1,y2,y3,area;
 FloatArray *b;
 
 b = new FloatArray (3);
 
   y1 = this -> giveNode(1) -> giveCoordinate(2);
   y2 = this -> giveNode(2) -> giveCoordinate(2);
   y3 = this -> giveNode(3) -> giveCoordinate(2);
 
 area = this -> giveArea();
 
   b->at(1)=(y2-y3)/2.0/area;
   b->at(2)=(y3-y1)/2.0/area;
   b->at(3)=(y1-y2)/2.0/area;
 
 return b;
}


FloatArray* LTrElementPPDE :: GivecCoeff ()
//
// Returns coefficients of partial derivatives of shape functions
// with respect to y
//
{
   double  x1,x2,x3,area;
 FloatArray *c;
 
 c = new FloatArray (3);

   x1 = this -> giveNode(1) -> giveCoordinate(1);
   x2 = this -> giveNode(2) -> giveCoordinate(1);
   x3 = this -> giveNode(3) -> giveCoordinate(1);
 
 area = this -> giveArea();
 
   c->at(1)=(x3-x2)/2.0/area;
   c->at(2)=(x1-x3)/2.0/area;
   c->at(3)=(x2-x1)/2.0/area;
 
 return c;
}


void
LTrElementPPDE :: computeBmatrixAt (GaussPoint *aGaussPoint, FloatMatrix& answer)
   // Returns the [2x3]  matrix {B} of the receiver, eva-
   // luated at aGaussPoint.
{
 // FloatMatrix *answer ;
 // Node *node1, *node2,*node3;
 // double x1,x2,x3,y1,y2,y3;
  // double area;
 FloatArray   *b,*c;
 
 b = this -> GivebCoeff();
 c = this -> GivecCoeff();
 // area = this->giveArea();

 // answer = new FloatMatrix(2,3) ;
 answer.resize (2,3);
 answer.zero();

 answer.at(1,1) = b->at(1);
 answer.at(1,2) = b->at(2);
 answer.at(1,3) = b->at(3);

 answer.at(2,1) = c->at(1);
 answer.at(2,2) = c->at(2);
 answer.at(2,3) = c->at(3);

 return  ;
}

   
void  LTrElementPPDE :: computeGaussPoints ()
   // Sets up the array containing the four Gauss points of the receiver.
{
   numberOfIntegrationRules = 1 ;
  integrationRulesArray = new IntegrationRule*[1];
  integrationRulesArray[0] = new GaussIntegrationRule (1,domain, 1, 2);
  integrationRulesArray[0]->setUpIntegrationPoints (_Triangle, numberOfGaussPoints, this, _2dHeat);

}


void
LTrElementPPDE :: computeNmatrixAt (GaussPoint* aGaussPoint, FloatMatrix& answer) 
   // Returns the displacement interpolation matrix {N} of the receiver, eva-
   // luated at aGaussPoint.
{

   Node *node1, *node2,*node3;
   double x1,x2,x3,y1,y2,y3,area;
   double l1,l2,l3;
   double ksi,eta,x,y;
   // FloatMatrix* answer ;

   ksi = aGaussPoint -> giveCoordinate(1) ;
   eta = aGaussPoint -> giveCoordinate(2) ;

   node1 = this -> giveNode(1) ;
   node2 = this -> giveNode(2) ;
   node3 = this -> giveNode(3) ;

   x1 = node1 -> giveCoordinate(1) ;
   x2 = node2 -> giveCoordinate(1) ;
   x3 = node3 -> giveCoordinate(1) ;

   y1 = node1 -> giveCoordinate(2) ;
   y2 = node2 -> giveCoordinate(2) ;
   y3 = node3 -> giveCoordinate(2) ;

   // transform "right-angled triangle" coordinates ksi, eta to real one

   x = x1 + ksi*(x2-x1)+eta*(x3-x1) ;
   y = y1 + ksi*(y2-y1)+eta*(y3-y1) ;
   area = 0.5*(x2*y3+x1*y2+y1*x3-x2*y1-x3*y2-x1*y3) ;

   l1 = ((x2*y3-x3*y2)+(y2-y3)*x+(x3-x2)*y)/2./area ;
   l2 = ((x3*y1-x1*y3)+(y3-y1)*x+(x1-x3)*y)/2./area ;
   l3 = ((x1*y2-x2*y1)+(y1-y2)*x+(x2-x1)*y)/2./area ;

   answer.resize(1,3) ;
  answer.zero();
  
   answer.at(1,1) = l1 ;
   answer.at(1,2) = l2 ;
   answer.at(1,3) = l3 ;

   return  ;
}


double  LTrElementPPDE :: computeVolumeAround (GaussPoint* aGaussPoint)
   // Returns the portion of the receiver which is attached to aGaussPoint.
{

  double area,weight ;

  weight  = aGaussPoint -> giveWeight() ;
  area    = this -> giveArea();

  return 2.0*area * weight * this -> giveCrossSection()->give('t');
}


IRResultType
LTrElementPPDE :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

  this->PPdeElement :: initializeFrom (ir);
  numberOfGaussPoints = 1;
 IR_GIVE_OPTIONAL_FIELD (ir, numberOfGaussPoints, IFT_LTrElementPPDE_nip, "nip"); // Macro

  if ( numberOfGaussPoints != 1) numberOfGaussPoints = 1;
  // set - up Gaussian integration points
  this -> computeGaussPoints();

  return IRRT_OK;
} 

void   LTrElementPPDE :: printOutputAt (FILE * file, TimeStep* stepN)
   // Performs end-of-step operations.
{
}


void
LTrElementPPDE ::   giveDofManDofIDMask  (int inode, EquationID, IntArray& answer) const {
// returns DofId mask array for inode element node.
// DofId mask array determines the dof ordering requsted from node.
// DofId mask array contains the DofID constants (defined in cltypes.h)
// describing physical meaning of particular DOFs.
 //IntArray* answer = new IntArray (1);
 answer.resize (1);

 answer.at(1) = T_f;

 return ;
}





#ifdef __OOFEG
#include "rcm2.h"
#define TR_LENGHT_REDUCT 0.3333

void LTrElementPPDE :: drawRawGeometry (oofegGraphicContext &gc)
{
  WCRec p[3];
  GraphicObj *go;

 if (!gc.testElementGraphicActivity(this)) return; 

  EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
  EASValsSetColor(gc.getElementColor());
  EASValsSetEdgeColor(gc.getElementEdgeColor());
  EASValsSetEdgeFlag(TRUE);
  EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
  p[0].x = (FPNum) this->giveNode(1)->giveCoordinate(1);
  p[0].y = (FPNum) this->giveNode(1)->giveCoordinate(2);
  p[0].z = 0.;
  p[1].x = (FPNum) this->giveNode(2)->giveCoordinate(1);
  p[1].y = (FPNum) this->giveNode(2)->giveCoordinate(2);
  p[1].z = 0.;
  p[2].x = (FPNum) this->giveNode(3)->giveCoordinate(1);
  p[2].y = (FPNum) this->giveNode(3)->giveCoordinate(2);
  p[2].z = 0.;
   
  go =  CreateTriangle3D(p);
  EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
  EGAttachObject(go, (EObjectP) this);
  EMAddGraphicsToModel(ESIModel(), go);
}

/*
void LTrElementPPDE :: drawInternalState (DrawMode mode)
//
// Draws internal state graphics representation
//
{
}
*/

#endif
 
