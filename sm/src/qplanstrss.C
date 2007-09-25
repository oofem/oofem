/* $Header: /home/cvs/bp/oofem/sm/src/qplanstrss.C,v 1.3.4.1 2004/04/05 15:19:47 bp Exp $ */
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

//   file QPLANSTRSS.CC

#include "qplanstrss.h"
#include "node.h"
#include "material.h"
#include "structuralcrosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "engngm.h"
#ifndef __MAKEDEPEND
#include <math.h>
#include <stdio.h>
#endif

QPlaneStress2d :: QPlaneStress2d (int n, Domain* aDomain)
                : StructuralElement (n,aDomain)
// Constructor.
{
   numberOfDofMans  = 8 ;
}

FloatArray*
QPlaneStress2d :: GiveDerivativeKsi (double ksi,double eta)
{
 FloatArray *n;

 n = new FloatArray(8);

   n->at(1) =  0.25*(1.+eta)*( 2.0*ksi+eta);
   n->at(2) = -0.25*(1.+eta)*(-2.0*ksi+eta);
   n->at(3) = -0.25*(1.-eta)*(-2.0*ksi-eta);
   n->at(4) =  0.25*(1.-eta)*( 2.0*ksi-eta);
   n->at(5) = -ksi*(1.+eta);
   n->at(6) = -0.5*(1.-eta*eta);
   n->at(7) = -ksi*(1.-eta);
   n->at(8) =  0.5*(1.-eta*eta);

 return n;
}

FloatArray*
QPlaneStress2d :: GiveDerivativeEta (double ksi,double eta)
{
 FloatArray *n;

 n = new FloatArray(8);

   n->at(1) =  0.25*(1.+ksi)*( 2.0*eta+ksi);
   n->at(2) =  0.25*(1.-ksi)*( 2.0*eta-ksi);
   n->at(3) = -0.25*(1.-ksi)*(-2.0*eta-ksi);
   n->at(4) = -0.25*(1.+ksi)*(-2.0*eta+ksi);
   n->at(5) =  0.5*(1.-ksi*ksi);
   n->at(6) = -eta*(1.-ksi);
   n->at(7) = -0.5*(1.-ksi*ksi);
   n->at(8) = -eta*(1.+ksi);

 return n;
}

void
QPlaneStress2d :: computeBmatrixAt (GaussPoint *aGaussPoint, FloatMatrix& answer, int li, int ui)
// Returns the [3x16] strain-displacement matrix {B} of the receiver,
// evaluated at aGaussPoint.
{
 int         i;
 FloatMatrix jacMtrx, inv ;
 FloatArray  *nx,*ny;
 double      ksi,eta;
 
 ksi = aGaussPoint -> giveCoordinate(1);
 eta = aGaussPoint -> giveCoordinate(2);
 
 nx = GiveDerivativeKsi(ksi,eta);
 ny = GiveDerivativeEta(ksi,eta);
 
 this->computeJacobianMatrixAt(jacMtrx, aGaussPoint);
 inv.beInverseOf (jacMtrx);
 
 // gm = new FloatMatrix(3,16);
 answer.resize (3, 16);
 answer.zero();

 for (i=1;i<=8;i++){
  answer.at(1,2*i-1)=nx->at(i)*inv.at(1,1)+ny->at(i)*inv.at(1,2);
  answer.at(2,2*i-0)=nx->at(i)*inv.at(2,1)+ny->at(i)*inv.at(2,2);
  answer.at(3,2*i-1)=nx->at(i)*inv.at(2,1)+ny->at(i)*inv.at(2,2);
  answer.at(3,2*i-0)=nx->at(i)*inv.at(1,1)+ny->at(i)*inv.at(1,2);
 }
 delete nx;  delete ny;
 return ;
}

void
QPlaneStress2d :: computeNmatrixAt (GaussPoint* aGaussPoint, FloatMatrix& answer) 
// Returns the displacement interpolation matrix {N} of the receiver,
// evaluated at aGaussPoint.
{
 int          i;
 double       ksi,eta;
 // FloatArray   *n;
 // FloatMatrix  *answer;

 ksi = aGaussPoint -> giveCoordinate(1);
 eta = aGaussPoint -> giveCoordinate(2);
 
 //n = new FloatArray (8);
 FloatArray n (8);

   n.at(1) = (1. + ksi) * (1. + eta) * 0.25 * ( ksi + eta - 1.) ;
   n.at(2) = (1. - ksi) * (1. + eta) * 0.25 * (-ksi + eta - 1.) ;
   n.at(3) = (1. - ksi) * (1. - eta) * 0.25 * (-ksi - eta - 1.) ;
   n.at(4) = (1. + ksi) * (1. - eta) * 0.25 * ( ksi - eta - 1.) ;
   n.at(5) = 0.5 * (1.-ksi*ksi) * (1.+eta);
   n.at(6) = 0.5 * (1.-ksi) * (1.-eta*eta);
   n.at(7) = 0.5 * (1.-ksi*ksi) * (1.-eta);
   n.at(8) = 0.5 * (1.+ksi) * (1.-eta*eta);

   //answer = new FloatMatrix(2,16);
 answer.resize (2,16);
 answer.zero();

 for (i=1;i<=8;i++){
  answer.at(1,2*i-1)=n.at(i);
  answer.at(2,2*i-0)=n.at(i);
 }
 return ;
}

void
QPlaneStress2d :: computeJacobianMatrixAt (FloatMatrix& answer, GaussPoint* aGaussPoint)
// Returns the jacobian matrix  J (x,y)/(ksi,eta)  of the receiver.
// Computes it if it does not exist yet.
{
 int         i;
 double      ksi,eta,x,y;
 FloatArray  *nx,*ny;
 
 answer.resize(2,2); answer.zero();

   ksi = aGaussPoint -> giveCoordinate (1);
   eta = aGaussPoint -> giveCoordinate (2);

 nx = this -> GiveDerivativeKsi (ksi,eta);
 ny = this -> GiveDerivativeEta (ksi,eta);

 for (i=1;i<=8;i++){
  x = this -> giveNode(i) -> giveCoordinate (1);
  y = this -> giveNode(i) -> giveCoordinate (2);

  answer.at(1,1)+=nx->at(i)*x;
  answer.at(1,2)+=nx->at(i)*y;
  answer.at(2,1)+=ny->at(i)*x;
  answer.at(2,2)+=ny->at(i)*y;
 }
 delete nx;  delete ny;
}

IRResultType
QPlaneStress2d :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

 this->Element :: initializeFrom (ir);
 numberOfGaussPoints = 4;
 IR_GIVE_OPTIONAL_FIELD (ir, numberOfGaussPoints, IFT_QPlaneStress2d_nip, "nip"); // Macro
 
 if (!((numberOfGaussPoints == 1) ||
   (numberOfGaussPoints == 4) ||
   (numberOfGaussPoints == 9) ||
   (numberOfGaussPoints == 16)))
  numberOfGaussPoints = 4;
 
 this -> computeGaussPoints();
 return IRRT_OK;
}

void
QPlaneStress2d :: computeGaussPoints ()
// Sets up the array containing the four Gauss points of the receiver.
{
  numberOfIntegrationRules = 1 ;
  integrationRulesArray = new IntegrationRule*[1];
  integrationRulesArray[0] = new GaussIntegrationRule (1,this, 1, 3);
  integrationRulesArray[0]->setUpIntegrationPoints (_Square, numberOfGaussPoints, _PlaneStress);
}

double
QPlaneStress2d :: computeVolumeAround (GaussPoint* aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
   double       determinant,weight,thickness,volume;
  FloatMatrix jacMtrx;

  this->computeJacobianMatrixAt(jacMtrx, aGaussPoint);
   determinant = fabs (jacMtrx.giveDeterminant());
   weight      = aGaussPoint -> giveWeight();
   thickness   = this -> giveCrossSection() -> give('t');
   volume      = determinant * weight * thickness ;

   return volume;
}


void
QPlaneStress2d ::   giveDofManDofIDMask  (int inode, EquationID, IntArray& answer) const {
// returns DofId mask array for inode element node.
// DofId mask array determines the dof ordering requsted from node.
// DofId mask array contains the DofID constants (defined in cltypes.h)
// describing physical meaning of particular DOFs.
 //IntArray* answer = new IntArray (2);
 answer.resize (2);

 answer.at(1) = D_u;
 answer.at(2) = D_v;

 return ;
}


int
QPlaneStress2d :: computeGlobalCoordinates (FloatArray& answer, const FloatArray& lcoords) 
{
 int i;
 double ksi, eta;
 FloatArray n(8);

 ksi = lcoords.at(1) ;
 eta = lcoords.at(2) ;
 
   n.at(1) = (1. + ksi) * (1. + eta) * 0.25 * ( ksi + eta - 1.) ;
   n.at(2) = (1. - ksi) * (1. + eta) * 0.25 * (-ksi + eta - 1.) ;
   n.at(3) = (1. - ksi) * (1. - eta) * 0.25 * (-ksi - eta - 1.) ;
   n.at(4) = (1. + ksi) * (1. - eta) * 0.25 * ( ksi - eta - 1.) ;
   n.at(5) = 0.5 * (1.-ksi*ksi) * (1.+eta);
   n.at(6) = 0.5 * (1.-ksi) * (1.-eta*eta);
   n.at(7) = 0.5 * (1.-ksi*ksi) * (1.-eta);
   n.at(8) = 0.5 * (1.+ksi) * (1.-eta*eta);

 
 answer.resize(2);
 answer.zero();
 for (i=1; i<=8; i++) {
  answer.at(1) += n.at(i)*this->giveNode(i)->giveCoordinate(1);
  answer.at(2) += n.at(i)*this->giveNode(i)->giveCoordinate(2);
 }

 return 1;
}


#ifdef __OOFEG
void QPlaneStress2d :: drawRawGeometry (oofegGraphicContext &gc)
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


void QPlaneStress2d :: drawDeformedGeometry (oofegGraphicContext &gc, UnknownType type)
{
  WCRec p[4];
  GraphicObj *go;
  TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
  double defScale = gc.getDefScale();

  if (!gc.testElementGraphicActivity(this)) return; 

  EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
  EASValsSetColor(gc.getDeformedElementColor());
  EASValsSetEdgeColor(gc.getElementEdgeColor());
  EASValsSetEdgeFlag(TRUE);
  EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
  EASValsSetFillStyle (FILL_HOLLOW);
  p[0].x = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
  p[0].y = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
  p[0].z = 0.;
  p[1].x = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
  p[1].y = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
  p[1].z = 0.;
  p[2].x = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
  p[2].y = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
  p[2].z = 0.;
  p[3].x = (FPNum) this->giveNode(4)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
  p[3].y = (FPNum) this->giveNode(4)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
  p[3].z = 0.;
   
  go =  CreateQuad3D(p);
  EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
  EMAddGraphicsToModel(ESIModel(), go);
}

#endif
