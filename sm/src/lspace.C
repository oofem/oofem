/* $Header: /home/cvs/bp/oofem/sm/src/lspace.C,v 1.8.4.1 2004/04/05 15:19:47 bp Exp $ */
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

//   file LSPACE.CC

#include "lspace.h"
#include "node.h"
#include "material.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "cltypes.h"
#include "structuralms.h"
#include "mathfem.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

#ifdef __OOFEG
#include "engngm.h"
#include "oofeggraphiccontext.h"
#include "oofegutils.h"
#include "conTable.h"
#endif


FEI3dHexaLin LSpace :: interpolation;

LSpace :: LSpace (int n, Domain* aDomain) : StructuralElement (n,aDomain), ZZNodalRecoveryModelInterface(), 
SPRNodalRecoveryModelInterface(), SpatialLocalizerInterface()
     , EIPrimaryUnknownMapperInterface(), HuertaErrorEstimatorInterface(), HuertaRemeshingCriteriaInterface()
     // Constructor.
{
  numberOfDofMans  = 8 ;
  numberOfGaussPoints = 8;
  // this -> computeGaussPoints() ; => moved to instanciateYourself();
}

Interface* 
LSpace :: giveInterface (InterfaceType interface)
{
 if (interface == ZZNodalRecoveryModelInterfaceType) return (ZZNodalRecoveryModelInterface*) this;
 else if (interface == SPRNodalRecoveryModelInterfaceType) return (SPRNodalRecoveryModelInterface*) this;
 else if (interface == NodalAveragingRecoveryModelInterfaceType) return (NodalAveragingRecoveryModelInterface*) this;
 else if (interface == SpatialLocalizerInterfaceType) return (SpatialLocalizerInterface*) this;
 else if (interface == EIPrimaryUnknownMapperInterfaceType) return (EIPrimaryUnknownMapperInterface*) this;
 else if (interface == HuertaErrorEstimatorInterfaceType) return (HuertaErrorEstimatorInterface*)this;
 else if (interface == HuertaRemeshingCriteriaInterfaceType) return (HuertaRemeshingCriteriaInterface*) this;
 return NULL;
}

void
LSpace :: computeBmatrixAt (GaussPoint *aGaussPoint, FloatMatrix& answer, int li, int ui)
   // Returns the [6x24] strain-displacement matrix {B} of the receiver, eva-
   // luated at aGaussPoint.
   // B matrix  -  6 rows : epsilon-X, epsilon-Y, epsilon-Z, gamma-YZ, gamma-ZX, gamma-XY  :
{
 int i;
 FloatMatrix  dnx ;

 this->interpolation.evaldNdx (dnx, this->giveDomain(), dofManArray, *aGaussPoint->giveCoordinates(), 0.0);

 answer.resize(6,24);
 answer.zero();
 
 for (i=1;i<=8;i++){
  answer.at(1,3*i-2) = dnx.at(i,1);
  answer.at(2,3*i-1) = dnx.at(i,2);
  answer.at(3,3*i-0) = dnx.at(i,3);
 }

 for (i=1;i<=8;i++){
  answer.at(4,3*i-1) = dnx.at(i,3);
  answer.at(4,3*i-0) = dnx.at(i,2);

  answer.at(5,3*i-2) = dnx.at(i,3);
  answer.at(5,3*i-0) = dnx.at(i,1);

  answer.at(6,3*i-2) = dnx.at(i,2);
  answer.at(6,3*i-1) = dnx.at(i,1);
 }

 return ;
}

   

void  LSpace :: computeGaussPoints ()
   // Sets up the array containing the four Gauss points of the receiver.
{
   numberOfIntegrationRules = 1 ;
  integrationRulesArray = new IntegrationRule*;
  integrationRulesArray[0] = new GaussIntegrationRule (1,domain,1, 6);
  integrationRulesArray[0]->setUpIntegrationPoints (_Cube, numberOfGaussPoints, this,  _3dMat);

}


void
LSpace :: computeNmatrixAt (GaussPoint* aGaussPoint, FloatMatrix& answer) 
   // Returns the displacement interpolation matrix {N} of the receiver, eva-
   // luated at aGaussPoint.
{
 int i;
 FloatArray n(8);

 answer.resize(3,24);
 answer.zero();
 this->interpolation.evalN (n, *aGaussPoint->giveCoordinates(), 0.0); 

 for (i=1;i<=8;i++){
  answer.at(1,3*i-2) = n.at(i);
  answer.at(2,3*i-1) = n.at(i);
  answer.at(3,3*i-0) = n.at(i);
 }
 
 return ;
}


double  LSpace :: computeVolumeAround (GaussPoint* aGaussPoint)
   // Returns the portion of the receiver which is attached to aGaussPoint.
{
 double       determinant,weight,volume;
 determinant = fabs(this->interpolation.giveTransformationJacobian (domain, dofManArray, 
                                   *aGaussPoint->giveCoordinates(), 0.0)); 
 
 
 weight      = aGaussPoint -> giveWeight();
 volume      = determinant * weight;
 
 return volume;
}


IRResultType
LSpace :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

  this->StructuralElement :: initializeFrom (ir);
 numberOfGaussPoints = 8;
 IR_GIVE_OPTIONAL_FIELD (ir, numberOfGaussPoints, IFT_LSpace_nip, "nip"); // Macro

  if (!((numberOfGaussPoints == 8) || (numberOfGaussPoints == 27))) numberOfGaussPoints = 8; 
  // set - up Gaussian integration points
  this -> computeGaussPoints();

  return IRRT_OK;
} 

void
LSpace ::   giveDofManDofIDMask  (int inode, EquationID ut, IntArray& answer) const {
// returns DofId mask array for inode element node.
// DofId mask array determines the dof ordering requsted from node.
// DofId mask array contains the DofID constants (defined in cltypes.h)
// describing physical meaning of particular DOFs.
 //IntArray* answer = new IntArray (3);
 answer.resize (3);

 answer.at(1) = D_u;
 answer.at(2) = D_v;
 answer.at(3) = D_w;

 return ;
}


double
LSpace::giveCharacteristicLenght(GaussPoint* gp, const FloatArray &normalToCrackPlane) 
{
 double factor = __OOFEM_POW((double) this->numberOfGaussPoints, 1./3.);
 return this -> giveLenghtInDir (normalToCrackPlane)/factor;
}




int
LSpace :: computeGlobalCoordinates (FloatArray& answer, const FloatArray& lcoords) 
{
 this->interpolation.local2global (answer, domain, dofManArray, lcoords, 0.0); 
 return 1;
}


int 
LSpace :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
 if ((type == IST_StressTensor) || (type == IST_StrainTensor)) return 6;

 GaussPoint *gp = integrationRulesArray[0]-> getIntegrationPoint(0) ;
 return this->giveIPValueSize (type, gp);
}


void
LSpace :: ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx  
(FloatMatrix& answer, GaussPoint* aGaussPoint, InternalStateType type)
{
  // evaluates N matrix (interpolation estimated stress matrix)
  // according to Zienkiewicz & Zhu paper
  // N(nsigma, nsigma*nnodes)
  // Definition : sigmaVector = N * nodalSigmaVector
 FloatArray n;
 this->interpolation.evalN (n, *aGaussPoint->giveCoordinates(), 0.0); 

 ////if (type == StressVector) answer.resize(6,48) ;
 if (this->giveIPValueSize(type, aGaussPoint)) answer.resize(1,8) ;
 else return;

 answer.at(1,1)  = n.at(1);
 answer.at(1,2)  = n.at(2);
 answer.at(1,3)  = n.at(3);
 answer.at(1,4)  = n.at(4);
 answer.at(1,5)  = n.at(5);
 answer.at(1,6)  = n.at(6);
 answer.at(1,7)  = n.at(7);
 answer.at(1,8)  = n.at(8);

/******
 answer.zero();
 for (i=1; i<=6; i++) 
  for (j=1; j<=8; j++) 
   answer.at(i,i+(j-1)*6) = n.at(j);
*******/
 return ;
}

void 
LSpace::SPRNodalRecoveryMI_giveSPRAssemblyPoints (IntArray& pap)
{
 int i;

 pap.resize(8);
 for (i=1; i<=8; i++) pap.at(i) = this->giveNode(i)->giveNumber();
}

void
LSpace::SPRNodalRecoveryMI_giveDofMansDeterminedByPatch (IntArray& answer, int pap)
{
 int i, found = 0;
 answer.resize(1);

 for (i=1; i<=8; i++) if (this->giveNode(i)->giveNumber() == pap) found = 1;
 
 if (found) answer.at(1) = pap;
 else _error ("SPRNodalRecoveryMI_giveDofMansDeterminedByPatch: node unknown");
}

int
LSpace::SPRNodalRecoveryMI_giveNumberOfIP ()
{ return this->giveDefaultIntegrationRulePtr()->getNumberOfIntegrationPoints(); }


void
LSpace::SPRNodalRecoveryMI_computeIPGlobalCoordinates (FloatArray& coords, GaussPoint* gp)
{
 this->computeGlobalCoordinates(coords, *gp->giveCoordinates());
}

SPRPatchType
LSpace::SPRNodalRecoveryMI_givePatchType()
{
 return SPRPatchType_3dBiLin;
}


void 
LSpace::NodalAveragingRecoveryMI_computeNodalValue (FloatArray& answer, int node,
                          InternalStateType type, TimeStep* tStep)
{
  int i,j;
  double x1=0.0, x2=0.0, x3=0.0, y = 0.0;
  GaussPoint *gp;
  //StructuralMaterialStatus* status;
  IntegrationRule* iRule = integrationRulesArray[0];
  FloatMatrix A(4,4);
  FloatMatrix b, r;
  FloatArray  val, *coord ;
  double   u,v,w;

 int size = NodalAveragingRecoveryMI_giveDofManRecordSize(type);
 if (size) {
  answer.resize (size);
  b.resize(4,size); r.resize(4,size);
  A.zero(); r.zero();
  for (i=0 ; i< iRule->getNumberOfIntegrationPoints() ; i++) {
   gp = iRule-> getIntegrationPoint(i) ;
   giveIPValue (val,gp, type, tStep);
   
   coord = gp -> giveCoordinates() ;
   u = coord -> at(1) ;
   v = coord -> at(2) ;
   w = coord -> at(3) ;
   
   //status = (StructuralMaterialStatus*) this->giveMaterial()->giveStatus(gp);
   
   A.at(1,1) += 1;     A.at(1,2) += u;    A.at(1,3) += v;    A.at(1,4) += w;
   A.at(2,1) += u;     A.at(2,2) += u*u;  A.at(2,3) += u*v;  A.at(2,4) += u*w;
   A.at(3,1) += v;     A.at(3,2) += v*u;  A.at(3,3) += v*v;  A.at(3,4) += v*w;
   A.at(4,1) += w;     A.at(4,2) += w*u;  A.at(4,3) += w*v;  A.at(4,4) += w*w;
   
   for (j=1; j<=size; j++) {
    y = val.at(j) ;
    r.at(1,j) += y;       
    r.at(2,j) += y*u;    
    r.at(3,j) += y*v;    
    r.at(4,j) += y*w;
   }
  }
  
  A.solveForRhs (r,b);
  
  switch (node) {
  case 1:
   x1 =  1.0;
   x2 =  1.0;
   x3 =  1.0;
   break;
  case 2:
   x1 = -1.0;
   x2 =  1.0;
   x3 =  1.0;
   break;
  case 3:
   x1 = -1.0;
   x2 = -1.0;
   x3 =  1.0;
   break;
  case 4:
   x1 =  1.0;
   x2 = -1.0;
   x3 =  1.0;
   break;
  case 5:
   x1 =  1.0;
   x2 =  1.0;
   x3 = -1.0;
   break;
  case 6:
   x1 = -1.0;
   x2 =  1.0;
   x3 = -1.0;
   break;
  case 7:
   x1 = -1.0;
   x2 = -1.0;
   x3 = -1.0;
   break;
  case 8:
   x1 =  1.0;
   x2 = -1.0;
   x3 = -1.0;
   break;
  default:
   _error ("LTRSpace ::giveInternalStateAtNode: unsupported node");
  }
  
  for (j=1; j<=size; j++) {
   answer.at(j) = b.at(1,j) + x1*b.at(2,j) * x2*b.at(3,j) * x3*b.at(4,j);
  }
  
  return;
  
 } else answer.resize(0);
}

void 
LSpace::NodalAveragingRecoveryMI_computeSideValue (FloatArray& answer, int side,
                          InternalStateType type, TimeStep* tStep)
{
 answer.resize(0);
}


int
LSpace::computeLocalCoordinates (FloatArray& answer, const FloatArray& coords)
{
 return this->interpolation.global2local (answer, domain, dofManArray, coords, 0.0); 
}


int 
LSpace::SpatialLocalizerI_containsPoint (const FloatArray& coords) 
{
 FloatArray lcoords;
 return this->computeLocalCoordinates (lcoords, coords);
}

double 
LSpace::SpatialLocalizerI_giveDistanceFromParametricCenter (const FloatArray& coords)
{
 FloatArray lcoords(3), gcoords;
 double dist;
 int size, gsize;

 lcoords.at(1) = lcoords.at(2) = lcoords.at(3) = 0.0;
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



int
LSpace::EIPrimaryUnknownMI_computePrimaryUnknownVectorAt (ValueModeType mode,
                             TimeStep* stepN, const FloatArray& coords,
                             FloatArray& answer)
{
 FloatArray lcoords, u;
 FloatMatrix n;
 FloatArray ni (8);
 int result;

 result = this->computeLocalCoordinates (lcoords, coords);

 this->interpolation.evalN (ni, lcoords, 0.0); 

 n.resize (3,24);
 n.zero();
 
 n.at(1,1)  = n.at(2,2)  = n.at(3,3)  = ni.at(1) ;
 n.at(1,4)  = n.at(2,5)  = n.at(3,6)  = ni.at(2) ;
 n.at(1,7)  = n.at(2,8)  = n.at(3,9)  = ni.at(3) ;
 n.at(1,10) = n.at(2,11) = n.at(3,12) = ni.at(4) ;
 n.at(1,13) = n.at(2,14) = n.at(3,15) = ni.at(5) ;
 n.at(1,16) = n.at(2,17) = n.at(3,18) = ni.at(6) ;
 n.at(1,19) = n.at(2,20) = n.at(3,21) = ni.at(7) ;
 n.at(1,22) = n.at(2,23) = n.at(3,24) = ni.at(8) ;

 this->computeVectorOf (EID_MomentumBalance, mode, stepN, u);
 answer.beProductOf (n,u);

 return result;
}

void 
LSpace::EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID (IntArray& answer)
{
  giveDofManDofIDMask  (1, EID_MomentumBalance, answer);
}

void
LSpace::HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId, 
                             IntArray &localNodeIdArray, IntArray &globalNodeIdArray, 
                             HuertaErrorEstimatorInterface::SetupMode sMode, TimeStep* tStep,
                             int &localNodeId, int &localElemId, int &localBcId,
                             IntArray &controlNode, IntArray &controlDof,
                             HuertaErrorEstimator::AnalysisMode aMode)
{
 Element *element = this->HuertaErrorEstimatorI_giveElement();
 FloatArray *corner[8], midSide[12], midFace[6], midNode;
 double x = 0.0, y = 0.0, z = 0.0;
 int inode, nodes = 8, iside, sides = 12, iface, faces = 6, nd, nd1, nd2;

 static int sideNode[12][2] = {{1, 2}, {2, 3}, {3, 4}, {4, 1},         // bottom
                {5, 6}, {6, 7}, {7, 8}, {8, 5},         // top
                {1, 5}, {2, 6}, {3, 7}, {4, 8}};        // vertical
 static int faceNode[6][4] = {{1, 2, 3, 4}, {5, 6, 7, 8},              // bottom, top
                {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 4, 8, 7}, {4, 1, 5, 8}};    // side

/* ordering of side and face hexa nodes must be compatible with refinedElement connectivity ordering;
  generally the ordering is: corner side side face side face face center;
  however the concrete ordering is element dependent (see refineMeshGlobally source if in doubts) */

 static int hexaSideNode[8][3] = {{1, 4, 9}, {2, 1, 10}, {3, 2, 11}, {4, 3, 12}, 
                  {8, 5, 9}, {5, 6, 10}, {6, 7, 11}, {7, 8, 12}};
 static int hexaFaceNode[8][3] = {{1, 3, 6}, {1, 4, 3}, {1, 5, 4}, {1, 6, 5}, 
                  {2, 6, 3}, {2, 3, 4}, {2, 4, 5}, {2, 5, 6}};

 if(sMode == HuertaErrorEstimatorInterface::NodeMode || 
   (sMode == HuertaErrorEstimatorInterface::BCMode && aMode == HuertaErrorEstimator::HEE_linear)){
  for(inode = 0; inode < nodes; inode++){
   corner[inode] = element -> giveNode(inode + 1) -> giveCoordinates();

   x += corner[inode] -> at(1);
   y += corner[inode] -> at(2);
   z += corner[inode] -> at(3);
  }
  
  for(iside = 0; iside < sides; iside++){
   midSide[iside].resize(3);

   nd1 = sideNode[iside][0] - 1;
   nd2 = sideNode[iside][1] - 1;

   midSide[iside].at(1) = (corner[nd1] -> at(1) + corner[nd2] -> at(1)) / 2.0;
   midSide[iside].at(2) = (corner[nd1] -> at(2) + corner[nd2] -> at(2)) / 2.0;
   midSide[iside].at(3) = (corner[nd1] -> at(3) + corner[nd2] -> at(3)) / 2.0;
  }

  midNode.resize(3);

  midNode.at(1) = x / nodes;
  midNode.at(2) = y / nodes;
  midNode.at(3) = z / nodes;

  for(iface = 0; iface < faces; iface++){
   x = y = z = 0.0;
   for(inode = 0; inode < 4; inode++){
    nd = faceNode[iface][inode] - 1;
    x += corner[nd] -> at(1);
    y += corner[nd] -> at(2);
    z += corner[nd] -> at(3);
   }

   midFace[iface].resize(3);

   midFace[iface].at(1) = x / 4;
   midFace[iface].at(2) = y / 4;
   midFace[iface].at(3) = z / 4;
  }
 }

 this -> setupRefinedElementProblem3D (element, refinedElement, level, nodeId, localNodeIdArray, globalNodeIdArray, 
                    sMode, tStep, nodes, corner, midSide, midFace, midNode, 
                    localNodeId, localElemId, localBcId, hexaSideNode, hexaFaceNode,
                    controlNode, controlDof, aMode, "LSpace");
}


double
LSpace::HuertaRemeshingCriteriaI_giveCharacteristicSize () {

 int i;
 IntegrationRule* iRule;
 GaussPoint* gp ;
 double volume = 0.0;

 iRule = integrationRulesArray[giveDefaultIntegrationRule()];
 for (i=0 ; i < iRule->getNumberOfIntegrationPoints() ; i++) {
  gp  = iRule->getIntegrationPoint(i) ;
  volume += this -> computeVolumeAround(gp) ;
 } 

 return __OOFEM_POW(volume, 1./3.);
}




#ifdef __OOFEG
#include "rcm2.h"
#define TR_LENGHT_REDUCT 0.3333

void LSpace :: drawRawGeometry (oofegGraphicContext &gc)
{
  int i;
  WCRec p[8];
  GraphicObj *go;

 if (!gc.testElementGraphicActivity(this)) return; 

  EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
  EASValsSetColor(gc.getElementColor());
  EASValsSetEdgeColor(gc.getElementEdgeColor());
  EASValsSetEdgeFlag(TRUE);
  EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
 EASValsSetFillStyle (FILL_SOLID);
  for (i=0; i< 8; i++) {
    p[i].x = (FPNum) this->giveNode(i+1)->giveCoordinate(1);
    p[i].y = (FPNum) this->giveNode(i+1)->giveCoordinate(2);
    p[i].z = (FPNum) this->giveNode(i+1)->giveCoordinate(3);
  }
   
  go =  CreateHexahedron(p);
  EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
  EGAttachObject(go, (EObjectP) this);
  EMAddGraphicsToModel(ESIModel(), go);

  FloatArray c(3); 
  c.at(1) = -1.0; c.at(2) = 0.0; c.at(3) = 0.0;
  this->drawTriad (c,4);
  c.at(1) = 1.0; c.at(2) = 0.0; c.at(3) = 0.0;
  this->drawTriad (c,6);
  c.at(1) = 0.0; c.at(2) = -1.0; c.at(3) = 0.0;
  this->drawTriad (c,5);
  c.at(1) = 0.0; c.at(2) =  1.0; c.at(3) = 0.0;
  this->drawTriad (c,3);
  c.at(1) = 0.0; c.at(2) = 0.0; c.at(3) = -1.0;
  this->drawTriad (c,2);
  c.at(1) = 0.0; c.at(2) = 0.0; c.at(3) =  1.0;
  this->drawTriad (c,1);
}

void LSpace :: drawTriad (FloatArray &coords, int isurf)
{
  FloatMatrix jm(3,3);
  FloatArray gc(3);
  GraphicObj *go;

  WCRec p[2];   // point
  double coeff = 1.0;
  int i, succ;
/*
   // version I
  this->interpolation.giveJacobianMatrixAt (jm, domain, nodeArray, coords); 
  // determine origin
  this->interpolation.local2global (gc, domain, nodeArray, coords, 0.0);
  // draw triad

*/

  // version II
  // determine surface center
  IntArray snodes(4);
  FloatArray h1(3),h2(3),nn(3),n(3);
  int j;
  this->interpolation.computeSurfaceMapping (snodes, dofManArray, isurf);
  for (i=1; i<=4; i++) 
    gc.add(domain->giveNode(snodes.at(i))->giveCoordinates());
  gc.times(1./4.);
  // determine "average normal"
  nn.zero();
  for (i=1; i<=4; i++) {
    j = (i)%4+1;
    h1 = *domain->giveNode(snodes.at(i))->giveCoordinates();
    h1.substract(gc);
    h2 = *domain->giveNode(snodes.at(j))->giveCoordinates();
    h2.substract(gc);
    n.beVectorProductOf(h1,h2);
    if (dotProduct(n,n,3) > 1.e-6) n.normalize();
    nn.add(n);
  }
  nn.times(1./4.); 
  if (dotProduct(nn,nn,3) < 1.e-6) return;
  nn.normalize();
  for (i=1; i<=3; i++) jm.at(i,3) = nn.at(i);
  // determine lcs of surface
  // local x axis in xy plane
  double test = fabs(fabs(nn.at(3)) -1.0);
  if (test < 1.e-5) {
    h1.at(1) = jm.at(1,1) = 1.0;
    h1.at(2) = jm.at(2,1) = 0.0;
  } else {
    h1.at(1) = jm.at(1,1) = jm.at(2,3);
    h1.at(2) = jm.at(2,1) = -jm.at(1,3);
  }
  h1.at(3) = jm.at(3,1) = 0.0;
  // local y axis perpendicular to local x,z axes
  h2.beVectorProductOf (nn,h1);
  for (i=1; i<=3; i++) jm.at(i,2) = h2.at(i);
  

  p[0].x = gc.at(1);
  p[0].y = gc.at(2);
  p[0].z = gc.at(3);
  for (i =1; i<=3; i++) {
    p[1].x = p[0].x + coeff*jm.at(1,i);
    p[1].y = p[0].y + coeff*jm.at(2,i);
    p[1].z = p[0].z + coeff*jm.at(3,i);
    
    if (i==1) EASValsSetColor(ColorGetPixelFromString("red", &succ));
    else if (i==2)  EASValsSetColor(ColorGetPixelFromString("green", &succ));
    else EASValsSetColor(ColorGetPixelFromString("blue", &succ));
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
  }
}


void LSpace :: drawDeformedGeometry (oofegGraphicContext &gc, UnknownType type)
{
  int i;
  WCRec p[8];
  GraphicObj *go;
  TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
 double defScale = gc.getDefScale();

 if (!gc.testElementGraphicActivity(this)) return; 

  EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
  EASValsSetColor(gc.getDeformedElementColor());
  EASValsSetEdgeColor(gc.getElementEdgeColor());
  EASValsSetEdgeFlag(TRUE);
  EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
 EASValsSetFillStyle (FILL_SOLID);
  for (i=0; i< 8; i++) {
    p[i].x = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
    p[i].y = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
    p[i].z = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(3,tStep,EID_MomentumBalance,defScale);
  }
   
  go =  CreateHexahedron (p);
  EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
  EMAddGraphicsToModel(ESIModel(), go);
}


void LSpace :: drawScalar   (oofegGraphicContext& context)
{
  int i, indx, result = 0;
  WCRec p[8];
  GraphicObj *tr;
  TimeStep* tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
  FloatArray v[8];
  double s[8], defScale=0.0;
  IntArray map;
  
  if (!context.testElementGraphicActivity(this)) return; 
  if (context.giveIntVarMode() == ISM_recovered) {
    for (i=1; i<=8; i++) 
      result+= this->giveInternalStateAtNode (v[i-1], context.giveIntVarType(), context.giveIntVarMode(), i, tStep);
    if (result !=8) return;
  } else if (context.giveIntVarMode() == ISM_local) {
    return;
  }
  this->giveIntVarCompFullIndx (map, context.giveIntVarType());
  if ((indx = map.at(context.giveIntVarIndx())) == 0) return;
  for (i=1; i<=8; i++) s[i-1]=v[i-1].at(indx);
  
  EASValsSetEdgeColor(context.getElementEdgeColor());
  EASValsSetEdgeFlag(TRUE);
  EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
  if (context.getScalarAlgo() == SA_ISO_SURF) {
    for (i=0; i<8; i++) {
      if (context.getInternalVarsDefGeoFlag()) {
        // use deformed geometry
        defScale = context.getDefScale();
        p[i].x = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
        p[i].y = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
        p[i].z = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(3,tStep,EID_MomentumBalance,defScale);
        
      } else {
        p[i].x = (FPNum) this->giveNode(i+1)->giveCoordinate(1);
        p[i].y = (FPNum) this->giveNode(i+1)->giveCoordinate(2);
        p[i].z = (FPNum) this->giveNode(i+1)->giveCoordinate(3);
      }
    }
    
    context.updateFringeTableMinMax (s, 8);
    tr = CreateHexahedronWD (p,s);
    EGWithMaskChangeAttributes(LAYER_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK, tr);
    EMAddGraphicsToModel(ESIModel(), tr);
  } 
  
}


/*
void LSpace :: drawInternalState (oofegGraphicContext &gc)
//
// Draws internal state graphics representation
//
{
  WCRec p[8];
  GraphicObj *tr;
  StructuralMaterial *mat = (StructuralMaterial*) this->giveMaterial();
  GaussPoint* gp;
  double v[8],ratio;
  int i, nPlastGp;
  IntegrationRule* iRule = integrationRulesArray[0];
 DrawMode mode = gc.getDrawMode();
 TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
 double defScale = gc.getDefScale();
 
 if (!gc.testElementGraphicActivity(this)) return; 

  // check for yield mode
  if (mode == yieldState) {
   // loop over available GPs
   nPlastGp = 0; 
   for (i=1 ; i<= iRule->getNumberOfIntegrationPoints() ; i++) {
    gp = iRule-> getIntegrationPoint(i) ;
    nPlastGp += (mat->giveStatusCharFlag(gp,ms_yield_flag) != 0);
   }
   if (nPlastGp == 0) return;
   // nPlastGp should contain number of yielded gp in element
   // good way should be select color accordingly
   ratio = nPlastGp / numberOfGaussPoints;
   EASValsSetLayer(OOFEG_YIELD_PATTERN_LAYER);
   for (i=0; i<8; i++) {
  if (gc.getInternalVarsDefGeoFlag()) {
   // use deformed geometry
   p[i].x = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
   p[i].y = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
   p[i].z = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(3,tStep,EID_MomentumBalance,defScale);
   
  } else {
     p[i].x = (FPNum) this->giveNode(i+1)->giveCoordinate(1);
     p[i].y = (FPNum) this->giveNode(i+1)->giveCoordinate(2);
     p[i].z = (FPNum) this->giveNode(i+1)->giveCoordinate(3);
  }
   }
   
   EASValsSetColor(gc.getYieldPlotColor(ratio));
   tr =  CreateHexahedron (p);
   EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | FILL_MASK | LAYER_MASK, tr);
   EMAddGraphicsToModel(ESIModel(), tr);
  }
  // check for valid stress-strain mode
  if (!((mode == sxForce) || (mode == syForce) || (mode == szForce) || 
 (mode == sxyForce) || (mode == syzForce) || (mode == szxForce))) return ;

 int result = 0;
  
  EASValsSetLayer(OOFEG_STRESS_CONTOUR_LAYER);
  for (i=0; i<8; i++) {
  if (gc.getInternalVarsDefGeoFlag()) {
   // use deformed geometry
   p[i].x = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
   p[i].y = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
   p[i].z = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(3,tStep,EID_MomentumBalance,defScale);
   
  } else {
   p[i].x = (FPNum) this->giveNode(i+1)->giveCoordinate(1);
   p[i].y = (FPNum) this->giveNode(i+1)->giveCoordinate(2);
      p[i].z = (FPNum) this->giveNode(i+1)->giveCoordinate(3);
  }
  result+= this->giveInternalStateAtNode (gc, i+1, &v[i]);
  }  

 if (result==8) {
  tr = CreateHexahedronWD (p,v);
  EGWithMaskChangeAttributes(LAYER_MASK, tr);
  EMAddGraphicsToModel(ESIModel(), tr);
 }
}
*/ 

#endif



void
LSpace :: computeEgdeNMatrixAt (FloatMatrix& answer, GaussPoint* aGaussPoint) 
{
  /*
  
  computes interpolation matrix for element edge.
  we assemble locally this matrix for only nonzero 
  shape functions. 
  (for example only two nonzero shape functions for 2 dofs are
  necessary for linear plane stress quad edge).
  These nonzero shape functions are then mapped to 
  global element functions.
  
  Using mapping technique will allow to assemble shape functions 
  without regarding particular side
  */
 
 FloatArray n(2);
 this->interpolation.edgeEvalN (n, *aGaussPoint->giveCoordinates(), 0.0); 

  answer.resize (3,6);
 answer.zero();
  
  answer.at(1,1) = n.at(1) ;
  answer.at(1,4) = n.at(2) ;
  answer.at(2,2) = n.at(1) ;
  answer.at(2,5) = n.at(2) ;
  answer.at(3,3) = n.at(1) ;
  answer.at(3,6) = n.at(2) ;
  
  return  ;
}

void
LSpace :: giveEdgeDofMapping (IntArray& answer, int iEdge) const
{
  /* 
   provides dof mapping of local edge dofs (only nonzero are taken into account)
   to global element dofs
  */

  answer.resize (6);
  if (iEdge == 1)  { // edge between nodes 1 2
  answer.at(1) = 1;
  answer.at(2) = 2;
  answer.at(3) = 3;
  answer.at(4) = 4;
  answer.at(5) = 5;
  answer.at(6) = 6;
  } else if (iEdge == 2) { // edge between nodes 2 3
  answer.at(1) = 4;
  answer.at(2) = 5;
  answer.at(3) = 6;
  answer.at(4) = 7;
  answer.at(5) = 8;
  answer.at(6) = 9;
  } else if (iEdge == 3) {// edge between nodes 3 4
  answer.at(1) = 7;
  answer.at(2) = 8;
  answer.at(3) = 9;
  answer.at(4) = 10;
  answer.at(5) = 11;
  answer.at(6) = 12;
  } else if (iEdge == 4) {// edge between nodes 4 1
  answer.at(1) = 10;
  answer.at(2) = 11;
  answer.at(3) = 12;
  answer.at(4) = 1;
  answer.at(5) = 2;
  answer.at(6) = 3;
  } else if (iEdge == 5) {// edge between nodes 1 5
  answer.at(1) = 1;
  answer.at(2) = 2;
  answer.at(3) = 3;
  answer.at(4) = 13;
  answer.at(5) = 14;
  answer.at(6) = 15;
  } else if (iEdge == 6) {// edge between nodes 2 6
  answer.at(1) = 4;
  answer.at(2) = 5;
  answer.at(3) = 6;
  answer.at(4) = 16;
  answer.at(5) = 17;
  answer.at(6) = 18;
  } else if (iEdge == 7) {// edge between nodes 3 7
  answer.at(1) = 7;
  answer.at(2) = 8;
  answer.at(3) = 9;
  answer.at(4) = 19;
  answer.at(5) = 20;
  answer.at(6) = 21;
  } else if (iEdge == 8) {// edge between nodes 4 8
  answer.at(1) = 10;
  answer.at(2) = 11;
  answer.at(3) = 12;
  answer.at(4) = 22;
  answer.at(5) = 23;
  answer.at(6) = 24;
  } else if (iEdge == 9) {// edge between nodes 5 6
  answer.at(1) = 13;
  answer.at(2) = 14;
  answer.at(3) = 15;
  answer.at(4) = 16;
  answer.at(5) = 17;
  answer.at(6) = 18;
  } else if (iEdge == 10) {// edge between nodes 6 7
  answer.at(1) = 16;
  answer.at(2) = 17;
  answer.at(3) = 18;
  answer.at(4) = 19;
  answer.at(5) = 22;
  answer.at(6) = 21;
  } else if (iEdge == 11) {// edge between nodes 7 8
  answer.at(1) = 19;
  answer.at(2) = 20;
  answer.at(3) = 21;
  answer.at(4) = 22;
  answer.at(5) = 23;
  answer.at(6) = 24;
  } else if (iEdge == 12) {// edge between nodes 8 5
  answer.at(1) = 22;
  answer.at(2) = 23;
  answer.at(3) = 24;
  answer.at(4) = 13;
  answer.at(5) = 14;
  answer.at(6) = 15;
  } else {
  _error ("giveEdgeDofMapping: wrong edge number");
  }

  return ;
}

double
LSpace ::   computeEdgeVolumeAround (GaussPoint* aGaussPoint, int iEdge) 
{
 double result = this->interpolation.edgeGiveTransformationJacobian (iEdge, domain, dofManArray, 
                                   *aGaussPoint->giveCoordinates(), 0.0); 
  return result * aGaussPoint -> giveWeight() ;
}


void
LSpace ::   computeEdgeIpGlobalCoords (FloatArray& answer, GaussPoint* gp, int iEdge)
{
 this->interpolation.edgeLocal2global (answer, iEdge, domain, dofManArray, *gp->giveCoordinates(), 0.0); 
}


int
LSpace :: computeLoadLEToLRotationMatrix (FloatMatrix& answer, int iEdge, GaussPoint* gp) 
{
  // returns transformation matrix from
  // edge local coordinate system
  // to element local c.s 
  // (same as global c.s in this case) 
  //
  // i.e. f(element local) = T * f(edge local)
  //
  _error ("computeLoadLEToLRotationMatrix: egde local coordinate system not supported");
  return 1;
}

void 
LSpace :: computeSurfaceNMatrixAt (FloatMatrix& answer, GaussPoint* sgp) 
{
 FloatArray n(4);
 interpolation.surfaceEvalN (n, *sgp->giveCoordinates(), 0.0);

 answer.resize (3,12);
 answer.zero();

 answer.at(1,1) = n.at(1);
 answer.at(1,4) = n.at(2);
 answer.at(1,7) = n.at(3);
 answer.at(1,10) = n.at(4);

 answer.at(2,2) = n.at(1);
 answer.at(2,5) = n.at(2);
 answer.at(2,8) = n.at(3);
 answer.at(2,11) = n.at(4);

 answer.at(3,3) = n.at(1);
 answer.at(3,6) = n.at(2);
 answer.at(3,9) = n.at(3);
 answer.at(3,12) = n.at(4);
}

void
LSpace :: giveSurfaceDofMapping (IntArray& answer, int iSurf) const
{
 answer.resize(12);
 if (iSurf == 1) {
  answer.at(1) = 1; // node 1
  answer.at(2) = 2;
  answer.at(3) = 3;

  answer.at(4) = 10; // node 4
  answer.at(5) = 11;
  answer.at(6) = 12;

  answer.at(7) = 7; // node 3
  answer.at(8) = 8;
  answer.at(9) = 9;

  answer.at(10) = 4; // node 2
  answer.at(11) = 5;
  answer.at(12) = 6;
 } else if (iSurf == 2) {
  answer.at(1) = 13; // node 5
  answer.at(2) = 14;
  answer.at(3) = 15;

  answer.at(4) = 16; // node 6
  answer.at(5) = 17;
  answer.at(6) = 18;

  answer.at(7) = 19; // node 7
  answer.at(8) = 20;
  answer.at(9) = 21;

  answer.at(10) = 22; // node 8
  answer.at(11) = 23;
  answer.at(12) = 24;
 } else if (iSurf == 3) {
  answer.at(1) = 1; // node 1
  answer.at(2) = 2;
  answer.at(3) = 3;

  answer.at(4) = 4; // node 2
  answer.at(5) = 5;
  answer.at(6) = 6;

  answer.at(7) = 16; // node 6
  answer.at(8) = 17;
  answer.at(9) = 18;

  answer.at(10) = 13; // node 5
  answer.at(11) = 14;
  answer.at(12) = 15;
 } else if (iSurf == 4) {
  answer.at(1) = 4; // node 2
  answer.at(2) = 5;
  answer.at(3) = 6;

  answer.at(4) = 7; // node 3
  answer.at(5) = 8;
  answer.at(6) = 9;

  answer.at(7) = 19; // node 7
  answer.at(8) = 20;
  answer.at(9) = 21;

  answer.at(10) = 16; // node 7
  answer.at(11) = 17;
  answer.at(12) = 18;
 } else if (iSurf == 5) {
  answer.at(1) = 7; // node 3
  answer.at(2) = 8;
  answer.at(3) = 9;

  answer.at(4) = 10; // node 4
  answer.at(5) = 11;
  answer.at(6) = 12;

  answer.at(7) = 22; // node 8
  answer.at(8) = 23;
  answer.at(9) = 24;

  answer.at(10) = 19; // node 7
  answer.at(11) = 20;
  answer.at(12) = 21;
 } else if (iSurf == 6) {
  answer.at(1) = 10; // node 4
  answer.at(2) = 11;
  answer.at(3) = 12;

  answer.at(4) = 1; // node 1
  answer.at(5) = 2;
  answer.at(6) = 3;

  answer.at(7) = 13; // node 5
  answer.at(8) = 14;
  answer.at(9) = 15;

  answer.at(10) = 22; // node 8
  answer.at(11) = 23;
  answer.at(12) = 24;
 } else {  
  _error ("giveSurfaceDofMapping: wrong surface number");
  }
}

IntegrationRule* 
LSpace :: GetSurfaceIntegrationRule (int approxOrder)
{
 IntegrationRule* iRule = new GaussIntegrationRule (1,domain, 1, 1);
 int npoints = iRule -> getRequiredNumberOfIntegrationPoints (_Square, approxOrder);
 iRule ->setUpIntegrationPoints (_Square, npoints, this, _Unknown);
 return iRule;
}

double
LSpace :: computeSurfaceVolumeAround (GaussPoint* gp, int iSurf) 
{
 double       determinant,weight,volume;
 determinant = fabs(interpolation.surfaceGiveTransformationJacobian (iSurf, domain, dofManArray, *gp->giveCoordinates(), 0.0));

 weight      = gp -> giveWeight();
 volume      = determinant * weight;
 
 return volume;
}

void
LSpace :: computeSurfIpGlobalCoords (FloatArray& answer, GaussPoint* gp, int isurf) 
{
 interpolation.surfaceLocal2global (answer, isurf, domain, dofManArray, *gp->giveCoordinates(), 0.0);
}


int
LSpace :: computeLoadLSToLRotationMatrix (FloatMatrix& answer, int isurf, GaussPoint* gp)
{
  // returns transformation matrix from
  // surface local coordinate system
  // to element local c.s 
  // (same as global c.s in this case) 
  //
  // i.e. f(element local) = T * f(edge local)
  
  // definition of local c.s on surface:
  // local z axis - perpendicular to surface, pointing outwards from element 
  // local x axis - is in global xy plane (perpendicular to global z axis)
  // local y axis - completes the righ hand side cs.
  
  /*
     _error ("computeLoadLSToLRotationMatrix: surface local coordinate system not supported");
     return 1;
  */
  int i, j;
  FloatArray gc(3);
  FloatArray h1(3),h2(3),nn(3),n(3);
  IntArray snodes(4);

  answer.resize(3,3);

  this->interpolation.computeSurfaceMapping (snodes, dofManArray, isurf);
  for (i=1; i<=4; i++) gc.add(domain->giveNode(snodes.at(i))->giveCoordinates());
  gc.times(1./4.);
  // determine "average normal"
  for (i=1; i<=4; i++) {
    j = (i)%4+1;
    h1 = *domain->giveNode(snodes.at(i))->giveCoordinates();
    h1.substract(gc);
    h2 = *domain->giveNode(snodes.at(j))->giveCoordinates();
    h2.substract(gc);
    n.beVectorProductOf(h1,h2);
    if (dotProduct(n,n,3) > 1.e-6) n.normalize();
    nn.add(n);
  }
  nn.times(1./4.); 
  if (dotProduct(nn,nn,3) < 1.e-6) {
    answer.zero();
    return 1;
  }
  nn.normalize();
  for (i=1; i<=3; i++) answer.at(i,3) = nn.at(i);
  // determine lcs of surface
  // local x axis in xy plane
  double test = fabs(fabs(nn.at(3)) -1.0);
  if (test < 1.e-5) {
    h1.at(1) = answer.at(1,1) = 1.0;
    h1.at(2) = answer.at(2,1) = 0.0;
  } else {
    h1.at(1) = answer.at(1,1) = answer.at(2,3);
    h1.at(2) = answer.at(2,1) = -answer.at(1,3);
  }
  h1.at(3) = answer.at(3,1) = 0.0;
  // local y axis perpendicular to local x,z axes
  h2.beVectorProductOf (nn,h1);
  for (i=1; i<=3; i++) answer.at(i,2) = h2.at(i);

  return 1;
}
