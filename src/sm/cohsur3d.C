/*
class CohesiveSurface3d added by Milan Jirasek on 1 Feb 2010

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2010   Borek Patzak                                       

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

#include "cohsur3d.h"
#include "node.h"
#include "particle.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "flotarry.h"

#include <stdio.h>

/*
#include "dof.h"
#include "engngm.h"
#include "structuralcrosssection.h"
#include "structuralmaterial.h"
#include "mathfem.h"
*/

#ifndef __MAKEDEPEND
#include <math.h>
#endif

#ifdef __OOFEG 
#include "oofeggraphiccontext.h"
#include "conTable.h"
#include "oofegutils.h"
#endif

namespace oofem {

CohesiveSurface3d :: CohesiveSurface3d (int n, Domain* aDomain) : StructuralElement (n,aDomain) 
// Constructor.
{
  numberOfDofMans = 2;
  area   = -1.;
  length = -1.;
}

void
CohesiveSurface3d :: computeBmatrixAt (GaussPoint* aGaussPoint, FloatMatrix& answer, int li, int ui)
   // Returns the strain-displacement matrix of the receiver.
{
  double x01, y01, z01, x02, y02, z02;
  FloatMatrix Bloc(3,12);
  Node *nodeA, *nodeB;

  // Coordinate differences

  nodeA   = this->giveNode(1) ;
  nodeB   = this->giveNode(2) ;
  x01 = nodeA->giveCoordinate(1) - center.at(1) ;
  y01 = nodeA->giveCoordinate(2) - center.at(2) ;
  z01 = nodeA->giveCoordinate(3) - center.at(3) ;
  x02 = nodeB->giveCoordinate(1) - center.at(1) ;
  y02 = nodeB->giveCoordinate(2) - center.at(2) ;
  z02 = nodeB->giveCoordinate(3) - center.at(3) ;

 // B matrix in local coordinates (and without the term 1/length)

 Bloc.zero();

 Bloc.at(1,1) =  -1.;
 Bloc.at(2,2) =  -1.;
 Bloc.at(3,3) =  -1.;

 Bloc.at(1,5) =   z01;
 Bloc.at(1,6) =  -y01;
 Bloc.at(2,4) =  -z01;
 Bloc.at(2,6) =   x01;
 Bloc.at(3,4) =   y01;
 Bloc.at(3,5) =  -x01;

 Bloc.at(1,7) =   1.;
 Bloc.at(2,8) =   1.;
 Bloc.at(3,9) =   1.;

 Bloc.at(1,11) =  -z02;
 Bloc.at(1,12) =   y02;
 Bloc.at(2,10) =   z02;
 Bloc.at(2,12) =  -x02;
 Bloc.at(3,10) =  -y02;
 Bloc.at(3,11) =   x02;

 // Transformation to global coordinates

 answer.resize(3,12);
 answer.beProductOf(lcs,Bloc);

 // Division by the length

 answer.times(1./length);

 return ;
}

void  CohesiveSurface3d :: computeGaussPoints ()
  // Sets up the array of Gauss points of the receiver.
{
  // The Gauss point is used only when methods from crosssection and/or material
  // classes are requested.
  numberOfIntegrationRules = 1 ;
  integrationRulesArray = new IntegrationRule* [1];
  integrationRulesArray[0] = new GaussIntegrationRule (1,this);
  integrationRulesArray[0]->setUpIntegrationPoints (_Line, 1, _3dInterface);
}

 
double  
CohesiveSurface3d :: computeVolumeAround (GaussPoint* aGaussPoint)
{
  return area * length;
}


void
CohesiveSurface3d :: giveDofManDofIDMask  (int inode, EquationID, IntArray& answer) const {
// returns DofId mask array for inode element node.
// DofId mask array determines the dof ordering requsted from node.
// DofId mask array contains the DofID constants (defined in cltypes.h)
// describing physical meaning of particular DOFs.

  answer.resize (6);

  answer.at(1) = D_u;
  answer.at(2) = D_v;
  answer.at(3) = D_w;
  answer.at(4) = R_u;
  answer.at(5) = R_v;
  answer.at(6) = R_w;
}

double  
CohesiveSurface3d :: giveLength ()
   // Returns the length of the receiver.
{
   double dx,dy,dz ;
   Node   *nodeA,*nodeB ;

   if (length <= 0.) {
      nodeA   = this->giveNode(1) ;
      nodeB   = this->giveNode(2) ;
      dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1) ;
      dy      = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2) ;
      dz      = nodeB->giveCoordinate(3) - nodeA->giveCoordinate(3) ;
      length  = sqrt(dx*dx + dy*dy + dz*dz) ;}

   return length ;
}

void
CohesiveSurface3d :: evaluateCenter()
{
  Particle *nodeA,*nodeB;
  double RA, RB, L, aux;
  int i;

  nodeA  = (Particle*) this -> giveNode(1);
  nodeB  = (Particle*) this -> giveNode(2);
  RA = nodeA -> giveRadius();
  RB = nodeB -> giveRadius();
  L  = giveLength();
  aux = 0.5 + (RA-RB)/(2.*L);

  center.resize(3);
  for (i=1; i<= 3; i++) 
    center.at(i) = aux*(nodeB->giveCoordinate(i)) + (1.-aux)*(nodeA->giveCoordinate(i));
}

void
CohesiveSurface3d :: evaluateLocalCoordinateSystem()
//
// Computes unit vectors of local coordinate system, stored by rows.
//
{
  FloatArray lx(3), ly(3), lz(3);
  Node   *nodeA,*nodeB;
  int i;

  nodeA  = this -> giveNode(1) ;
  nodeB  = this -> giveNode(2) ;
 
  for (i=1; i<= 3; i++) 
    lx.at(i) = nodeB->giveCoordinate(i)-nodeA->giveCoordinate(i);
  lx.normalize(); 

  ly.zero();
  if (abs(lx.at(1))>abs(lx.at(2)))
    ly.at(2) = 1.;
  else
    ly.at(1) = 1.;

  lz.beVectorProductOf (lx, ly);
  lz.normalize();
  ly.beVectorProductOf (lz, lx);
  ly.normalize();

  lcs.resize (3,3);
  for (i=1; i<= 3; i++) {
    lcs.at(1,i) = lx.at(i);
    lcs.at(2,i) = ly.at(i);
    lcs.at(3,i) = lz.at(i);
  }
}

IRResultType
CohesiveSurface3d :: initializeFrom (InputRecord* ir)
{
  const char *__proc = "initializeFrom"; 
  IRResultType result;                   

  // first call parent 
  StructuralElement :: initializeFrom (ir);

  // read the area from the input file
  IR_GIVE_FIELD (ir, area, IFT_Beam3d_refnode, "area"); 
  if (area < 0.)
    _error ("CohesiveSurface3d::instanciateFrom: negative area specified");

  // initialize one Gauss point
  this -> computeGaussPoints();

  // evaluate the length
  giveLength();

  // evaluate the coordinates of the center
  evaluateCenter();

  // evaluate the local coordinate system
  evaluateLocalCoordinateSystem();

  return IRRT_OK;
}


int
CohesiveSurface3d :: computeGlobalCoordinates (FloatArray& answer, const FloatArray& lcoords) 
{
  answer = center;
  return 1;
}

void 
CohesiveSurface3d :: printOutputAt (FILE* File, TimeStep* stepN)
{
 // Performs end-of-step operations.
 
 int         i;
 FloatArray rg, rl, Fg, Fl;
 FloatMatrix T;

 fprintf (File,"element %d :\n",number) ;
 
 for (i=0 ; i < numberOfIntegrationRules ; i++) 
   integrationRulesArray[i]->printOutputAt(File,stepN);
 
 fprintf (File,"\n") ;
}

#ifdef __OOFEG
void CohesiveSurface3d :: drawRawGeometry (oofegGraphicContext& gc)
{
  if (!gc.testElementGraphicActivity(this)) return; 

  WCRec p[4];
  GraphicObj *go;

  Particle* nodeA = (Particle*) this->giveNode(1);
  Particle* nodeB = (Particle*) this->giveNode(2);
  double rA = nodeA -> giveRadius(); 
  double rB = nodeB -> giveRadius(); 
  double r = (rA+rB)/4.;

  EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
  EASValsSetColor(gc.getElementColor());
  EASValsSetEdgeColor(gc.getElementEdgeColor());
  EASValsSetEdgeFlag(TRUE);
  EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);

//  plot a line segment connecting the particles
  WCRec pl[2];   
  pl[0].x = (FPNum) nodeA->giveCoordinate(1);
  pl[0].y = (FPNum) nodeA->giveCoordinate(2);
  pl[0].z = (FPNum) nodeA->giveCoordinate(3);
  pl[1].x = (FPNum) nodeB->giveCoordinate(1);
  pl[1].y = (FPNum) nodeB->giveCoordinate(2);
  pl[1].z = (FPNum) nodeB->giveCoordinate(3);
  go = CreateLine3D(pl);
  EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
  EGAttachObject(go, (EObjectP) this);
  EMAddGraphicsToModel(ESIModel(), go);
}

void CohesiveSurface3d :: drawDeformedGeometry (oofegGraphicContext& gc, UnknownType type)
{
  GraphicObj *go, *go1, *go2;

 if (!gc.testElementGraphicActivity(this)) return; 

  TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
  double defScale = gc.getDefScale();
  WCRec p[2];   
  EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
  EASValsSetColor(gc.getDeformedElementColor());
  EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);


//  plot a line segment connecting the displaced particles
  p[0].x = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
  p[0].y = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
  p[0].z = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(3,tStep,EID_MomentumBalance,defScale);

  p[1].x = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
  p[1].y = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
  p[1].z = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(3,tStep,EID_MomentumBalance,defScale);

// plot the displaced particles
  
//  EASValsSetLayer(OOFEG_NODE_ANNOTATION_LAYER);
  EASValsSetMType(FILLED_CIRCLE_MARKER);
  EASValsSetColor(gc.getNodeColor());
  EASValsSetMSize(4);

  go1 = CreateMarker3D(p);
  EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK | MTYPE_MASK | MSIZE_MASK, go);
  EMAddGraphicsToModel(ESIModel(), go1);
} 


void 
CohesiveSurface3d :: drawScalar(oofegGraphicContext& context)
{
  if (!context.testElementGraphicActivity(this)) 
    return;
  
  FloatArray val; 
  TimeStep* tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
  GaussPoint* gp = integrationRulesArray[0]-> getIntegrationPoint(0);
  if (!giveIPValue (val, gp, context.giveIntVarType(), tStep))
    return;
 
  IntArray map; int indx;
  this->giveIntVarCompFullIndx (map, context.giveIntVarType());
  if ((indx = map.at(context.giveIntVarIndx())) == 0) 
    return;
  double s = val.at(indx);
   
  int i; WCRec p[3];
  for (i=0; i<2; i++) {
    if (context.getInternalVarsDefGeoFlag()) {
      // use deformed geometry
      double defScale = context.getDefScale();
      p[i].x = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
      p[i].y = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
      p[i].z = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(3,tStep,EID_MomentumBalance,defScale);
    } 
    else {
      p[i].x = (FPNum) this->giveNode(i+1)->giveCoordinate(1);
      p[i].y = (FPNum) this->giveNode(i+1)->giveCoordinate(2);
      p[i].z = (FPNum) this->giveNode(i+1)->giveCoordinate(3);
    }
  }

  EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
  EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
  if (s>0.)
    EASValsSetColor(context.getActiveCrackColor());
  else 
    EASValsSetColor(context.getCrackPatternColor());
  GraphicObj *go = CreateLine3D (p);
  EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
  EMAddGraphicsToModel(ESIModel(), go);
}
#endif

} // namespace oofem
