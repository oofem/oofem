/* $Header: /home/cvs/bp/oofem/oofemlib/src/element.h,v 1.27 2003/04/06 14:08:24 bp Exp $ */
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


/*
  oofem nodes - controll points (coordinates ) + dofs
  oofem elements - NURBS patches as well as integration elements
  
  
  NURBS PATCH:
  knot vector - store knot coordinates + multiplicity
  patch integration rule - keep list of elements
  
  integration element
  
  FEInterpolation:
  - need to be enriched, as one should (or can) pass knot span to evaluation routines 
  (in this way the patch by patch evaluation can be faster)
  
*/

#include "inputrecord.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "matresponsemode.h"
#include "crosssection.h"
#include "structuralcrosssection.h"
#include "mathfem.h"
#include "iga.h"

void BSplineInterpolation::evalN(FloatArray &answer, const FloatArray &lcoords, const IntArray *span,  double time) {
  int i, l,k,c;
  FloatArray N[nsd];
  
  for (i=0; i< nsd; i++) {
    //span(i) = this->findSpan (numberOfControllPoints[i], degree[i], lcoords(i), knotVector[i]);
    this->basisFuns (N[i], span->at(i+1), lcoords(i), degree[i], knotVector[i]);
  }
  
  if (nsd == 2) {
    c=0;
    for (l=0;l<=degree[0]; l++) {
      for (k=0; k<=degree[1]; k++) {
        answer.at(c++) = N[0](l)*N[1](k);
      }
    }
  }
  else {
    OOFEM_ERROR2 ("evalN not implemented for nsd = %d", nsd);
  }
}


void BSplineInterpolation::evalN(FloatArray &answer, const FloatArray &lcoords, double time) {
  int i;
  IntArray span(nsd);
  
  for (i=0; i< nsd; i++) {
    span(i) = this->findSpan (numberOfControllPoints[i], degree[i], lcoords(i), knotVector[i]);
  }
  this->evalN (answer, lcoords, &span, time);
}

void BSplineInterpolation::evaldNdx(FloatMatrix &answer, const FloatArray **coords, const FloatArray &lcoords, const IntArray* knotSpan, double time) {
  OOFEM_ERROR ("Not yet inplemented, contact lazy dr for implementation");
}

void BSplineInterpolation::evaldNdx(FloatMatrix &answer, const FloatArray **coords, const FloatArray &lcoords, double time) {
  OOFEM_ERROR ("Not yet inplemented, contact lazy dr for implementation");
}

void BSplineInterpolation::evaldNdx(FloatMatrix &answer, Domain *d, IntArray &nodes, const FloatArray &lcoords, const IntArray* span, double time) {
  FloatMatrix jacobian(nsd,nsd);
  FloatMatrix ders[nsd];
  FloatArray temp(nsd);
  double Jacob;
  int count = 1, cnt, i, l, k, uind, vind;
  for (i=0; i< nsd; i++) {
    this->DersBasisFuns(1, lcoords(i), span->at(i+1), degree[i], knotVector[i], ders[i]);
    count *=giveNumberOfKnotBasisFunctions ();
  }
  jacobian.zero();
  answer.resize(count, nsd);
  answer.zero();
  if (nsd == 2) {
    uind = span->at(1)-degree[0];
    for (l=0;l<=degree[1]; l++) {
      temp.zero();
      vind = span->at(2)-degree[1]+1;
      for (k=0; k<=degree[0]; k++) {
        temp(0) += ders[0](1,k)*d->giveNode(nodes.at((uind+k)*numberOfControllPoints[1]+vind+1))->giveCoordinate(1);  // HAHA see local2global
        temp(1) += ders[0](1,k)*d->giveNode(nodes.at((uind+k)*numberOfControllPoints[1]+vind+1))->giveCoordinate(2);  // HAHA coords
      }
      jacobian(0,0) += ders[1](0,l)*temp(0);
      jacobian(0,1) += ders[1](0,l)*temp(1);
    }
    for (l=0;l<=degree[1]; l++) {
      temp.zero();
      vind = span->at(2)-degree[1]+1;
      for (k=0; k<=degree[0]; k++) {
        temp(0) += ders[0](0,k)*d->giveNode(nodes.at((uind+k)*numberOfControllPoints[1]+vind+1))->giveCoordinate(1);  // HAHA see local2global
        temp(1) += ders[0](0,k)*d->giveNode(nodes.at((uind+k)*numberOfControllPoints[1]+vind+1))->giveCoordinate(2);  // HAHA coords
      }
      jacobian(1,0) += ders[1](1,l)*temp(0);
      jacobian(1,1) += ders[1](1,l)*temp(1);
    }
    Jacob = jacobian.giveDeterminant();
    
    if(fabs(Jacob) < 1.0e-10){
      OOFEM_ERROR ("evaldNdx - zero Jacobian");
    }
		
    cnt=0;
    // HAHA ordering ??? 
    for (k=0; k<=degree[0]; k++) {
      for (l=0;l<=degree[1]; l++) {
        temp(0) = ders[0](1,k)*ders[1](0,l);
        temp(1) = ders[0](0,k)*ders[1](1,0);
        answer(cnt,0) += (jacobian(1,1)*temp(0)-jacobian(0,1)*temp(1)) / Jacob;
        answer(cnt,1) += (-jacobian(1,0)*temp(0)+jacobian(0,0)*temp(1)) / Jacob;
      }
    }
  }
  else {
    OOFEM_ERROR2 ("evaldNdx not implemented for nsd = %d", nsd);
  }
}

void BSplineInterpolation::evaldNdx(FloatMatrix &answer, Domain *d, IntArray &nodes, const FloatArray &lcoords, double time) {
  IntArray span(nsd);
  int i;
  for (i=0; i< nsd; i++) {
    span(i) = this->findSpan (numberOfControllPoints[i], degree[i], lcoords(i), knotVector[i]);
  } 
  evaldNdx (answer, d, nodes, lcoords, &span, time);
}

void BSplineInterpolation::local2global(FloatArray &answer, Domain *d, IntArray &nodes, const FloatArray &lcoords, double time) {
  /* Based on SurfacePoint A3.5 implementation*/
  int i,l,k;
  IntArray span(nsd);
  FloatArray temp(nsd);
  FloatArray N[nsd];
  answer.resize(nsd);
  
  for (i=0; i< nsd; i++) {
    span(i) = this->findSpan (numberOfControllPoints[i], degree[i], lcoords(i), knotVector[i]);
    this->basisFuns (N[i], span(i), lcoords(i), degree[i], knotVector[i]); 
  }
  if (nsd == 2) {
    int uind, vind;
    
    uind = span(0)-degree[0]; 
    answer.zero();
    for (l=0;l<=degree[1]; l++) {
      temp.zero();
      vind = span(1)-degree[1]+1;  
      for (k=0; k<=degree[0]; k++) {
        //temp(0) = temp(0) + N[0][k]*d->giveNode(nodes.at((uind+k)*numberOfControllPoints[0]+vind))->giveCoordinate(1);
        //temp(1) = temp(1) + N[0][k]*d->giveNode(nodes.at((uind+k)*numberOfControllPoints[0]+vind))->giveCoordinate(2);
        temp(0) += N[0](k)*d->giveNode(nodes.at((uind+k)*numberOfControllPoints[0]+vind+1))->giveCoordinate(1);
        temp(1) += N[0](k)*d->giveNode(nodes.at((uind+k)*numberOfControllPoints[0]+vind+1))->giveCoordinate(2);
      }
      answer(0) += N[1](l)*temp(0);
      answer(1) += N[1](l)*temp(1);
    }
  } else {
    OOFEM_ERROR2 ("local2global not implemented for nsd = %d", nsd);
  }
}

double BSplineInterpolation::giveTransformationJacobian(Domain *d, IntArray &nodes, const FloatArray &lcoords, const IntArray* span, double time) {
  FloatMatrix jacobian(nsd,nsd);
  FloatMatrix ders[nsd];
  FloatArray temp(nsd);
  double Jacob;
  int count = 1, i, l, k, uind, vind;
  for (i=0; i< nsd; i++) {
    this->DersBasisFuns(1, lcoords(i), span->at(i+1), degree[i], knotVector[i], ders[i]);
    count *=giveNumberOfKnotBasisFunctions ();
  }
  jacobian.zero();
  if (nsd == 2) {
    uind = span->at(1)-degree[0];
    for (l=0;l<=degree[1]; l++) {
      temp.zero();
      vind = span->at(2)-degree[1]+1;
      for (k=0; k<=degree[0]; k++) {
        temp(0) += ders[0](1,k)*d->giveNode(nodes.at((uind+k)*numberOfControllPoints[1]+vind+1))->giveCoordinate(1);  // HAHA see local2global
        temp(1) += ders[0](1,k)*d->giveNode(nodes.at((uind+k)*numberOfControllPoints[1]+vind+1))->giveCoordinate(2);  // HAHA coords
      }
      jacobian(0,0) += ders[1](0,l)*temp(0);
      jacobian(0,1) += ders[1](0,l)*temp(1);
    }
    for (l=0;l<=degree[1]; l++) {
      temp.zero();
      vind = span->at(2)-degree[1]+1;
      for (k=0; k<=degree[0]; k++) {
        temp(0) += ders[0](0,k)*d->giveNode(nodes.at((uind+k)*numberOfControllPoints[1]+vind+1))->giveCoordinate(1);  // HAHA see local2global
        temp(1) += ders[0](0,k)*d->giveNode(nodes.at((uind+k)*numberOfControllPoints[1]+vind+1))->giveCoordinate(2);  // HAHA coords
      }
      jacobian(1,0) += ders[1](1,l)*temp(0);
      jacobian(1,1) += ders[1](1,l)*temp(1);
    }
    Jacob = jacobian.giveDeterminant();
    
    if(fabs(Jacob) < 1.0e-10){
      OOFEM_ERROR ("evaldNdx - zero Jacobian");
    }
		
  }
  else {
    OOFEM_ERROR2 ("evaldNdx not implemented for nsd = %d", nsd);
  }
  return Jacob;
}

double BSplineInterpolation::giveTransformationJacobian(Domain *d, IntArray &nodes, const FloatArray &lcoords, double time) {
  IntArray span(nsd);
  int i;
  for (i=0; i< nsd; i++) {
    span(i) = this->findSpan (numberOfControllPoints[i], degree[i], lcoords(i), knotVector[i]);
  } 
  return giveTransformationJacobian (d, nodes, lcoords, &span, time);
}

int BSplineInterpolation::giveKnotBasisFuncMask (const IntArray& knotSpan, IntArray& mask) {
  int c, _i, _j, _iindx, _jindx;
  
  if (nsd == 2) {
    mask.resize((degree[0]+1)*(degree[1]+1));
    c=1;
    for (_i=0; _i<=degree[0]; _i++) {
      _iindx = (_i+knotSpan(0)-degree[0]);
      for (_j=0; _j<=degree[1]; _j++) {
        _jindx = (_j+knotSpan(1)-degree[1]);
        mask.at(c++) = _iindx+(_jindx)*numberOfControllPoints[0];
      }
    }          
  } else {
    OOFEM_ERROR2 ("giveKnotBasisFunctMask not implemented for nsd = %d", nsd);
  }
  return 1;
}


int  BSplineInterpolation::giveNumberOfKnotBasisFunctions () { 
  int i, answer  = 1;
  for (i=0; i<nsd; i++) answer*=(degree[i]+1);
  return answer;
}


void BSplineInterpolation::basisFuns (FloatArray &N, int i, double u, int p, const double* U) 
{
  int j, r;
  FloatArray right(p+1);
  FloatArray left(p+1);
  double saved, temp;
  
  N.resize(p);
  N(0) = 1.0;
  for (j=1; j<=p; j++) {
    left(j)=u-U[i+1-j];
    right(j)=U[i+j]-u;
    saved = 0.0;
    for (r=0; r<j; r++) {
      temp=N(r)/(right(r+1)+left(j-r));
      N(r)=saved+right(r+1)*temp;
      saved=left(j-r)*temp;
    }
    N(j)=saved;
  }
}


void BSplineInterpolation::DersBasisFuns(int n, double u, int span, int deg, double* const U, FloatMatrix& ders) {
  
  FloatArray left(deg+1);
  FloatArray right(deg+1);
  
  FloatMatrix ndu(deg+1,deg+1) ;
  double saved,temp ;
  int j,r ;
  
  ders.resize(n+1,deg+1) ;
  
  ndu(0,0) = 1.0 ;
  for(j=1; j<= deg ;j++){
    left(j) = u-U[span+1-j] ;
    right(j) = U[span+j]-u ;
    saved = 0.0 ;
    
    for(r=0;r<j ; r++){
      // Lower triangle
      ndu(j,r) = right(r+1)+left(j-r) ;
      temp = ndu(r,j-1)/ndu(j,r) ;
      // Upper triangle
      ndu(r,j) = saved+right(r+1) * temp ;
      saved = left(j-r) * temp ;
    }
    
    ndu(j,j) = saved ;
  }
  
  for(j=0;j<=deg;j++)
    ders(0,j) = ndu(j,deg) ;
  
  // Compute the derivatives
  FloatMatrix a(2,deg+1) ;
  for(r=0;r<=deg;r++){
    int s1,s2 ;
    s1 = 0 ; s2 = 1 ; // alternate rows in array a
    a(0,0) = 1.0 ;
    // Compute the kth derivative
    for(int k=1;k<=n;k++){
      double d ;
      int rk,pk,j1,j2 ;
      d = 0.0 ;
      rk = r-k ; pk = deg-k ;
      
      if(r>=k){
        a(s2,0) = a(s1,0)/ndu(pk+1,rk) ;
        d = a(s2,0)*ndu(rk,pk) ;
      }
      
      if(rk>=-1){
        j1 = 1 ;
      } else{
        j1 = -rk ;
      }
      
      if(r-1 <= pk){
        j2 = k-1 ;
      }else{
        j2 = deg-r ;
      }
      
      for(j=j1;j<=j2;j++){
        a(s2,j) = (a(s1,j)-a(s1,j-1))/ndu(pk+1,rk+j) ;
        d += a(s2,j)*ndu(rk+j,pk) ;
      }
      
      if(r<=pk){
        a(s2,k) = -a(s1,k-1)/ndu(pk+1,r) ;
        d += a(s2,k)*ndu(r,pk) ;
      }
      ders(k,r) = d ;
      j = s1 ; s1 = s2 ; s2 = j ; // Switch rows
    }
  }
  
  // Multiply through by the correct factors
  r = deg ;
  for(int k=1;k<=n;k++){
    for(j=0;j<=deg;j++)
      ders(k,j) *= r ;
    r *= deg-k ;
  }
}

int BSplineInterpolation::findSpan(int n, int p, double u, const double* U) const {
  
  if(u == U[n+1]) return n; 
  
  int low  = p ;
  int high = n+1 ; 
  int mid = (low+high)/2 ;
  
  while(u<U[mid] || u>= U[mid+1]){
    if(u<U[mid])
      high = mid ;
    else
        low = mid ;
    mid = (low+high)/2 ;
  }
  return mid ;
  
}

IRResultType IGAElement::initializeFrom(InputRecord *ir) {

  Element::initializeFrom (ir); // read nodes , material, cross section
  this->giveInterpolation()->initializeFrom (ir); // read geometry
  
  int indx, ui,vi, i, numberOfGaussPoints; // HUHU
  int **knotMultiplicity; // HUHU
  double du,dv;
  double ** const knotVector = this->giveInterpolation()->giveKnotVector();
  const FloatArray *gpcoords;
  FloatArray newgpcoords;
  
  // generate individual IntegrationElements; one for each nonzero knot span
  if (this->giveNsd() == 2) {
    // HUHU mame pristup na numberOfKnotSpans[]?
    this->numberOfIntegrationRules = this->giveInterpolation()->giveNumberOfKnotSpans(1)*this->giveInterpolation()->giveNumberOfKnotSpans(2);
    integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];
    indx = 0; 
    IntArray knotSpan(2); 
    knotSpan.at(1) = -1; 
    for (ui=0; ui<this->giveInterpolation()->giveNumberOfKnotSpans(1); ui++) {
      knotSpan.at(1) += knotMultiplicity[0][ui];
      knotSpan.at(2) = -1;
      for (vi=0; vi<this->giveInterpolation()->giveNumberOfKnotSpans(2); vi++) {
        indx++;
        knotSpan.at(2)+=knotMultiplicity[1][vi];
        integrationRulesArray [ indx ] = new IGA_IntegrationElement (indx, this, knotSpan);
        integrationRulesArray [ indx ] ->setUpIntegrationPoints(_Square, numberOfGaussPoints, _PlaneStress); // HUHU _PlaneStress
        // remap local subelement gp coordinates into knot span coordinates and update integration weight 
        for (i=0; i<numberOfGaussPoints; i++) {
          gpcoords = integrationRulesArray [ indx ]->getIntegrationPoint(i)->giveCoordinates();
          du = knotVector[0][knotSpan.at(1)+1]-knotVector[0][knotSpan.at(1)];
          dv = knotVector[1][knotSpan.at(2)+1]-knotVector[1][knotSpan.at(2)];
          newgpcoords.at(1) = knotVector[0][knotSpan.at(1)]+du*gpcoords->at(1)/2.0;
          newgpcoords.at(2) = knotVector[1][knotSpan.at(1)]+dv*gpcoords->at(2)/2.0;
          integrationRulesArray [ indx ]->getIntegrationPoint(i)->setCoordinates(newgpcoords);
          integrationRulesArray [ indx ]->getIntegrationPoint(i)->setWeight(integrationRulesArray [ indx ]->getIntegrationPoint(i)->giveWeight()/4.0*du*dv);
        }
      }
    }
  } else {
    OOFEM_ERROR2 ("unsupported number of spatial dimensions (nsd = %d)", this->giveNsd());
  }
  return IRRT_OK; 
}

StructuralElementEvaluator::StructuralElementEvaluator() {
  this->rotationMatrix=NULL;
}
         
int StructuralElementEvaluator::giveIntegrationElementCodeNumbers (IntArray& answer, Element* elem, 
                                                                   IntegrationRule* ie, EquationID ut) {
  int i;
  IntArray mask, nodeDofIDMask, nodalArray;
  // first evaluate nonzero basis function mask
  if (elem->giveInterpolation()->hasSubPatchFormulation()) {
    IGA_IntegrationElement *ee = (IGA_IntegrationElement *)ie;
    elem->giveInterpolation()->giveKnotBasisFuncMask (*ee->giveKnotSpan(), mask) ;
    // loop over nonzero shape functions and assemble localization array
    answer.resize(0);
    for (i=1; i<=mask.giveSize(); i++) {
      elem->giveDofManDofIDMask (mask.at(i), ut, nodeDofIDMask);
      elem->giveDofManager(i)->giveLocationArray(nodeDofIDMask, nodalArray);
      answer.followedBy(nodalArray);
    }
    return 1;
  } else {
    return 0;
  }
}


void PlaneStressStructuralElementEvaluator::computeNMatrixAt (FloatMatrix& answer, GaussPoint* gp) {
  
  int i, nDofMan;
  FloatArray N;
  FEInterpolation* interp = gp->giveElement()->giveInterpolation();
  
  interp->evalN (N, *gp->giveCoordinates(), gp->giveIntegrationRule()->giveKnotSpan(), 0.0);
  
  if ((nDofMan = interp->giveNumberOfKnotBasisFunctions()) == 0)
    nDofMan = gp->giveElement()->giveNumberOfDofManagers();
  
  answer.resize(2, nDofMan*2);
  answer.zero();
  
  for (i=1; i <= nDofMan; i++) {
    answer.at(1, i*2-1) = N.at(i);
    answer.at(2, i*2)   = N.at(i);
  }
}

void PlaneStressStructuralElementEvaluator::computeBMatrixAt (FloatMatrix& answer, GaussPoint* gp) {
  
  int i, nDofMan;
  IntArray dofmanSubElementMask;
  FloatMatrix d;
  
  FEInterpolation* interp = gp->giveElement()->giveInterpolation();
  // this uses FEIInterpolation::nodes2coords - quite inefficient in this case (large num of dofmans)
  interp->evaldNdx (d, gp->giveElement()->giveDomain(), gp->giveElement()->giveDofManArray(), 
                    *gp->giveCoordinates(), gp->giveIntegrationRule()->giveKnotSpan(), 0.0);
  
  if ((nDofMan = interp->giveNumberOfKnotBasisFunctions ()) == 0)
    nDofMan = gp->giveElement()->giveNumberOfDofManagers();
  
  answer.resize(3, nDofMan*2);
  answer.zero();
  
  for (i=1; i <= nDofMan; i++) {
    answer.at(1, i*2-1) = d.at(i, 1);
    answer.at(2, i*2)   = d.at(i, 2);
    
    answer.at(3, 2*i-1) = d.at(i, 2);
    answer.at(3, 2*i-0) = d.at(i, 1);
  }
}
  

void
StructuralElementEvaluator ::  giveCharacteristicMatrix(FloatMatrix &answer,
                                                        CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver accordind to mtrx
//
{
  if ( mtrx == StiffnessMatrix ) {
    this->computeStiffnessMatrix(answer, TangentStiffness, tStep);
  } else {
    OOFEM_ERROR2( "giveCharacteristicMatrix: Unknown Type of characteristic mtrx (%s)", __CharTypeToString(mtrx) );
  }
  
  return;
}
 
void StructuralElementEvaluator::computeStiffnessMatrix (FloatMatrix& answer, MatResponseMode rMode, TimeStep *tStep) {
  
  int ir, j, numberOfIntegrationRules;
  FloatMatrix temp, bj, d, dbj;
  IntegrationRule* iRule;
  GaussPoint* gp;
  Element *elem = this->giveElement();
  int ndofs = elem->computeNumberOfDofs(EID_MomentumBalance);
  bool matStiffSymmFlag = elem->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode, elem->giveMaterial()->giveNumber());
  IntArray irlocnum;
  double dV;
  
  answer.resize( ndofs, ndofs );
  answer.zero();	  
  
  FloatMatrix &m = answer;
  if (elem->giveInterpolation()->hasSubPatchFormulation()) m = temp;
  
  numberOfIntegrationRules = elem->giveNumberOfIntegrationRules();
  // loop over individual integration rules
  for (ir=0; ir < numberOfIntegrationRules; ir++) {
    m.resize(0,0);
    iRule = elem->giveIntegrationRule(ir);
    // loop over individual integration points
    for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
      gp = iRule->getIntegrationPoint(j);
      this->computeBMatrixAt(bj, gp);
      //elem->computeConstitutiveMatrixAt(d, rMode, gp, tStep);
      ( ( StructuralCrossSection * ) elem->giveCrossSection() )
        ->giveCharMaterialStiffnessMatrix(d, rMode, gp, tStep);
      
      dV = this->computeVolumeAround(gp);
      dbj.beProductOf(d, bj);
      if ( matStiffSymmFlag ) {
        m.plusProductSymmUpper(bj, dbj, dV);
      } else {
        m.plusProductUnsym(bj, dbj, dV);
      }
    }
    
    if ( matStiffSymmFlag ) {
      m.symmetrized();
    }
    // localize irule contribution into element matrix
    if (this->giveIntegrationElementCodeNumbers (irlocnum, elem, iRule, EID_MomentumBalance)) 
      answer.assemble (m, irlocnum);
    
  } // end loop over irules
  
  if ( this->updateRotationMatrix() ) {
    answer.rotatedWith(* this->rotationMatrix);
  }
  return;
}
  
double PlaneStressStructuralElementEvaluator::computeVolumeAround(GaussPoint *gp) { 
  double determinant, weight, thickness, volume;
  determinant = fabs( this->giveElement()->giveInterpolation()
                      ->giveTransformationJacobian(gp->giveElement()->giveDomain(), gp->giveElement()->giveDofManArray(), 
                                                   *gp->giveCoordinates(), 0.0) );
  weight      = gp->giveWeight();
  thickness   = this->giveElement()->giveCrossSection()->give('t');
  volume      = determinant * weight * thickness;
    
  return volume;
}

BsplinePlaneStressElement::BsplinePlaneStressElement (int n, Domain *aDomain) : IGAElement (n, aDomain), PlaneStressStructuralElementEvaluator(), interpolation() {}
