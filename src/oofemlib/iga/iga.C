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

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#include "oofegutils.h"
#endif

BSplineInterpolation::~BSplineInterpolation() {
  delete [] degree;
  delete [] numberOfControllPoints;
  delete [] numberOfKnotSpans;
  // delete also knotVector and knotMultiplicity (huhu)
}

IRResultType 
BSplineInterpolation::initializeFrom(InputRecord *ir) {
  const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
  IRResultType result;                            // Required by IR_GIVE_FIELD macro

  IntArray degree_tmp, knotMultiplicity_tmp;
  FloatArray knotVector_tmp;
  double *knotVec, knotVal;
  int i, j, n, sum, pos, *knotMul;
  const char *IFT_knotVectorString[3] = {"knotvectoru", "knotvectorv", "knotvectorw"};
  InputFieldType IFT_knotVectorType[3] = {IFT_BSplineInterpolation_knotVectorU, 
                                          IFT_BSplineInterpolation_knotVectorV, 
                                          IFT_BSplineInterpolation_knotVectorW};
  const char *IFT_knotMultiplicityString[3] = {"knotmultiplicityu", "knotmultiplicityv", "knotmultiplicityw"};
  InputFieldType IFT_knotMultiplicityType[3] = {IFT_BSplineInterpolation_knotMultiplicityU, 
                                                IFT_BSplineInterpolation_knotMultiplicityV, 
                                                IFT_BSplineInterpolation_knotMultiplicityW};

  degree = new int [ nsd ];
  knotVector = new double * [ nsd ];
  numberOfKnotSpans = new int [ nsd ];
  knotMultiplicity = new int * [ nsd ];
  numberOfControllPoints = new int  [ nsd ];

  IR_GIVE_FIELD(ir, degree_tmp, IFT_BSplineInterpolation_degree, "degree"); // Macro
  if(degree_tmp.giveSize() != nsd){
    OOFEM_ERROR("BSplineInterpolation::initializeFrom - degree size mismatch");
  }
  for(i=0; i<nsd; i++)degree[i] = degree_tmp.at(i+1);

  for(n=0; n<nsd; n++){
    knotVector_tmp.resize(0);
    IR_GIVE_FIELD(ir, knotVector_tmp, IFT_knotVectorType[n], IFT_knotVectorString[n]); // Macro
    if(knotVector_tmp.giveSize() < 2){
      OOFEM_ERROR2("BSplineInterpolation::initializeFrom - invalid size of knot vector %s", IFT_knotVectorString[n]);
    }

    // check for monotonicity of knot vector without multiplicity
    knotVal = knotVector_tmp.at(1);
    for(i=1; i<knotVector_tmp.giveSize(); i++){
      if(knotVector_tmp.at(i+1) <= knotVal){
        OOFEM_ERROR2("BSplineInterpolation::initializeFrom - knot vector %s is not monotonic", IFT_knotVectorString[n]);
      }
      knotVal = knotVector_tmp.at(i+1);
    }

    numberOfKnotSpans[n] = knotVector_tmp.giveSize() - 1;

    knotMul = knotMultiplicity[n] = new int [ knotVector_tmp.giveSize() ];

    knotMultiplicity_tmp.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, knotMultiplicity_tmp, IFT_knotMultiplicityType[n], IFT_knotMultiplicityString[n]); // Macro
    if(knotMultiplicity_tmp.giveSize() == 0){
      // default multiplicity
      for(i=1; i<knotVector_tmp.giveSize()-1; i++)knotMul[i] = 1;
      knotMultiplicity_tmp.resize(knotVector_tmp.giveSize());
    }
    else{
      if(knotMultiplicity_tmp.giveSize() != knotVector_tmp.giveSize()){
        OOFEM_ERROR2("BSplineInterpolation::initializeFrom - knot multiplicity %s size mismatch", IFT_knotMultiplicityString[n]);
      }
      knotMul[0] = knotMultiplicity_tmp.at(1);
      knotMul[knotVector_tmp.giveSize()-1] = knotMultiplicity_tmp.at(knotVector_tmp.giveSize());
      for(i=1; i<knotMultiplicity_tmp.giveSize()-1; i++){
        knotMul[i] = knotMultiplicity_tmp.at(i+1);
        // check for multiplicity range
        if(knotMul[i] < 1 || knotMul[i] > degree[n]){
          OOFEM_ERROR3("BSplineInterpolation::initializeFrom - knot multiplicity %s out of range - value %d", 
                       IFT_knotMultiplicityString[n], knotMul[i]);
        }
      }
      if(knotMul[0] != degree[n] + 1){
        OOFEM_LOG_RELEVANT("Multiplicity of the first knot in knot vector %s changed to %d\n", IFT_knotVectorString[n], degree[n] + 1);
      }
      if(knotMul[knotVector_tmp.giveSize()-1] != degree[n] + 1){
        OOFEM_LOG_RELEVANT("Multiplicity of the last knot in knot vector %s changed to %d\n", IFT_knotVectorString[n], degree[n] + 1);
      }
    }

    // multiplicity of the 1st and last knot set to degree + 1
    knotMul[0] = knotMul[knotVector_tmp.giveSize()-1] = degree[n]+1;

    // sum the size of knot vector with multiplicity values
    sum = 0;
    for(i=0; i<knotVector_tmp.giveSize(); i++)sum+=knotMul[i];
			
    knotVec = knotVector[n] = new double [ sum ];

    // fill knot vector including multiplicity values
    pos = 0;
    for(i=0; i<knotVector_tmp.giveSize(); i++){
      for(j=0; j<knotMul[i]; j++){
        knotVec[pos++] = knotVector_tmp.at(i+1);
      }
    }

    numberOfControllPoints[n] = sum - degree[n] - 1;
  }
  return IRRT_OK;
}



void BSplineInterpolation::evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time){
  FEIIGAElementGeometryWrapper* gw = (FEIIGAElementGeometryWrapper*) &cellgeo;

  int i, l,k,c;
  FloatArray N[nsd];
  IntArray span(nsd);

  if (gw->knotSpan) {
    span = *gw->knotSpan;
  } else {
    for (i=0; i< nsd; i++) {
      span(i) = this->findSpan (numberOfControllPoints[i], degree[i], lcoords(i), knotVector[i]);
    }
  } 
  
  for (i=0; i< nsd; i++) {
    //span(i) = this->findSpan (numberOfControllPoints[i], degree[i], lcoords(i), knotVector[i]);
    this->basisFuns (N[i], span.at(i+1), lcoords(i), degree[i], knotVector[i]);
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

 
void BSplineInterpolation::evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) {

  FEIIGAElementGeometryWrapper* gw = (FEIIGAElementGeometryWrapper*) &cellgeo;

  FloatMatrix jacobian(nsd,nsd);
  FloatMatrix ders[nsd];
  FloatArray temp(nsd);
  double Jacob;
  int count = 1, cnt, i, l, k, uind, vind;
  IntArray span(nsd);

  if (gw->knotSpan) {
    span = *gw->knotSpan;
  } else {
    for (i=0; i< nsd; i++) {
      span(i) = this->findSpan (numberOfControllPoints[i], degree[i], lcoords(i), knotVector[i]);
    }
  } 


  for (i=0; i< nsd; i++) {
    this->DersBasisFuns(1, lcoords(i), span.at(i+1), degree[i], knotVector[i], ders[i]);
  }
  count =giveNumberOfKnotBasisFunctions ();

  jacobian.zero();
  answer.resize(count, nsd);
  //answer.zero();
  if (nsd == 2) {
    uind = span.at(1)-degree[0];
    for (l=0;l<=degree[1]; l++) {
      temp.zero();
      vind = span.at(2)-degree[1]+l;
      for (k=0; k<=degree[0]; k++) {
        temp(0) += ders[0](1,k)*cellgeo.giveVertexCoordinates((vind)*numberOfControllPoints[0]+uind+k+1)->at(1);  // HAHA see local2global
        temp(1) += ders[0](1,k)*cellgeo.giveVertexCoordinates((vind)*numberOfControllPoints[0]+uind+k+1)->at(2);  // HAHA coords
      }
      jacobian(0,0) += ders[1](0,l)*temp(0);
      jacobian(0,1) += ders[1](0,l)*temp(1);
    }
    for (l=0;l<=degree[1]; l++) {
      temp.zero();
      vind = span.at(2)-degree[1]+l;
      for (k=0; k<=degree[0]; k++) {
        temp(0) += ders[0](0,k)*cellgeo.giveVertexCoordinates((vind)*numberOfControllPoints[0]+uind+k+1)->at(1);  // HAHA see local2global
        temp(1) += ders[0](0,k)*cellgeo.giveVertexCoordinates((vind)*numberOfControllPoints[0]+uind+k+1)->at(2);  // HAHA coords
      }
      jacobian(1,0) += ders[1](1,l)*temp(0);
      jacobian(1,1) += ders[1](1,l)*temp(1);
    }
    Jacob = jacobian.giveDeterminant();
    
    if(fabs(Jacob) < 1.0e-10){
      OOFEM_ERROR ("evaldNdx - zero Jacobian");
    }
		
    cnt=0;
    for (l=0;l<=degree[1]; l++) {
      for (k=0; k<=degree[0]; k++) {
        temp(0) = ders[0](1,k)*ders[1](0,l);
        temp(1) = ders[0](0,k)*ders[1](1,l);
        answer(cnt,0) = (jacobian(1,1)*temp(0)-jacobian(0,1)*temp(1)) / Jacob;
        answer(cnt,1) = (-jacobian(1,0)*temp(0)+jacobian(0,0)*temp(1)) / Jacob;
        cnt++;
      }
    }
  }
  else {
    OOFEM_ERROR2 ("evaldNdx not implemented for nsd = %d", nsd);
  }
}

void BSplineInterpolation::local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) {
  /* Based on SurfacePoint A3.5 implementation*/
  FEIIGAElementGeometryWrapper* gw = (FEIIGAElementGeometryWrapper*) &cellgeo;
  int i,l,k;
  FloatArray temp(nsd);
  FloatArray N[nsd];
  answer.resize(nsd);
  
  IntArray span(nsd);

  if (gw->knotSpan) {
    span = *gw->knotSpan;
  } else {
    for (i=0; i< nsd; i++) {
      span(i) = this->findSpan (numberOfControllPoints[i], degree[i], lcoords(i), knotVector[i]);
    }
  } 


  for (i=0; i< nsd; i++) {
    this->basisFuns (N[i], span(i), lcoords(i), degree[i], knotVector[i]); 
  }
  if (nsd == 2) {
    int uind, vind;
    
    uind = span(0)-degree[0]; 
    answer.zero();
    for (l=0;l<=degree[1]; l++) {
      temp.zero();
      vind = span(1)-degree[1]+l;  
      for (k=0; k<=degree[0]; k++) {
        //temp(0) = temp(0) + N[0][k]*d->giveNode(nodes.at((uind+k)*numberOfControllPoints[0]+vind))->giveCoordinate(1);
        //temp(1) = temp(1) + N[0][k]*d->giveNode(nodes.at((uind+k)*numberOfControllPoints[0]+vind))->giveCoordinate(2);
        temp(0) += N[0](k)*cellgeo.giveVertexCoordinates((vind)*numberOfControllPoints[0]+uind+k+1)->at(1);
        temp(1) += N[0](k)*cellgeo.giveVertexCoordinates((vind)*numberOfControllPoints[0]+uind+k+1)->at(2);
      }
      answer(0) += N[1](l)*temp(0);
      answer(1) += N[1](l)*temp(1);
    }
  } else {
    OOFEM_ERROR2 ("local2global not implemented for nsd = %d", nsd);
  }
}

double BSplineInterpolation::giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) {

  FEIIGAElementGeometryWrapper* gw = (FEIIGAElementGeometryWrapper*) &cellgeo;
  FloatMatrix jacobian(nsd,nsd);
  FloatMatrix ders[nsd];
  FloatArray temp(nsd);
  double Jacob;
  int i, l, k, uind, vind;

  IntArray span(nsd);

  if (gw->knotSpan) {
    span = *gw->knotSpan;
  } else {
    for (i=0; i< nsd; i++) {
      span(i) = this->findSpan (numberOfControllPoints[i], degree[i], lcoords(i), knotVector[i]);
    }
  } 

  for (i=0; i< nsd; i++) {
    this->DersBasisFuns(1, lcoords(i), span.at(i+1), degree[i], knotVector[i], ders[i]);
  }

  jacobian.zero();
  if (nsd == 2) {
    uind = span.at(1)-degree[0];
    for (l=0;l<=degree[1]; l++) {
      temp.zero();
      vind = span.at(2)-degree[1]+l;
      for (k=0; k<=degree[0]; k++) {
        temp(0) += ders[0](1,k)*cellgeo.giveVertexCoordinates((vind)*numberOfControllPoints[0]+uind+k+1)->at(1);  // HAHA see local2global
        temp(1) += ders[0](1,k)*cellgeo.giveVertexCoordinates((vind)*numberOfControllPoints[0]+uind+k+1)->at(2);  // HAHA coords
      }
      jacobian(0,0) += ders[1](0,l)*temp(0);
      jacobian(0,1) += ders[1](0,l)*temp(1);
    }
    for (l=0;l<=degree[1]; l++) {
      temp.zero();
      vind = span.at(2)-degree[1]+l;
      for (k=0; k<=degree[0]; k++) {
        temp(0) += ders[0](0,k)*cellgeo.giveVertexCoordinates((vind)*numberOfControllPoints[0]+uind+k+1)->at(1);  // HAHA see local2global
        temp(1) += ders[0](0,k)*cellgeo.giveVertexCoordinates((vind)*numberOfControllPoints[0]+uind+k+1)->at(2);  // HAHA coords
      }
      jacobian(1,0) += ders[1](1,l)*temp(0);
      jacobian(1,1) += ders[1](1,l)*temp(1);
    }
    Jacob = jacobian.giveDeterminant();
    
    if(fabs(Jacob) < 1.0e-10){
      OOFEM_ERROR ("evaldNdx - zero Jacobian");
    }
    
  } else {
    OOFEM_ERROR2 ("evaldNdx not implemented for nsd = %d", nsd);
  }
  return Jacob;
}

int BSplineInterpolation::giveKnotBasisFuncMask (const IntArray& knotSpan, IntArray& mask) {
  int c, _i, _j, _iindx, _jindx;
  
  if (nsd == 2) {
    mask.resize((degree[0]+1)*(degree[1]+1));
    c=1;
    for (_i=0; _i<=degree[1]; _i++) {
      _iindx = (_i+knotSpan(1)-degree[1]);
      for (_j=0; _j<=degree[0]; _j++) {
        _jindx = (_j+knotSpan(0)-degree[0]);
        mask.at(c++) = _jindx+(_iindx)*numberOfControllPoints[0]+1;
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
  //
  // Based on Algorithm A2.2 (p. 70)
  //
  int j, r;
  FloatArray right(p+1);
  FloatArray left(p+1);
  double saved, temp;
  
  N.resize(p+1);
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
  const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
  IRResultType result;                            // Required by IR_GIVE_FIELD macro


  int indx, ui,vi, i, numberOfGaussPoints=1;
  double du,dv;
  const FloatArray *gpcoords;
  FloatArray newgpcoords;

  Element::initializeFrom (ir); // read nodes , material, cross section
  // set number of dofmanagers
  this->numberOfDofMans = dofManArray.giveSize();
  this->giveInterpolation()->initializeFrom (ir); // read geometry

  int ** const knotMultiplicity = this->giveInterpolation()->giveKnotMultiplicity(); 
  double ** const knotVector = this->giveInterpolation()->giveKnotVector();
  IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_IGAElement_NIP, "nip"); // Macro
  
  // generate individual IntegrationElements; one for each nonzero knot span
  if (this->giveNsd() == 2) {
    // HUHU mame pristup na numberOfKnotSpans[]?
    this->numberOfIntegrationRules = this->giveInterpolation()->giveNumberOfKnotSpans(1)*this->giveInterpolation()->giveNumberOfKnotSpans(2);
    integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];
    newgpcoords.resize(2);
    indx = -1; 
    IntArray knotSpan(2); 
    knotSpan.at(2) = -1; 
    for (vi=0; vi<this->giveInterpolation()->giveNumberOfKnotSpans(2); vi++) {
      knotSpan.at(2) += knotMultiplicity[1][vi];
      knotSpan.at(1) = -1;
      for (ui=0; ui<this->giveInterpolation()->giveNumberOfKnotSpans(1); ui++) {
        indx++;
        knotSpan.at(1)+=knotMultiplicity[0][ui];
        integrationRulesArray [ indx ] = new IGA_IntegrationElement (indx, this, knotSpan);
        integrationRulesArray [ indx ] ->setUpIntegrationPoints(_Square, numberOfGaussPoints, _PlaneStress); // HUHU _PlaneStress
        // remap local subelement gp coordinates into knot span coordinates and update integration weight 
        for (i=0; i<numberOfGaussPoints; i++) {
          gpcoords = integrationRulesArray [ indx ]->getIntegrationPoint(i)->giveCoordinates();
          du = knotVector[0][knotSpan.at(1)+1]-knotVector[0][knotSpan.at(1)];
          dv = knotVector[1][knotSpan.at(2)+1]-knotVector[1][knotSpan.at(2)];
          newgpcoords.at(1) = knotVector[0][knotSpan.at(1)]+du*(gpcoords->at(1)/2.0+0.5);
          newgpcoords.at(2) = knotVector[1][knotSpan.at(2)]+dv*(gpcoords->at(2)/2.0+0.5);
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




void
IGAElement::drawRawGeometry(oofegGraphicContext &gc)
{

    WCRec p [ 4 ];
    GraphicObj *go;
    FEInterpolation* interp = this->giveInterpolation();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(TRUE);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    EASValsSetFillStyle(FILL_HOLLOW);

    numberOfIntegrationRules = this->giveNumberOfIntegrationRules();
    double** const  knotVector = interp->giveKnotVector();
    const IntArray *span;
    IntegrationRule* iRule;
    int ir, j, nsd = this->giveNsd();
    FloatArray c[4], cg[4];

    if (nsd == 2) {
      for (j=0; j<4; j++) {
        c[j].resize(2);
        cg[j].resize(2);
      }
    } else {
      OOFEM_ERROR2 ("drawRawGeometry: not implemented for nsd = %d", nsd); 
    }


    // loop over individual integration rules (i.e., knot spans)
    for (ir=0; ir < numberOfIntegrationRules; ir++) {
      iRule = this->giveIntegrationRule(ir);
      span = iRule->giveKnotSpan();
      if (nsd == 2) {

        // divide span locally to get finer geometry rep.
        int i,j,k, nseg = 4;
        double du = (knotVector[0][span->at(1)+1] - knotVector[0][span->at(1)])/ nseg;
        double dv = (knotVector[1][span->at(2)+1] - knotVector[1][span->at(2)])/ nseg;
        for (i=1; i<=4; i++) {
          for (j=1; j<=4; j++) {

            c[0].at(1) = knotVector[0][span->at(1)] + du*(i-1);
            c[0].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
            c[1].at(1) = knotVector[0][span->at(1)] + du*i;
            c[1].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
            c[2].at(1) = knotVector[0][span->at(1)] + du*i;
            c[2].at(2) = knotVector[1][span->at(2)] + dv*j;
            c[3].at(1) = knotVector[0][span->at(1)] + du*(i-1);
            c[3].at(2) = knotVector[1][span->at(2)] + dv*j;
            
            for (k=0;k<4; k++) {
              interp->local2global (cg[k], c[k], FEIIGAElementGeometryWrapper(this, iRule->giveKnotSpan()), 0.0);
              p [ k ].x = ( FPNum ) cg[k].at(1);
              p [ k ].y = ( FPNum ) cg[k].at(2);
              p [ k ].z = 0.;
            }
            go =  CreateQuad3D(p);
            EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
            EGAttachObject(go, ( EObjectP ) this);
            EMAddGraphicsToModel(ESIModel(), go);
          }
        }
      }

    } // end loop over knot spans (irules)
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
      elem->giveDofManager(mask.at(i))->giveLocationArray(nodeDofIDMask, nodalArray);
      answer.followedBy(nodalArray);
    }
    return 1;
  } else {
    return 0;
  }
}

int StructuralElementEvaluator::giveIntegrationElementLocalCodeNumbers (IntArray& answer, Element* elem, 
                                                                        IntegrationRule* ie, EquationID ut) {
  int i;
  IntArray mask, nodeDofIDMask, nodalArray;
  int dofmandof ;

  // get number of dofs in node
  elem->giveDofManDofIDMask (1, ut, nodeDofIDMask);
  dofmandof=nodeDofIDMask.giveSize();

  // first evaluate nonzero basis function mask


  if (elem->giveInterpolation()->hasSubPatchFormulation()) {
    IGA_IntegrationElement *ee = (IGA_IntegrationElement *)ie;
    elem->giveInterpolation()->giveKnotBasisFuncMask (*ee->giveKnotSpan(), mask) ;
    // loop over nonzero shape functions and assemble localization array
    answer.resize(0);
    for (i=1; i<=mask.giveSize(); i++) {
      nodalArray.resize(nodeDofIDMask.giveSize());
      nodalArray.at(1) = dofmandof*(mask.at(i)-1)+1;
      nodalArray.at(2) = dofmandof*(mask.at(i)-1)+2;
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
  
  interp->evalN (N, *gp->giveCoordinates(), FEIIGAElementGeometryWrapper(gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan()), 0.0);
  
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
  interp->evaldNdx (d, *gp->giveCoordinates(), 
		    FEIIGAElementGeometryWrapper(gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan()), 0.0);
  
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

void
StructuralElementEvaluator :: computeBcLoadVectorAt(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// Computes the load vector due to the boundary conditions acting on the
// receiver's nodes, at stepN. Returns NULL if this array contains only
// zeroes.
{
    FloatArray d, dp;
    FloatMatrix s;
    Element* elem = this->giveElement();
    int numberOfDofMans = elem->giveNumberOfDofManagers();
    /*
     * this -> computeVectorOfPrescribed(DisplacementVector,TotalMode,stepN, d) ;
     * if ((stepN->giveLoadResponseMode()==IncrementOfLoad) && (!stepN->isTheFirstStep())) {
     * this -> computeVectorOfPrescribed(DisplacementVector,TotalMode,stepN->givePreviousStep(), dp);
     * d.substract (dp);
     * //delete dp;
     * }
     */
    this->computeVectorOfPrescribed(EID_MomentumBalance, mode, stepN, d);
    //this -> computeVectorOfPrescribed(DisplacementVector,umode,stepN, d) ;

    if ( d.containsOnlyZeroes() ) {
        answer.resize(0);
    } else {
        this->computeStiffnessMatrix(s, TangentStiffness, stepN);
        answer.beProductOf(s, d);
        answer.negated();
    }

    // delete d ;

    // if engngmodel supports dynamic change of static system
    // we must test if element has not been removed in previous step 
    // if not, we must also test if there was previous BC on some DOF and now it is released.
    // if was, it is necessary to load it by reaction force.
    if ( elem->giveDomain()->giveEngngModel()->requiresUnknownsDictionaryUpdate() ) {

        FloatArray prevInternalForces;
        IntArray elementNodeMask, dofMask;
        DofManager *nodeI;
        Dof *dofJ;
        int nDofs, i, j, k = 0;

        if ( ( mode == VM_Incremental ) && ( !stepN->isTheFirstStep() ) ) {

          for ( i = 1; i <= numberOfDofMans; i++ ) {
            nodeI = elem->giveDofManager(i);
            elem->giveDofManDofIDMask(i, EID_MomentumBalance, elementNodeMask);
            nodeI->giveDofArray(elementNodeMask, dofMask);
            nDofs = dofMask.giveSize();
            for ( j = 1; j <= nDofs; j++ ) {
              dofJ = nodeI->giveDof( dofMask.at(j) );
              k++;
              if ( !dofJ->hasBc(stepN) && dofJ->hasBc( stepN->givePreviousStep() ) ) {
                if ( prevInternalForces.giveSize() == 0 ) {
                  // allocate and compute only if needed
                  // use updated gp record
                  this->giveInternalForcesVector(prevInternalForces,
                                                 stepN->givePreviousStep(), 1);
                }
                
                // check for allocated answer
                if ( answer.giveSize() == 0 ) {
                  answer.resize( elem->computeNumberOfDofs(EID_MomentumBalance) );
                  answer.zero();
                }
                
                // add element part of reaction  to load vector
                answer.at(k) -= prevInternalForces.at(k);
              }
            }
            
            //delete elementNodeMask;
            // delete dofMask;
          }
        }
    }
    
    return;
}


void
StructuralElementEvaluator :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step stepN. The nature of these strains depends
// on the element's type.
{
  int i;
    FloatMatrix b;
    FloatArray u, ur;
    Element* elem=this->giveElement();

    if (!this->isActivated(stepN)) {
      answer.resize(elem->giveCrossSection()->giveIPValueSize(IST_StrainTensor, gp));
      answer.zero();
      return;}

    this->computeBMatrixAt(b, gp);
    elem->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);

    /*
    // substract initial displacements, if defined
    if (initialDisplacements) u.substract(initialDisplacements);
    */
    if ( this->updateRotationMatrix() ) {
        u.rotatedWith(this->rotationMatrix, 'n');
    }

    // get local code numbers corresponding to ir
    IntArray lc;
    this->giveIntegrationElementLocalCodeNumbers (lc, elem, gp->giveIntegrationRule(), EID_MomentumBalance);
    ur.resize(b.giveNumberOfColumns());
    for (i=1; i<=lc.giveSize(); i++) ur.at(i) = u.at(lc.at(i));

    answer.beProductOf(b, ur);
    return;
}


void
StructuralElementEvaluator :: computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
// Computes the vector containing the stresses at the Gauss point gp of
// the receiver, at time step stepN. The nature of these stresses depends
// on the element's type.
// this version assumes TOTAL LAGRANGE APPROACH
{
    /*
     * StructuralCrossSection* cs = (StructuralCrossSection*) this->giveCrossSection();
     * FloatArray totalEpsilon;
     * // FloatArray *help;
     *
     *
     * this->computeStrainVector(totalEpsilon, gp,stepN) ;
     * cs->giveRealStresses (answer, ReducedForm, gp,totalEpsilon,stepN);
     *
     * return ;
     */
    FloatArray Epsilon;
    Element* elem=this->giveElement();
    StructuralCrossSection *cs = ( StructuralCrossSection * ) elem->giveCrossSection();

    this->computeStrainVector(Epsilon, gp, stepN);
    cs->giveRealStresses(answer, ReducedForm, gp, Epsilon, stepN);

    return;
}

void
StructuralElementEvaluator :: updateInternalState(TimeStep *stepN)
// Updates the receiver at end of step.
{
    int i, j;
    IntegrationRule *iRule;
    FloatArray stress;
    Element *elem = this->giveElement();
    
    // force updating strains & stresses
    for ( i = 0; i < elem->giveNumberOfIntegrationRules(); i++ ) {
      iRule = elem->giveIntegrationRule(i);
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
          computeStressVector(stress, iRule->getIntegrationPoint(j), stepN);
        }
    }
}




void
StructuralElementEvaluator :: computeNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// Computes the load vector of the receiver, at stepN.
{
    FloatArray helpLoadVector;

    answer.resize(0);

    // test for deactivation of receiver
    if ( ( mode == VM_Incremental ) && ( !stepN->isTheFirstStep() ) ) {
      if (isActivated(stepN->givePreviousStep()) && !isActivated(stepN)) {
        // use updated gp record
        this->giveInternalForcesVector(answer, stepN->givePreviousStep(), 1);
      }
    }
    if (!this->isActivated(stepN)) return;
    
    /*
    this->computePrescribedStrainLoadVectorAt(helpLoadVector, stepN, mode);
    if ( helpLoadVector.giveSize() ) {
      answer.add(helpLoadVector);
    }
    */
    
    this->computeBcLoadVectorAt(helpLoadVector, stepN, mode);
    if ( helpLoadVector.giveSize() ) {
      answer.add(helpLoadVector);
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
  
  FloatMatrix *m = &answer;
  if (elem->giveInterpolation()->hasSubPatchFormulation()) m = &temp;
  
  numberOfIntegrationRules = elem->giveNumberOfIntegrationRules();
  // loop over individual integration rules
  for (ir=0; ir < numberOfIntegrationRules; ir++) {
    m->resize(0,0);
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
        m->plusProductSymmUpper(bj, dbj, dV);
      } else {
        m->plusProductUnsym(bj, dbj, dV);
      }
    }
    
    if ( matStiffSymmFlag ) {
      m->symmetrized();
    }
    // localize irule contribution into element matrix
    if (this->giveIntegrationElementLocalCodeNumbers (irlocnum, elem, iRule, EID_MomentumBalance)) 
      answer.assemble (*m, irlocnum);
    
  } // end loop over irules
  
  if ( this->updateRotationMatrix() ) {
    answer.rotatedWith(* this->rotationMatrix);
  }
  return;
}
  
double PlaneStressStructuralElementEvaluator::computeVolumeAround(GaussPoint *gp) { 
  double determinant, weight, thickness, volume;
  determinant = fabs( this->giveElement()->giveInterpolation()
                      ->giveTransformationJacobian(*gp->giveCoordinates(), 
						   FEIIGAElementGeometryWrapper(this->giveElement(), 
										gp->giveIntegrationRule()->giveKnotSpan()),
						   0.0) );
  weight      = gp->giveWeight();
  thickness   = this->giveElement()->giveCrossSection()->give('t');
  volume      = determinant * weight * thickness;
    
  return volume;
}

BsplinePlaneStressElement::BsplinePlaneStressElement (int n, Domain *aDomain) : IGAElement (n, aDomain), PlaneStressStructuralElementEvaluator(), interpolation(2) {}


void BsplinePlaneStressElement :: drawScalar(oofegGraphicContext &context)
{
    int indx;
    WCRec p [ 4 ];
    GraphicObj *go;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    double s [ 4 ];
    IntArray map;
    FEInterpolation* interp = this->giveInterpolation();
    FloatArray val;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    numberOfIntegrationRules = this->giveNumberOfIntegrationRules();
    double** const  knotVector = interp->giveKnotVector();
    const IntArray *span;
    IntegrationRule* iRule;
    int ir, j, nsd = this->giveNsd();
    FloatArray c[4], cg[4];

    if (nsd == 2) {
      for (j=0; j<4; j++) {
        c[j].resize(2);
        cg[j].resize(2);
      }
    } else {
      OOFEM_ERROR2 ("drawRawGeometry: not implemented for nsd = %d", nsd); 
    }

    this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
    if ( ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
      return;
    }


    // loop over individual integration rules (i.e., knot spans)
    for (ir=0; ir < numberOfIntegrationRules; ir++) {
      iRule = this->giveIntegrationRule(ir);
      span = iRule->giveKnotSpan();
      if (nsd == 2) {

        // divide span locally to get finer geometry rep.
        int i,j,k, nseg = 4;
        double du = (knotVector[0][span->at(1)+1] - knotVector[0][span->at(1)])/ nseg;
        double dv = (knotVector[1][span->at(2)+1] - knotVector[1][span->at(2)])/ nseg;
        for (i=1; i<=4; i++) {
          for (j=1; j<=4; j++) {

            c[0].at(1) = knotVector[0][span->at(1)] + du*(i-1);
            c[0].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
            c[1].at(1) = knotVector[0][span->at(1)] + du*i;
            c[1].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
            c[2].at(1) = knotVector[0][span->at(1)] + du*i;
            c[2].at(2) = knotVector[1][span->at(2)] + dv*j;
            c[3].at(1) = knotVector[0][span->at(1)] + du*(i-1);
            c[3].at(2) = knotVector[1][span->at(2)] + dv*j;
            
            for (k=0;k<4; k++) {
              interp->local2global (cg[k], c[k], FEIIGAElementGeometryWrapper(this, iRule->giveKnotSpan()), 0.0);
              p [ k ].x = ( FPNum ) cg[k].at(1);
              p [ k ].y = ( FPNum ) cg[k].at(2);
              p [ k ].z = 0.;
              // create a dummy ip's
              FloatArray *cc = new FloatArray;
              cc->beCopyOf (&c[k]);
              GaussPoint gp(iRule, 999, cc, 1.0, _PlaneStress);            
              this->computeStrainVector (val, &gp, tStep);
              s[k]=val.at(indx);
            }
            go =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
            EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
            EGAttachObject(go, ( EObjectP ) this);
            EMAddGraphicsToModel(ESIModel(), go);
          }
        }
      }

    } // end loop over knot spans (irules)
}
