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
  oofem nodes - control points (coordinates ) + dofs
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


// optimized version of A4.4 for d=1
#define OPTIMIZED_VERSION_A4dot4


// jaky je rozdil mezi FloatArray ders[nsd] a FloatArray ders(nsd) ==> prvni je pole (o velikosti nsd) FloatArrays o nulove velikosti
//                                                                     druhe je FloatArray inicializovane na velikost nsd
// jaky je rozdil mezi span(0) a spat.at(1) ==> zadny
// jaky je rozdil mezi span(0) a span[0] ==> prvni implementovano jako C++ objekt Array, druhe jako ceckove pole


// there is a duplication of code for evaluation of jacobian in methods evaldNdx and giveTransformationJacobian;
// this is to avoid duplication of calculation of basis functions and their derivatives;
// modification is desirable;



// BSpline

BSplineInterpolation::~BSplineInterpolation() {
	int i;

  delete [] degree;
  delete [] numberOfControlPoints;
  delete [] numberOfKnotSpans;

	for (i=0;i<nsd;i++){
		delete [] knotVector[i];
	}
	delete [] knotValues;
	delete [] knotMultiplicity;
	delete [] knotVector;
}


IRResultType 
BSplineInterpolation::initializeFrom(InputRecord *ir) {
  const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
  IRResultType result;                   // Required by IR_GIVE_FIELD macro

  IntArray degree_tmp;
  double *knotVec, knotVal;
  int i, j, n, sum, pos, size;
  const char *IFT_knotVectorString[3] = {"knotvectoru", "knotvectorv", "knotvectorw"};
  InputFieldType IFT_knotVectorType[3] = {IFT_BSplineInterpolation_knotVectorU, 
                                          IFT_BSplineInterpolation_knotVectorV, 
                                          IFT_BSplineInterpolation_knotVectorW};
  const char *IFT_knotMultiplicityString[3] = {"knotmultiplicityu", "knotmultiplicityv", "knotmultiplicityw"};
  InputFieldType IFT_knotMultiplicityType[3] = {IFT_BSplineInterpolation_knotMultiplicityU, 
                                                IFT_BSplineInterpolation_knotMultiplicityV, 
                                                IFT_BSplineInterpolation_knotMultiplicityW};

	knotValues = new FloatArray [ nsd ];
	knotMultiplicity = new IntArray [ nsd ];
  degree = new int [ nsd ];
  knotVector = new double * [ nsd ];
  numberOfKnotSpans = new int [ nsd ];
  numberOfControlPoints = new int  [ nsd ];

  IR_GIVE_FIELD(ir, degree_tmp, IFT_BSplineInterpolation_degree, "degree"); // Macro
  if(degree_tmp.giveSize() != nsd){
    OOFEM_ERROR("BSplineInterpolation::initializeFrom - degree size mismatch");
  }
  for(i=0; i<nsd; i++)degree[i] = degree_tmp.at(i+1);

  for(n=0; n<nsd; n++){
    IR_GIVE_FIELD(ir, knotValues[n], IFT_knotVectorType[n], IFT_knotVectorString[n]); // Macro
		size = knotValues[n].giveSize();
    if(size < 2){
      OOFEM_ERROR2("BSplineInterpolation::initializeFrom - invalid size of knot vector %s", IFT_knotVectorString[n]);
    }

    // check for monotonicity of knot vector without multiplicity
    knotVal = knotValues[n].at(1);
    for(i=1; i<size; i++){
      if(knotValues[n].at(i+1) <= knotVal){
        OOFEM_ERROR2("BSplineInterpolation::initializeFrom - knot vector %s is not monotonic", IFT_knotVectorString[n]);
      }
      knotVal = knotValues[n].at(i+1);
    }
		
    IR_GIVE_OPTIONAL_FIELD(ir, knotMultiplicity[n], IFT_knotMultiplicityType[n], IFT_knotMultiplicityString[n]); // Macro
    if(knotMultiplicity[n].giveSize() == 0){
      // default multiplicity
      knotMultiplicity[n].resize(size);
			// skip the first and last one
      for(i=1; i<size-1; i++)knotMultiplicity[n].at(i+1) = 1;
    }
    else{
      if(knotMultiplicity[n].giveSize() != size){
        OOFEM_ERROR2("BSplineInterpolation::initializeFrom - knot multiplicity %s size mismatch", IFT_knotMultiplicityString[n]);
      }
			// check for multiplicity range (skip the first and last one)
			for(i=1; i<size-1; i++){
        if(knotMultiplicity[n].at(i+1) < 1 || knotMultiplicity[n].at(i+1) > degree[n]){
          OOFEM_ERROR3("BSplineInterpolation::initializeFrom - knot multiplicity %s out of range - value %d", 
                       IFT_knotMultiplicityString[n], knotMultiplicity[n].at(i+1));
        }
      }
			// check for multiplicity of the first and last one
      if(knotMultiplicity[n].at(1) != degree[n] + 1){
        OOFEM_LOG_RELEVANT("Multiplicity of the first knot in knot vector %s changed to %d\n", IFT_knotVectorString[n], degree[n] + 1);
      }
      if(knotMultiplicity[n].at(size) != degree[n] + 1){
        OOFEM_LOG_RELEVANT("Multiplicity of the last knot in knot vector %s changed to %d\n", IFT_knotVectorString[n], degree[n] + 1);
      }
    }

    // multiplicity of the 1st and last knot set to degree + 1
		knotMultiplicity[n].at(1) = knotMultiplicity[n].at(size) = degree[n]+1;

    // sum the size of knot vector with multiplicity values
    sum = 0;
    for(i=0; i<size; i++)sum+=knotMultiplicity[n].at(i+1);
			
    knotVec = knotVector[n] = new double [ sum ];

    // fill knot vector including multiplicity values
    pos = 0;
    for(i=0; i<size; i++){
      for(j=0; j<knotMultiplicity[n].at(i+1); j++){
        knotVec[pos++] = knotValues[n].at(i+1);
      }
    }

    numberOfKnotSpans[n] = size - 1;
    numberOfControlPoints[n] = sum - degree[n] - 1;
  }
  return IRRT_OK;
}



void BSplineInterpolation::evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) {
  FEIIGAElementGeometryWrapper* gw = (FEIIGAElementGeometryWrapper*) &cellgeo;
  FloatArray N[nsd];
  IntArray span(nsd);
  int i,l,k,m,c=1,count;

  if (gw->knotSpan) {
    span = *gw->knotSpan;
  } else {
    for (i=0; i< nsd; i++) {
      span(i) = this->findSpan (numberOfControlPoints[i], degree[i], lcoords(i), knotVector[i]);
    }
  } 
  
  for (i=0; i< nsd; i++) {
    this->basisFuns (N[i], span(i), lcoords(i), degree[i], knotVector[i]);
  }
  
  count =giveNumberOfKnotSpanBasisFunctions (span);
  answer.resize(count);

	if(nsd == 1){
    for (k=0; k<=degree[0]; k++) {
			answer.at(c++) = N[0](k);
    }
	}
	else if(nsd == 2){
    for (l=0; l<=degree[1]; l++) {
      for (k=0; k<=degree[0]; k++) {
        answer.at(c++) = N[0](k)*N[1](l);
      }
    }
	}
	else if(nsd == 3){
		for (m=0; m<=degree[2]; k++) {
			for (l=0; l<=degree[1]; l++) {
				for (k=0; k<=degree[0]; k++) {
					answer.at(c++) = N[0](k)*N[1](l)*N[2](m);
				}
			}
    }
	}
	else{
    OOFEM_ERROR2 ("evalN not implemented for nsd = %d", nsd);
  }
}

 
void BSplineInterpolation::evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) {
  FEIIGAElementGeometryWrapper* gw = (FEIIGAElementGeometryWrapper*) &cellgeo;
  const FloatArray *vertexCoordsPtr;
  FloatMatrix jacobian(nsd,nsd);
  FloatMatrix ders[nsd];
  IntArray span(nsd);
  double Jacob;
  int count,cnt,i,l,k,m,ind,indx,uind,vind,tind;

  if (gw->knotSpan) {
    span = *gw->knotSpan;
  } else {
    for (i=0; i< nsd; i++) {
      span(i) = this->findSpan (numberOfControlPoints[i], degree[i], lcoords(i), knotVector[i]);
    }
  } 

  for (i=0; i< nsd; i++) {
    this->dersBasisFuns(1, lcoords(i), span(i), degree[i], knotVector[i], ders[i]);
  }

  count =giveNumberOfKnotSpanBasisFunctions (span);
  answer.resize(count, nsd);
  jacobian.zero();

	if(nsd == 1){
    uind = span(0)-degree[0];
		ind = uind+1;
		for (k=0; k<=degree[0]; k++) {
			vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);
			jacobian(0,0) += ders[0](1,k)*vertexCoordsPtr->at(1);  // dx/du=sum(dNu/du*x)
		}

    Jacob = jacobian.giveDeterminant();
    
    if(fabs(Jacob) < 1.0e-10){
      OOFEM_ERROR ("evaldNdx - zero Jacobian");
    }

    cnt=0;
		for (k=0; k<=degree[0]; k++) {
			answer(cnt,0) = ders[0](1,k) / Jacob;  // dN/dx=dN/du / dx/du
			cnt++;
		}
	}
	else if(nsd == 2){
		FloatArray tmp1(nsd),tmp2(nsd);

    uind = span(0)-degree[0];
		vind = span(1)-degree[1];
		ind = vind*numberOfControlPoints[0]+uind+1;
    for (l=0; l<=degree[1]; l++) {
      tmp1.zero();
      tmp2.zero();
      for (k=0; k<=degree[0]; k++) {
				vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);
				
				tmp1(0) += ders[0](1,k)*vertexCoordsPtr->at(1);  // sum(dNu/du*x)
        tmp1(1) += ders[0](1,k)*vertexCoordsPtr->at(2);  // sum(dNu/du*y)
				
        tmp2(0) += ders[0](0,k)*vertexCoordsPtr->at(1);  // sum(Nu*x)
        tmp2(1) += ders[0](0,k)*vertexCoordsPtr->at(2);  // sum(Nu*y)
      }
			ind += numberOfControlPoints[0];
			
      jacobian(0,0) += ders[1](0,l)*tmp1(0);  // dx/du=sum(Nv*sum(dNu/du*x))
      jacobian(0,1) += ders[1](0,l)*tmp1(1);  // dy/du=sum(Nv*sum(dNu/du*y))
			
      jacobian(1,0) += ders[1](1,l)*tmp2(0);  // dx/dv=sum(dNv/dv*sum(Nu*x))
      jacobian(1,1) += ders[1](1,l)*tmp2(1);  // dy/dv=sum(dNv/dv*sum(Nu*y))
    }
		
    Jacob = jacobian.giveDeterminant();
    
    if(fabs(Jacob) < 1.0e-10){
      OOFEM_ERROR ("evaldNdx - zero Jacobian");
    }
		
    cnt=0;
    for (l=0; l<=degree[1]; l++) {
      for (k=0; k<=degree[0]; k++) {
        tmp1(0) = ders[0](1,k)*ders[1](0,l);   // dN/du=dNu/du*Nv
        tmp1(1) = ders[0](0,k)*ders[1](1,l);   // dN/dv=Nu*dNv/dv
        answer(cnt,0) = (+jacobian(1,1)*tmp1(0)-jacobian(0,1)*tmp1(1)) / Jacob; // dN/dx
        answer(cnt,1) = (-jacobian(1,0)*tmp1(0)+jacobian(0,0)*tmp1(1)) / Jacob; // dN/dy
        cnt++;
      }
    }
	}
	else if(nsd == 3){
		FloatArray tmp1(nsd),tmp2(nsd);
		FloatArray temp1(nsd),temp2(nsd),temp3(nsd);

    uind = span(0)-degree[0];
		vind = span(1)-degree[1];
		tind = span(2)-degree[2];
		ind = tind*numberOfControlPoints[0]*numberOfControlPoints[1]+vind*numberOfControlPoints[0]+uind+1;
		for (m=0; m<=degree[2]; m++) {
			temp1.zero();
			temp2.zero();
			temp3.zero();
			indx = ind;
			for (l=0; l<=degree[1]; l++) {
				tmp1.zero();
				tmp2.zero();
				for (k=0; k<=degree[0]; k++) {
					vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);
					
					tmp1(0) += ders[0](1,k)*vertexCoordsPtr->at(1);  // sum(dNu/du*x)
					tmp1(1) += ders[0](1,k)*vertexCoordsPtr->at(2);  // sum(dNu/du*y)
					tmp1(2) += ders[0](1,k)*vertexCoordsPtr->at(3);  // sum(dNu/du*z)
					
					tmp2(0) += ders[0](0,k)*vertexCoordsPtr->at(1);  // sum(Nu*x)
					tmp2(1) += ders[0](0,k)*vertexCoordsPtr->at(2);  // sum(Nu*y)
					tmp2(2) += ders[0](0,k)*vertexCoordsPtr->at(3);  // sum(Nu*y)
				}
				ind += numberOfControlPoints[0];
				
				temp1(0) += ders[1](0,l)*tmp1(0);  // sum(Nv*sum(dNu/du*x))
				temp1(1) += ders[1](0,l)*tmp1(1);  // sum(Nv*sum(dNu/du*y))
				temp1(2) += ders[1](0,l)*tmp1(2);  // sum(Nv*sum(dNu/du*z))
				
				temp2(0) += ders[1](1,l)*tmp2(0);  // sum(dNv/dv*sum(Nu*x))
				temp2(1) += ders[1](1,l)*tmp2(1);  // sum(dNv/dv*sum(Nu*y))
				temp2(2) += ders[1](1,l)*tmp2(1);  // sum(dNv/dv*sum(Nu*z))
				
				temp3(0) += ders[1](0,l)*tmp2(0);  // sum(Nv*sum(Nu*x))
				temp3(1) += ders[1](0,l)*tmp2(1);  // sum(Nv*sum(Nu*y))
				temp3(2) += ders[1](0,l)*tmp2(1);  // sum(Nv*sum(Nu*z))
			}
			ind = indx+numberOfControlPoints[0]*numberOfControlPoints[1];
			
      jacobian(0,0) += ders[2](0,m)*temp1(0);  // dx/du=sum(Nt*sum(Nv*sum(dNu/du*x)))
      jacobian(0,1) += ders[2](0,m)*temp1(1);  // dy/du=sum(Nt*sum(Nv*sum(dNu/du*y)))
      jacobian(0,2) += ders[2](0,m)*temp1(2);  // dz/du=sum(Nt*sum(Nv*sum(dNu/du*z)))
			
      jacobian(1,0) += ders[2](0,m)*temp2(0);  // dx/dv=sum(Nt*sum(dNv/dv*sum(Nu*x)))
      jacobian(1,1) += ders[2](0,m)*temp2(1);  // dy/dv=sum(Nt*sum(dNv/dv*sum(Nu*y)))
      jacobian(1,2) += ders[2](0,m)*temp2(2);  // dz/dv=sum(Nt*sum(dNv/dv*sum(Nu*z)))
			
      jacobian(2,0) += ders[2](1,m)*temp3(0);  // dx/dt=sum(dNt/dt*sum(Nv*sum(Nu*x)))
      jacobian(2,1) += ders[2](1,m)*temp3(1);  // dy/dt=sum(dNt/dt*sum(Nv*sum(Nu*y)))
      jacobian(2,2) += ders[2](1,m)*temp3(2);  // dz/dt=sum(dNt/dt*sum(Nv*sum(Nu*z)))
    }
		
    Jacob = jacobian.giveDeterminant();
    
    if(fabs(Jacob) < 1.0e-10){
      OOFEM_ERROR ("evaldNdx - zero Jacobian");
    }
		
    cnt=0;
    for (m=0; m<=degree[2]; m++) {
			for (l=0; l<=degree[1]; l++) {
				for (k=0; k<=degree[0]; k++) {
					tmp1(0) = ders[0](1,k)*ders[1](0,l)*ders[2](0,m);   // dN/du=dNu/du*Nv*Nt
					tmp1(1) = ders[0](0,k)*ders[1](1,l)*ders[2](0,m);   // dN/dv=Nu*dNv/dv*Nt
					tmp1(2) = ders[0](0,k)*ders[1](0,l)*ders[2](1,m);   // dN/dt=Nu*Nv*dNt/dt
					answer(cnt,0) = ((jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1))*tmp1(0) +
													 (jacobian(0,2)*jacobian(2,1)-jacobian(0,1)*jacobian(2,2))*tmp1(1) +
													 (jacobian(0,1)*jacobian(1,2)-jacobian(0,2)*jacobian(1,1))*tmp1(2)) / Jacob;  // dN/dx
					answer(cnt,1) = ((jacobian(1,2)*jacobian(2,0)-jacobian(1,0)*jacobian(2,2))*tmp1(0) +
													 (jacobian(0,0)*jacobian(2,2)-jacobian(0,2)*jacobian(2,0))*tmp1(1) +
													 (jacobian(0,2)*jacobian(1,0)-jacobian(0,0)*jacobian(1,2))*tmp1(2)) / Jacob;  // dN/dy
					answer(cnt,2) = ((jacobian(1,0)*jacobian(2,1)-jacobian(1,1)*jacobian(2,0))*tmp1(0) +
													 (jacobian(0,1)*jacobian(2,0)-jacobian(0,0)*jacobian(2,1))*tmp1(1) +
													 (jacobian(0,0)*jacobian(1,1)-jacobian(0,1)*jacobian(1,0))*tmp1(2)) / Jacob;  // dN/dz
					cnt++;
				}
			}
		}
	}
	else{
    OOFEM_ERROR2 ("evaldNdx not implemented for nsd = %d", nsd);
  }
}


void BSplineInterpolation::local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) {
  /* Based on SurfacePoint A3.5 implementation*/
  FEIIGAElementGeometryWrapper* gw = (FEIIGAElementGeometryWrapper*) &cellgeo;
  const FloatArray *vertexCoordsPtr;
  FloatArray N[nsd];
  IntArray span(nsd);
  int i,l,k,m,ind,indx,uind,vind,tind;

  if (gw->knotSpan) {
    span = *gw->knotSpan;
  } else {
    for (i=0; i< nsd; i++) {
      span(i) = this->findSpan (numberOfControlPoints[i], degree[i], lcoords(i), knotVector[i]);
    }
  } 

  for (i=0; i< nsd; i++) {
    this->basisFuns (N[i], span(i), lcoords(i), degree[i], knotVector[i]); 
  }

  answer.resize(nsd);
	answer.zero();

	if(nsd == 1){
    uind = span(0)-degree[0];
		ind = uind+1;
		for (k=0; k<=degree[0]; k++) {
			vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);
			answer(0) += N[0](k)*vertexCoordsPtr->at(1);
		}
	}
	else if(nsd == 2){
		FloatArray tmp(nsd);

    uind = span(0)-degree[0]; 
		vind = span(1)-degree[1];  
		ind = vind*numberOfControlPoints[0]+uind+1;
    for (l=0; l<=degree[1]; l++) {
      tmp.zero();
      for (k=0; k<=degree[0]; k++) {
				vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);

        tmp(0) += N[0](k)*vertexCoordsPtr->at(1);
        tmp(1) += N[0](k)*vertexCoordsPtr->at(2);
      }
			ind += numberOfControlPoints[0];

      answer(0) += N[1](l)*tmp(0);
      answer(1) += N[1](l)*tmp(1);
    }
	}
	else if(nsd == 3){
		FloatArray tmp(nsd),temp(nsd);

    uind = span(0)-degree[0]; 
		vind = span(1)-degree[1];  
		tind = span(2)-degree[2];  
		ind = tind*numberOfControlPoints[0]*numberOfControlPoints[1]+vind*numberOfControlPoints[0]+uind+1;
		for (m=0; m<=degree[2]; m++) {
			temp.zero();
			indx = ind;
			for (l=0; l<=degree[1]; l++) {
				tmp.zero();
				for (k=0; k<=degree[0]; k++) {
					vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);
					tmp(0) += N[0](k)*vertexCoordsPtr->at(1);
					tmp(1) += N[0](k)*vertexCoordsPtr->at(2);
					tmp(2) += N[0](k)*vertexCoordsPtr->at(3);
				}
				ind += numberOfControlPoints[0];

				temp(0) += N[1](l)*tmp(0);
				temp(1) += N[1](l)*tmp(1);
				temp(2) += N[1](l)*tmp(2);
			}
			ind = indx+numberOfControlPoints[0]*numberOfControlPoints[1];

      answer(0) += N[2](m)*temp(0);
      answer(1) += N[2](m)*temp(1);
      answer(2) += N[2](m)*temp(2);
    }
	}
	else{
    OOFEM_ERROR2 ("local2global not implemented for nsd = %d", nsd);
  }
}


double BSplineInterpolation::giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) {
  FEIIGAElementGeometryWrapper* gw = (FEIIGAElementGeometryWrapper*) &cellgeo;
  const FloatArray *vertexCoordsPtr;
  FloatMatrix jacobian(nsd,nsd);
  FloatMatrix ders[nsd];
  IntArray span(nsd);
  double Jacob;
  int i,l,k,m,indx,ind,uind,vind,tind;

  if (gw->knotSpan) {
    span = *gw->knotSpan;
  } else {
    for (i=0; i< nsd; i++) {
      span(i) = this->findSpan (numberOfControlPoints[i], degree[i], lcoords(i), knotVector[i]);
    }
  } 

  for (i=0; i< nsd; i++) {
    this->dersBasisFuns(1, lcoords(i), span(i), degree[i], knotVector[i], ders[i]);
  }

  jacobian.zero();

	if(nsd == 1){
    uind = span(0)-degree[0];
		ind = uind+1;
		for (k=0; k<=degree[0]; k++) {
			vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);
			jacobian(0,0) += ders[0](1,k)*vertexCoordsPtr->at(1);  // dx/du=sum(dNu/du*x)
		}
	}
	else if(nsd == 2){
		FloatArray tmp1(nsd),tmp2(nsd);

    uind = span(0)-degree[0];
		vind = span(1)-degree[1];
		ind = vind*numberOfControlPoints[0]+uind+1;
    for (l=0; l<=degree[1]; l++) {
      tmp1.zero();
      tmp2.zero();
      for (k=0; k<=degree[0]; k++) {
				vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);

        tmp1(0) += ders[0](1,k)*vertexCoordsPtr->at(1);  // sum(dNu/du*x)
        tmp1(1) += ders[0](1,k)*vertexCoordsPtr->at(2);  // sum(dNu/du*y)

        tmp2(0) += ders[0](0,k)*vertexCoordsPtr->at(1);  // sum(Nu*x)
        tmp2(1) += ders[0](0,k)*vertexCoordsPtr->at(2);  // sum(Nu*y)
      }
			ind += numberOfControlPoints[0];

      jacobian(0,0) += ders[1](0,l)*tmp1(0);  // dx/du=sum(Nv*sum(dNu/du*x))
      jacobian(0,1) += ders[1](0,l)*tmp1(1);  // dy/du=sum(Nv*sum(dNu/du*y))

      jacobian(1,0) += ders[1](1,l)*tmp2(0);  // dx/dv=sum(dNv/dv*sum(Nu*x))
      jacobian(1,1) += ders[1](1,l)*tmp2(1);  // dy/dv=sum(dNv/dv*sum(Nu*y))
    }
	}
	else if(nsd == 3){
		FloatArray tmp1(nsd),tmp2(nsd);
		FloatArray temp1(nsd),temp2(nsd),temp3(nsd);

    uind = span(0)-degree[0];
		vind = span(1)-degree[1];
		tind = span(2)-degree[2];
		ind = tind*numberOfControlPoints[0]*numberOfControlPoints[1]+vind*numberOfControlPoints[0]+uind+1;
		for (m=0; m<=degree[2]; m++) {
			temp1.zero();
			temp2.zero();
			temp3.zero();
			indx = ind;
			for (l=0; l<=degree[1]; l++) {
				tmp1.zero();
				tmp2.zero();
				for (k=0; k<=degree[0]; k++) {
					vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);

					tmp1(0) += ders[0](1,k)*vertexCoordsPtr->at(1);  // sum(dNu/du*x)
					tmp1(1) += ders[0](1,k)*vertexCoordsPtr->at(2);  // sum(dNu/du*y)
					tmp1(2) += ders[0](1,k)*vertexCoordsPtr->at(3);  // sum(dNu/du*z)

					tmp2(0) += ders[0](0,k)*vertexCoordsPtr->at(1);  // sum(Nu*x)
					tmp2(1) += ders[0](0,k)*vertexCoordsPtr->at(2);  // sum(Nu*y)
					tmp2(2) += ders[0](0,k)*vertexCoordsPtr->at(3);  // sum(Nu*y)
				}
				ind += numberOfControlPoints[0];

				temp1(0) += ders[1](0,l)*tmp1(0);  // sum(Nv*sum(dNu/du*x)
				temp1(1) += ders[1](0,l)*tmp1(1);  // sum(Nv*sum(dNu/du*y)
				temp1(2) += ders[1](0,l)*tmp1(2);  // sum(Nv*sum(dNu/du*z)

				temp2(0) += ders[1](1,l)*tmp2(0);  // sum(dNv/dv*sum(Nu*x)
				temp2(1) += ders[1](1,l)*tmp2(1);  // sum(dNv/dv*sum(Nu*y)
				temp2(2) += ders[1](1,l)*tmp2(1);  // sum(dNv/dv*sum(Nu*z)

				temp3(0) += ders[1](0,l)*tmp2(0);  // sum(Nv*sum(Nu*x)
				temp3(1) += ders[1](0,l)*tmp2(1);  // sum(Nv*sum(Nu*y)
				temp3(2) += ders[1](0,l)*tmp2(1);  // sum(Nv*sum(Nu*z)
			}
			ind = indx+numberOfControlPoints[0]*numberOfControlPoints[1];

      jacobian(0,0) += ders[2](0,m)*temp1(0);  // dx/du=sum(Nt*sum(Nv*sum(dNu/du*x)))
      jacobian(0,1) += ders[2](0,m)*temp1(1);  // dy/du=sum(Nt*sum(Nv*sum(dNu/du*y)))
      jacobian(0,2) += ders[2](0,m)*temp1(2);  // dz/du=sum(Nt*sum(Nv*sum(dNu/du*z)))

      jacobian(1,0) += ders[2](0,m)*temp2(0);  // dx/dv=sum(Nt*sum(dNv/dv*sum(Nu*x)))
      jacobian(1,1) += ders[2](0,m)*temp2(1);  // dy/dv=sum(Nt*sum(dNv/dv*sum(Nu*y)))
      jacobian(1,2) += ders[2](0,m)*temp2(2);  // dz/dv=sum(Nt*sum(dNv/dv*sum(Nu*z)))

      jacobian(2,0) += ders[2](1,m)*temp3(0);  // dx/dt=sum(dNt/dt*sum(Nv*sum(Nu*x)))
      jacobian(2,1) += ders[2](1,m)*temp3(1);  // dy/dt=sum(dNt/dt*sum(Nv*sum(Nu*y)))
      jacobian(2,2) += ders[2](1,m)*temp3(2);  // dz/dt=sum(dNt/dt*sum(Nv*sum(Nu*z)))
    }
	}
	else{
    OOFEM_ERROR2 ("giveTransformationJacobian not implemented for nsd = %d", nsd);
  }

	Jacob = jacobian.giveDeterminant();
    
	if(fabs(Jacob) < 1.0e-10){
		OOFEM_ERROR ("giveTransformationJacobian - zero Jacobian");
	}
    
  return Jacob;
}


int BSplineInterpolation::giveKnotSpanBasisFuncMask (const IntArray& knotSpan, IntArray& mask) {
  int size,c=1,i,j,k,iindx,jindx,kindx;

	size = giveNumberOfKnotSpanBasisFunctions(knotSpan);
	mask.resize(size);

	if(nsd == 1){
		for (i=0; i<=degree[0]; i++) {
			iindx = (i+knotSpan(0)-degree[0]);
			mask.at(c++) = iindx+1;
    }          
	}
	else if(nsd == 2){
    for (j=0; j<=degree[1]; j++) {
      jindx = (j+knotSpan(1)-degree[1]);
      for (i=0; i<=degree[0]; i++) {
        iindx = (i+knotSpan(0)-degree[0]);
        mask.at(c++) = jindx*numberOfControlPoints[0]+iindx+1;
      }
    }          
	}
	else if(nsd == 3){
    for (k=0; k<=degree[2]; k++) {
      kindx = (k+knotSpan(2)-degree[2]);
			for (j=0; j<=degree[1]; j++) {
				jindx = (j+knotSpan(1)-degree[1]);
				for (i=0; i<=degree[0]; i++) {
					iindx = (i+knotSpan(0)-degree[0]);
					mask.at(c++) = kindx*numberOfControlPoints[0]*numberOfControlPoints[1]+jindx*numberOfControlPoints[0]+iindx+1;
				}
			}
    }          
	}
	else{
    OOFEM_ERROR2 ("giveKnotSpanBasisFunctMask not implemented for nsd = %d", nsd);
  }
  return 1;
}


// for pure Bspline the number of nonzero basis functions is the same for each knot span
int BSplineInterpolation::giveNumberOfKnotSpanBasisFunctions (const IntArray& knotSpan) { 
  int i, answer=1;
	// there are always degree+1 nonzero basis functions on each knot span
  for (i=0; i<nsd; i++) answer*=(degree[i]+1);
  return answer;
}


// generally it is redundant to pass p and U as these data are part of BSplineInterpolation 
// and can be retrieved for given spatial dimension;
// however in such a case this function could not be used for calculation on local knot vector of TSpline;
// it is also redundant to pass the span which can be calculated
// but we want to profit from knowing the span appriori

void BSplineInterpolation::basisFuns (FloatArray &N, int span, double u, int p, const double* U) {
  //
  // Based on Algorithm A2.2 (p. 70)
  //
  FloatArray right(p+1);
  FloatArray left(p+1);
  double saved,temp;
  int j,r;
  
  N.resize(p+1);
  N(0) = 1.0;
  for (j=1; j<=p; j++) {
    left(j)=u-U[span+1-j];
    right(j)=U[span+j]-u;
    saved = 0.0;
    for (r=0; r<j; r++) {
      temp=N(r)/(right(r+1)+left(j-r));
      N(r)=saved+right(r+1)*temp;
      saved=left(j-r)*temp;
    }
    N(j)=saved;
  }
}


// generally it is redundant to pass p and U as these data are part of BSplineInterpolation 
// and can be retrieved for given spatial dimension;
// however in such a case this function could not be used for calculation on local knot vector of TSpline;
// it is also redundant to pass the span which can be calculated
// but we want to profit from knowing the span appriori

void BSplineInterpolation::dersBasisFuns(int n, double u, int span, int p, double* const U, FloatMatrix& ders) {
  //
  // Based on Algorithm A2.3 (p. 72)
  //
  FloatArray left(p+1);
  FloatArray right(p+1);
  FloatMatrix ndu(p+1,p+1) ;
  double saved,temp ;
  int j,r ;
  
  ders.resize(n+1,p+1) ;
  
  ndu(0,0) = 1.0 ;
  for(j=1; j<= p ;j++){
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
  
  for(j=0;j<=p;j++)
    ders(0,j) = ndu(j,p) ;
  
  // Compute the derivatives
  FloatMatrix a(2,p+1) ;
  for(r=0;r<=p;r++){
    int s1,s2 ;
    s1 = 0 ; s2 = 1 ; // alternate rows in array a
    a(0,0) = 1.0 ;
    // Compute the kth derivative
    for(int k=1;k<=n;k++){
      double d ;
      int rk,pk,j1,j2 ;
      d = 0.0 ;
      rk = r-k ; pk = p-k ;
      
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
        j2 = p-r ;
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
  r = p ;
  for(int k=1;k<=n;k++){
    for(j=0;j<=p;j++)
      ders(k,j) *= r ;
    r *= p-k ;
  }
}



// generally it is redundant to pass n, p and U as these data are part of BSplineInterpolation 
// and can be retrieved for given spatial dimension;
// however in such a case this function could not be used for span localization in local knot vector of TSpline

// jaky ma vyznam const = ve funkci se objekt nesmi zmenit

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


// NURBS

NURBSInterpolation::~NURBSInterpolation() {
}


void NURBSInterpolation::evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) {
  FEIIGAElementGeometryWrapper* gw = (FEIIGAElementGeometryWrapper*) &cellgeo;
  FloatArray N[nsd];
  IntArray span(nsd);
	double sum=0.0,val;
  int count,c=1,i,l,k,m,ind,indx,uind,vind,tind;

  if (gw->knotSpan) {
    span = *gw->knotSpan;
  } else {
    for (i=0; i< nsd; i++) {
      span(i) = this->findSpan (numberOfControlPoints[i], degree[i], lcoords(i), knotVector[i]);
    }
  } 
  
  for (i=0; i< nsd; i++) {
    this->basisFuns (N[i], span(i), lcoords(i), degree[i], knotVector[i]);
  }

  count =giveNumberOfKnotSpanBasisFunctions (span);
  answer.resize(count);

	if(nsd == 1){
    uind = span(0)-degree[0];
		ind = uind+1;
		for (k=0; k<=degree[0]; k++) {
			answer.at(c++) = val = N[0](k)*cellgeo.giveVertexCoordinates(ind+k)->at(2);  // Nu*w
			sum += val;
		}
	}
	else if(nsd == 2){
    uind = span(0)-degree[0];
		vind = span(1)-degree[1];
		ind = vind*numberOfControlPoints[0]+uind+1;
    for (l=0; l<=degree[1]; l++) {
			for (k=0; k<=degree[0]; k++) {
        answer.at(c++) = val = N[0](k)*N[1](l)*cellgeo.giveVertexCoordinates(ind+k)->at(3);  // Nu*Nv*w
				sum += val;
			}
			ind += numberOfControlPoints[0];
		}
	}
	else if(nsd == 3){
		uind = span(0)-degree[0];
		vind = span(1)-degree[1];
		tind = span(2)-degree[2];
		ind = tind*numberOfControlPoints[0]*numberOfControlPoints[1]+vind*numberOfControlPoints[0]+uind+1;
		for (m=0; m<=degree[2]; m++) {
			indx = ind;
			for (l=0; l<=degree[1]; l++) {
				for (k=0; k<=degree[0]; k++) {
					answer.at(c++) = val = N[0](k)*N[1](l)*N[2](m)*cellgeo.giveVertexCoordinates(ind+k)->at(4);  // Nu*Nv*Nt*w
					sum += val;
				}
				ind += numberOfControlPoints[0];
			}
			ind = indx+numberOfControlPoints[0]*numberOfControlPoints[1];
		}
	}
	else{
    OOFEM_ERROR2 ("evalN not implemented for nsd = %d", nsd);
  }

	while(count)answer.at(count--) /= sum;
}



void NURBSInterpolation::evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) {

  FEIIGAElementGeometryWrapper* gw = (FEIIGAElementGeometryWrapper*) &cellgeo;
  const FloatArray *vertexCoordsPtr;
  FloatMatrix jacobian(nsd,nsd);
  FloatMatrix ders[nsd];
	IntArray span(nsd);
  double Jacob,product,w,weight;
  int count,cnt,i,l,k,m,ind,indx,uind,vind,tind;

  if (gw->knotSpan) {
    span = *gw->knotSpan;
  } else {
    for (i=0; i< nsd; i++) {
      span(i) = this->findSpan (numberOfControlPoints[i], degree[i], lcoords(i), knotVector[i]);
    }
  } 

  for (i=0; i< nsd; i++) {
    this->dersBasisFuns(1, lcoords(i), span(i), degree[i], knotVector[i], ders[i]);
  }

  count =giveNumberOfKnotSpanBasisFunctions (span);
  answer.resize(count, nsd);

#if 0                       // code according NURBS book (too general allowing higher derivatives)
	if(nsd == 2){
		FloatArray tmp1(nsd+1),tmp2(nsd+1); // allow for weight

		FloatMatrix Aders[nsd]; // derivatives in each coordinate direction on BSpline
#ifndef OPTIMIZED_VERSION_A4dot4
		FloatMatrix Sders[nsd]; // derivatives in each coordinate direction on NURBS
#endif
		FloatMatrix wders;      // derivatives in w direction on BSpline
		/*
		IntArray Bin(2,2);      // binomial coefficients from 0 to d=1 
			                      // Bin(n,k)=(n above k)=n!/k!(n-k)! for n>=k>=0
														// lower triangle corresponds to Pascal triangle
														// according to A4.4 it seems that only coefficients in lower triangle except the first column are used
														*/
		// resizing to (2,2) has nothing common with nsd
		// it is related to the fact that 0th and 1st derivatives are computed in each direction
		for(i=0;i<nsd;i++){
			Aders[i].resize(2,2);
			Aders[i].zero();
#ifndef OPTIMIZED_VERSION_A4dot4
			Sders[i].resize(2,2);
#endif
		}
		wders.resize(2,2);
		wders.zero();

		// calculation of jacobian matrix according to A4.4
		// calculate values and derivatives of nonrational Bspline surface with weights at first (Aders, wders)
    uind = span(0)-degree[0];
		vind = span(1)-degree[1];
		ind = vind*numberOfControlPoints[0]+uind+1;
    for (l=0;l<=degree[1]; l++) {
      tmp1.zero();
      tmp2.zero();
      for (k=0; k<=degree[0]; k++) {
				vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);
				w = vertexCoordsPtr->at(3);

        tmp1(0) += ders[0](0,k)*vertexCoordsPtr->at(1)*w;  // sum(Nu*x*w)
        tmp1(1) += ders[0](0,k)*vertexCoordsPtr->at(2)*w;  // sum(Nu*y*w)
        tmp1(2) += ders[0](0,k)*w;                         // sum(Nu*w)

        tmp2(0) += ders[0](1,k)*vertexCoordsPtr->at(1)*w;  // sum(dNu/du*x*w)
        tmp2(1) += ders[0](1,k)*vertexCoordsPtr->at(2)*w;  // sum(dNu/du*y*w)
        tmp2(2) += ders[0](1,k)*w;                         // sum(dNu/du*w)
      }
			ind += numberOfControlPoints[0];

      Aders[0](0,0) += ders[1](0,l)*tmp1(0);  // xw=sum(Nv*sum(Nu*x*w))
      Aders[1](0,0) += ders[1](0,l)*tmp1(1);  // yw=sum(Nv*sum(Nu*y*w))
      wders(0,0)    += ders[1](0,l)*tmp1(2);  // w=sum(Nv*sum(Nu*w))

      Aders[0](0,1) += ders[1](1,l)*tmp1(0);  // dxw/dv=sum(dNv/dv*sum(Nu*x*w))
      Aders[1](0,1) += ders[1](1,l)*tmp1(1);  // dyw/dv=sum(dNv/dv*sum(Nu*y*w))
      wders(0,1)    += ders[1](1,l)*tmp1(2);  // dw/dv=sum(dNv/dv*sum(Nu*w))
			
      Aders[0](1,0) += ders[1](0,l)*tmp2(0);  // dxw/du=sum(Nv*sum(dNu/du*x*w))
      Aders[1](1,0) += ders[1](0,l)*tmp2(1);  // dyw/du=sum(Nv*sum(dNu/du*y*w))
			wders(1,0)    += ders[1](0,l)*tmp2(2);  // dw/du=sum(Nv*sum(dNu/du*w))
    }

		weight = wders(0,0);
		
#ifndef OPTIMIZED_VERSION_A4dot4
    int j;
    const int d=1;
		// calculate values and derivatives of NURBS surface (A4.4)
		// since all entries in Pascal triangle up to d=1 are 1, binomial coefficients are ignored
		for(k=0;k<=d;k++){
			for(l=0;l<=d-k;l++){
				tmp1(0) = Aders[0](k,l);
				tmp1(1) = Aders[1](k,l);
				for(j=1;j<=l;j++){
					tmp1(0) -= wders(0,j)*Sders[0](k,l-j); // *Bin(l,j)
					tmp1(1) -= wders(0,j)*Sders[1](k,l-j); // *Bin(l,j)
				}
				for(i=1;i<=k;i++){
					tmp1(0) -= wders(i,0)*Sders[0](k-i,l); // *Bin(k,i)
					tmp1(1) -= wders(i,0)*Sders[1](k-i,l); // *Bin(k,i)
					tmp2.zero();
					for(j=1;j<=l;j++){
						tmp2(0) += wders(i,j)*Sders[0](k-i,l-j); // *Bin(l,j)
						tmp2(1) += wders(i,j)*Sders[1](k-i,l-j); // *Bin(l,j)
					}
					tmp1(0) -= tmp2(0); // *Bin(k,i)
					tmp1(1) -= tmp2(1); // *Bin(k,i)
				}
				Sders[0](k,l) = tmp1(0) / weight;
				Sders[1](k,l) = tmp1(1) / weight;
			}
		}

		jacobian(0,0) = Sders[0](1,0);   // dx/du
		jacobian(0,1) = Sders[1](1,0);   // dy/du
		jacobian(1,0) = Sders[0](0,1);   // dx/dv
		jacobian(1,1) = Sders[1](0,1);   // dy/dv
#else
		// optimized version of A4.4 for d=1, binomial coefficients ignored

		/*
		// k=0 l=0 loop
		Sders[0](0,0) = Aders[0](0,0) / weight;
		Sders[1](0,0) = Aders[1](0,0) / weight;
		// k=1 l=0 loop
		Sders[0](1,0) = (Aders[0](1,0)-wders(1,0)*Sders[0](0,0)) / weight;
		Sders[1](1,0) = (Aders[1](1,0)-wders(1,0)*Sders[1](0,0)) / weight;
		// k=0 l=1 loop
		Sders[0](0,1) = (Aders[0](0,1)-wders(0,1)*Sders[0](0,0)) / weight;
		Sders[1](0,1) = (Aders[1](0,1)-wders(0,1)*Sders[1](0,0)) / weight;

		jacobian(0,0) = Sders[0](1,0);   // dx/du
		jacobian(0,1) = Sders[1](1,0);   // dy/du
		jacobian(1,0) = Sders[0](0,1);   // dx/dv
		jacobian(1,1) = Sders[1](0,1);   // dy/dv
		*/

		// k=0 l=0 loop
		tmp1(0) = Aders[0](0,0) / weight;
		tmp1(1) = Aders[1](0,0) / weight;
		// k=1 l=0 loop
		jacobian(0,0) = (Aders[0](1,0)-wders(1,0)*tmp1(0)) / weight;   // dx/du
		jacobian(0,1) = (Aders[1](1,0)-wders(1,0)*tmp1(1)) / weight;   // dy/du
		// k=0 l=1 loop
		jacobian(1,0) = (Aders[0](0,1)-wders(0,1)*tmp1(0)) / weight;   // dx/dv
		jacobian(1,1) = (Aders[1](0,1)-wders(0,1)*tmp1(1)) / weight;   // dy/dv
#endif

		Jacob = jacobian.giveDeterminant();

		//calculation of derivatives of NURBS basis functions with respect to local parameters is not covered by NURBS book
		product = Jacob * weight * weight;
    cnt=0;
		ind = vind*numberOfControlPoints[0]+uind+1;
    for (l=0;l<=degree[1]; l++) {
      for (k=0; k<=degree[0]; k++) {
				w = cellgeo.giveVertexCoordinates(ind+k)->at(3);
        // dNu/du*Nv*w*sum(Nv*Nu*w) - Nu*Nv*w*sum(dNu/du*Nv*w)
        tmp1(0) = ders[0](1,k)*ders[1](0,l)*w*weight - ders[0](0,k)*ders[1](0,l)*w*wders(1,0);   
        // Nu*dNv/dv*w*sum(Nv*Nu*w) - Nu*Nv*w*sum(Nu*dNv/dv*w)
        tmp1(1) = ders[0](0,k)*ders[1](1,l)*w*weight - ders[0](0,k)*ders[1](0,l)*w*wders(0,1);   

        answer(cnt,0) = (+jacobian(1,1)*tmp1(0)-jacobian(0,1)*tmp1(1)) / product;
        answer(cnt,1) = (-jacobian(1,0)*tmp1(0)+jacobian(0,0)*tmp1(1)) / product;
        cnt++;
      }
			ind += numberOfControlPoints[0];
    }
	}
	else{
    OOFEM_ERROR2 ("evaldNdx not implemented for nsd = %d", nsd);
	}
#else
	FloatArray Aders[nsd]; // 0th and 1st derivatives in each coordinate direction on BSpline
	FloatArray wders;      // 0th and 1st derivatives in w direction on BSpline

	for(i=0;i<nsd;i++){
		Aders[i].resize(nsd+1);
		Aders[i].zero();
	}
	wders.resize(nsd+1);
	wders.zero();

	if(nsd == 1){
		// calculate values and derivatives of nonrational Bspline curve with weights at first (Aders, wders)
    uind = span(0)-degree[0];
		ind = uind+1;
		for (k=0; k<=degree[0]; k++) {
			vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);
			w = vertexCoordsPtr->at(2);

			Aders[0](0) += ders[0](0,k)*vertexCoordsPtr->at(1)*w;  // xw=sum(Nu*x*w)
			wders(0)    += ders[0](0,k)*w;                         // w=sum(Nu*w)

			Aders[0](1) += ders[0](1,k)*vertexCoordsPtr->at(1)*w;  // dxw/du=sum(dNu/du*x*w)
			wders(1)    += ders[0](1,k)*w;                         // dw/du=sum(dNu/du*w)
		}

		weight = wders(0);

		// calculation of jacobian matrix according to Eq 4.7
		jacobian(0,0) = (Aders[0](1)-wders(1)*Aders[0](0)/weight) / weight;   // dx/du

		Jacob = jacobian.giveDeterminant();

		//calculation of derivatives of NURBS basis functions with respect to local parameters is not covered by NURBS book
		product = Jacob * weight * weight;
    cnt=0;
		ind = uind+1;
		for (k=0; k<=degree[0]; k++) {
			w = cellgeo.giveVertexCoordinates(ind+k)->at(2);
			// [dNu/du*w*sum(Nu*w) - Nu*w*sum(dNu/du*w)] / [J*sum(Nu*w)^2]
			answer(cnt,0) = ders[0](1,k)*w*weight - ders[0](0,k)*w*wders(1) / product;
			cnt++;
		}
	}
	else if(nsd == 2){
		FloatArray tmp1(nsd+1),tmp2(nsd+1); // allow for weight

		// calculate values and derivatives of nonrational Bspline surface with weights at first (Aders, wders)
    uind = span(0)-degree[0];
		vind = span(1)-degree[1];
		ind = vind*numberOfControlPoints[0]+uind+1;
    for (l=0; l<=degree[1]; l++) {
      tmp1.zero();
      tmp2.zero();
      for (k=0; k<=degree[0]; k++) {
				vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);
				w = vertexCoordsPtr->at(3);

        tmp1(0) += ders[0](0,k)*vertexCoordsPtr->at(1)*w;  // sum(Nu*x*w)
        tmp1(1) += ders[0](0,k)*vertexCoordsPtr->at(2)*w;  // sum(Nu*y*w)
        tmp1(2) += ders[0](0,k)*w;                         // sum(Nu*w)

        tmp2(0) += ders[0](1,k)*vertexCoordsPtr->at(1)*w;  // sum(dNu/du*x*w)
        tmp2(1) += ders[0](1,k)*vertexCoordsPtr->at(2)*w;  // sum(dNu/du*y*w)
        tmp2(2) += ders[0](1,k)*w;                         // sum(dNu/du*w)
      }
			ind += numberOfControlPoints[0];

      Aders[0](0) += ders[1](0,l)*tmp1(0);  // xw=sum(Nv*sum(Nu*x*w)
      Aders[1](0) += ders[1](0,l)*tmp1(1);  // yw=sum(Nv*sum(Nu*y*w)
      wders(0)    += ders[1](0,l)*tmp1(2);  // w=sum(Nv*sum(Nu*w)

      Aders[0](1) += ders[1](0,l)*tmp2(0);  // dxw/du=sum(Nv*sum(dNu/du*x*w)
      Aders[1](1) += ders[1](0,l)*tmp2(1);  // dyw/du=sum(Nv*sum(dNu/du*y*w)
			wders(1)    += ders[1](0,l)*tmp2(2);  // dw/du=sum(Nv*sum(dNu/du*w)

      Aders[0](2) += ders[1](1,l)*tmp1(0);  // dxw/dv=sum(dNv/dv*sum(Nu*x*w)
      Aders[1](2) += ders[1](1,l)*tmp1(1);  // dyw/dv=sum(dNv/dv*sum(Nu*y*w)
      wders(2)    += ders[1](1,l)*tmp1(2);  // dw/dv=sum(dNv/dv*sum(Nu*w)
    }

		weight = wders(0);

		// calculation of jacobian matrix according to Eq 4.19
		tmp1(0) = Aders[0](0) / weight;
		tmp1(1) = Aders[1](0) / weight;
		jacobian(0,0) = (Aders[0](1)-wders(1)*tmp1(0)) / weight;   // dx/du
		jacobian(0,1) = (Aders[1](1)-wders(1)*tmp1(1)) / weight;   // dy/du
		jacobian(1,0) = (Aders[0](2)-wders(2)*tmp1(0)) / weight;   // dx/dv
		jacobian(1,1) = (Aders[1](2)-wders(2)*tmp1(1)) / weight;   // dy/dv

		Jacob = jacobian.giveDeterminant();

		//calculation of derivatives of NURBS basis functions with respect to local parameters is not covered by NURBS book
		product = Jacob * weight * weight;
    cnt=0;
		ind = vind*numberOfControlPoints[0]+uind+1;
    for (l=0; l<=degree[1]; l++) {
      for (k=0; k<=degree[0]; k++) {
				w = cellgeo.giveVertexCoordinates(ind+k)->at(3);
        // dNu/du*Nv*w*sum(Nu*Nv*w) - Nu*Nv*w*sum(dNu/du*Nv*w)
        tmp1(0) = ders[0](1,k)*ders[1](0,l)*w*weight - ders[0](0,k)*ders[1](0,l)*w*wders(1);   
        // Nu*dNv/dv*w*sum(Nu*Nv*w) - Nu*Nv*w*sum(Nu*dNv/dv*w)
        tmp1(1) = ders[0](0,k)*ders[1](1,l)*w*weight - ders[0](0,k)*ders[1](0,l)*w*wders(2);   

        answer(cnt,0) = (+jacobian(1,1)*tmp1(0)-jacobian(0,1)*tmp1(1)) / product;
        answer(cnt,1) = (-jacobian(1,0)*tmp1(0)+jacobian(0,0)*tmp1(1)) / product;
        cnt++;
      }
			ind += numberOfControlPoints[0];
    }
	}
	else if(nsd == 3){
		FloatArray tmp1(nsd+1),tmp2(nsd+1); // allow for weight
		FloatArray temp1(nsd+1),temp2(nsd+1),temp3(nsd+1); // allow for weight

		// calculate values and derivatives of nonrational Bspline solid with weights at first (Aders, wders)
    uind = span(0)-degree[0];
		vind = span(1)-degree[1];
		tind = span(2)-degree[2];
		ind = tind*numberOfControlPoints[0]*numberOfControlPoints[1]+vind*numberOfControlPoints[0]+uind+1;
		for (m=0; m<=degree[2]; m++) {
			temp1.zero();
			temp2.zero();
			temp3.zero();
			indx = ind;
			for (l=0; l<=degree[1]; l++) {
				tmp1.zero();
				tmp2.zero();
				for (k=0; k<=degree[0]; k++) {
					vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);
					w = vertexCoordsPtr->at(4);

					tmp1(0) += ders[0](0,k)*vertexCoordsPtr->at(1)*w;  // sum(Nu*x*w)
					tmp1(1) += ders[0](0,k)*vertexCoordsPtr->at(2)*w;  // sum(Nu*y*w)
					tmp1(2) += ders[0](0,k)*vertexCoordsPtr->at(3)*w;  // sum(Nu*z*w)
					tmp1(3) += ders[0](0,k)*w;                         // sum(Nu*w)
					
					tmp2(0) += ders[0](1,k)*vertexCoordsPtr->at(1)*w;  // sum(dNu/du*x*w)
					tmp2(1) += ders[0](1,k)*vertexCoordsPtr->at(2)*w;  // sum(dNu/du*y*w)
					tmp2(2) += ders[0](1,k)*vertexCoordsPtr->at(3)*w;  // sum(dNu/du*z*w)
					tmp2(3) += ders[0](1,k)*w;                         // sum(dNu/du*w)
				}
				ind += numberOfControlPoints[0];
				
				temp1(0) += ders[1](0,l)*tmp1(0);  // sum(Nv*sum(Nu*x*w))
				temp1(1) += ders[1](0,l)*tmp1(1);  // sum(Nv*sum(Nu*y*w))
				temp1(2) += ders[1](0,l)*tmp1(2);  // sum(Nv*sum(Nu*z*w))
				temp1(3) += ders[1](0,l)*tmp1(3);  // sum(Nv*sum(Nu*w))
				
				temp2(0) += ders[1](0,l)*tmp2(0);  // sum(Nv*sum(dNu/du*x*w))
				temp2(1) += ders[1](0,l)*tmp2(1);  // sum(Nv*sum(dNu/du*y*w))
				temp2(2) += ders[1](0,l)*tmp2(2);  // sum(Nv*sum(dNu/du*z*w))
				temp2(3) += ders[1](0,l)*tmp2(3);  // sum(Nv*sum(dNu/du*w))

				temp3(0) += ders[1](1,l)*tmp1(0);  // sum(dNv/dv*sum(Nu*x*w))
				temp3(1) += ders[1](1,l)*tmp1(1);  // sum(dNv/dv*sum(Nu*y*w))
				temp3(2) += ders[1](1,l)*tmp1(2);  // sum(dNv/dv*sum(Nu*z*w))
				temp3(3) += ders[1](1,l)*tmp1(3);  // sum(dNv/dv*sum(Nu*w))
				
			}
			ind = indx+numberOfControlPoints[0]*numberOfControlPoints[1];
				
			Aders[0](0) += ders[2](0,m)*temp1(0);  // x=sum(Nt*sum(Nv*sum(Nu*x*w)))
			Aders[1](0) += ders[2](0,m)*temp1(1);  // y=sum(Nt*sum(Nv*sum(Nu*y*w)))
			Aders[2](0) += ders[2](0,m)*temp1(2);  // y=sum(Nt*sum(Nv*sum(Nu*y*w)))
			wders(0)    += ders[2](0,m)*temp1(3);  // w=sum(Nt*sum(Nv*sum(Nu*w)))
			
			Aders[0](1) += ders[2](0,m)*temp2(0);  // dx/du=sum(Nt*sum(Nv*sum(dNu/du*x*w)))
			Aders[1](1) += ders[2](0,m)*temp2(1);  // dy/du=sum(Nt*sum(Nv*sum(dNu/du*y*w)))
			Aders[2](1) += ders[2](0,m)*temp2(2);  // dy/du=sum(Nt*sum(Nv*sum(dNu/du*y*w)))
			wders(1)    += ders[2](0,m)*temp2(3);  // dw/du=sum(Nt*sum(Nv*sum(dNu/du*w)))
			
			Aders[0](2) += ders[2](0,m)*temp3(0);  // dx/dv=sum(Nt*sum(dNv/dv*sum(Nu*x*w)))
			Aders[1](2) += ders[2](0,m)*temp3(1);  // dy/dv=sum(Nt*sum(dNv/dv*sum(Nu*y*w)))
			Aders[2](2) += ders[2](0,m)*temp3(2);  // dy/dv=sum(Nt*sum(dNv/dv*sum(Nu*y*w)))
			wders(2)    += ders[2](0,m)*temp3(3);  // dw/dv=sum(Nt*sum(dNv/dv*sum(Nu*w)))
			
			Aders[0](3) += ders[2](1,m)*temp1(0);  // dx/dt=sum(dNt/dt*sum(Nv*sum(Nu*x*w)))
			Aders[1](3) += ders[2](1,m)*temp1(1);  // dy/dt=sum(dNt/dt*sum(Nv*sum(Nu*y*w)))
			Aders[2](3) += ders[2](1,m)*temp1(2);  // dy/dt=sum(dNt/dt*sum(Nv*sum(Nu*y*w)))
			wders(3)    += ders[2](1,m)*temp1(3);  // dw/dt=sum(dNt/dt*sum(Nv*sum(Nu*w)))
		}

		weight = wders(0);

		// calculation of jacobian matrix
		tmp1(0) = Aders[0](0) / weight;
		tmp1(1) = Aders[1](0) / weight;
		tmp1(2) = Aders[2](0) / weight;
		jacobian(0,0) = (Aders[0](1)-wders(1)*tmp1(0)) / weight;   // dx/du
		jacobian(0,1) = (Aders[1](1)-wders(1)*tmp1(1)) / weight;   // dy/du
		jacobian(0,2) = (Aders[2](1)-wders(1)*tmp1(2)) / weight;   // dz/du
		jacobian(1,0) = (Aders[0](2)-wders(2)*tmp1(0)) / weight;   // dx/dv
		jacobian(1,1) = (Aders[1](2)-wders(2)*tmp1(1)) / weight;   // dy/dv
		jacobian(1,2) = (Aders[2](2)-wders(2)*tmp1(2)) / weight;   // dz/dv
		jacobian(2,0) = (Aders[0](3)-wders(3)*tmp1(0)) / weight;   // dx/dt
		jacobian(2,1) = (Aders[1](3)-wders(3)*tmp1(1)) / weight;   // dy/dt
		jacobian(2,2) = (Aders[2](3)-wders(3)*tmp1(2)) / weight;   // dz/dt

		Jacob = jacobian.giveDeterminant();

		//calculation of derivatives of NURBS basis functions with respect to local parameters is not covered by NURBS book
		product = Jacob * weight * weight;
    cnt=0;
		ind = tind*numberOfControlPoints[0]*numberOfControlPoints[1]+vind*numberOfControlPoints[0]+uind+1;
		for (m=0; m<=degree[2]; m++) {
			indx = ind;
			for (l=0;l<=degree[1]; l++) {
				for (k=0; k<=degree[0]; k++) {
					w = cellgeo.giveVertexCoordinates(ind+k)->at(4);
					// dNu/du*Nv*Nt*w*sum(Nu*Nv*Nt*w) - Nu*Nv*Nt*w*sum(dNu/du*Nv*Nt*w)
					tmp1(0) = ders[0](1,k)*ders[1](0,l)*ders[2](0,m)*w*weight - ders[0](0,k)*ders[1](0,l)*ders[2](0,m)*w*wders(1);   
					// Nu*dNv/dv*Nt*w*sum(Nu*Nv*Nt*w) - Nu*Nv*Nt*w*sum(Nu*dNv/dv*Nt*w)
					tmp1(1) = ders[0](0,k)*ders[1](1,l)*ders[2](0,m)*w*weight - ders[0](0,k)*ders[1](0,l)*ders[2](0,m)*w*wders(2);   
					// Nu*Nv*dNt/dt*w*sum(Nu*Nv*Nt*w) - Nu*Nv*Nt*w*sum(Nu*Nv*dNt/dt*w)
					tmp1(2) = ders[0](0,k)*ders[1](0,l)*ders[2](1,m)*w*weight - ders[0](0,k)*ders[1](0,l)*ders[2](0,m)*w*wders(3);   

					answer(cnt,0) = ((jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1))*tmp1(0) +
													 (jacobian(0,2)*jacobian(2,1)-jacobian(0,1)*jacobian(2,2))*tmp1(1) +
													 (jacobian(0,1)*jacobian(1,2)-jacobian(0,2)*jacobian(1,1))*tmp1(2)) / product;  // dN/dx
					answer(cnt,1) = ((jacobian(1,2)*jacobian(2,0)-jacobian(1,0)*jacobian(2,2))*tmp1(0) +
													 (jacobian(0,0)*jacobian(2,2)-jacobian(0,2)*jacobian(2,0))*tmp1(1) +
													 (jacobian(0,2)*jacobian(1,0)-jacobian(0,0)*jacobian(1,2))*tmp1(2)) / product;  // dN/dy
					answer(cnt,2) = ((jacobian(1,0)*jacobian(2,1)-jacobian(1,1)*jacobian(2,0))*tmp1(0) +
													 (jacobian(0,1)*jacobian(2,0)-jacobian(0,0)*jacobian(2,1))*tmp1(1) +
													 (jacobian(0,0)*jacobian(1,1)-jacobian(0,1)*jacobian(1,0))*tmp1(2)) / product;  // dN/dz
					cnt++;
				}
				ind += numberOfControlPoints[0];
			}
			ind = indx+numberOfControlPoints[0]*numberOfControlPoints[1];
    }
	}
	else{
    OOFEM_ERROR2 ("evaldNdx not implemented for nsd = %d", nsd);
  }
#endif

}	


void NURBSInterpolation::local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) {
  /* Based on SurfacePoint A4.3 implementation*/
  FEIIGAElementGeometryWrapper* gw = (FEIIGAElementGeometryWrapper*) &cellgeo;
  const FloatArray *vertexCoordsPtr;
  FloatArray N[nsd];
  IntArray span(nsd);
	double w,weight = 0.0;
  int i,l,k,m,ind,indx,uind,vind,tind;

  if (gw->knotSpan) {
    span = *gw->knotSpan;
  } else {
    for (i=0; i< nsd; i++) {
      span(i) = this->findSpan (numberOfControlPoints[i], degree[i], lcoords(i), knotVector[i]);
    }
  } 

  for (i=0; i< nsd; i++) {
    this->basisFuns (N[i], span(i), lcoords(i), degree[i], knotVector[i]);
  }

  answer.resize(nsd);
	answer.zero();

	if(nsd == 1){
    uind = span(0)-degree[0]; 
		ind = uind+1;
		for (k=0; k<=degree[0]; k++) {
			vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);
			w = vertexCoordsPtr->at(2);
			answer(0) += N[0](k)*vertexCoordsPtr->at(1)*w;  // xw=sum(Nu*x*w)
			weight    += N[0](k)*w;                         // w=sum(Nu*w)
		}
	}
	else if(nsd == 2){
		FloatArray tmp(nsd+1); // allow for weight

    uind = span(0)-degree[0]; 
		vind = span(1)-degree[1];  
		ind = vind*numberOfControlPoints[0]+uind+1;
		for (l=0; l<=degree[1]; l++) {
      tmp.zero();
      for (k=0; k<=degree[0]; k++) {
        vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);
        w = vertexCoordsPtr->at(3);
        tmp(0) += N[0](k)*vertexCoordsPtr->at(1)*w;  // sum(Nu*x*w)
        tmp(1) += N[0](k)*vertexCoordsPtr->at(2)*w;  // sum(Nu*y*w)
        tmp(2) += N[0](k)*w;                         // sum(Nu*w)
      }
			ind += numberOfControlPoints[0];
      answer(0) += N[1](l)*tmp(0);  // xw=sum(Nv*Nu*x*w)
      answer(1) += N[1](l)*tmp(1);  // yw=sum(Nv*Nu*y*w)
      weight    += N[1](l)*tmp(2);  // w=sum(Nv*Nu*w)
    }
	}
	else if(nsd == 3){
		FloatArray tmp(nsd+1),temp(nsd+1); // allow for weight

    uind = span(0)-degree[0]; 
		vind = span(1)-degree[1];  
		tind = span(2)-degree[2];  
		ind = tind*numberOfControlPoints[0]*numberOfControlPoints[1]+vind*numberOfControlPoints[0]+uind+1;
		for (m=0; m<=degree[2]; m++) {
			temp.zero();
			indx = ind;
			for (l=0; l<=degree[1]; l++) {
				tmp.zero();
				for (k=0; k<=degree[0]; k++) {
					vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);
					w = vertexCoordsPtr->at(4);
					tmp(0) += N[0](k)*vertexCoordsPtr->at(1)*w;  // sum(Nu*x*w)
					tmp(1) += N[0](k)*vertexCoordsPtr->at(2)*w;  // sum(Nu*y*w)
					tmp(2) += N[0](k)*vertexCoordsPtr->at(3)*w;  // sum(Nu*z*w)
					tmp(3) += N[0](k)*w;                         // sum(Nu*w)
				}
				ind += numberOfControlPoints[0];
				
				temp(0) += N[1](l)*tmp(0);  // sum(Nv*Nu*x*w)
				temp(1) += N[1](l)*tmp(1);  // sum(Nv*Nu*y*w)
				temp(2) += N[1](l)*tmp(2);  // sum(Nv*Nu*z*w)
				temp(3) += N[1](l)*tmp(3);  // sum(Nv*Nu*w)
			}
			ind = indx+numberOfControlPoints[0]*numberOfControlPoints[1];
			
      answer(0) += N[2](m)*temp(0);  // xw=sum(Nv*Nu*Nt*x*w)
      answer(1) += N[2](m)*temp(1);  // yw=sum(Nv*Nu*Nt*y*w)
      answer(2) += N[2](m)*temp(2);  // zw=sum(Nv*Nu*Nt*z*w)
      weight    += N[2](m)*temp(3);  // w=sum(Nv*Nu*Nt*w)
		}
	}
	else{
    OOFEM_ERROR2 ("local2global not implemented for nsd = %d", nsd);
  }

	answer.times(1.0/weight);
}


double NURBSInterpolation::giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) {
  //
  // Based on Algorithm A4.4 (p. 137) for d=1
  //
  FEIIGAElementGeometryWrapper* gw = (FEIIGAElementGeometryWrapper*) &cellgeo;
  const FloatArray *vertexCoordsPtr;
  FloatMatrix jacobian(nsd,nsd);
  FloatMatrix ders[nsd];
  IntArray span(nsd);
  double Jacob,w,weight;
  int i,l,k,m,ind,indx,uind,vind,tind;

  if (gw->knotSpan) {
    span = *gw->knotSpan;
  } else {
    for (i=0; i< nsd; i++) {
      span(i) = this->findSpan (numberOfControlPoints[i], degree[i], lcoords(i), knotVector[i]);
    }
  } 

  for (i=0; i< nsd; i++) {
    this->dersBasisFuns(1, lcoords(i), span(i), degree[i], knotVector[i], ders[i]);
  }

#if 0                       // code according NURBS book (too general allowing higher derivatives)
  if (nsd == 2) {
		FloatArray tmp1(nsd+1),tmp2(nsd+1); // allow for weight

		FloatMatrix Aders[nsd]; // derivatives in each coordinate direction on BSpline
#ifndef OPTIMIZED_VERSION_A4dot4
		FloatMatrix Sders[nsd]; // derivatives in each coordinate direction on NURBS
#endif
		FloatMatrix wders;      // derivatives in w direction on BSpline
		/*
		IntArray Bin(2,2);      // binomial coefficients from 0 to d=1 
                            // Bin(n,k)=(n above k)=n!/k!(n-k)! for n>=k>=0
														// lower triangle corresponds to Pascal triangle
														// according to A4.4 it seems that only coefficients in lower triangle except the first column are used
														*/
		// resizing to (2,2) has nothing common with nsd
		// it is related to the fact that 0th and 1st derivatives are computed
		for(i=0;i<nsd;i++){
			Aders[i].resize(2,2);
			Aders[i].zero();
#ifndef OPTIMIZED_VERSION_A4dot4
			Sders[i].resize(2,2);
#endif
		}
		wders.resize(2,2);
		wders.zero();

		// calculation of jacobian matrix according to A4.4
		// calculate values and derivatives of nonrational Bspline surface with weights at first (Aders, wders)
    uind = span(0)-degree[0];
		vind = span(1)-degree[1];
		ind = vind*numberOfControlPoints[0]+uind+1;
    for (l=0;l<=degree[1]; l++) {
      tmp1.zero();
      tmp2.zero();
      for (k=0; k<=degree[0]; k++) {
				vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);
				w = vertexCoordsPtr->at(3);

        tmp1(0) += ders[0](0,k)*vertexCoordsPtr->at(1)*w;  // sum(Nu*x*w)
        tmp1(1) += ders[0](0,k)*vertexCoordsPtr->at(2)*w;  // sum(Nu*y*w)
        tmp1(2) += ders[0](0,k)*w;                         // sum(Nu*w)

        tmp2(0) += ders[0](1,k)*vertexCoordsPtr->at(1)*w;  // sum(dNu/du*x*w)
        tmp2(1) += ders[0](1,k)*vertexCoordsPtr->at(2)*w;  // sum(dNu/du*y*w)
        tmp2(2) += ders[0](1,k)*w;                         // sum(dNu/du*w)
      }
			ind += numberOfControlPoints[0];

      Aders[0](0,0) += ders[1](0,l)*tmp1(0);  // xw=sum(Nv*sum(Nu*x*w))
      Aders[1](0,0) += ders[1](0,l)*tmp1(1);  // yw=sum(Nv*sum(Nu*y*w))
      wders(0,0)    += ders[1](0,l)*tmp1(2);  // w=sum(Nv*sum(Nu*w))

      Aders[0](0,1) += ders[1](1,l)*tmp1(0);  // dxw/dv=sum(dNv/dv*sum(Nu*x*w))
      Aders[1](0,1) += ders[1](1,l)*tmp1(1);  // dyw/dv=sum(dNv/dv*sum(Nu*y*w))
      wders(0,1)    += ders[1](1,l)*tmp1(2);  // dw/dv=sum(dNv/dv*sum(Nu*w))
			
      Aders[0](1,0) += ders[1](0,l)*tmp2(0);  // dxw/du=sum(Nv*sum(dNu/du*x*w))
      Aders[1](1,0) += ders[1](0,l)*tmp2(1);  // dyw/du=sum(Nv*sum(dNu/du*y*w))
			wders(1,0)    += ders[1](0,l)*tmp2(2);  // dw/du=sum(Nv*sum(dNu/du*w))
    }

		weight = wders(0,0);
		
#ifndef OPTIMIZED_VERSION_A4dot4
    int j;
    const int d=1;
		// calculate values and derivatives of NURBS surface (A4.4)
		// since all entries in Pascal triangle up to d=1 are 1, binomial coefficients are ignored
		for(k=0;k<=d;k++){
			for(l=0;l<=d-k;l++){
				tmp1(0) = Aders[0](k,l);
				tmp1(1) = Aders[1](k,l);
				for(j=1;j<=l;j++){
					tmp1(0) -= wders(0,j)*Sders[0](k,l-j); // *Bin(l,j)
					tmp1(1) -= wders(0,j)*Sders[1](k,l-j); // *Bin(l,j)
				}
				for(i=1;i<=k;i++){
					tmp1(0) -= wders(i,0)*Sders[0](k-i,l); // *Bin(k,i)
					tmp1(1) -= wders(i,0)*Sders[1](k-i,l); // *Bin(k,i)
					tmp2.zero();
					for(j=1;j<=l;j++){
						tmp2(0) += wders(i,j)*Sders[0](k-i,l-j); // *Bin(l,j)
						tmp2(1) += wders(i,j)*Sders[1](k-i,l-j); // *Bin(l,j)
					}
					tmp1(0) -= tmp2(0); // *Bin(k,i)
					tmp1(1) -= tmp2(1); // *Bin(k,i)
				}
				Sders[0](k,l) = tmp1(0) / weight;
				Sders[1](k,l) = tmp1(1) / weight;
			}
		}

		jacobian(0,0) = Sders[0](1,0);   // dx/du
		jacobian(0,1) = Sders[1](1,0);   // dy/du
		jacobian(1,0) = Sders[0](0,1);   // dx/dv
		jacobian(1,1) = Sders[1](0,1);   // dy/dv
#else
		// optimized version of A4.4 for d=1, binomial coefficients ignored

		/*
		// k=0 l=0 loop
		Sders[0](0,0) = Aders[0](0,0) / weight;
		Sders[1](0,0) = Aders[1](0,0) / weight;
		// k=1 l=0 loop
		Sders[0](1,0) = (Aders[0](1,0)-wders(1,0)*Sders[0](0,0)) / weight;
		Sders[1](1,0) = (Aders[1](1,0)-wders(1,0)*Sders[1](0,0)) / weight;
		// k=0 l=1 loop
		Sders[0](0,1) = (Aders[0](0,1)-wders(0,1)*Sders[0](0,0)) / weight;
		Sders[1](0,1) = (Aders[1](0,1)-wders(0,1)*Sders[1](0,0)) / weight;

		jacobian(0,0) = Sders[0](1,0);   // dx/du
		jacobian(0,1) = Sders[1](1,0);   // dy/du
		jacobian(1,0) = Sders[0](0,1);   // dx/dv
		jacobian(1,1) = Sders[1](0,1);   // dy/dv
		*/

		// k=0 l=0 loop
		tmp1(0) = Aders[0](0,0) / weight;
		tmp1(1) = Aders[1](0,0) / weight;
		// k=1 l=0 loop
		jacobian(0,0) = (Aders[0](1,0)-wders(1,0)*tmp1(0)) / weight;   // dx/du
		jacobian(0,1) = (Aders[1](1,0)-wders(1,0)*tmp1(1)) / weight;   // dy/du
		// k=0 l=1 loop
		jacobian(1,0) = (Aders[0](0,1)-wders(0,1)*tmp1(0)) / weight;   // dx/dv
		jacobian(1,1) = (Aders[1](0,1)-wders(0,1)*tmp1(1)) / weight;   // dy/dv
#endif

  }
	else{
    OOFEM_ERROR2 ("giveTransformationJacobian not implemented for nsd = %d", nsd);
	}
#else
	FloatArray Aders[nsd]; // 0th and 1st derivatives in each coordinate direction on BSpline
	FloatArray wders;      // 0th and 1st derivatives in w direction on BSpline

	for(i=0;i<nsd;i++){
		Aders[i].resize(nsd+1);
		Aders[i].zero();
	}
	wders.resize(nsd+1);
	wders.zero();

	if(nsd == 1){
		// calculate values and derivatives of nonrational Bspline curve with weights at first (Aders, wders)
    uind = span(0)-degree[0];
		ind = uind+1;
		for (k=0; k<=degree[0]; k++) {
			vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);
			w = vertexCoordsPtr->at(2);

			Aders[0](0) += ders[0](0,k)*vertexCoordsPtr->at(1)*w;  // xw=sum(Nu*x*w)
			wders(0)    += ders[0](0,k)*w;                         // w=sum(Nu*w)

			Aders[0](1) += ders[0](1,k)*vertexCoordsPtr->at(1)*w;  // dxw/du=sum(dNu/du*x*w)
			wders(1)    += ders[0](1,k)*w;                         // dw/du=sum(dNu/du*w)
		}

		weight = wders(0);

		// calculation of jacobian matrix according to Eq 4.7
		jacobian(0,0) = (Aders[0](1)-wders(1)*Aders[0](0)/weight) / weight;   // dx/du
	}
	else if(nsd == 2){
		FloatArray tmp1(nsd+1),tmp2(nsd+1); // allow for weight

		// calculate values and derivatives of nonrational Bspline surface with weights at first (Aders, wders)
    uind = span(0)-degree[0];
		vind = span(1)-degree[1];
		ind = vind*numberOfControlPoints[0]+uind+1;
    for (l=0; l<=degree[1]; l++) {
      tmp1.zero();
      tmp2.zero();
      for (k=0; k<=degree[0]; k++) {
				vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);
				w = vertexCoordsPtr->at(3);

        tmp1(0) += ders[0](0,k)*vertexCoordsPtr->at(1)*w;  // sum(Nu*x*w)
        tmp1(1) += ders[0](0,k)*vertexCoordsPtr->at(2)*w;  // sum(Nu*y*w)
        tmp1(2) += ders[0](0,k)*w;                         // sum(Nu*w)

        tmp2(0) += ders[0](1,k)*vertexCoordsPtr->at(1)*w;  // sum(dNu/du*x*w)
        tmp2(1) += ders[0](1,k)*vertexCoordsPtr->at(2)*w;  // sum(dNu/du*y*w)
        tmp2(2) += ders[0](1,k)*w;                         // sum(dNu/du*w)
      }
			ind += numberOfControlPoints[0];

      Aders[0](0) += ders[1](0,l)*tmp1(0);  // xw=sum(Nv*sum(Nu*x*w)
      Aders[1](0) += ders[1](0,l)*tmp1(1);  // yw=sum(Nv*sum(Nu*y*w)
      wders(0)    += ders[1](0,l)*tmp1(2);  // w=sum(Nv*sum(Nu*w)

      Aders[0](1) += ders[1](0,l)*tmp2(0);  // dxw/du=sum(Nv*sum(dNu/du*x*w)
      Aders[1](1) += ders[1](0,l)*tmp2(1);  // dyw/du=sum(Nv*sum(dNu/du*y*w)
			wders(1)    += ders[1](0,l)*tmp2(2);  // dw/du=sum(Nv*sum(dNu/du*w)

      Aders[0](2) += ders[1](1,l)*tmp1(0);  // dxw/dv=sum(dNv/dv*sum(Nu*x*w)
      Aders[1](2) += ders[1](1,l)*tmp1(1);  // dyw/dv=sum(dNv/dv*sum(Nu*y*w)
      wders(2)    += ders[1](1,l)*tmp1(2);  // dw/dv=sum(dNv/dv*sum(Nu*w)
    }

		weight = wders(0);

		// calculation of jacobian matrix according to Eq 4.19
		tmp1(0) = Aders[0](0) / weight;
		tmp1(1) = Aders[1](0) / weight;
		jacobian(0,0) = (Aders[0](1)-wders(1)*tmp1(0)) / weight;   // dx/du
		jacobian(0,1) = (Aders[1](1)-wders(1)*tmp1(1)) / weight;   // dy/du
		jacobian(1,0) = (Aders[0](2)-wders(2)*tmp1(0)) / weight;   // dx/dv
		jacobian(1,1) = (Aders[1](2)-wders(2)*tmp1(1)) / weight;   // dy/dv
	}
	else if(nsd == 3){
		FloatArray tmp1(nsd+1),tmp2(nsd+1); // allow for weight
		FloatArray temp1(nsd+1),temp2(nsd+1),temp3(nsd+1); // allow for weight

		// calculate values and derivatives of nonrational Bspline solid with weights at first (Aders, wders)
    uind = span(0)-degree[0];
		vind = span(1)-degree[1];
		tind = span(2)-degree[2];
		ind = tind*numberOfControlPoints[0]*numberOfControlPoints[1]+vind*numberOfControlPoints[0]+uind+1;
		for (m=0; m<=degree[2]; m++) {
			temp1.zero();
			temp2.zero();
			temp3.zero();
			indx = ind;
			for (l=0; l<=degree[1]; l++) {
				tmp1.zero();
				tmp2.zero();
				for (k=0; k<=degree[0]; k++) {
					vertexCoordsPtr = cellgeo.giveVertexCoordinates(ind+k);
					w = vertexCoordsPtr->at(4);

					tmp1(0) += ders[0](0,k)*vertexCoordsPtr->at(1)*w;  // sum(Nu*x*w)
					tmp1(1) += ders[0](0,k)*vertexCoordsPtr->at(2)*w;  // sum(Nu*y*w)
					tmp1(2) += ders[0](0,k)*vertexCoordsPtr->at(3)*w;  // sum(Nu*z*w)
					tmp1(3) += ders[0](0,k)*w;                         // sum(Nu*w)
					
					tmp2(0) += ders[0](1,k)*vertexCoordsPtr->at(1)*w;  // sum(dNu/du*x*w)
					tmp2(1) += ders[0](1,k)*vertexCoordsPtr->at(2)*w;  // sum(dNu/du*y*w)
					tmp2(2) += ders[0](1,k)*vertexCoordsPtr->at(3)*w;  // sum(dNu/du*z*w)
					tmp2(3) += ders[0](1,k)*w;                         // sum(dNu/du*w)
				}
				ind += numberOfControlPoints[0];
				
				temp1(0) += ders[1](0,l)*tmp1(0);  // sum(Nv*sum(Nu*x*w))
				temp1(1) += ders[1](0,l)*tmp1(1);  // sum(Nv*sum(Nu*y*w))
				temp1(2) += ders[1](0,l)*tmp1(2);  // sum(Nv*sum(Nu*z*w))
				temp1(3) += ders[1](0,l)*tmp1(3);  // sum(Nv*sum(Nu*w))
				
				temp2(0) += ders[1](0,l)*tmp2(0);  // sum(Nv*sum(dNu/du*x*w))
				temp2(1) += ders[1](0,l)*tmp2(1);  // sum(Nv*sum(dNu/du*y*w))
				temp2(2) += ders[1](0,l)*tmp2(2);  // sum(Nv*sum(dNu/du*z*w))
				temp2(3) += ders[1](0,l)*tmp2(3);  // sum(Nv*sum(dNu/du*w))

				temp3(0) += ders[1](1,l)*tmp1(0);  // sum(dNv/dv*sum(Nu*x*w))
				temp3(1) += ders[1](1,l)*tmp1(1);  // sum(dNv/dv*sum(Nu*y*w))
				temp3(2) += ders[1](1,l)*tmp1(2);  // sum(dNv/dv*sum(Nu*z*w))
				temp3(3) += ders[1](1,l)*tmp1(3);  // sum(dNv/dv*sum(Nu*w))
				
			}
			ind = indx+numberOfControlPoints[0]*numberOfControlPoints[1];
				
			Aders[0](0) += ders[2](0,m)*temp1(0);  // x=sum(Nt*sum(Nv*sum(Nu*x*w)))
			Aders[1](0) += ders[2](0,m)*temp1(1);  // y=sum(Nt*sum(Nv*sum(Nu*y*w)))
			Aders[2](0) += ders[2](0,m)*temp1(2);  // y=sum(Nt*sum(Nv*sum(Nu*y*w)))
			wders(0)    += ders[2](0,m)*temp1(3);  // w=sum(Nt*sum(Nv*sum(Nu*w)))
			
			Aders[0](1) += ders[2](0,m)*temp2(0);  // dx/du=sum(Nt*sum(Nv*sum(dNu/du*x*w)))
			Aders[1](1) += ders[2](0,m)*temp2(1);  // dy/du=sum(Nt*sum(Nv*sum(dNu/du*y*w)))
			Aders[2](1) += ders[2](0,m)*temp2(2);  // dy/du=sum(Nt*sum(Nv*sum(dNu/du*y*w)))
			wders(1)    += ders[2](0,m)*temp2(3);  // dw/du=sum(Nt*sum(Nv*sum(dNu/du*w)))
			
			Aders[0](2) += ders[2](0,m)*temp3(0);  // dx/dv=sum(Nt*sum(dNv/dv*sum(Nu*x*w)))
			Aders[1](2) += ders[2](0,m)*temp3(1);  // dy/dv=sum(Nt*sum(dNv/dv*sum(Nu*y*w)))
			Aders[2](2) += ders[2](0,m)*temp3(2);  // dy/dv=sum(Nt*sum(dNv/dv*sum(Nu*y*w)))
			wders(2)    += ders[2](0,m)*temp3(3);  // dw/dv=sum(Nt*sum(dNv/dv*sum(Nu*w)))
			
			Aders[0](3) += ders[2](1,m)*temp1(0);  // dx/dt=sum(dNt/dt*sum(Nv*sum(Nu*x*w)))
			Aders[1](3) += ders[2](1,m)*temp1(1);  // dy/dt=sum(dNt/dt*sum(Nv*sum(Nu*y*w)))
			Aders[2](3) += ders[2](1,m)*temp1(2);  // dy/dt=sum(dNt/dt*sum(Nv*sum(Nu*y*w)))
			wders(3)    += ders[2](1,m)*temp1(3);  // dw/dt=sum(dNt/dt*sum(Nv*sum(Nu*w)))
		}

		weight = wders(0);

		// calculation of jacobian matrix
		tmp1(0) = Aders[0](0) / weight;
		tmp1(1) = Aders[1](0) / weight;
		tmp1(2) = Aders[2](0) / weight;
		jacobian(0,0) = (Aders[0](1)-wders(1)*tmp1(0)) / weight;   // dx/du
		jacobian(0,1) = (Aders[1](1)-wders(1)*tmp1(1)) / weight;   // dy/du
		jacobian(0,2) = (Aders[2](1)-wders(1)*tmp1(2)) / weight;   // dz/du
		jacobian(1,0) = (Aders[0](2)-wders(2)*tmp1(0)) / weight;   // dx/dv
		jacobian(1,1) = (Aders[1](2)-wders(2)*tmp1(1)) / weight;   // dy/dv
		jacobian(1,2) = (Aders[2](2)-wders(2)*tmp1(2)) / weight;   // dz/dv
		jacobian(2,0) = (Aders[0](3)-wders(3)*tmp1(0)) / weight;   // dx/dt
		jacobian(2,1) = (Aders[1](3)-wders(3)*tmp1(1)) / weight;   // dy/dt
		jacobian(2,2) = (Aders[2](3)-wders(3)*tmp1(2)) / weight;   // dz/dt
	}
	else{
    OOFEM_ERROR2 ("giveTransformationJacobianMatrix not implemented for nsd = %d", nsd);
  }
#endif

	Jacob = jacobian.giveDeterminant();
    
	if(fabs(Jacob) < 1.0e-10){
		OOFEM_ERROR ("giveTransformationJacobianMatrix - zero Jacobian");
	}
    
  return Jacob;
}



// TSpline

TSplineInterpolation::~TSplineInterpolation() {
	int i, j;
	
	for(i=0;i<=numberOfControlPoints[0]; i++){
		for(j=0;j<nsd;j++){
			delete [] localIndexKnotVector[i][j];
		}
		delete [] localIndexKnotVector[i];
	}
	delete [] localIndexKnotVector;

	delete [] openLocalKnotVector;
}



IRResultType 
TSplineInterpolation::initializeFrom(InputRecord *ir) {
  const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
  IRResultType result;                   // Required by IR_GIVE_FIELD macro

	BSplineInterpolation::initializeFrom (ir);

  IntArray localIndexKnotVector_tmp;
  int *indexKnotVec, indexKnotVal;
  int i, j, n, pos, p;

  const char *IFT_localIndexKnotVectorString[3] = {"localindexknotvectoru", "localindexknotvectorv", "localindexknotvectorw"};
  InputFieldType IFT_localIndexKnotVectorType[3] = {IFT_TSplineInterpolation_localIndexKnotVectorU, 
																										IFT_TSplineInterpolation_localIndexKnotVectorV, 
																										IFT_TSplineInterpolation_localIndexKnotVectorW};
	int max_deg = 0;
	for(i=0;i<nsd;i++){
		if(degree[i]>max_deg)max_deg=degree[i];
	}
	openLocalKnotVector = new double [ 3*max_deg + 2 ];

  localIndexKnotVector = new int ** [ totalNumberOfControlPoints ];
	for(i=0;i<totalNumberOfControlPoints;i++)localIndexKnotVector[i] = new int * [ nsd ];

	for(n=0;n<nsd;n++){
		localIndexKnotVector_tmp.resize(0);
    IR_GIVE_FIELD(ir, localIndexKnotVector_tmp, IFT_localIndexKnotVectorType[n], IFT_localIndexKnotVectorString[n]); // Macro
    if(localIndexKnotVector_tmp.giveSize() != totalNumberOfControlPoints*(degree[n]+2)){
      OOFEM_ERROR2("BSplineInterpolation::initializeFrom - invalid size of knot vector %s", IFT_localIndexKnotVectorString[n]);
    }

		pos = 0;
		for(i=0;i<totalNumberOfControlPoints;i++){
			indexKnotVec = localIndexKnotVector[i][n] = new int [ degree[n]+2 ];
			
			p = 0;
			for(j=0;j<degree[n]+2;j++)indexKnotVec[p++]=localIndexKnotVector_tmp(pos++);

			// check for monotonicity of local index knot vector with multiplicity
			indexKnotVal = indexKnotVec[0];
			for(j=1;j<degree[n]+2;j++){
				if(indexKnotVal > indexKnotVec[j])
					OOFEM_ERROR3("TSplineInterpolation::initializeFrom - local index knot vector %s of control point %d is not monotonic", 
											 IFT_localIndexKnotVectorString[n], i+1);
				/* this is only for the case when TSpline = NURBS
				if(indexKnotVal+1 < indexKnotVec[j])
					OOFEM_ERROR3("TSplineInterpolation::initializeFrom - local index knot vector %s of control point %d is not continuous", 
											 IFT_localIndexKnotVectorString[n], i+1);
				*/
				indexKnotVal=indexKnotVec[j];
			}
			// check for nondegeneracy of local index knot vector
			if(indexKnotVal == indexKnotVec[0])
				OOFEM_ERROR3("TSplineInterpolation::initializeFrom - local index knot vector %s of control point %d is degenerated", 
										 IFT_localIndexKnotVectorString[n], i+1);
			// check for range of local index knot vector
			if(indexKnotVec[0] <= 0 || indexKnotVal > knotValues[n].giveSize())
				OOFEM_ERROR3("TSplineInterpolation::initializeFrom - local index knot vector %s of control point %d out of range", 
										 IFT_localIndexKnotVectorString[n], i+1);
		}
	}

  return IRRT_OK;
}



void TSplineInterpolation::evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) {
  FEIIGAElementGeometryWrapper* gw = (FEIIGAElementGeometryWrapper*) &cellgeo;
  FloatArray N(nsd);
  IntArray span(nsd);
	IntArray mask;
	double sum=0.0,val;
  int count,i,k,uind,vind;
	
	if(nsd != 2){
		OOFEM_ERROR2 ("evalN not implemented for nsd = %d", nsd);
	}
	
  if (gw->knotSpan) {
    span = *gw->knotSpan;
  } else {
    for (i=0; i< nsd; i++) {
      span(i) = this->findSpan (numberOfControlPoints[i], degree[i], lcoords(i), knotVector[i]);
    }
  } 
	
	// identify which basis functions are nonzero
	giveKnotSpanBasisFuncMask (span, mask) ;
	count = mask.giveSize();
  answer.resize(count);
	
	if(nsd == 2){
		for(k=0;k<count;k++){
			for(i=0;i<nsd;i++){
				N(i)=this->basisFunction(lcoords(i), degree[i], *giveKnotValues(i+1), localIndexKnotVector[mask(k)-1][i]);
			}
			
			answer(k) = val = N(0)*N(1)*cellgeo.giveVertexCoordinates(mask(k))->at(3);  // Nu*Nv*w
			sum += val;
		}
	}
	
	while(count)answer.at(count--) /= sum;
}



void TSplineInterpolation::evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) {
	
  FEIIGAElementGeometryWrapper* gw = (FEIIGAElementGeometryWrapper*) &cellgeo;
  const FloatArray *vertexCoordsPtr;
  FloatMatrix jacobian(nsd,nsd);
  FloatArray tmp_ders[nsd];
	FloatMatrix ders[nsd];
  FloatArray temp(nsd);
	IntArray span(nsd);
	IntArray mask;
  double Jacob,product,w,xw,yw,weight;
  int count,i,k;
 	/*
	IntArray Bin(2,2);      // binomial coefficients from 0 to d=1 
                          // Bin(n,k)=(n above k)=n!/k!(n-k)! for n>=k>=0
                          // lower triangle corresponds to Pascal triangle
													// according to A4.4 it seems that only coefficients in lower triangle except the first column are used
													*/
	if(nsd != 2){
		OOFEM_ERROR2 ("evaldNdx not implemented for nsd = %d", nsd);
	}

  if (gw->knotSpan) {
    span = *gw->knotSpan;
  } else {
    for (i=0; i< nsd; i++) {
      span(i) = this->findSpan (numberOfControlPoints[i], degree[i], lcoords(i), knotVector[i]);
    }
  } 

	// identify which basis functions are nonzero
	giveKnotSpanBasisFuncMask (span, mask) ;
	count = mask.giveSize();
  answer.resize(count, nsd);

	for(i=0;i<nsd;i++){
		ders[i].resize(2,count);
	}

	if(nsd == 2){
		FloatMatrix Aders[nsd]; // derivatives in each coordinate direction on BSpline
		//FloatMatrix Sders[nsd]; // derivatives in each coordinate direction on TSpline
		FloatMatrix wders;      // derivatives in w direction on BSpline
	
		// resizing to (2,2) has nothing common with nsd
		// it is related to the fact that 0th and 1st derivatives are computed
		for(i=0;i<nsd;i++){
			Aders[i].resize(2,2);
			Aders[i].zero();
			//Sders[i].resize(2,2);
		}
		wders.resize(2,2);
		wders.zero();

		for(k=0;k<count;k++){
			for(i=0;i<nsd;i++){
				// it would be simpler if I could pass k-th column of ders[i] directly to dersBasisFunction HUHU array
				this->dersBasisFunction(1, lcoords(i), degree[i], *giveKnotValues(i+1), localIndexKnotVector[mask(k)-1][i], tmp_ders[i]);
				ders[i](0,k)=tmp_ders[i](0);
				ders[i](1,k)=tmp_ders[i](1);
			}

			// calculation of jacobian matrix in similar fashion as A4.4
			// calculate values and derivatives of nonrational Bspline surface with weights at first (Aders, wders)
			vertexCoordsPtr = cellgeo.giveVertexCoordinates(mask(k));
			w = vertexCoordsPtr->at(3);
			xw = vertexCoordsPtr->at(1)*w;
			yw = vertexCoordsPtr->at(2)*w;

			product = tmp_ders[0](0)*tmp_ders[1](0);   // Nu*Nv
			Aders[0](0,0) += product*xw;       // x=sum Nu*Nv*x*w
			Aders[1](0,0) += product*yw;       // x=sum Nu*Nv*y*w
			wders(0,0)    += product*w;        // w=sum Nu*Nv*w
			
			product = tmp_ders[0](1)*tmp_ders[1](0);   // dNu/du*Nv
      Aders[0](1,0) += product*xw;       // dx/du=sum dNu/du*Nv*x*w
      Aders[1](1,0) += product*yw;       // dy/du=sum dNu/du*Nv*y*w
			wders(1,0)    += product*w;        // dw/du=sum dNu/du*Nv*w
			
			product = tmp_ders[0](0)*tmp_ders[1](1);   // Nu*dNv/dv
      Aders[0](0,1) += product*xw;       // dx/dv=sum Nu*dNv/dv*x*w
      Aders[1](0,1) += product*yw;       // dy/dv=sum Nu*dNv/dv*y*w
      wders(0,1)    += product*w;        // dw/dv=sum Nu*dNv/dv*w
		} 

		weight = wders(0,0);

		// optimized version of A4.4 for d=1, binomial coefficients ignored
		/*
			Sders[0](0,0) = Aders[0](0,0) / weight;
			Sders[1](0,0) = Aders[1](0,0) / weight;
			Sders[0](0,1) = (Aders[0](0,1)-wders(0,1)*Sders[0](0,0)) / weight;
			Sders[1](0,1) = (Aders[1](0,1)-wders(0,1)*Sders[1](0,0)) / weight;
			Sders[0](1,0) = (Aders[0](1,0)-wders(1,0)*Sders[0](0,0)) / weight;
			Sders[1](1,0) = (Aders[1](1,0)-wders(1,0)*Sders[1](0,0)) / weight;
			
			jacobian(0,0) = Sders[0](1,0);   // dx/du
			jacobian(0,1) = Sders[1](1,0);   // dy/du
			jacobian(1,0) = Sders[0](0,1);   // dx/dv
			jacobian(1,1) = Sders[1](0,1);   // dy/dv
		*/
		
		temp(0) = Aders[0](0,0) / weight;
		temp(1) = Aders[1](0,0) / weight;
		jacobian(1,0) = (Aders[0](0,1)-wders(0,1)*temp(0)) / weight;   // dx/dv
		jacobian(1,1) = (Aders[1](0,1)-wders(0,1)*temp(1)) / weight;   // dy/dv
		jacobian(0,0) = (Aders[0](1,0)-wders(1,0)*temp(0)) / weight;   // dx/du
		jacobian(0,1) = (Aders[1](1,0)-wders(1,0)*temp(1)) / weight;   // dy/du
		
		Jacob = jacobian.giveDeterminant();

		//calculation of derivatives of TSpline basis functions with respect to local parameters
		product = Jacob * weight * weight;

		for(k=0;k<count;k++){
			w = cellgeo.giveVertexCoordinates(mask(k))->at(3);
			// dNu/du*Nv*w*sum(Nv*Nu*w) - Nu*Nv*w*sum(dNu/du*Nv*w)
			temp(0) = ders[0](1,k)*ders[1](0,k)*w*weight - ders[0](0,k)*ders[1](0,k)*w*wders(1,0);   
			// Nu*dNv/dv*w*sum(Nv*Nu*w) - Nu*Nv*w*sum(Nu*dNv/dv*w)
			temp(1) = ders[0](0,k)*ders[1](1,k)*w*weight - ders[0](0,k)*ders[1](0,k)*w*wders(0,1);   
				
			answer(k,0) = (jacobian(1,1)*temp(0)-jacobian(0,1)*temp(1)) / product;
			answer(k,1) = (-jacobian(1,0)*temp(0)+jacobian(0,0)*temp(1)) / product;
		}
	}
}	


void TSplineInterpolation::local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) {
  /* Based on SurfacePoint A4.3 implementation*/
  FEIIGAElementGeometryWrapper* gw = (FEIIGAElementGeometryWrapper*) &cellgeo;
  const FloatArray *vertexCoordsPtr;
  FloatArray N(nsd);
  IntArray span(nsd);
	IntArray mask;
	double w,xw,yw,product,weight = 0.0;
  int i,k,count;

	if(nsd != 2){
    OOFEM_ERROR2 ("local2global not implemented for nsd = %d", nsd);
  }

  if (gw->knotSpan) {
    span = *gw->knotSpan;
  } else {
    for (i=0; i< nsd; i++) {
      span(i) = this->findSpan (numberOfControlPoints[i], degree[i], lcoords(i), knotVector[i]);
    }
  } 

	// identify which basis functions are nonzero
	giveKnotSpanBasisFuncMask (span, mask) ;
	count = mask.giveSize();

  answer.resize(nsd);
	answer.zero();
  
	if(nsd == 2){
		for(k=0;k<count;k++){
			for(i=0;i<nsd;i++){
				N(i)=this->basisFunction(lcoords(i), degree[i], *giveKnotValues(i+1), localIndexKnotVector[mask(k)-1][i]);
			}

			vertexCoordsPtr = cellgeo.giveVertexCoordinates(mask(k));
			w = vertexCoordsPtr->at(3);
			xw = vertexCoordsPtr->at(1)*w;
			yw = vertexCoordsPtr->at(2)*w;

			product = N(0)*N(1);      // Nu*Nv
      answer(0) += product*xw;  // x=sum Nu*Nv*x*w
      answer(1) += product*yw;  // y=sum Nu*Nv*y*w
      weight    += product*w;   // w=sum Nu*Nv*w
		}
	}

	answer.times(1.0/weight);
}


double TSplineInterpolation::giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) {
  //
  // Based on Algorithm A4.4 (p. 137) for d=1
  //
  FEIIGAElementGeometryWrapper* gw = (FEIIGAElementGeometryWrapper*) &cellgeo;
  const FloatArray *vertexCoordsPtr;
  FloatMatrix jacobian(nsd,nsd);
  FloatArray ders[nsd];
  FloatArray temp(nsd);
  IntArray span(nsd);
	IntArray mask;
  double Jacob,w,xw,yw,product,weight;
  int i,k,count;
	/*
	IntArray Bin(2,2);      // binomial coefficients from 0 to d=1 
                          // Bin(n,k)=(n above k)=n!/k!(n-k)! for n>=k>=0
                          // lower triangle corresponds to Pascal triangle
													// according to A4.4 it seems that only coefficients in lower triangle except the first column are used
													*/
	if(nsd != 2){
    OOFEM_ERROR2 ("giveTransformationJacobianMatrix not implemented for nsd = %d", nsd);
  }

  if (gw->knotSpan) {
    span = *gw->knotSpan;
  } else {
    for (i=0; i< nsd; i++) {
      span(i) = this->findSpan (numberOfControlPoints[i], degree[i], lcoords(i), knotVector[i]);
    }
  } 

	// identify which basis functions are nonzero
	giveKnotSpanBasisFuncMask (span, mask) ;
	count = mask.giveSize();

	if(nsd == 2){
		FloatMatrix Aders[nsd]; // derivatives in each coordinate direction on BSpline
		//FloatMatrix Sders[nsd]; // derivatives in each coordinate direction on TSpline
		FloatMatrix wders;      // derivatives in w direction on BSpline
	
		// resizing to (2,2) has nothing common with nsd
		// it is related to the fact that 0th and 1st derivatives are computed
		for(i=0;i<nsd;i++){
			Aders[i].resize(2,2);
			Aders[i].zero();
			//Sders[i].resize(2,2);
		}
		wders.resize(2,2);
		wders.zero();

		for(k=0;k<count;k++){
			for(i=0;i<nsd;i++){
				this->dersBasisFunction(1, lcoords(i), degree[i], *giveKnotValues(i+1), localIndexKnotVector[mask(k)-1][i], ders[i]);
			}

			// calculation of jacobian matrix in similar fashion as A4.4
			// calculate values and derivatives of nonrational Bspline surface with weights at first (Aders, wders)
			vertexCoordsPtr = cellgeo.giveVertexCoordinates(mask(k));
			w = vertexCoordsPtr->at(3);
			xw = vertexCoordsPtr->at(1)*w;
			yw = vertexCoordsPtr->at(2)*w;

			product = ders[0](0)*ders[1](0);   // Nu*Nv
			Aders[0](0,0) += product*xw;       // x=sum Nu*Nv*x*w
			Aders[1](0,0) += product*yw;       // x=sum Nu*Nv*y*w
			wders(0,0)    += product*w;        // w=sum Nu*Nv*w
			
			product = ders[0](1)*ders[1](0);   // dNu/du*Nv
      Aders[0](1,0) += product*xw;       // dx/du=sum dNu/du*Nv*x*w
      Aders[1](1,0) += product*yw;       // dy/du=sum dNu/du*Nv*y*w
			wders(1,0)    += product*w;        // dw/du=sum dNu/du*Nv*w
			
			product = ders[0](0)*ders[1](1);   // Nu*dNv/dv
      Aders[0](0,1) += product*xw;       // dx/dv=sum Nu*dNv/dv*x*w
      Aders[1](0,1) += product*yw;       // dy/dv=sum Nu*dNv/dv*y*w
      wders(0,1)    += product*w;        // dw/dv=sum Nu*dNv/dv*w
		} 

		weight = wders(0,0);

		// optimized version of A4.4 for d=1, binomial coefficients ignored
		/*
			Sders[0](0,0) = Aders[0](0,0) / weight;
			Sders[1](0,0) = Aders[1](0,0) / weight;
			Sders[0](0,1) = (Aders[0](0,1)-wders(0,1)*Sders[0](0,0)) / weight;
			Sders[1](0,1) = (Aders[1](0,1)-wders(0,1)*Sders[1](0,0)) / weight;
			Sders[0](1,0) = (Aders[0](1,0)-wders(1,0)*Sders[0](0,0)) / weight;
			Sders[1](1,0) = (Aders[1](1,0)-wders(1,0)*Sders[1](0,0)) / weight;
			
			jacobian(0,0) = Sders[0](1,0);   // dx/du
			jacobian(0,1) = Sders[1](1,0);   // dy/du
			jacobian(1,0) = Sders[0](0,1);   // dx/dv
			jacobian(1,1) = Sders[1](0,1);   // dy/dv
		*/
		
		temp(0) = Aders[0](0,0) / weight;
		temp(1) = Aders[1](0,0) / weight;
		jacobian(1,0) = (Aders[0](0,1)-wders(0,1)*temp(0)) / weight;   // dx/dv
		jacobian(1,1) = (Aders[1](0,1)-wders(0,1)*temp(1)) / weight;   // dy/dv
		jacobian(0,0) = (Aders[0](1,0)-wders(1,0)*temp(0)) / weight;   // dx/du
		jacobian(0,1) = (Aders[1](1,0)-wders(1,0)*temp(1)) / weight;   // dy/du
	}
		
	Jacob = jacobian.giveDeterminant();

	if(fabs(Jacob) < 1.0e-10){
		OOFEM_ERROR ("giveTransformationJacobianMatrix - zero Jacobian");
	}
    
  return Jacob;
}


// knotSpan corresponds to knot span in terms of BSpline;
// it should not matter which part of IGAIntegrationElement if covering more spans is addressed;
// this implementation relies on the fact that IGAIntegrationElements are those subsets
// of T-mesh cells on which there are fully (this means not only partially) nonzero relevant basis functions

int TSplineInterpolation::giveKnotSpanBasisFuncMask (const IntArray& knotSpan, IntArray& mask) {
  int i,j,nonzero;
	FloatArray knotStart(nsd), knotEnd(nsd);
  
	// resize the mask initially to the size corresponding to BSpline case
	// but there may be more nonzero basis functions 
  if (nsd == 2) {
    mask.preallocate((degree[0]+1)*(degree[1]+1));
  } else {
    OOFEM_ERROR2 ("giveKnotSpanBasisFunctMask not implemented for nsd = %d", nsd);
  }

	// get starting and ending knots
	for(j=0; j<nsd; j++){
		knotStart(j)=knotVector[j][knotSpan(j)];
		knotEnd(j)=knotVector[j][knotSpan(j)+1];
	}

	// for each control point check
	for (i=0; i<totalNumberOfControlPoints; i++){
		// whether local knot vector overlaps the given knot span
		nonzero=1;
		for(j=0; j<nsd; j++){
			if((knotEnd(j) <= knotValues[j].at(localIndexKnotVector[i][j][0])) ||
				 (knotStart(j) >= knotValues[j].at(localIndexKnotVector[i][j][degree[j]+1]))){
				nonzero=0;
				break;
			}
		}
		if(nonzero)mask.followedBy(i+1, 4);
	}

	return 1;
}



// knotSpan corresponds to knot span in terms of BSpline;
// it should not matter which part of IGAIntegrationElement if covering more spans is addressed;
// this implementation relies on the fact that IGAIntegrationElements are those subsets
// of T-mesh cells on which there are fully (this means not only partially) nonzero relevant basis functions

int TSplineInterpolation::giveNumberOfKnotSpanBasisFunctions (const IntArray& knotSpan) { 
  int i, j, answer=0;
	FloatArray knotStart(nsd), knotEnd(nsd);

	// get starting and ending knots
	for(j=0; j<nsd; j++){
		knotStart(j)=knotVector[j][knotSpan(j)];
		knotEnd(j)=knotVector[j][knotSpan(j)+1];
	}

	// for each control point check
	for (i=0; i<totalNumberOfControlPoints; i++){
		answer++;
		// whether local knot vector overlaps the given knot span
		for(j=0; j<nsd; j++){
			if((knotEnd(j) <= knotValues[j].at(localIndexKnotVector[i][j][0])) ||
				 (knotStart(j) >= knotValues[j].at(localIndexKnotVector[i][j][degree[j]+1]))){
				answer--;
				break;
			}
		}
	}
  return answer;
}



// starKnotSpan and endKnotSpan correspond to knot span in terms of BSpline;
// should the number of non-zero basis function be calculated for single knot span
// starKnotSpan and endKnotSpan are equal;
// some of the basis function may not cover the whole knot span interval !!!

int TSplineInterpolation::giveKnotSpanBasisFuncMask (const IntArray& startKnotSpan, const IntArray& endKnotSpan, IntArray& mask) {
  int i,j,nonzero;
	FloatArray knotStart(nsd), knotEnd(nsd);
  
	// resize the mask initially to the size corresponding to BSpline case
	// but there may be more nonzero basis functions 
  if (nsd == 2) {
    mask.preallocate((degree[0]+1)*(degree[1]+1));
  } else {
    OOFEM_ERROR2 ("giveKnotSpanBasisFunctMask not implemented for nsd = %d", nsd);
  }

	// get starting and ending knots
	for(j=0; j<nsd; j++){
		knotStart(j)=knotVector[j][startKnotSpan(j)];
		knotEnd(j)=knotVector[j][endKnotSpan(j)+1];
	}

	// for each control point check
	for (i=0; i<totalNumberOfControlPoints; i++){
		// whether local knot vector overlaps at least partially the knot span interval
		nonzero=1;
		for(j=0; j<nsd; j++){
			if((knotEnd(j) <= knotValues[j].at(localIndexKnotVector[i][j][0])) ||
				 (knotStart(j) >= knotValues[j].at(localIndexKnotVector[i][j][degree[j]+1]))){
				nonzero=0;
				break;
			}
		}
		if(nonzero)mask.followedBy(i+1,4);
	}

	return 1;
}


// starKnotSpan and endKnotSpan correspond to knot apan in terms of BSpline;
// should the number of non-zero basis function be calculated for single knot span
// starKnotSpan and endKnotSpan are equal
// some of the basis function may not cover the whole knot span interval !!!

int TSplineInterpolation::giveNumberOfKnotSpanBasisFunctions (const IntArray& startKnotSpan, const IntArray& endKnotSpan) { 
  int i, j, answer=0;
	FloatArray knotStart(nsd), knotEnd(nsd);

	// get starting and ending knots
	for(j=0; j<nsd; j++){
		knotStart(j)=knotVector[j][startKnotSpan(j)];
		knotEnd(j)=knotVector[j][endKnotSpan(j)+1];
	}

	// for each control point check
	for (i=0; i<totalNumberOfControlPoints; i++){
		answer++;
		// whether local knot vector overlaps at least partially the knot span interval
		for(j=0; j<nsd; j++){
			if((knotEnd(j) <= knotValues[j].at(localIndexKnotVector[i][j][0])) ||
				 (knotStart(j) >= knotValues[j].at(localIndexKnotVector[i][j][degree[j]+1]))){
				answer--;
				break;
			}
		}
	}
  return answer;
}



// call corresponding BSpline methods for open local knot vector

double TSplineInterpolation::basisFunction (double u, int p, const FloatArray& U, const int* I) {
	int span, prepend, append;
	FloatArray N;

	createLocalKnotVector(p, U, I, &prepend, &append);
	span=BSplineInterpolation::findSpan(prepend+append, p, u, openLocalKnotVector);
	BSplineInterpolation::basisFuns(N, span, u, p, openLocalKnotVector);
	
	// extract the middle basis function
	// this corresponds to index p-span, however prepended knotspans must be considered
	return N(p-span+prepend);
}


// call corresponding BSpline methods for open local knot vector

void TSplineInterpolation::dersBasisFunction(int n, double u, int p, const FloatArray& U, const int* I, FloatArray& ders) {
	int i, span, prepend, append;
	FloatMatrix Ders;

	createLocalKnotVector(p, U, I, &prepend, &append);
	span=BSplineInterpolation::findSpan(prepend+append, p, u, openLocalKnotVector);
	BSplineInterpolation::dersBasisFuns(n, u, span, p, openLocalKnotVector, Ders);

	// extract the middle basis function and its derivatives
	// this corresponds to index p-span, however prepended knotspans must be considered
	ders.resize(n+1);
	for(i=0;i<=n;i++)ders(i)=Ders(i,p-span+prepend);
}


void TSplineInterpolation::createLocalKnotVector(int p, const FloatArray& U, const int *I, int *prepend, int *append) {
	int i,j=0,index_first=I[0],index_last=I[p+1],mult_first=1,mult_last=1;
	double first=U.at(index_first), last=U.at(index_last);

	for(i=1;i<p+1;i++){
		if(I[i] != index_first)break;
		mult_first++;
	}
	for(i=p;i>0;i--){
		if(I[i] != index_last)break;
		mult_last++;
	}
	*prepend=p+1-mult_first;
	*append=p+1-mult_last;

	// prepend first knot (once more)
	for(i=0;i<=*prepend;i++)openLocalKnotVector[j++]=first;
	// copy middle of knot vector (without first and last)
	for(i=1;i<=p;i++)openLocalKnotVector[j++]=U.at(I[i]);
	// append last knot (once more)
	for(i=0;i<=*append;i++)openLocalKnotVector[j++]=last;
}


IRResultType IGAElement::initializeFrom(InputRecord *ir) {
  const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
  IRResultType result;                   // Required by IR_GIVE_FIELD macro

  int indx=0, ui,vi,wi, i, nsd, numberOfGaussPoints=1;
  double du,dv,dw;
  const FloatArray *gpcoords;
  FloatArray newgpcoords;
	IntArray knotSpan; 

  Element::initializeFrom (ir); // read nodes , material, cross section
  // set number of dofmanagers
  this->numberOfDofMans = dofManArray.giveSize();
  this->giveInterpolation()->initializeFrom (ir); // read geometry

  IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_IGAElement_NIP, "nip"); // Macro
  
  // generate individual IntegrationElements; one for each nonzero knot span
	nsd = this->giveNsd();
  if (nsd == 1) {
		//HUHU

	} else if (nsd == 2) {
		int numberOfKnotSpansU = this->giveInterpolation()->giveNumberOfKnotSpans(1);
		int numberOfKnotSpansV = this->giveInterpolation()->giveNumberOfKnotSpans(2);
		IntArray* const knotMultiplicityU = this->giveInterpolation()->giveKnotMultiplicity(1); 
		IntArray* const knotMultiplicityV = this->giveInterpolation()->giveKnotMultiplicity(2);
		FloatArray* const knotValuesU = this->giveInterpolation()->giveKnotValues(1);
		FloatArray* const knotValuesV = this->giveInterpolation()->giveKnotValues(2);

    newgpcoords.resize(2);
		knotSpan.resize(2);

    this->numberOfIntegrationRules = numberOfKnotSpansU * numberOfKnotSpansV;
    integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];

    knotSpan.at(2) = -1; 
    for (vi=1; vi<=numberOfKnotSpansV; vi++) {
			dv = knotValuesV->at(vi+1)-knotValuesV->at(vi);
      knotSpan.at(2) += knotMultiplicityV->at(vi);

      knotSpan.at(1) = -1;
      for (ui=1; ui<=numberOfKnotSpansU; ui++) {
				du = knotValuesU->at(ui+1)-knotValuesU->at(ui);
        knotSpan.at(1)+=knotMultiplicityU->at(ui);

        integrationRulesArray [ indx ] = new IGAIntegrationElement (indx, this, knotSpan);
        integrationRulesArray [ indx ] -> setUpIntegrationPoints(_Square, numberOfGaussPoints, _PlaneStress); // HUHU _PlaneStress, rectangle

        // remap local subelement gp coordinates into knot span coordinates and update integration weight 
        for (i=0; i<numberOfGaussPoints; i++) {
          gpcoords = integrationRulesArray [ indx ]->getIntegrationPoint(i)->giveCoordinates();
      
          newgpcoords.at(1) = knotValuesU->at(ui)+du*(gpcoords->at(1)/2.0+0.5);
          newgpcoords.at(2) = knotValuesV->at(vi)+dv*(gpcoords->at(2)/2.0+0.5);
          integrationRulesArray [ indx ]->getIntegrationPoint(i)->setCoordinates(newgpcoords);
          integrationRulesArray [ indx ]->getIntegrationPoint(i)->setWeight(integrationRulesArray [ indx ]->getIntegrationPoint(i)->giveWeight()/4.0*du*dv);
        }
				indx++;
      }
    }
  } else if (nsd == 3) {
		int numberOfKnotSpansU = this->giveInterpolation()->giveNumberOfKnotSpans(1);
		int numberOfKnotSpansV = this->giveInterpolation()->giveNumberOfKnotSpans(2);
		int numberOfKnotSpansW = this->giveInterpolation()->giveNumberOfKnotSpans(3);
		IntArray* const knotMultiplicityU = this->giveInterpolation()->giveKnotMultiplicity(1); 
		IntArray* const knotMultiplicityV = this->giveInterpolation()->giveKnotMultiplicity(2);
		IntArray* const knotMultiplicityW = this->giveInterpolation()->giveKnotMultiplicity(3);
		FloatArray* const knotValuesU = this->giveInterpolation()->giveKnotValues(1);
		FloatArray* const knotValuesV = this->giveInterpolation()->giveKnotValues(2);
		FloatArray* const knotValuesW = this->giveInterpolation()->giveKnotValues(3);

    newgpcoords.resize(3);
		knotSpan.resize(3);

    this->numberOfIntegrationRules = numberOfKnotSpansU * numberOfKnotSpansV * numberOfKnotSpansW;
    integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];

    knotSpan.at(3) = -1; 
    for (wi=1; wi<=numberOfKnotSpansW; wi++) {
			dw = knotValuesW->at(wi+1)-knotValuesW->at(wi);
      knotSpan.at(3) += knotMultiplicityW->at(wi);

			knotSpan.at(2) = -1; 
			for (vi=1; vi<=numberOfKnotSpansV; vi++) {
				dv = knotValuesV->at(vi+1)-knotValuesV->at(vi);
				knotSpan.at(2) += knotMultiplicityV->at(vi);
				
				knotSpan.at(1) = -1;
				for (ui=1; ui<=numberOfKnotSpansU; ui++) {
					du = knotValuesU->at(ui+1)-knotValuesU->at(ui);
					knotSpan.at(1)+=knotMultiplicityU->at(ui);

					integrationRulesArray [ indx ] = new IGAIntegrationElement (indx, this, knotSpan);
					integrationRulesArray [ indx ] -> setUpIntegrationPoints(_Cube, numberOfGaussPoints, _3dMat);

					// remap local subelement gp coordinates into knot span coordinates and update integration weight 
					for (i=0; i<numberOfGaussPoints; i++) {
						gpcoords = integrationRulesArray [ indx ]->getIntegrationPoint(i)->giveCoordinates();
      
						newgpcoords.at(1) = knotValuesU->at(ui)+du*(gpcoords->at(1)/2.0+0.5);
						newgpcoords.at(2) = knotValuesV->at(vi)+dv*(gpcoords->at(2)/2.0+0.5);
						newgpcoords.at(3) = knotValuesW->at(wi)+dw*(gpcoords->at(3)/2.0+0.5);
						integrationRulesArray [ indx ]->getIntegrationPoint(i)->setCoordinates(newgpcoords);
						integrationRulesArray [ indx ]->getIntegrationPoint(i)->setWeight(integrationRulesArray [ indx ]->getIntegrationPoint(i)->giveWeight()/8.0*du*dv*dw);
					}
					indx++;
				}
      }
    }
	} else {
    OOFEM_ERROR2 ("unsupported number of spatial dimensions (nsd = %d)", nsd);
  }
  return IRRT_OK; 
}



// integration elements are setup in the same way as for IGAElement for now HUHU

IRResultType IGATSplineElement::initializeFrom(InputRecord *ir) {
  const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
  IRResultType result;                   // Required by IR_GIVE_FIELD macro
	TSplineInterpolation *interpol = (TSplineInterpolation *)this->giveInterpolation();

  int indx=0, ui,vi, i, nsd,numberOfGaussPoints=1;
  double du,dv;
  const FloatArray *gpcoords;
  FloatArray newgpcoords;
	IntArray knotSpan; 

  Element::initializeFrom (ir); // read nodes , material, cross section
  // set number of dofmanagers
  this->numberOfDofMans = dofManArray.giveSize();
	// set number of control points before initialization HUHU HAHA
	interpol->setNumberOfControlPoints(this->numberOfDofMans);
  this->giveInterpolation()->initializeFrom (ir); // read geometry

  IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_IGAElement_NIP, "nip"); // Macro
  
  // generate individual IntegrationElements; one for each nonzero knot span
	nsd = giveNsd();
  if (nsd == 2) {
		int numberOfKnotSpansU = this->giveInterpolation()->giveNumberOfKnotSpans(1);
		int numberOfKnotSpansV = this->giveInterpolation()->giveNumberOfKnotSpans(2);
		IntArray* const knotMultiplicityU = this->giveInterpolation()->giveKnotMultiplicity(1); 
		IntArray* const knotMultiplicityV = this->giveInterpolation()->giveKnotMultiplicity(2);
		FloatArray* const knotValuesU = this->giveInterpolation()->giveKnotValues(1);
		FloatArray* const knotValuesV = this->giveInterpolation()->giveKnotValues(2);

    newgpcoords.resize(2);
		knotSpan.resize(2);

    this->numberOfIntegrationRules = numberOfKnotSpansU * numberOfKnotSpansV;
    integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];

    knotSpan.at(2) = -1; 
    for (vi=1; vi<=numberOfKnotSpansV; vi++) {
			dv = knotValuesV->at(vi+1)-knotValuesV->at(vi);
      knotSpan.at(2) += knotMultiplicityV->at(vi);

      knotSpan.at(1) = -1;
      for (ui=1; ui<=numberOfKnotSpansU; ui++) {
				du = knotValuesU->at(ui+1)-knotValuesU->at(ui);
        knotSpan.at(1)+=knotMultiplicityU->at(ui);

        integrationRulesArray [ indx ] = new IGAIntegrationElement (indx, this, knotSpan);
        integrationRulesArray [ indx ] -> setUpIntegrationPoints(_Square, numberOfGaussPoints, _PlaneStress); // HUHU _PlaneStress, rectangle

        // remap local subelement gp coordinates into knot span coordinates and update integration weight 
        for (i=0; i<numberOfGaussPoints; i++) {
          gpcoords = integrationRulesArray [ indx ]->getIntegrationPoint(i)->giveCoordinates();
      
          newgpcoords.at(1) = knotValuesU->at(ui)+du*(gpcoords->at(1)/2.0+0.5);
          newgpcoords.at(2) = knotValuesV->at(vi)+dv*(gpcoords->at(2)/2.0+0.5);
          integrationRulesArray [ indx ]->getIntegrationPoint(i)->setCoordinates(newgpcoords);
          integrationRulesArray [ indx ]->getIntegrationPoint(i)->setWeight(integrationRulesArray [ indx ]->getIntegrationPoint(i)->giveWeight()/4.0*du*dv);
        }
				indx++;
			}
    }
  } else {
    OOFEM_ERROR2 ("unsupported number of spatial dimensions (nsd = %d)", nsd);
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
    IGAIntegrationElement *ee = (IGAIntegrationElement *)ie;
    elem->giveInterpolation()->giveKnotSpanBasisFuncMask (*ee->giveKnotSpan(), mask) ;
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
  int i,j, nsd;
  IntArray mask, nodeDofIDMask, nodalArray;
  int dofmandof;

  // get number of dofs in node
  elem->giveDofManDofIDMask (1, ut, nodeDofIDMask);
  dofmandof=nodeDofIDMask.giveSize();

	nsd = elem->giveInterpolation()->giveNsd();

  // first evaluate nonzero basis function mask
  if (elem->giveInterpolation()->hasSubPatchFormulation()) {
		IGAIntegrationElement *ee = (IGAIntegrationElement *)ie;
		elem->giveInterpolation()->giveKnotSpanBasisFuncMask (*ee->giveKnotSpan(), mask) ;
		// loop over nonzero shape functions and assemble localization array
		answer.resize(0);
		for (i=1; i<=mask.giveSize(); i++) {
			nodalArray.resize(nodeDofIDMask.giveSize());
			for(j=1; j<=nsd; j++)nodalArray.at(j) = dofmandof*(mask.at(i)-1)+j;
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
  
  if ((nDofMan = interp->giveNumberOfKnotSpanBasisFunctions(*(gp->giveIntegrationRule()->giveKnotSpan()))) == 0)  // HUHU
    nDofMan = gp->giveElement()->giveNumberOfDofManagers();
  
  answer.resize(2, nDofMan*2);
  answer.zero();
  
  for (i=1; i <= nDofMan; i++) {
    answer.at(1, i*2-1) = N.at(i);
    answer.at(2, i*2-0)   = N.at(i);
  }
}

void PlaneStressStructuralElementEvaluator::computeBMatrixAt (FloatMatrix& answer, GaussPoint* gp) {
  int i, nDofMan;
  //IntArray dofmanSubElementMask;
  FloatMatrix d;
  
  FEInterpolation* interp = gp->giveElement()->giveInterpolation();
  // this uses FEIInterpolation::nodes2coords - quite inefficient in this case (large num of dofmans)
  interp->evaldNdx (d, *gp->giveCoordinates(), 
		    FEIIGAElementGeometryWrapper(gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan()), 0.0);
  
  if ((nDofMan = interp->giveNumberOfKnotSpanBasisFunctions (*(gp->giveIntegrationRule()->giveKnotSpan()))) == 0)  // HUHU
    nDofMan = gp->giveElement()->giveNumberOfDofManagers();
  
  answer.resize(3, nDofMan*2);
  answer.zero();
  
  for (i=1; i <= nDofMan; i++) {
    answer.at(1, i*2-1) = d.at(i, 1);
    answer.at(2, i*2-0)   = d.at(i, 2);
    
    answer.at(3, 2*i-1) = d.at(i, 2);
    answer.at(3, 2*i-0) = d.at(i, 1);
  }
}
  

void StructuralElementEvaluator::giveCharacteristicMatrix(FloatMatrix &answer,
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

void StructuralElementEvaluator::computeBcLoadVectorAt(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
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
            //delete dofMask;
          }
        }
    }
    
    return;
}




void StructuralElementEvaluator::computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN) 
{
  FloatArray u;
  Element* elem=this->giveElement();


  elem->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);
    
  /*
  // substract initial displacements, if defined
  if (initialDisplacements) u.substract(initialDisplacements);
  */
  if ( this->updateRotationMatrix() ) {
    u.rotatedWith(this->rotationMatrix, 'n');
  }
  this->computeStrainVector(answer, gp, stepN, u); 
}

void StructuralElementEvaluator::computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN, FloatArray& u)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step stepN. The nature of these strains depends
// on the element's type.
{
  int i;
    FloatMatrix b;
    FloatArray ur;
    Element* elem=this->giveElement();

    if (!this->isActivated(stepN)) {
      answer.resize(elem->giveCrossSection()->giveIPValueSize(IST_StrainTensor, gp));
      answer.zero();
      return;}

    this->computeBMatrixAt(b, gp);

    // get local code numbers corresponding to ir
    IntArray lc;
    this->giveIntegrationElementLocalCodeNumbers (lc, elem, gp->giveIntegrationRule(), EID_MomentumBalance);
    ur.resize(b.giveNumberOfColumns());
    for (i=1; i<=lc.giveSize(); i++) ur.at(i) = u.at(lc.at(i));

    answer.beProductOf(b, ur);

    return;
}


void StructuralElementEvaluator::computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
// Computes the vector containing the stresses at the Gauss point gp of
// the receiver, at time step stepN. The nature of these stresses depends
// on the element's type.
// this version assumes TOTAL LAGRANGE APPROACH
{
  FloatArray Epsilon;
  Element* elem=this->giveElement();
  StructuralCrossSection *cs = ( StructuralCrossSection * ) elem->giveCrossSection();
  
  this->computeStrainVector(Epsilon, gp, stepN);
  cs->giveRealStresses(answer, ReducedForm, gp, Epsilon, stepN);
  
  return;
} 


void StructuralElementEvaluator::computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN, FloatArray& u)
// Computes the vector containing the stresses at the Gauss point gp of
// the receiver, at time step stepN. The nature of these stresses depends
// on the element's type.
// this version assumes TOTAL LAGRANGE APPROACH
{
  FloatArray Epsilon;
  Element* elem=this->giveElement();
  StructuralCrossSection *cs = ( StructuralCrossSection * ) elem->giveCrossSection();
  
  this->computeStrainVector(Epsilon, gp, stepN, u);
  cs->giveRealStresses(answer, ReducedForm, gp, Epsilon, stepN);
  
  return;
}

void StructuralElementEvaluator::updateInternalState(TimeStep *stepN)
// Updates the receiver at end of step.
{
  FloatArray u;
  Element* elem=this->giveElement();

  elem->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);
    
  /*
  // substract initial displacements, if defined
  if (initialDisplacements) u.substract(initialDisplacements);
  */
  if ( this->updateRotationMatrix() ) {
    u.rotatedWith(this->rotationMatrix, 'n');
  }
  
  int i, j;
  IntegrationRule *iRule;
  FloatArray stress;
  
  // force updating strains & stresses
  for ( i = 0; i < elem->giveNumberOfIntegrationRules(); i++ ) {
    iRule = elem->giveIntegrationRule(i);
    for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
      computeStressVector(stress, iRule->getIntegrationPoint(j), stepN, u);
    }
  }

  /* 
     // Original unoptimized version
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
  */
}


void StructuralElementEvaluator::computeNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
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


IRResultType BsplinePlaneStressElement::initializeFrom(InputRecord *ir) {
	BSplineInterpolation *interpol = (BSplineInterpolation *)this->giveInterpolation();
	IGAElement::initializeFrom(ir);
	//PlaneStressStructuralElementEvaluator::initializeFrom(ir);

	// HUHU checkConsistency()
	if(giveNumberOfDofManagers() != interpol->giveNumberOfControlPoints(1)*interpol->giveNumberOfControlPoints(2)){
    OOFEM_ERROR("BsplinePlaneStressElement::initializeFrom - number of control points mismatch");
	}

	return IRRT_OK;
}



NURBSPlaneStressElement::NURBSPlaneStressElement (int n, Domain *aDomain) : IGAElement (n, aDomain), PlaneStressStructuralElementEvaluator(), interpolation(2) {}


IRResultType NURBSPlaneStressElement::initializeFrom(InputRecord *ir) {
	NURBSInterpolation *interpol = (NURBSInterpolation *)this->giveInterpolation();
	IGAElement::initializeFrom(ir);
	//PlaneStressStructuralElementEvaluator::initializeFrom(ir);

	// HUHU
	if(giveNumberOfDofManagers() != interpol->giveNumberOfControlPoints(1)*interpol->giveNumberOfControlPoints(2)){
    OOFEM_ERROR("NURBSPlaneStressElement::initializeFrom - number of control points mismatch");
	}

	return IRRT_OK;
}


TSplinePlaneStressElement::TSplinePlaneStressElement (int n, Domain *aDomain) : IGATSplineElement (n, aDomain), PlaneStressStructuralElementEvaluator(), interpolation(2) {}



/* 3D Space Elements */
void Space3dStructuralElementEvaluator::computeNMatrixAt (FloatMatrix& answer, GaussPoint* gp) {
  int i, nDofMan;
  FloatArray N;
  FEInterpolation* interp = gp->giveElement()->giveInterpolation();
  
  interp->evalN (N, *gp->giveCoordinates(), FEIIGAElementGeometryWrapper(gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan()), 0.0);
  
  if ((nDofMan = interp->giveNumberOfKnotSpanBasisFunctions(*(gp->giveIntegrationRule()->giveKnotSpan()))) == 0)  // HUHU
    nDofMan = gp->giveElement()->giveNumberOfDofManagers();
  
  answer.resize(3, nDofMan*3);
  answer.zero();
  
  for (i=1; i <= nDofMan; i++) {
    answer.at(1, i*3-2) = N.at(i);
    answer.at(2, i*3-1) = N.at(i);
    answer.at(3, i*3-0) = N.at(i);
  }
}

void Space3dStructuralElementEvaluator::computeBMatrixAt (FloatMatrix& answer, GaussPoint* gp) {
  int i, nDofMan;
  //IntArray dofmanSubElementMask;
  FloatMatrix d;
  
  FEInterpolation* interp = gp->giveElement()->giveInterpolation();
  // this uses FEIInterpolation::nodes2coords - quite inefficient in this case (large num of dofmans)
  interp->evaldNdx (d, *gp->giveCoordinates(), 
		    FEIIGAElementGeometryWrapper(gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan()), 0.0);
  
  if ((nDofMan = interp->giveNumberOfKnotSpanBasisFunctions (*(gp->giveIntegrationRule()->giveKnotSpan()))) == 0)  // HUHU
    nDofMan = gp->giveElement()->giveNumberOfDofManagers();
  
  answer.resize(6, nDofMan*3);
  answer.zero();
  
  for (i=1; i <= nDofMan; i++) {
    answer.at(1, i*3-2) = d.at(i, 1);
    answer.at(2, i*3-1) = d.at(i, 2);
		answer.at(3, i*3-0) = d.at(i, 3);
    
    answer.at(4, 3*i-1) = d.at(i, 3);
    answer.at(4, 3*i-0) = d.at(i, 2);

    answer.at(5, 3*i-2) = d.at(i, 3);
    answer.at(5, 3*i-0) = d.at(i, 1);

    answer.at(6, 3*i-2) = d.at(i, 2);
    answer.at(6, 3*i-1) = d.at(i, 1);
  }
}

double Space3dStructuralElementEvaluator::computeVolumeAround(GaussPoint *gp) { 
  double determinant, weight, volume;
  determinant = fabs( this->giveElement()->giveInterpolation()
                      ->giveTransformationJacobian(*gp->giveCoordinates(), 
						   FEIIGAElementGeometryWrapper(this->giveElement(), 
										gp->giveIntegrationRule()->giveKnotSpan()),
						   0.0) );
  weight      = gp->giveWeight();
  volume      = determinant * weight;
    
  return volume;
}


NURBSSpace3dElement::NURBSSpace3dElement (int n, Domain *aDomain) : IGAElement (n, aDomain), Space3dStructuralElementEvaluator(), interpolation(3) {}


IRResultType NURBSSpace3dElement::initializeFrom(InputRecord *ir) {
	NURBSInterpolation *interpol = (NURBSInterpolation *)this->giveInterpolation();
	IGAElement::initializeFrom(ir);
	//PlaneStressStructuralElementEvaluator::initializeFrom(ir);

	// HUHU
	if(giveNumberOfDofManagers() != interpol->giveNumberOfControlPoints(1)*interpol->giveNumberOfControlPoints(2)*interpol->giveNumberOfControlPoints(3)){
    OOFEM_ERROR("NURBSSpace3dElement::initializeFrom - number of control points mismatch");
	}

	return IRRT_OK;
}



#ifdef __OOFEG

#define DRAW_MESH

// if DRAW_MESH is defined only boundary of integration elements are drawn;
// currently mesh is not properly drawn for tsplines
// because integration elements (does not matter whether single span or multi span)
// are generaly finer than T-mesh;

void IGAElement::drawRawGeometry(oofegGraphicContext &gc) {
    WCRec p [ 8 ];
    GraphicObj *go;
    FEInterpolation* interp = this->giveInterpolation();
    int i, j, k, m, nseg = 4;

#ifdef DRAW_MESH
    WCRec pp [ 2 ];
#endif

    if ( !gc.testElementGraphicActivity(this) ) {
			return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(TRUE);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);

#ifdef DRAW_MESH
		EASValsSetLineWidth(0);
		EASValsSetLineStyle(SOLID_STYLE);
		nseg = 8;
#else
    EASValsSetFillStyle(FILL_HOLLOW);
#endif

    int numberOfIntegrationRules = this->giveNumberOfIntegrationRules();
    double** const  knotVector = interp->giveKnotVector();
    const IntArray *span;
    IntegrationRule* iRule;
    int ir, nsd = this->giveNsd();;
		
    if (nsd == 1) {
			FloatArray c[2], cg[2];
			double du;
			
      for (j=0; j<2; j++) {
        c[j].resize(1);
        cg[j].resize(1);
      }
			
			// loop over individual integration rules (i.e., knot spans)
			for (ir=0; ir < numberOfIntegrationRules; ir++) {
				iRule = this->giveIntegrationRule(ir);
				span = iRule->giveKnotSpan();
				// divide span locally to get finer geometry rep.
				du = (knotVector[0][span->at(1)+1] - knotVector[0][span->at(1)])/ nseg;
				for (i=1; i<=nseg; i++) {
					c[0].at(1) = knotVector[0][span->at(1)] + du*(i-1);
					c[1].at(1) = knotVector[0][span->at(1)] + du*i;
					
					for (k=0; k<2; k++) {
						interp->local2global (cg[k], c[k], FEIIGAElementGeometryWrapper(this, iRule->giveKnotSpan()), 0.0);
						p [ k ].x = ( FPNum ) cg[k].at(1);
						p [ k ].y = 0.;
						p [ k ].z = 0.;
					}
					
					go =  CreateLine3D(p);
					EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
					EGAttachObject(go, ( EObjectP ) this);
					EMAddGraphicsToModel(ESIModel(), go);
        }
			} // end loop over knot spans (irules)
		}
		else if(nsd == 2){
			FloatArray c[4], cg[4];
			double du, dv;
			
      for (j=0; j<4; j++) {
        c[j].resize(2);
        cg[j].resize(2);
      }
			
			// loop over individual integration rules (i.e., knot spans)
			for (ir=0; ir < numberOfIntegrationRules; ir++) {
				iRule = this->giveIntegrationRule(ir);
				span = iRule->giveKnotSpan();
				// divide span locally to get finer geometry rep.
				du = (knotVector[0][span->at(1)+1] - knotVector[0][span->at(1)])/ nseg;
				dv = (knotVector[1][span->at(2)+1] - knotVector[1][span->at(2)])/ nseg;
				for (i=1; i<=nseg; i++) {
          for (j=1; j<=nseg; j++) {
            c[0].at(1) = knotVector[0][span->at(1)] + du*(i-1);
            c[0].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
            c[1].at(1) = knotVector[0][span->at(1)] + du*i;
            c[1].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
            c[2].at(1) = knotVector[0][span->at(1)] + du*i;
            c[2].at(2) = knotVector[1][span->at(2)] + dv*j;
            c[3].at(1) = knotVector[0][span->at(1)] + du*(i-1);
            c[3].at(2) = knotVector[1][span->at(2)] + dv*j;
            
            for (k=0; k<4; k++) {
              interp->local2global (cg[k], c[k], FEIIGAElementGeometryWrapper(this, iRule->giveKnotSpan()), 0.0);
              p [ k ].x = ( FPNum ) cg[k].at(1);
              p [ k ].y = ( FPNum ) cg[k].at(2);
              p [ k ].z = 0.;
            }
						
#ifdef DRAW_MESH
						if(i == 1){
							pp[0].x = p[0].x;
							pp[0].y = p[0].y;
							pp[0].z = p[0].z;
							
							pp[1].x = p[3].x;
							pp[1].y = p[3].y;
							pp[1].z = p[3].z;

						  go =  CreateLine3D(pp);
							EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
							EGAttachObject(go, ( EObjectP ) this);
							EMAddGraphicsToModel(ESIModel(), go);
						}
						if(j == 1){
							pp[0].x = p[0].x;
							pp[0].y = p[0].y;
							pp[0].z = p[0].z;
							
							pp[1].x = p[1].x;
							pp[1].y = p[1].y;
							pp[1].z = p[1].z;

						  go =  CreateLine3D(pp);
							EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
							EGAttachObject(go, ( EObjectP ) this);
							EMAddGraphicsToModel(ESIModel(), go);
						}
						if(i == nseg){
							pp[0].x = p[1].x;
							pp[0].y = p[1].y;
							pp[0].z = p[1].z;
							
							pp[1].x = p[2].x;
							pp[1].y = p[2].y;
							pp[1].z = p[2].z;

						  go =  CreateLine3D(pp);
							EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
							EGAttachObject(go, ( EObjectP ) this);
							EMAddGraphicsToModel(ESIModel(), go);
						}
						if(j == nseg){
							pp[0].x = p[2].x;
							pp[0].y = p[2].y;
							pp[0].z = p[2].z;
							
							pp[1].x = p[3].x;
							pp[1].y = p[3].y;
							pp[1].z = p[3].z;

						  go =  CreateLine3D(pp);
							EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
							EGAttachObject(go, ( EObjectP ) this);
							EMAddGraphicsToModel(ESIModel(), go);
						}
#else
            go =  CreateQuad3D(p);
            EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
            EGAttachObject(go, ( EObjectP ) this);
            EMAddGraphicsToModel(ESIModel(), go);
#endif
          }
        }
			} // end loop over knot spans (irules)
		}
		else if(nsd == 3){
			FloatArray c[8], cg[8];
			double du, dv, dt;
			
      for (j=0; j<8; j++) {
        c[j].resize(3);
        cg[j].resize(3);
      }
			
			// loop over individual integration rules (i.e., knot spans)
			for (ir=0; ir < numberOfIntegrationRules; ir++) {
				iRule = this->giveIntegrationRule(ir);
				span = iRule->giveKnotSpan();
				// divide span locally to get finer geometry rep.
				du = (knotVector[0][span->at(1)+1] - knotVector[0][span->at(1)])/ nseg;
				dv = (knotVector[1][span->at(2)+1] - knotVector[1][span->at(2)])/ nseg;
				dt = (knotVector[2][span->at(3)+1] - knotVector[2][span->at(3)])/ nseg;
				for (i=1; i<=nseg; i++) {
          for (j=1; j<=nseg; j++) {
						for (k=1; k<=nseg; k++) {
							c[0].at(1) = knotVector[0][span->at(1)] + du*(i-1);
							c[0].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
							c[0].at(3) = knotVector[2][span->at(3)] + dt*(k-1);
							c[1].at(1) = knotVector[0][span->at(1)] + du*i;
							c[1].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
							c[1].at(3) = knotVector[2][span->at(3)] + dt*(k-1);
							c[2].at(1) = knotVector[0][span->at(1)] + du*i;
							c[2].at(2) = knotVector[1][span->at(2)] + dv*j;
							c[2].at(3) = knotVector[2][span->at(3)] + dt*(k-1);
							c[3].at(1) = knotVector[0][span->at(1)] + du*(i-1);
							c[3].at(2) = knotVector[1][span->at(2)] + dv*j;
							c[3].at(3) = knotVector[2][span->at(3)] + dt*(k-1);
							c[4].at(1) = knotVector[0][span->at(1)] + du*(i-1);
							c[4].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
							c[4].at(3) = knotVector[2][span->at(3)] + dt*k;
							c[5].at(1) = knotVector[0][span->at(1)] + du*i;
							c[5].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
							c[5].at(3) = knotVector[2][span->at(3)] + dt*k;
							c[6].at(1) = knotVector[0][span->at(1)] + du*i;
							c[6].at(2) = knotVector[1][span->at(2)] + dv*j;
							c[6].at(3) = knotVector[2][span->at(3)] + dt*k;
							c[7].at(1) = knotVector[0][span->at(1)] + du*(i-1);
							c[7].at(2) = knotVector[1][span->at(2)] + dv*j;
							c[7].at(3) = knotVector[2][span->at(3)] + dt*k;
							
							for (m=0; m<8; m++) {
								interp->local2global (cg[m], c[m], FEIIGAElementGeometryWrapper(this, iRule->giveKnotSpan()), 0.0);
								p [ m ].x = ( FPNum ) cg[m].at(1);
								p [ m ].y = ( FPNum ) cg[m].at(2);
								p [ m ].z = ( FPNum ) cg[m].at(3);
							}

#ifdef DRAW_MESH
							if(i == 1 && j == 1){
								pp[0].x = p[0].x;
								pp[0].y = p[0].y;
								pp[0].z = p[0].z;
								
								pp[1].x = p[4].x;
								pp[1].y = p[4].y;
								pp[1].z = p[4].z;

								go =  CreateLine3D(pp);
								EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
								EGAttachObject(go, ( EObjectP ) this);
								EMAddGraphicsToModel(ESIModel(), go);
							}
							if(i == 1 && j == nseg){
								pp[0].x = p[3].x;
								pp[0].y = p[3].y;
								pp[0].z = p[3].z;
								
								pp[1].x = p[7].x;
								pp[1].y = p[7].y;
								pp[1].z = p[7].z;

								go =  CreateLine3D(pp);
								EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
								EGAttachObject(go, ( EObjectP ) this);
								EMAddGraphicsToModel(ESIModel(), go);
							}
							if(i == nseg && j == 1){
								pp[0].x = p[1].x;
								pp[0].y = p[1].y;
								pp[0].z = p[1].z;
								
								pp[1].x = p[5].x;
								pp[1].y = p[5].y;
								pp[1].z = p[5].z;

								go =  CreateLine3D(pp);
								EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
								EGAttachObject(go, ( EObjectP ) this);
								EMAddGraphicsToModel(ESIModel(), go);
							}
							if(i == nseg && j == nseg){
								pp[0].x = p[2].x;
								pp[0].y = p[2].y;
								pp[0].z = p[2].z;
								
								pp[1].x = p[6].x;
								pp[1].y = p[6].y;
								pp[1].z = p[6].z;

								go =  CreateLine3D(pp);
								EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
								EGAttachObject(go, ( EObjectP ) this);
								EMAddGraphicsToModel(ESIModel(), go);
							}
							
							if(j == 1 && k == 1){
								pp[0].x = p[0].x;
								pp[0].y = p[0].y;
								pp[0].z = p[0].z;
								
								pp[1].x = p[1].x;
								pp[1].y = p[1].y;
								pp[1].z = p[1].z;

								go =  CreateLine3D(pp);
								EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
								EGAttachObject(go, ( EObjectP ) this);
								EMAddGraphicsToModel(ESIModel(), go);
							}
							if(j == 1 && k == nseg){
								pp[0].x = p[4].x;
								pp[0].y = p[4].y;
								pp[0].z = p[4].z;
								
								pp[1].x = p[5].x;
								pp[1].y = p[5].y;
								pp[1].z = p[5].z;

								go =  CreateLine3D(pp);
								EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
								EGAttachObject(go, ( EObjectP ) this);
								EMAddGraphicsToModel(ESIModel(), go);
							}
							if(j == nseg && k == 1){
								pp[0].x = p[3].x;
								pp[0].y = p[3].y;
								pp[0].z = p[3].z;
								
								pp[1].x = p[2].x;
								pp[1].y = p[2].y;
								pp[1].z = p[2].z;

								go =  CreateLine3D(pp);
								EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
								EGAttachObject(go, ( EObjectP ) this);
								EMAddGraphicsToModel(ESIModel(), go);
							}
							if(j == nseg && k == nseg){
								pp[0].x = p[7].x;
								pp[0].y = p[7].y;
								pp[0].z = p[7].z;
								
								pp[1].x = p[6].x;
								pp[1].y = p[6].y;
								pp[1].z = p[6].z;

								go =  CreateLine3D(pp);
								EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
								EGAttachObject(go, ( EObjectP ) this);
								EMAddGraphicsToModel(ESIModel(), go);
							}
							
							if(k == 1 && i == 1){
								pp[0].x = p[0].x;
								pp[0].y = p[0].y;
								pp[0].z = p[0].z;
								
								pp[1].x = p[3].x;
								pp[1].y = p[3].y;
								pp[1].z = p[3].z;

								go =  CreateLine3D(pp);
								EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
								EGAttachObject(go, ( EObjectP ) this);
								EMAddGraphicsToModel(ESIModel(), go);
							}
							if(k == 1 && i == nseg){
								pp[0].x = p[1].x;
								pp[0].y = p[1].y;
								pp[0].z = p[1].z;
								
								pp[1].x = p[2].x;
								pp[1].y = p[2].y;
								pp[1].z = p[2].z;

								go =  CreateLine3D(pp);
								EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
								EGAttachObject(go, ( EObjectP ) this);
								EMAddGraphicsToModel(ESIModel(), go);
							}
							if(k == nseg && i == 1){
								pp[0].x = p[4].x;
								pp[0].y = p[4].y;
								pp[0].z = p[4].z;
								
								pp[1].x = p[7].x;
								pp[1].y = p[7].y;
								pp[1].z = p[7].z;

								go =  CreateLine3D(pp);
								EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
								EGAttachObject(go, ( EObjectP ) this);
								EMAddGraphicsToModel(ESIModel(), go);
							}
							if(k == nseg && i == nseg){
								pp[0].x = p[5].x;
								pp[0].y = p[5].y;
								pp[0].z = p[5].z;
								
								pp[1].x = p[6].x;
								pp[1].y = p[6].y;
								pp[1].z = p[6].z;

								go =  CreateLine3D(pp);
								EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
								EGAttachObject(go, ( EObjectP ) this);
								EMAddGraphicsToModel(ESIModel(), go);
							}
#else
							go =  CreateHexahedron(p);
							EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
							EGAttachObject(go, ( EObjectP ) this);
							EMAddGraphicsToModel(ESIModel(), go);
#endif
						}
					}
				}
			} // end loop over knot spans (irules)
		}
		else{
      OOFEM_ERROR2 ("drawRawGeometry: not implemented for nsd = %d", nsd); 
		}
}


void drawIGAPatchDeformedGeometry(Element* elem, StructuralElementEvaluator* se, oofegGraphicContext &gc, UnknownType) {
    WCRec p [ 8 ];
    GraphicObj *go;
    int i, j, k, m, n, nseg = 4;
		FloatArray u;
    FloatMatrix N;
    FloatArray ur, d;
    IntArray lc;
    FEInterpolation* interp = elem->giveInterpolation();
    TimeStep *stepN = elem->giveDomain()->giveEngngModel()->giveCurrentStep();
		double defScale = gc.getDefScale();


#ifdef DRAW_MESH
    WCRec pp [ 2 ];
#endif

    if ( !gc.testElementGraphicActivity(elem) ) {
			return;
    }

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(TRUE);
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);

#ifdef DRAW_MESH
		EASValsSetLineWidth(0);
		EASValsSetLineStyle(SOLID_STYLE);
    EASValsSetFillStyle(FILL_HOLLOW);
		nseg = 8;
#else
    EASValsSetFillStyle(FILL_SOLID);
#endif

    int numberOfIntegrationRules = elem->giveNumberOfIntegrationRules();
    double** const  knotVector = interp->giveKnotVector();
    const IntArray *span;
    IntegrationRule* iRule;
    int ir, nsd = interp->giveNsd();;
		
		se->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);

		if ( se->updateRotationMatrix() ) {
			u.rotatedWith(se->rotationMatrix, 'n');
		}

    if (nsd == 1) {
			FloatArray c[2], cg[2];
			double du;
			
      for (j=0; j<2; j++) {
        c[j].resize(1);
        cg[j].resize(1);
      }
			
			// loop over individual integration rules (i.e., knot spans)
			for (ir=0; ir < numberOfIntegrationRules; ir++) {
				iRule = elem->giveIntegrationRule(ir);
				span = iRule->giveKnotSpan();
				// divide span locally to get finer geometry rep.
				du = (knotVector[0][span->at(1)+1] - knotVector[0][span->at(1)])/ nseg;
				for (i=1; i<=nseg; i++) {
					c[0].at(1) = knotVector[0][span->at(1)] + du*(i-1);
					c[1].at(1) = knotVector[0][span->at(1)] + du*i;
					
					for (k=0; k<2; k++) {
					// create a dummy ip's
						FloatArray *cc = new FloatArray;
						cc->beCopyOf (&c[k]);  // constructor of gp does not make its own copy  
						GaussPoint gp(iRule, 999, cc, 1.0, _PlaneStress);            
						
						// compute displacements at gp
						se->computeNMatrixAt(N, &gp);
						
						// get local code numbers corresponding to ir
						se->giveIntegrationElementLocalCodeNumbers (lc, elem, gp.giveIntegrationRule(), EID_MomentumBalance);
						ur.resize(N.giveNumberOfColumns());
						for (n=1; n<=lc.giveSize(); n++) ur.at(n) = u.at(lc.at(n));
						
						// interpolate displacements
						d.beProductOf(N, ur);
						
						interp->local2global (cg[k], c[k], FEIIGAElementGeometryWrapper(elem, iRule->giveKnotSpan()), 0.0);
						p [ k ].x = ( FPNum ) (cg[k].at(1) + d.at(1)*defScale);
						p [ k ].y = 0.;
						p [ k ].z = 0.;
					}
					
					go =  CreateLine3D(p);
					EGWithMaskChangeAttributes(WIDTH_MASK | STYLE_MASK | COLOR_MASK | LAYER_MASK, go);
					EGAttachObject(go, ( EObjectP ) elem);
					EMAddGraphicsToModel(ESIModel(), go);
        }
			} // end loop over knot spans (irules)
		}
		else if(nsd == 2){
			FloatArray c[4], cg[4];
			double du, dv;
			
      for (j=0; j<4; j++) {
        c[j].resize(2);
        cg[j].resize(2);
      }
			
			// loop over individual integration rules (i.e., knot spans)
			for (ir=0; ir < numberOfIntegrationRules; ir++) {
				iRule = elem->giveIntegrationRule(ir);
				span = iRule->giveKnotSpan();
				// divide span locally to get finer geometry rep.
				du = (knotVector[0][span->at(1)+1] - knotVector[0][span->at(1)])/ nseg;
				dv = (knotVector[1][span->at(2)+1] - knotVector[1][span->at(2)])/ nseg;
				for (i=1; i<=nseg; i++) {
          for (j=1; j<=nseg; j++) {
            c[0].at(1) = knotVector[0][span->at(1)] + du*(i-1);
            c[0].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
            c[1].at(1) = knotVector[0][span->at(1)] + du*i;
            c[1].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
            c[2].at(1) = knotVector[0][span->at(1)] + du*i;
            c[2].at(2) = knotVector[1][span->at(2)] + dv*j;
            c[3].at(1) = knotVector[0][span->at(1)] + du*(i-1);
            c[3].at(2) = knotVector[1][span->at(2)] + dv*j;
            
						for (k=0; k<4; k++) {
							// create a dummy ip's
							FloatArray *cc = new FloatArray;
							cc->beCopyOf (&c[k]);  // constructor of gp does not make its own copy  
							GaussPoint gp(iRule, 999, cc, 1.0, _PlaneStress);            
							
							// compute displacements at gp
							se->computeNMatrixAt(N, &gp);
							
							// get local code numbers corresponding to ir
							se->giveIntegrationElementLocalCodeNumbers (lc, elem, gp.giveIntegrationRule(), EID_MomentumBalance);
							ur.resize(N.giveNumberOfColumns());
							for (n=1; n<=lc.giveSize(); n++) ur.at(n) = u.at(lc.at(n));
							
							// interpolate displacements
							d.beProductOf(N, ur);
							
							interp->local2global (cg[k], c[k], FEIIGAElementGeometryWrapper(elem, iRule->giveKnotSpan()), 0.0);
							p [ k ].x = ( FPNum ) (cg[k].at(1) + d.at(1)*defScale);
							p [ k ].y = ( FPNum ) (cg[k].at(2) + d.at(2)*defScale);
							p [ k ].z = 0.;
						}

            go =  CreateQuad3D(p);
            EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
            EGAttachObject(go, ( EObjectP ) elem);
            EMAddGraphicsToModel(ESIModel(), go);
          }
        }
			} // end loop over knot spans (irules)
		}
		else if(nsd == 3){
			FloatArray c[8], cg[8];
			double du, dv, dt;
			
      for (j=0; j<8; j++) {
        c[j].resize(3);
        cg[j].resize(3);
      }
			
			// loop over individual integration rules (i.e., knot spans)
			for (ir=0; ir < numberOfIntegrationRules; ir++) {
				iRule = elem->giveIntegrationRule(ir);
				span = iRule->giveKnotSpan();
				// divide span locally to get finer geometry rep.
				du = (knotVector[0][span->at(1)+1] - knotVector[0][span->at(1)])/ nseg;
				dv = (knotVector[1][span->at(2)+1] - knotVector[1][span->at(2)])/ nseg;
				dt = (knotVector[2][span->at(3)+1] - knotVector[2][span->at(3)])/ nseg;
				for (i=1; i<=nseg; i++) {
          for (j=1; j<=nseg; j++) {
						for (k=1; k<=nseg; k++) {
							c[0].at(1) = knotVector[0][span->at(1)] + du*(i-1);
							c[0].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
							c[0].at(3) = knotVector[2][span->at(3)] + dt*(k-1);
							c[1].at(1) = knotVector[0][span->at(1)] + du*i;
							c[1].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
							c[1].at(3) = knotVector[2][span->at(3)] + dt*(k-1);
							c[2].at(1) = knotVector[0][span->at(1)] + du*i;
							c[2].at(2) = knotVector[1][span->at(2)] + dv*j;
							c[2].at(3) = knotVector[2][span->at(3)] + dt*(k-1);
							c[3].at(1) = knotVector[0][span->at(1)] + du*(i-1);
							c[3].at(2) = knotVector[1][span->at(2)] + dv*j;
							c[3].at(3) = knotVector[2][span->at(3)] + dt*(k-1);
							c[4].at(1) = knotVector[0][span->at(1)] + du*(i-1);
							c[4].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
							c[4].at(3) = knotVector[2][span->at(3)] + dt*k;
							c[5].at(1) = knotVector[0][span->at(1)] + du*i;
							c[5].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
							c[5].at(3) = knotVector[2][span->at(3)] + dt*k;
							c[6].at(1) = knotVector[0][span->at(1)] + du*i;
							c[6].at(2) = knotVector[1][span->at(2)] + dv*j;
							c[6].at(3) = knotVector[2][span->at(3)] + dt*k;
							c[7].at(1) = knotVector[0][span->at(1)] + du*(i-1);
							c[7].at(2) = knotVector[1][span->at(2)] + dv*j;
							c[7].at(3) = knotVector[2][span->at(3)] + dt*k;
							
							for (m=0; m<8; m++) {
								// create a dummy ip's
								FloatArray *cc = new FloatArray;
								cc->beCopyOf (&c[m]);  // constructor of gp does not make its own copy  
								GaussPoint gp(iRule, 999, cc, 1.0, _3dMat);            
								
								// compute displacements at gp
								se->computeNMatrixAt(N, &gp);
								
								// get local code numbers corresponding to ir
								se->giveIntegrationElementLocalCodeNumbers (lc, elem, gp.giveIntegrationRule(), EID_MomentumBalance);
								ur.resize(N.giveNumberOfColumns());
								for (n=1; n<=lc.giveSize(); n++) ur.at(n) = u.at(lc.at(n));
								
								// interpolate displacements
								d.beProductOf(N, ur);
								
								interp->local2global (cg[m], c[m], FEIIGAElementGeometryWrapper(elem, iRule->giveKnotSpan()), 0.0);
								p [ m ].x = ( FPNum ) (cg[m].at(1) + d.at(1)*defScale);
								p [ m ].y = ( FPNum ) (cg[m].at(2) + d.at(2)*defScale);
								p [ m ].z = ( FPNum ) (cg[m].at(3) + d.at(3)*defScale);
							}

							go =  CreateHexahedron(p);
							EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
							EGAttachObject(go, ( EObjectP ) elem);
							EMAddGraphicsToModel(ESIModel(), go);
						}
					}
				}
			} // end loop over knot spans (irules)
		}
		else{
      OOFEM_ERROR2 ("drawDeformedGeometry: not implemented for nsd = %d", nsd); 
		}
}





// HUHU should be implemented by IGA element (it is the same for Bspline NURBS and TSpline)
// however in such a case it should be generalized in terms of appropriately multiplying 
// nseq for those integration elements which span more tham just a single knot span
// the reason is to ensure compatible division to quads over which scalar quantity is interpolated
// bilinearly !!!

void BsplinePlaneStressElement::drawScalar(oofegGraphicContext &context) {
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
    EASValsSetFillStyle(FILL_SOLID);
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
        for (i=1; i<=nseg; i++) {
          for (j=1; j<=nseg; j++) {
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

							/*
							//move sampling gp out of boundary to overcome degeneracy on quarter plate with hole modelled by single patch
							if(c[k].at(1) > 0.99999 && c[k].at(2) > 0.495 && c[k].at(2) < 0.505){
								c[k].at(1) += sign[k].at(1)*du/10.0;
								c[k].at(2) += sign[k].at(2)*dv/10.0;
								c[k].at(3) += sign[k].at(3)*dw/10.0;
							}
							*/

              // create a dummy ip's
              FloatArray *cc = new FloatArray;
              cc->beCopyOf (&c[k]);  // constructor of gp does not make its own copy  
              GaussPoint gp(iRule, 999, cc, 1.0, _PlaneStress);            
              // this->computeStrainVector (val, &gp, tStep);
              this->computeStressVector (val, &gp, tStep);
              s[k]=val.at(indx);
            }

            if ((isnan(s[0])) || (isnan(s[1])) || (isnan(s[2])) || (isnan(s[3]))) continue;
            if ((fabs(s[0])>1.e5) || (fabs(s[1])>1.e5) || (fabs(s[2])>1.e5) || (fabs(s[3])>1.e5)) continue;
            //printf ("QWD: %e %e %e %e\n", s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);

            go =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
            EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
            EGAttachObject(go, ( EObjectP ) this);
            EMAddGraphicsToModel(ESIModel(), go);
          }
        }
      }
    } // end loop over knot spans (irules)
}


void NURBSPlaneStressElement::drawScalar(oofegGraphicContext &context) {
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
    EASValsSetFillStyle(FILL_SOLID);
    EASValsSetEdgeFlag(FALSE);
    numberOfIntegrationRules = this->giveNumberOfIntegrationRules();
    double** const  knotVector = interp->giveKnotVector();
    const IntArray *span;
    IntegrationRule* iRule;
    int ir, j, nsd = this->giveNsd();
    int numberOfKnotSpans_u = interp->giveNumberOfKnotSpans(1);
    int numberOfKnotSpans_v = interp->giveNumberOfKnotSpans(2);
    FloatArray c[4], cg[4];
		IntArray sign[4];

    if (nsd == 2) {
      for (j=0; j<4; j++) {
        c[j].resize(2);
        cg[j].resize(2);
				sign[j].resize(2);
      }
			
			sign[0].at(1) = 1;
			sign[0].at(2) = 1;
			sign[1].at(1) = -1;
			sign[1].at(2) = 1;
			sign[2].at(1) = -1;
			sign[2].at(2) = -1;
			sign[3].at(1) = 1;
			sign[3].at(2) = -1;
    } else {
      OOFEM_ERROR2 ("drawRawGeometry: not implemented for nsd = %d", nsd); 
    }

    this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
    if ( ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
      return;
    }

		double maxs=-100, mins=100.0;

    // loop over individual integration rules (i.e., knot spans)
    for (ir=0; ir < numberOfIntegrationRules; ir++) {
      iRule = this->giveIntegrationRule(ir);
      span = iRule->giveKnotSpan();
      if (nsd == 2) {

        // divide span locally to get finer geometry rep.
        int i,j,k, nseg = 8;
        double du = (knotVector[0][span->at(1)+1] - knotVector[0][span->at(1)])/ nseg;
        double dv = (knotVector[1][span->at(2)+1] - knotVector[1][span->at(2)])/ nseg;
        for (i=1; i<=nseg; i++) {
          for (j=1; j<=nseg; j++) {
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

							//move sampling gp out of boundary to overcome degeneracy on quarter plate with hole modelled by single patch
							if(c[k].at(1) > 0.99999 && c[k].at(2) > 0.495 && c[k].at(2) < 0.505){
								c[k].at(1) += sign[k].at(1)*du/10.0;
								c[k].at(2) += sign[k].at(2)*dv/10.0;
							}

              // create a dummy ip's
              FloatArray *cc = new FloatArray;
              cc->beCopyOf (&c[k]);  // constructor of gp does not make its own copy     
              GaussPoint gp(iRule, 999, cc, 1.0, _PlaneStress);            
              // this->computeStrainVector (val, &gp, tStep);
              this->computeStressVector (val, &gp, tStep);
              s[k]=val.at(indx);

							/*
							double x, y, r, phi, rate, E, G, kap, ny;
							x = cg[k].at(1);
							y = cg[k].at(2);
							if(x < 1.0e-10){
								phi = M_PI / 2.0;
								r = y;
							}
							else{
								phi = atan(y/x);
								r = x/cos(phi);
							}

							// exact stresses quarter plate with hole s0=1 a=1
							rate=1.0/r;
							rate *= rate;
							if(indx == 1)s[k]=0.5*(2.0+3.0*rate*rate*cos(4.0*phi)-rate*(3*cos(2.0*phi)+2.0*cos(4.0*phi)));
							if(indx == 2)s[k]=0.5*(-3.0*rate*rate*cos(4.0*phi)-rate*(cos(2.0*phi)-2.0*cos(4.0*phi)));
							if(indx == 3)s[k]=0.5*(3.0*rate*rate*sin(4.0*phi)-rate*(sin(2.0*phi)+2.0*sin(4.0*phi)));

							if(indx ==2){
								if(cg[k].at(1) <= 1.0 + 1.0e-10 && cg[k].at(2) <= 1.0e-10){
									fprintf(stderr, "A: syy = %e\n", s[k]);
								}
							}
							if(indx ==1){
								if(cg[k].at(1) <= 1.0e-10 && cg[k].at(2) <= 1.0 + 1.0e-10){
									fprintf(stderr, "B: sxx = %e\n", s[k]);
								}
							}

							// exact displ quarter plate with hole s0=1 a=1
							E = 15000.0;
							ny = 0.25;
							G = E/(2.0*(1.0+ny));
							kap = (3.0 - ny)/(1.0+ny);
							rate=1.0/r;
							if(indx == 1)s[k]=(r*(kap+1.0)*cos(phi)+2.0*rate*((1.0+kap)*cos(phi)+cos(3.0*phi))-2.0*rate*rate*rate*cos(3.0*phi))/(8.0*G);
							if(indx == 2)s[k]=(r*(kap-3.0)*sin(phi)+2.0*rate*((1.0-kap)*sin(phi)+sin(3.0*phi))-2.0*rate*rate*rate*sin(3.0*phi))/(8.0*G);

							if(indx ==1){
								if(cg[k].at(1) <= 1.0 + 1.0e-10 && cg[k].at(2) <= 1.0e-10){
									fprintf(stderr, "A: ux = %e\n", s[k]);
								}
							}
							if(indx ==2){
								if(cg[k].at(1) <= 1.0e-10 && cg[k].at(2) <= 1.0 + 1.0e-10){
									fprintf(stderr, "B: uy = %e\n", s[k]);
								}
							}
 
							if(s[k] < mins)mins=s[k];
							if(s[k] > maxs)maxs=s[k];
							*/

							if(indx ==2){
								if(cg[k].at(1) <= 1.0 + 1.0e-10 && cg[k].at(2) <= 1.0e-10){
									fprintf(stderr, "A: syy = %e\n", s[k]);
								}
							}
							if(indx ==1){
								if(cg[k].at(1) <= 1.0e-10 && cg[k].at(2) <= 1.0 + 1.0e-10){
									fprintf(stderr, "B: sxx = %e\n", s[k]);
								}
							}
            }

            if ((isnan(s[0])) || (isnan(s[1])) || (isnan(s[2])) || (isnan(s[3]))) continue;
            if ((fabs(s[0])>1.e5) || (fabs(s[1])>1.e5) || (fabs(s[2])>1.e5) || (fabs(s[3])>1.e5)) continue;
            //printf ("QWD: %e %e %e %e\n", s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);

            go =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
            EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
            EGAttachObject(go, ( EObjectP ) this);
            EMAddGraphicsToModel(ESIModel(), go);
          }
        }
      }
    } // end loop over knot spans (irules)

		//fprintf(stderr, "%d %e %e %e %e\n", indx, mins, maxs, (10.0*mins+maxs)/11.0, (10.0*maxs+mins)/11.0);
}


// refinement of integration elements should be generalized in terms of appropriately multiplying 
// nseq for those integration elements which span more tham just a single knot span
// the reason is to ensure compatible division to quads over which scalar quantity is interpolated
// bilinearly !!!

void TSplinePlaneStressElement::drawScalar(oofegGraphicContext &context) {
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
    EASValsSetFillStyle(FILL_SOLID);
    EASValsSetEdgeFlag(FALSE);
    numberOfIntegrationRules = this->giveNumberOfIntegrationRules();
    double** const  knotVector = interp->giveKnotVector();
    const IntArray *span;
    IntegrationRule* iRule;
    int ir, j, nsd = this->giveNsd();
    int numberOfKnotSpans_u = interp->giveNumberOfKnotSpans(1);
    int numberOfKnotSpans_v = interp->giveNumberOfKnotSpans(2);
    FloatArray c[4], cg[4];
		IntArray sign[4];

    if (nsd == 2) {
      for (j=0; j<4; j++) {
        c[j].resize(2);
        cg[j].resize(2);
				sign[j].resize(2);
      }
			
			sign[0].at(1) = 1;
			sign[0].at(2) = 1;
			sign[1].at(1) = -1;
			sign[1].at(2) = 1;
			sign[2].at(1) = -1;
			sign[2].at(2) = -1;
			sign[3].at(1) = 1;
			sign[3].at(2) = -1;
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
        for (i=1; i<=nseg; i++) {
          for (j=1; j<=nseg; j++) {
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

							/*
							//move sampling gp out of boundary to overcome degeneracy on quarter plate with hole modelled by single patch
							if(c[k].at(1) > 0.99999 && c[k].at(2) > 0.495 && c[k].at(2) < 0.505){
								c[k].at(1) += sign[k].at(1)*du/10.0;
								c[k].at(2) += sign[k].at(2)*dv/10.0;
							}
							*/

              // create a dummy ip's
              FloatArray *cc = new FloatArray;
              cc->beCopyOf (&c[k]);  // constructor of gp does not make its own copy
              GaussPoint gp(iRule, 999, cc, 1.0, _PlaneStress);            
              // this->computeStrainVector (val, &gp, tStep);
              this->computeStressVector (val, &gp, tStep);
              s[k]=val.at(indx);
            }

            if ((isnan(s[0])) || (isnan(s[1])) || (isnan(s[2])) || (isnan(s[3]))) continue;
            if ((fabs(s[0])>1.e5) || (fabs(s[1])>1.e5) || (fabs(s[2])>1.e5) || (fabs(s[3])>1.e5)) continue;
            //printf ("QWD: %e %e %e %e\n", s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);

            go =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
            EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
            EGAttachObject(go, ( EObjectP ) this);
            EMAddGraphicsToModel(ESIModel(), go);
          }
        }
      }
    } // end loop over knot spans (irules)
}



void NURBSSpace3dElement::drawScalar(oofegGraphicContext &context) {
    int indx;
    WCRec p [ 8 ];
    GraphicObj *go;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    double s [ 8 ];
    IntArray map;
    FEInterpolation* interp = this->giveInterpolation();
    FloatArray val;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    EASValsSetFillStyle(FILL_SOLID);
    EASValsSetEdgeFlag(FALSE);
    numberOfIntegrationRules = this->giveNumberOfIntegrationRules();
    double** const  knotVector = interp->giveKnotVector();
    const IntArray *span;
    IntegrationRule* iRule;
    int ir, j, nsd = this->giveNsd();
    int numberOfKnotSpans_u = interp->giveNumberOfKnotSpans(1);
    int numberOfKnotSpans_v = interp->giveNumberOfKnotSpans(2);
    int numberOfKnotSpans_w = interp->giveNumberOfKnotSpans(3);
    FloatArray c[8], cg[8];
		IntArray sign[8];

    if (nsd == 3) {
      for (j=0; j<8; j++) {
        c[j].resize(3);
        cg[j].resize(3);
				sign[j].resize(3);
      }
			
			sign[0].at(1) = 1;
			sign[0].at(2) = 1;
			sign[0].at(3) = 1;
			sign[1].at(1) = -1;
			sign[1].at(2) = 1;
			sign[1].at(3) = 1;
			sign[2].at(1) = -1;
			sign[2].at(2) = -1;
			sign[2].at(3) = 1;
			sign[3].at(1) = 1;
			sign[3].at(2) = -1;
			sign[3].at(3) = 1;
			sign[4].at(1) = 1;
			sign[4].at(2) = 1;
			sign[4].at(3) = -1;
			sign[5].at(1) = -1;
			sign[5].at(2) = 1;
			sign[5].at(3) = -1;
			sign[6].at(1) = -1;
			sign[6].at(2) = -1;
			sign[6].at(3) = -1;
			sign[7].at(1) = 1;
			sign[7].at(2) = -1;
			sign[7].at(3) = -1;
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
      if (nsd == 3) {

        // divide span locally to get finer geometry rep.
        int i,j,k,m, nseg = 8;
        double du = (knotVector[0][span->at(1)+1] - knotVector[0][span->at(1)])/ nseg;
        double dv = (knotVector[1][span->at(2)+1] - knotVector[1][span->at(2)])/ nseg;
        double dw = (knotVector[2][span->at(3)+1] - knotVector[2][span->at(3)])/ nseg;
        for (i=1; i<=nseg; i++) {
          for (j=1; j<=nseg; j++) {
						for (m=1; m<=nseg; m++) {
							c[0].at(1) = knotVector[0][span->at(1)] + du*(i-1);
							c[0].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
							c[0].at(3) = knotVector[2][span->at(3)] + dw*(m-1);
							c[1].at(1) = knotVector[0][span->at(1)] + du*i;
							c[1].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
							c[1].at(3) = knotVector[2][span->at(3)] + dw*(m-1);
							c[2].at(1) = knotVector[0][span->at(1)] + du*i;
							c[2].at(2) = knotVector[1][span->at(2)] + dv*j;
							c[2].at(3) = knotVector[2][span->at(3)] + dw*(m-1);
							c[3].at(1) = knotVector[0][span->at(1)] + du*(i-1);
							c[3].at(2) = knotVector[1][span->at(2)] + dv*j;
							c[3].at(3) = knotVector[2][span->at(3)] + dw*(m-1);
							c[4].at(1) = knotVector[0][span->at(1)] + du*(i-1);
							c[4].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
							c[4].at(3) = knotVector[2][span->at(3)] + dw*m;
							c[5].at(1) = knotVector[0][span->at(1)] + du*i;
							c[5].at(2) = knotVector[1][span->at(2)] + dv*(j-1);
							c[5].at(3) = knotVector[2][span->at(3)] + dw*m;
							c[6].at(1) = knotVector[0][span->at(1)] + du*i;
							c[6].at(2) = knotVector[1][span->at(2)] + dv*j;
							c[6].at(3) = knotVector[2][span->at(3)] + dw*m;
							c[7].at(1) = knotVector[0][span->at(1)] + du*(i-1);
							c[7].at(2) = knotVector[1][span->at(2)] + dv*j;
							c[7].at(3) = knotVector[2][span->at(3)] + dw*m;
            
							for (k=0;k<8; k++) {
								interp->local2global (cg[k], c[k], FEIIGAElementGeometryWrapper(this, iRule->giveKnotSpan()), 0.0);
								p [ k ].x = ( FPNum ) cg[k].at(1);
								p [ k ].y = ( FPNum ) cg[k].at(2);
								p [ k ].z = ( FPNum ) cg[k].at(3);

								// create a dummy ip's
								FloatArray *cc = new FloatArray;
								cc->beCopyOf (&c[k]);  // constructor of gp does not make its own copy     
								GaussPoint gp(iRule, 999, cc, 1.0, _3dMat);            
								// this->computeStrainVector (val, &gp, tStep);
								this->computeStressVector (val, &gp, tStep);
								s[k]=val.at(indx);
							}

							if ((isnan(s[0])) || (isnan(s[1])) || (isnan(s[2])) || (isnan(s[3]))) continue;
							if ((isnan(s[4])) || (isnan(s[5])) || (isnan(s[6])) || (isnan(s[7]))) continue;
							if ((fabs(s[0])>1.e5) || (fabs(s[1])>1.e5) || (fabs(s[2])>1.e5) || (fabs(s[3])>1.e5)) continue;
							if ((fabs(s[4])>1.e5) || (fabs(s[5])>1.e5) || (fabs(s[6])>1.e5) || (fabs(s[7])>1.e5)) continue;
							//printf ("HWD: %e %e %e %e %e %e %e %e\n", s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ], s [ 4 ], s [ 5 ], s [ 6 ], s [ 7 ]);
							
							go =  CreateHexahedronWD(p, s);
							EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
							EGAttachObject(go, ( EObjectP ) this);
							EMAddGraphicsToModel(ESIModel(), go);
						}
					}
        }
      }
    } // end loop over knot spans (irules)
}


#endif

