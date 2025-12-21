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
 *               Copyright (C) 1993 - 2025   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


#include "fecontactsurface.h"
#include "contactpoint.h"
#include "datastream.h"
#include "contextioerr.h"
#include "dynamicinputrecord.h"
#include "set.h"
#include "feinterpol2d.h"
#include "feinterpol3d.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"

namespace oofem {


FEContactSurface :: FEContactSurface(int n, Domain *d) : ContactSurface(n, d)
{
}


  
void FEContactSurface :: initializeFrom(InputRecord &ir)
{
    ContactSurface :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, this->contactElementSetNumber, _IFT_FEContactSurface_contactElementSetNumber);
}


void
FEContactSurface :: postInitialize()
{
  this->contactElementSet = domain->giveSet(this->contactElementSetNumber)->giveElementList();
}


ContactElement*
FEContactSurface :: giveContactElement(int i)
{
  if(i == -1) {
    return nullptr;
  } else {
    auto e = dynamic_cast<ContactElement*>( this->giveDomain()->giveElement(i));
    if(e) {
      return e;
    } else {
      OOFEM_ERROR("Element %d is not a contact element", i);
    }
  }
}

ContactElement*
FEContactSurface :: giveContactElement_InSet(int i)
{
  auto e = dynamic_cast<ContactElement*>( this->giveDomain()->giveElement(this->contactElementSet.at(i)));
  if(e) {
    return e;
  } else {
    OOFEM_ERROR("Element is not a contact element");
  }
}

  
  std::tuple<bool, FloatArrayF<2>, double, FloatArrayF<3>, FloatArrayF<3>, FloatArrayF<3>>
FEContactSurface :: findContactPointInElement_3d(ContactPoint *cp, ContactElement *contactElement, TimeStep *tStep)
{
    // detailed search
    return computeContactPointLocalCoordinates_3d(cp, contactElement, tStep);
}



  std::tuple<bool, FloatArrayF<1>, double, FloatArrayF<2>, FloatArrayF<2>>
FEContactSurface :: findContactPointInElement_2d(ContactPoint *cp, ContactElement *contactElement, TimeStep *tStep)
{
    // detailed search
    return computeContactPointLocalCoordinates_2d(cp, contactElement, tStep);
}
  
  
  
  std::tuple <bool, FloatArrayF<2>, double, FloatArrayF<3>,FloatArrayF<3>,FloatArrayF<3>> 
FEContactSurface :: computeContactPointLocalCoordinates_3d(ContactPoint *cp, ContactElement *contactElement, TimeStep *tStep)
{

  
  FEInterpolation3d *interpolation = dynamic_cast<FEInterpolation3d *>(contactElement->giveInterpolation());
  if (interpolation == nullptr) {
    OOFEM_ERROR("Wrong contact element");
  }
  //now apply the Newton-Raphson closest point projection
  FEIElementDeformedGeometryWrapper cellgeo(contactElement, tStep);
  FloatArray r, ro;
  FloatArrayF<3> gapVector;
  FloatArrayF<2> contactPointLocalCoords, dKsi, dFdXi;
  FloatArrayF<3> t1, t2, unitNormal;
  FloatMatrixF<2,3> dRodXi;
  FloatMatrixF<2,2> kappa, m,  G, invG;
  
  cp->giveUpdatedCoordinates(r, tStep);
  double error = 1.;
  int iter = 0;
  int maxIter = 100;
  double tol = 1.e-12, point_tol = 1.e-6;
  // test if inside
  bool inside = true;  
  //lets minimize F = 0.5 * (r-ro)*(r-ro)
  //initial guess is zero local coordinate
  do {
    std::tie(t1,t2) = interpolation->surfaceEvalBaseVectorsAt(0, contactPointLocalCoords, cellgeo);
    interpolation->local2global(ro, contactPointLocalCoords, cellgeo);

    dRodXi.setRow(t1, 0);
    dRodXi.setRow(t2, 1);

    
    gapVector = FloatArrayF<3>(r-ro);
    //dFdXi = - dRodXi * (r - ro)
    dFdXi = -dot(dRodXi, gapVector);
    error = norm(dFdXi);
    
    //get unit normal
    unitNormal = cross(t1,t2);
    //    unitNormal /= norm(unitNormal); //the orientation of the unit normal is not relevant here, since we only ever use its square
    if (error <= tol) {
      break;
    }
    
    //get curvature
    // kappa = [dro/dksidksi * normal, dro/dksideta * normal;
    //          dro/detadksi * normal, dro/detadeta * normal]
    FloatMatrix d2Ndxi2;
    interpolation->surfaceEvald2Ndxi2(d2Ndxi2, contactPointLocalCoords);
    
    FloatArray dRodXidXi, dRodEtadEta, dRodXidEta;
    for (int i = 1; i <= contactElement->giveNumberOfNodes(); ++i) {
      dRodXidXi.add(d2Ndxi2.at(i, 1), cellgeo.giveVertexCoordinates(i));
      dRodEtadEta.add(d2Ndxi2.at(i, 2), cellgeo.giveVertexCoordinates(i));
      dRodXidEta.add(d2Ndxi2.at(i, 3), cellgeo.giveVertexCoordinates(i));
    }
    kappa = {
        dRodXidXi.dotProduct(unitNormal),  dRodXidEta.dotProduct(unitNormal),
	dRodXidEta.dotProduct(unitNormal), dRodEtadEta.dotProduct(unitNormal)
    };
    //kappa *= 1./norm(unitNormal);
    //get metric tensor
    m= {
        dot(t1,t1), dot(t1,t2),
        dot(t2,t1), dot(t2,t2)
    };
    //construct G = m - gap*kappa
    //G = m - dot(gapVector, unitNormal)/norm(unitNormal) * kappa;
    G = m - dot(gapVector, unitNormal) * kappa;   
    //construct ksi increment
    //dksi = -G^-1 * dFdXi
    invG = inv(G);
    dKsi = dot(invG, dFdXi);
    //    
    contactPointLocalCoords -= dKsi;
    
    iter++;
  } while (iter < maxIter);
  
  if (iter == maxIter) {
    //OOFEM_WARNING("Closest point projection on contact surface failed to converge in %i iterations (error = %d)",iter, error);
    inside = false;
  }

    
  for (int i = 1; i <= 2; i++) {
    if (contactPointLocalCoords.at(i) < (-1. - point_tol)) {
      contactPointLocalCoords.at(i) = -1.;
      inside = false;
      //@todo:verify
      //inside = true;
    }
    else if (contactPointLocalCoords.at(i) > (1. + point_tol)) {
      contactPointLocalCoords.at(i) = 1.;
      inside = false;
      //@todo:verify
      //inside = true;
    }
  }
  //
  auto gap = dot(gapVector,unitNormal)/norm(unitNormal);
  //
  return std::make_tuple(inside, contactPointLocalCoords, gap, unitNormal, t1, t2);
}





  
  std::tuple <bool, FloatArrayF<1>, double, FloatArrayF<2>,FloatArrayF<2>> 
FEContactSurface :: computeContactPointLocalCoordinates_2d(ContactPoint *cp, ContactElement *contactElement, TimeStep *tStep)
{

  
  FEInterpolation2d *interpolation = dynamic_cast<FEInterpolation2d *>(contactElement->giveInterpolation());
  if (interpolation == nullptr) {
    OOFEM_ERROR("Wrong contact element");
  }
  //now apply the Newton-Raphson closest point projection
  FEIElementDeformedGeometryWrapper cellgeo(contactElement, tStep);
  FloatArray r, ro;
  FloatArrayF<2> gapVector;
  double dFdXi;
  FloatArrayF<1> contactPointLocalCoords;
  FloatArrayF<2> t1, normal;
  FloatMatrixF<2,3> dRodXi;
  
  cp->giveUpdatedCoordinates(r, tStep);
  double error = 1.;
  int iter = 0;
  int maxIter = 100;;
  double tol = 1.e-8, point_tol = 5.e-2;
  // test if inside
  bool inside = true;
  //lets minimize F = 0.5 * (r-ro)*(r-ro)
  //initial guess is zero local coordinate
  do {
    t1 = interpolation->surfaceEvalBaseVectorsAt(1, contactPointLocalCoords, cellgeo);
    interpolation->local2global(ro, contactPointLocalCoords, cellgeo);
    gapVector = FloatArrayF<2>(r-ro);
    //dFdXi = - dRodXi * (r - ro)
    dFdXi = -dot(t1, gapVector);
    error = fabs(dFdXi);
    
    //get normal. note: not normalized
    normal = {+t1.at(2), -t1.at(1)};
    if (error <= tol) {
      break;
    }
    
    //get curvature
    // kappa = [dro/dksidksi * normal, dro/dksideta * normal;
    //          dro/detadksi * normal, dro/detadeta * normal]
    FloatMatrix d2Ndxi2;
    interpolation->surfaceEvald2Ndxi2(d2Ndxi2, contactPointLocalCoords);
    
    FloatArray dRodXidXi;
    for (int i = 1; i <= contactElement->giveNumberOfNodes(); ++i) {
      dRodXidXi.add(d2Ndxi2.at(i, 1), cellgeo.giveVertexCoordinates(i));
    }
    auto G = dot(t1,t1) - dot(gapVector, normal) * dRodXidXi.dotProduct(normal);
    //@todo: the factor 2 here is because the element goes from -1 to 1 so the length is 2 - we have to somehow introduce it into the interpolation or change somethig as it shouldn't be there for element parametrized on the interval 0,1
    contactPointLocalCoords.at(1) -= dFdXi/G ;
    iter++;
  } while (iter < maxIter);
  
  if (iter == maxIter) {
    //OOFEM_WARNING("Closest point projection on contact surface failed to converge in %i iterations (error = %d)",iter, error);
  }
  

  if (contactPointLocalCoords.at(1) < (-1. - point_tol)) {
    contactPointLocalCoords.at(1) = -1.;
    inside = false;
  } else if (contactPointLocalCoords.at(1) > (1. + point_tol)) {
    contactPointLocalCoords.at(1) = 1.;
    inside = false;
  }
  //
  auto gap = dot(gapVector,normal)/norm(normal);
  //
  return std::make_tuple(inside, contactPointLocalCoords, gap, normal, t1);
}

  
  
  

}
