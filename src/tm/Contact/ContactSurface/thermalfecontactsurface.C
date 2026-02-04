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
 *               Copyright (C) 1993 - 2023   Borek Patzak
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


#include "thermalfecontactsurface.h"
#include "classfactory.h"
#include "gausspoint.h"
#include "fei3dquadlin.h"
namespace oofem {
REGISTER_ContactSurface(ThermalFEContactSurface);

  /*
void ThermalFEContactSurface :: initializeFrom(InputRecord &ir)
{
    ContactSurface :: initializeFrom(ir);
}
  */
  /*void
ThermalFEContactSurface :: postInitialize()
{
  this->createContactElements();
}
  */

  /*
bool
ThermalFEContactSurface :: computeContactPointLocalCoordinates(FloatArray &cPointLocalCoords, Node *node, ContactElement *contactElement, TimeStep *tStep)
{

    auto interpolation = contactElement->giveInterpolation();
    //now apply the Newton-Raphson closest point projection
    FEIElementDeformedGeometryWrapper cellgeo(elem, tStep);
    FloatArray r, ro, dKsi, t1, t2, gapVector, dFdXi, unitNormal;
    FloatMatrix dRodXi, kappa, m, G, invG;
    
    node->giveUpdatedCoordinates(r, tStep);
    double error = 1.;
    int iter = 0;
    int maxIter = 100;
    double tol = 1.e-6, point_tol = 1.e-4;
    
    //lets minimize F = 0.5 * (r-ro)*(r-ro)
    ksi.resize(2); //initial guess
    do {
      interpolation->surfaceEvalBaseVectorsAt(t1, t2, isurf, ksi, cellgeo);
      interpolation->surfaceLocal2global(ro, isurf, ksi, cellgeo);
      dRodXi.resize(2, 3);
      dRodXi.addSubVectorRow(t1, 1, 1);
      dRodXi.addSubVectorRow(t2, 2, 1);
      gapVector.beDifferenceOf(r, ro);
      //dFdXi = - dRodXi * (r - ro)
      dFdXi.beProductOf(dRodXi, gapVector);
      dFdXi.times(-1.);
      error = dFdXi.computeNorm();
      if (error <= tol) {
	break;
      }
      
      //get unit normal
      unitNormal.beVectorProductOf(t1, t2);
      unitNormal.normalize(); //the orientation of the unit normal is not relevant here, since we only ever use its square
      
      //get curvature
      // kappa = [dro/dksidksi * normal, dro/dksideta * normal;
      //          dro/detadksi * normal, dro/detadeta * normal]
      FloatMatrix d2Ndxidxi;
      interpolation->surfaceEvald2Ndxidxi(d2Ndxidxi, isurf, ksi, cellgeo);
      
      IntArray nodeIndices;
      nodeIndices = interpolation->computeLocalSurfaceMapping(isurf);
      FloatArray dRodXidXi, dRodEtadEta, dRodXidEta;
      for (int i = 1; i <= nodeIndices.giveSize(); ++i) {
	dRodXidXi.add(d2Ndxidxi.at(i, 1), cellgeo.giveVertexCoordinates(nodeIndices.at(i)));
	dRodEtadEta.add(d2Ndxidxi.at(i, 2), cellgeo.giveVertexCoordinates(nodeIndices.at(i)));
	dRodXidEta.add(d2Ndxidxi.at(i, 3), cellgeo.giveVertexCoordinates(nodeIndices.at(i)));
      }
      
      kappa.resize(2, 2);
      kappa.at(1, 1) = dRodXidXi.dotProduct(unitNormal);
      kappa.at(1, 2) = dRodXidEta.dotProduct(unitNormal);
      kappa.at(2, 1) = dRodXidEta.dotProduct(unitNormal);
      kappa.at(2, 2) = dRodEtadEta.dotProduct(unitNormal);
      
      //get metric tensor
      m.resize(2, 2);
      m.at(1, 1) = t1.dotProduct(t1);
      m.at(1, 2) = t1.dotProduct(t2);
      m.at(2, 1) = t2.dotProduct(t1);
      m.at(2, 2) = t2.dotProduct(t2);
      
      //construct G = m - gap*kappa
      G = kappa;
      G.times(-1.*gapVector.dotProduct(unitNormal));
      G.add(m);
      
      //construct ksi increment
      //dksi = -G^-1 * dFdXi
      invG.beInverseOf(G);
      dKsi.beProductOf(invG, dFdXi);
      dKsi.times(-1.);
      
      ksi.add(dKsi);
      
      iter++;
    } while (iter < maxIter);
    
    if (iter == maxIter) {
      OOFEM_WARNING("Closest point projection on contact surface failed to converge in %i iterations (error = %d)",iter, error);
    }
    
    // test if inside
    bool inside = true;
    for (int i = 1; i <= 2; i++) {
      if (ksi.at(i) < (-1. - point_tol)) {
	ksi.at(i) = -1.;
	inside = false;
      }
      else if (ksi.at(i) > (1. + point_tol)) {
	ksi.at(i) = 1.;
	inside = false;
      }
    }
    
    
    
    if (!inside) {
      //outside of surface, return zeros and false
      ksi.resize(2);
      return false;
    }
    
    return true;
}
  
  return false;
}
*/
  


  
  // generalization to a point with given coordinates not just a node, so we can easily extend it to surface-to-surface case???
  // send there node, or gp, or just coords ?
  /*std::list<int, double>
ThermalFEContactSurface :: findClosestContactElement(const FloatArray &coords)
{
    int i;
    double gap;
    for(i = 0; i < this->numberOfContactElements, i++) {
      auto iter = 0;
      if(this->contactElements.at(i)->containPointInBoundingBox(coords)) {
	iter++;
	gap_new = this->contactElements.at(i)->computeGap(coords);
	if(iter == 1) {
	  gap_old = gap = gap_new;
	} else {
	  gap = (gap_new < gap_old) ? gap_new : gap_old; 
	}
      }
    }
    return {i, gap};
  
}
  */
  

} //end namespace oofem
