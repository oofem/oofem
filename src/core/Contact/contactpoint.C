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

#include "contactpoint.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "feinterpol2d.h"
#include "feinterpol3d.h"

namespace oofem {



FEInterpolation*
FEContactPoint:: giveInterpolation()
{
  int spatial_dimension = this->surface_dimension+1;
  if(spatial_dimension == 3) {
    auto fe= dynamic_cast<FEInterpolation3d*>(contactSurface->giveContactElement(contactElementId)->giveInterpolation());
    return fe;
  } else if(spatial_dimension == 2) {
    auto fe= dynamic_cast<FEInterpolation2d*>(contactSurface->giveContactElement(contactElementId)->giveInterpolation());
  return fe;
  } else {
    OOFEM_ERROR("Incorrect spatial dimension");
  }
}

void FEContactPoint:: computeNmatrix(FloatMatrix &answer)
{
  contactSurface->giveContactElement(contactElementId)->computeNmatrixAt(this->giveLocalCoordinates(), answer);
}


  
void
FEContactPoint :: compute_dNdxi_matrix(FloatMatrix &dNdxi, int index)
{

  FloatMatrix dN;
  this->giveInterpolation()->surfaceEvaldNdxi(dN, this->giveLocalCoordinates());
  //
  int spatial_dimension = this->surface_dimension+1;
  dNdxi.resize(spatial_dimension, spatial_dimension * dN.giveNumberOfRows());
  for (int i = 1; i <= dN.giveNumberOfRows(); i++) {
    FloatMatrix dn(spatial_dimension,spatial_dimension), dNdi;
    dn.beUnitMatrix();
    //
    dn.times(dN.at(i, index));
    //
    dNdxi.setSubMatrix(dn, 1, 1 + (i - 1) * spatial_dimension);
  }
 
}


FloatArray
FEContactPoint :: giveNormalVector()
{
  return contactSurface->giveContactElement(contactElementId)->computeNormalVectorAt(this->giveLocalCoordinates());
}

  
void  
FEContactPoint :: computeCurvature(FloatMatrix &kappa, const FloatArray &normal, TimeStep *tStep)
{
    //get curvature
    // kappa = [dro/dksidksi * normal, dro/dksideta * normal;
    //          dro/detadksi * normal, dro/detadeta * normal]
  auto ce = contactSurface->giveContactElement(contactElementId);
  FEIElementDeformedGeometryWrapper cellgeo(ce, tStep);
  auto fe = this->giveInterpolation();
  FloatMatrix d2Ndxi2;
  fe->surfaceEvald2Ndxi2(d2Ndxi2, this->giveLocalCoordinates());
  
  FloatArray d2Ndxidxi, d2Ndetadeta, d2Ndxideta;
  for (int i = 1; i <= ce->giveNumberOfNodes(); ++i) {
    d2Ndxidxi.add(d2Ndxi2.at(i, 1), cellgeo.giveVertexCoordinates(i));
    if(surface_dimension == 2) {
      d2Ndetadeta.add(d2Ndxi2.at(i, 2), cellgeo.giveVertexCoordinates(i));
      d2Ndxideta.add(d2Ndxi2.at(i, 3), cellgeo.giveVertexCoordinates(i));
    }
  }
  kappa.resize(surface_dimension,surface_dimension);
  kappa.at(1,1) = d2Ndxidxi.dotProduct(normal);
  if(surface_dimension == 2) {
    kappa.at(1,2) = d2Ndxideta.dotProduct(normal);
    kappa.at(2,1) = kappa.at(1,2);
    kappa.at(2,2) = d2Ndetadeta.dotProduct(normal);
  }
}
  

void
FEContactPoint :: computeVectorOf(ValueModeType mode, TimeStep *tStep, FloatArray &answer)
{
  contactSurface->giveContactElement(contactElementId)->computeVectorOf(mode, tStep, answer);
}

  

void
FEContactPoint :: giveUpdatedCoordinates(FloatArray &coords, TimeStep* tStep)
{
  this->giveInterpolation()->local2global(coords, this->giveLocalCoordinates(), FEIElementDeformedGeometryWrapper(contactSurface->giveContactElement(contactElementId), tStep));
}



  
bool
FEContactPoint :: giveLocationArray(IntArray &locationArray, const IntArray &dofIDArry, const UnknownNumberingScheme &s) const
{
  if(contactElementId != -1) {
    contactSurface->giveContactElement(contactElementId)->giveLocationArray(locationArray,dofIDArry, s);
    return true;
  } else {
    return false;
  }
}


void
FEContactPoint ::giveUnknownVector(FloatArray &answer, const IntArray &dofMask, ValueModeType mode, TimeStep *tStep, bool padding)
{
  this->contactSurface->giveContactElement(contactElementId)->computeVectorOf(dofMask, mode, tStep, answer);
}


FloatArray
FEContactPoint_Master :: giveGlobalCoordinates()
{
  FloatArray ret;
  this->contactSurface->giveContactElement(contactElementId)->computeGlobalCoordinates(ret, this->localCoordinates);
  return ret;
}
  

} // end namespace oofem
