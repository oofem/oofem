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
#include "floatarrayf.h"

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

ContactElement *
FEContactPoint::giveContactElement()
{
    return contactSurface->giveContactElement(contactElementId);
}

ContactElement *
FEContactPoint::giveContactElementOnSurface(int elementId) const
{
    return contactSurface->giveContactElement(elementId);
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


double
FEContactPoint :: giveSurfaceMeasure(TimeStep *tStep)
{
  auto ce = contactSurface->giveContactElement(contactElementId);
  FEIElementDeformedGeometryWrapper cellgeo(ce, tStep);
  if (surface_dimension == 2) {
    auto fe = dynamic_cast<FEInterpolation3d *>(this->giveInterpolation());
    if (fe == nullptr) {
      OOFEM_ERROR("Wrong 3d contact interpolation");
    }
    auto [g1, g2] = fe->surfaceEvalBaseVectorsAt(0, this->giveLocalCoordinates(), cellgeo);
    return norm(cross(g1, g2));
  } else if (surface_dimension == 1) {
    auto fe = dynamic_cast<FEInterpolation2d *>(this->giveInterpolation());
    if (fe == nullptr) {
      OOFEM_ERROR("Wrong 2d contact interpolation");
    }
    return norm(fe->surfaceEvalBaseVectorsAt(1, this->giveLocalCoordinates(), cellgeo));
  } else {
    OOFEM_ERROR("Incorrect contact surface dimension");
  }
}

double
FEContactPoint :: giveReferenceSurfaceMeasure()
{
  auto ce = contactSurface->giveContactElement(contactElementId);
  FEIElementGeometryWrapper cellgeo(ce);
  if (surface_dimension == 2) {
    auto fe = dynamic_cast<FEInterpolation3d *>(this->giveInterpolation());
    if (fe == nullptr) {
      OOFEM_ERROR("Wrong 3d contact interpolation");
    }
    auto [g1, g2] = fe->surfaceEvalBaseVectorsAt(0, this->giveLocalCoordinates(), cellgeo);
    return norm(cross(g1, g2));
  } else if (surface_dimension == 1) {
    auto fe = dynamic_cast<FEInterpolation2d *>(this->giveInterpolation());
    if (fe == nullptr) {
      OOFEM_ERROR("Wrong 2d contact interpolation");
    }
    return norm(fe->surfaceEvalBaseVectorsAt(1, this->giveLocalCoordinates(), cellgeo));
  }
  OOFEM_ERROR("Incorrect contact surface dimension");
}

  
void
FEContactPoint :: computeSecondBaseVectors(std::vector<std::vector<FloatArray>> &answer, TimeStep *tStep)
{
  if (surface_dimension < 1 || surface_dimension > 2) {
    OOFEM_ERROR("FEContactPoint: unsupported surface dimension for second base vectors");
  }

  ContactElement *ce = contactSurface->giveContactElement(contactElementId);
  FEIElementDeformedGeometryWrapper cellgeo(ce, tStep);
  FloatMatrix d2Ndxi2;
  this->giveInterpolation()->surfaceEvald2Ndxi2(d2Ndxi2, this->giveLocalCoordinates());

  const int spatialDimension = surface_dimension + 1;
  answer.assign(surface_dimension, std::vector<FloatArray>(surface_dimension));
  for (int node = 1; node <= ce->giveNumberOfNodes(); ++node) {
    const Coordinates &fullCoordinates = cellgeo.giveVertexCoordinates(node);
    FloatArray coordinates(spatialDimension);
    for (int k = 1; k <= spatialDimension; ++k) {
      coordinates.at(k) = fullCoordinates.at(k);
    }
    answer[0][0].add(d2Ndxi2.at(node, 1), coordinates);
    if (surface_dimension == 2) {
      answer[1][1].add(d2Ndxi2.at(node, 2), coordinates);
      answer[0][1].add(d2Ndxi2.at(node, 3), coordinates);
    }
  }
  if (surface_dimension == 2) {
    answer[1][0] = answer[0][1];
  }
}

void
FEContactPoint :: computeCurvature(FloatMatrix &kappa, const FloatArray &normal, TimeStep *tStep)
{
  std::vector<std::vector<FloatArray>> secondBaseVectors;
  this->computeSecondBaseVectors(secondBaseVectors, tStep);

  kappa.resize(surface_dimension, surface_dimension);
  for (int i = 0; i < surface_dimension; ++i) {
    for (int j = 0; j < surface_dimension; ++j) {
      kappa(i, j) = secondBaseVectors[i][j].dotProduct(normal);
    }
  }
}

void
FEContactPoint :: computeVectorOf(ValueModeType mode, TimeStep *tStep, FloatArray &answer)
{
  contactSurface->giveContactElement(contactElementId)->computeVectorOf(mode, tStep, answer);
}

  

void
FEContactPoint :: giveUpdatedCoordinates(Coordinates &coords, TimeStep* tStep)
{
  this->giveInterpolation()->local2global(coords, this->giveLocalCoordinates(), FEIElementDeformedGeometryWrapper(contactSurface->giveContactElement(contactElementId), tStep));
}

void
FEContactPoint :: giveUpdatedCoordinatesOnElement(Coordinates &coords, int elementId,
                                                   const FloatArray &localCoords, TimeStep *tStep) const
{
  ContactElement *element = contactSurface->giveContactElement(elementId);
  if (element == nullptr) {
    OOFEM_ERROR("FEContactPoint: invalid previous master contact element");
  }
  element->giveInterpolation()->local2global(coords, localCoords,
                                              FEIElementDeformedGeometryWrapper(element, tStep));
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


Coordinates
FEContactPoint_Master :: giveGlobalCoordinates()
{
  Coordinates ret;
  this->contactSurface->giveContactElement(contactElementId)->computeGlobalCoordinates(ret, this->localCoordinates);
  return ret;
}
  

} // end namespace oofem
