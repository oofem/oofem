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
#include "node.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace oofem {

namespace {

constexpr double cppRelativeTolerance = 1.e-10;
constexpr double cppAlgebraTolerance = 64.0 * std::numeric_limits<double>::epsilon();

bool isFinite(const FloatArrayF<1> &value)
{
  return std::isfinite(value.at(1));
}

bool isFinite(const FloatArrayF<2> &value)
{
  return std::isfinite(value.at(1)) && std::isfinite(value.at(2));
}

bool isFinite(const FloatArrayF<3> &value)
{
  return std::isfinite(value.at(1)) && std::isfinite(value.at(2))
      && std::isfinite(value.at(3));
}

bool isPositiveDefinite2x2(double h11, double h12, double h22)
{
  if (!std::isfinite(h11) || !std::isfinite(h12) || !std::isfinite(h22)) {
    return false;
  }

  const double scale = std::max({std::abs(h11), std::abs(h12), std::abs(h22)});
  if (!(scale > 0.0) || !std::isfinite(scale)) {
    return false;
  }

  const double determinant = h11 * h22 - h12 * h12;
  return h11 > cppAlgebraTolerance * scale
      && h22 > cppAlgebraTolerance * scale
      && determinant > cppAlgebraTolerance * scale * scale;
}

double stationarityScale2d(double tangentSquared, double gapNorm)
{
  return std::max(tangentSquared, std::sqrt(tangentSquared) * gapNorm);
}

double stationarityScale3d(double a11, double a22, double gapNorm)
{
  const double metricScale = std::max(a11, a22);
  return std::max(metricScale, std::sqrt(metricScale) * gapNorm);
}

std::vector<FloatArray> giveFeatureLocalVertices(ContactElement *element)
{
  if (element->giveGeometryType() == EGT_quad_1_interface) {
    return {
      Vec2(-1.0, -1.0), Vec2(1.0, -1.0),
      Vec2(1.0, 1.0), Vec2(-1.0, 1.0)
    };
  }
  if (element->giveGeometryType() == EGT_triangle_1) {
    return {
      Vec2(0.0, 0.0), Vec2(1.0, 0.0),
      Vec2(0.0, 1.0)
    };
  }
  if (element->giveGeometryType() == EGT_line_1) {
    return { Vec1(-1.0), Vec1(1.0) };
  }
  return {};
}

std::pair<int, int> giveLocalEdgeVertices(ContactElement *element, int edgeIndex)
{
  const int nVertices = static_cast<int>(giveFeatureLocalVertices(element).size());
  if ((nVertices != 3 && nVertices != 4) || edgeIndex < 1 || edgeIndex > nVertices) {
    return { 0, 0 };
  }
  return { edgeIndex, edgeIndex == nVertices ? 1 : edgeIndex + 1 };
}

FloatArrayF<3> giveUpdatedNodeCoordinates3d(Node *node, TimeStep *tStep)
{
  if (tStep == nullptr) {
    return {
      node->giveCoordinate(1),
      node->giveCoordinate(2),
      node->giveCoordinate(3)
    };
  }
  return {
    node->giveUpdatedCoordinate(1, tStep),
    node->giveUpdatedCoordinate(2, tStep),
    node->giveUpdatedCoordinate(3, tStep)
  };
}

FloatArrayF<2> giveUpdatedNodeCoordinates2d(Node *node, TimeStep *tStep)
{
  if (tStep == nullptr) {
    return {
      node->giveCoordinate(1),
      node->giveCoordinate(2)
    };
  }
  return {
    node->giveUpdatedCoordinate(1, tStep),
    node->giveUpdatedCoordinate(2, tStep)
  };
}

} // namespace


FEContactSurface :: FEContactSurface(int n, Domain *d) : ContactSurface(n, d)
{
}



void FEContactSurface :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    ContactSurface :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, this->contactElementSetNumber, _IFT_FEContactSurface_contactElementSetNumber);
}


void
FEContactSurface :: postInitialize()
{
  this->contactElementSet = domain->giveSet(this->contactElementSetNumber)->giveElementList();
}

void
FEContactSurface :: saveContext(DataStream &stream, ContextMode mode)
{
  ContactSurface :: saveContext(stream, mode);
  if ((mode & CM_Definition) && !stream.write(contactElementSetNumber)) {
    THROW_CIOERR(CIO_IOERR);
  }
}

void
FEContactSurface :: restoreContext(DataStream &stream, ContextMode mode)
{
  ContactSurface :: restoreContext(stream, mode);
  if (mode & CM_Definition) {
    if (!stream.read(contactElementSetNumber)) {
      THROW_CIOERR(CIO_IOERR);
    }
    this->contactElementSet = domain->giveSet(contactElementSetNumber)->giveElementList();
  }
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


  std::tuple<bool, FloatArrayF<2>, double, double, FloatArrayF<3>, FloatArrayF<3>, FloatArrayF<3>>
FEContactSurface :: findContactPointInElement_3d(ContactPoint *cp, ContactElement *contactElement, TimeStep *tStep)
{
    // detailed search
    return computeContactPointLocalCoordinates_3d(cp, contactElement, tStep);
}


FloatArrayF<3>
FEContactSurface :: computeContactPointNormal_3d(ContactPoint *cp, TimeStep *tStep) const
{
  ContactElement *element = cp ? cp->giveContactElement() : nullptr;
  auto *interpolation = element
    ? dynamic_cast<FEInterpolation3d *>(element->giveInterpolation()) : nullptr;
  if (interpolation == nullptr) {
    return {};
  }

  FEIElementDeformedGeometryWrapper geometry(element, tStep);
  const auto [tangent1, tangent2] = interpolation->surfaceEvalBaseVectorsAt(
    0, cp->giveLocalCoordinates(), geometry);
  FloatArrayF<3> normal = cross(tangent1, tangent2);
  const double normalNorm = norm(normal);
  if (!isFinite(normal) || !std::isfinite(normalNorm) || !(normalNorm > 0.0)) {
    return {};
  }
  return normal / normalNorm;
}


std::tuple<bool, FloatArrayF<2>, double, double, FloatArrayF<3>, FloatArrayF<3>, FloatArrayF<3>>
FEContactSurface :: findContactPointAlongDirectionInElement_3d(
    ContactPoint *cp, ContactElement *contactElement,
    const FloatArrayF<3> &inputDirection, TimeStep *tStep)
{
  auto *interpolation = dynamic_cast<FEInterpolation3d *>(contactElement->giveInterpolation());
  FloatArrayF<2> localCoordinates;
  const auto invalidProjection = [&localCoordinates]() {
    const double infinity = std::numeric_limits<double>::infinity();
    return std::make_tuple(false, localCoordinates, infinity, infinity,
                           FloatArrayF<3> {}, FloatArrayF<3> {}, FloatArrayF<3> {});
  };
  if (interpolation == nullptr) {
    return invalidProjection();
  }

  const auto geometryType = contactElement->giveGeometryType();
  const bool isTriangle = geometryType == EGT_triangle_1;
  if (!isTriangle && geometryType != EGT_quad_1_interface) {
    return invalidProjection();
  }

  const double directionNorm = norm(inputDirection);
  if (!isFinite(inputDirection) || !std::isfinite(directionNorm)
      || !(directionNorm > 0.0)) {
    return invalidProjection();
  }
  const FloatArrayF<3> direction = inputDirection / directionNorm;

  Coordinates slaveCoordinatesDynamic;
  cp->giveUpdatedCoordinates(slaveCoordinatesDynamic, tStep);
  if (slaveCoordinatesDynamic.giveSize() != 3) {
    return invalidProjection();
  }
  const FloatArrayF<3> slaveCoordinates(slaveCoordinatesDynamic);
  FEIElementDeformedGeometryWrapper geometry(contactElement, tStep);

  double rayCoordinate = 0.0;
  Coordinates masterCoordinatesDynamic;
  FloatArrayF<3> tangent1, tangent2, residual;
  bool converged = false;
  constexpr int maxIterations = 50;
  for (int iteration = 0; iteration < maxIterations; ++iteration) {
    std::tie(tangent1, tangent2) = interpolation->surfaceEvalBaseVectorsAt(
      0, localCoordinates, geometry);
    interpolation->local2global(masterCoordinatesDynamic, localCoordinates, geometry);
    if (masterCoordinatesDynamic.giveSize() != 3) {
      return invalidProjection();
    }
    const FloatArrayF<3> masterCoordinates(masterCoordinatesDynamic);
    residual = masterCoordinates + rayCoordinate * direction - slaveCoordinates;
    const double residualNorm = norm(residual);
    const double lengthScale = std::max({1.0, norm(masterCoordinates), norm(slaveCoordinates),
                                         std::abs(rayCoordinate)});
    if (!isFinite(tangent1) || !isFinite(tangent2) || !isFinite(residual)
        || !std::isfinite(residualNorm)) {
      return invalidProjection();
    }
    if (residualNorm <= cppRelativeTolerance * lengthScale) {
      converged = true;
      break;
    }

    FloatMatrixF<3,3> jacobian;
    jacobian.setColumn(tangent1, 0);
    jacobian.setColumn(tangent2, 1);
    jacobian.setColumn(direction, 2);
    const double jacobianDeterminant = det(jacobian);
    const double jacobianScale = std::max(1.0, norm(tangent1) * norm(tangent2));
    if (!std::isfinite(jacobianDeterminant)
        || std::abs(jacobianDeterminant) <= cppAlgebraTolerance * jacobianScale) {
      return invalidProjection();
    }
    const FloatArrayF<3> increment = dot(inv_33(jacobian, 0.0), residual);
    if (!isFinite(increment)) {
      return invalidProjection();
    }
    localCoordinates.at(1) -= increment.at(1);
    localCoordinates.at(2) -= increment.at(2);
    rayCoordinate -= increment.at(3);
  }
  if (!converged) {
    return invalidProjection();
  }

  double xi = localCoordinates.at(1);
  double eta = localCoordinates.at(2);
  if (isTriangle) {
    if (xi < -domainTolerance || eta < -domainTolerance
        || xi + eta > 1.0 + domainTolerance) {
      return invalidProjection();
    }
  } else {
    if (xi < -1.0 - domainTolerance || xi > 1.0 + domainTolerance
        || eta < -1.0 - domainTolerance || eta > 1.0 + domainTolerance) {
      return invalidProjection();
    }
  }
  // domainTolerance deliberately extends the admissible facet parameter
  // domain. Keep the converged ray intersection in that narrow extension:
  // clamping it back to the topological edge moves the master point off the
  // ray, so the transverse-gap validation below rejects every genuinely
  // exterior intersection and reduces searchtol to a roundoff tolerance.
  localCoordinates = Vec2(xi, eta);

  std::tie(tangent1, tangent2) = interpolation->surfaceEvalBaseVectorsAt(
    0, localCoordinates, geometry);
  interpolation->local2global(masterCoordinatesDynamic, localCoordinates, geometry);
  const FloatArrayF<3> masterCoordinates(masterCoordinatesDynamic);
  const FloatArrayF<3> gapVector = slaveCoordinates - masterCoordinates;
  const double rayGap = dot(gapVector, direction);
  const FloatArrayF<3> transverseGap = gapVector - rayGap * direction;
  const double validationScale = std::max({1.0, norm(masterCoordinates), norm(slaveCoordinates)});
  if (!isFinite(tangent1) || !isFinite(tangent2) || !isFinite(gapVector)
      || norm(transverseGap) > cppRelativeTolerance * validationScale
      || !std::isfinite(rayGap)) {
    return invalidProjection();
  }

  // The slave normal defines only the projection ray.  OOFEM's normal-gap
  // convention is positive in separation and negative in penetration, so the
  // contact normal must face opposite the slave ray (the facing master-side
  // orientation).  Returning the ray direction itself reverses the active set:
  // it activates across a free gap and releases after actual interpenetration.
  FloatArrayF<3> contactNormal = cross(tangent1, tangent2);
  const double contactNormalNorm = norm(contactNormal);
  if (!isFinite(contactNormal) || !std::isfinite(contactNormalNorm)
      || !(contactNormalNorm > 0.0)) {
    return invalidProjection();
  }
  contactNormal /= contactNormalNorm;
  if (dot(contactNormal, direction) > 0.0) {
    contactNormal *= -1.0;
  }
  const double gap = dot(gapVector, contactNormal);
  if (!std::isfinite(gap)) {
    return invalidProjection();
  }

  return std::make_tuple(true, localCoordinates, gap, dot(gapVector, gapVector),
                         contactNormal, tangent1, tangent2);
}



  std::tuple<bool, FloatArrayF<1>, double, double, FloatArrayF<2>, FloatArrayF<2>>
FEContactSurface :: findContactPointInElement_2d(ContactPoint *cp, ContactElement *contactElement, TimeStep *tStep)
{
    // detailed search
    return computeContactPointLocalCoordinates_2d(cp, contactElement, tStep);
}



  std::tuple <bool, FloatArrayF<2>, double, double, FloatArrayF<3>,FloatArrayF<3>,FloatArrayF<3>>
FEContactSurface :: computeContactPointLocalCoordinates_3d(ContactPoint *cp, ContactElement *contactElement, TimeStep *tStep)
{
  FEInterpolation3d *interpolation = dynamic_cast<FEInterpolation3d *>(contactElement->giveInterpolation());
  FloatArrayF<2> contactPointLocalCoords;
  const auto invalidProjection = [&contactPointLocalCoords]() {
    const double infinity = std::numeric_limits<double>::infinity();
    return std::make_tuple(false, contactPointLocalCoords, infinity, infinity,
                           FloatArrayF<3> {}, FloatArrayF<3> {}, FloatArrayF<3> {});
  };
  if (interpolation == nullptr) {
    return invalidProjection();
  }

  const auto geometryType = contactElement->giveGeometryType();
  const bool isTriangle = geometryType == EGT_triangle_1;
  if (!isTriangle && geometryType != EGT_quad_1_interface) {
    return invalidProjection();
  }

  // Apply Newton's method to the unconstrained smooth-surface CPP. A root
  // outside the element domain is not clamped to a lower-dimensional feature:
  // edge and corner contact require different kinematics and tangents.
  FEIElementDeformedGeometryWrapper cellgeo(contactElement, tStep);
  Coordinates r, ro;
  FloatArrayF<3> gapVector;
  FloatArrayF<2> dKsi, dFdXi;
  FloatArrayF<3> t1, t2, unitNormal;
  FloatMatrixF<2,3> dRodXi;
  FloatMatrixF<2,2> G, invG;

  cp->giveUpdatedCoordinates(r, tStep);
  double error = std::numeric_limits<double>::infinity();
  int iter = 0;
  constexpr int maxIter = 100;
  bool converged = false;

  auto evaluateHessian = [&](const FloatArrayF<3> &currentGap,
                             double &h11, double &h12, double &h22) {
    FloatMatrix d2Ndxi2;
    interpolation->surfaceEvald2Ndxi2(d2Ndxi2, contactPointLocalCoords);

    FloatArray dRodXidXi, dRodEtadEta, dRodXidEta;
    for (int i = 1; i <= contactElement->giveNumberOfNodes(); ++i) {
      dRodXidXi.add(d2Ndxi2.at(i, 1), cellgeo.giveVertexCoordinates(i));
      dRodEtadEta.add(d2Ndxi2.at(i, 2), cellgeo.giveVertexCoordinates(i));
      dRodXidEta.add(d2Ndxi2.at(i, 3), cellgeo.giveVertexCoordinates(i));
    }

    h11 = dot(t1, t1) - dRodXidXi.dotProduct(currentGap);
    h12 = dot(t1, t2) - dRodXidEta.dotProduct(currentGap);
    h22 = dot(t2, t2) - dRodEtadEta.dotProduct(currentGap);
    return isPositiveDefinite2x2(h11, h12, h22);
  };

  while (iter < maxIter) {
    std::tie(t1,t2) = interpolation->surfaceEvalBaseVectorsAt(0, contactPointLocalCoords, cellgeo);
    interpolation->local2global(ro, contactPointLocalCoords, cellgeo);

    if (!isFinite(contactPointLocalCoords) || !isFinite(t1) || !isFinite(t2)) {
      return invalidProjection();
    }

    dRodXi.setRow(t1, 0);
    dRodXi.setRow(t2, 1);
    gapVector = FloatArrayF<3>(r-ro);
    if (!isFinite(gapVector)) {
      return invalidProjection();
    }

    const double a11 = dot(t1, t1);
    const double a12 = dot(t1, t2);
    const double a22 = dot(t2, t2);
    if (!isPositiveDefinite2x2(a11, a12, a22)) {
      return invalidProjection();
    }

    // F_,i = -rho_i . (r-rho), for F = 0.5 |r-rho|^2.
    dFdXi = -dot(dRodXi, gapVector);
    error = norm(dFdXi);
    const double gapNorm = norm(gapVector);
    const double errorScale = stationarityScale3d(a11, a22, gapNorm);
    if (!std::isfinite(error) || !std::isfinite(errorScale) || !(errorScale > 0.0)) {
      return invalidProjection();
    }

    unitNormal = cross(t1,t2);
    if (!isFinite(unitNormal)) {
      return invalidProjection();
    }
    if (error <= cppRelativeTolerance * errorScale) {
      converged = true;
      break;
    }

    double h11, h12, h22;
    if (!evaluateHessian(gapVector, h11, h12, h22)) {
      return invalidProjection();
    }
    G = { h11, h12, h12, h22 };
    invG = inv(G);
    dKsi = dot(invG, dFdXi);
    if (!isFinite(dKsi)) {
      return invalidProjection();
    }
    contactPointLocalCoords -= dKsi;
    iter++;
  }

  if (!converged) {
    return invalidProjection();
  }

  // Admit a stationary root on the closed element boundary as a one-sided
  // smooth-surface limit. Only roundoff-sized coordinate drift is snapped;
  // a genuinely exterior root is rejected rather than constrained onto an
  // edge or corner.
  double xi = contactPointLocalCoords.at(1);
  double eta = contactPointLocalCoords.at(2);
  if (isTriangle) {
    if (xi < -domainTolerance || eta < -domainTolerance
        || xi + eta > 1.0 + domainTolerance) {
      return invalidProjection();
    }
    xi = std::max(0.0, xi);
    eta = std::max(0.0, eta);
    const double sum = xi + eta;
    if (sum > 1.0) {
      xi /= sum;
      eta /= sum;
    }
  } else {
    if (xi < -1.0 - domainTolerance || xi > 1.0 + domainTolerance
        || eta < -1.0 - domainTolerance || eta > 1.0 + domainTolerance) {
      return invalidProjection();
    }
    xi = std::clamp(xi, -1.0, 1.0);
    eta = std::clamp(eta, -1.0, 1.0);
  }
  contactPointLocalCoords.at(1) = xi;
  contactPointLocalCoords.at(2) = eta;

  // Re-evaluate after the possible roundoff snap. This also validates a root
  // that converged before Newton needed to form a Hessian.
  std::tie(t1,t2) = interpolation->surfaceEvalBaseVectorsAt(0, contactPointLocalCoords, cellgeo);
  interpolation->local2global(ro, contactPointLocalCoords, cellgeo);
  gapVector = FloatArrayF<3>(r-ro);
  if (!isFinite(contactPointLocalCoords) || !isFinite(t1) || !isFinite(t2)
      || !isFinite(gapVector)) {
    return invalidProjection();
  }

  const double a11 = dot(t1, t1);
  const double a12 = dot(t1, t2);
  const double a22 = dot(t2, t2);
  if (!isPositiveDefinite2x2(a11, a12, a22)) {
    return invalidProjection();
  }
  dRodXi.setRow(t1, 0);
  dRodXi.setRow(t2, 1);
  dFdXi = -dot(dRodXi, gapVector);
  const double gapNorm = norm(gapVector);
  const double finalErrorScale = stationarityScale3d(a11, a22, gapNorm);
  if (!std::isfinite(finalErrorScale) || !(finalErrorScale > 0.0)
      || norm(dFdXi) > cppRelativeTolerance * finalErrorScale) {
    return invalidProjection();
  }

  double h11, h12, h22;
  if (!evaluateHessian(gapVector, h11, h12, h22)) {
    return invalidProjection();
  }

  unitNormal = cross(t1,t2);
  const double unitNormalNorm = norm(unitNormal);
  if (!isFinite(unitNormal) || !std::isfinite(unitNormalNorm) || !(unitNormalNorm > 0.0)) {
    return invalidProjection();
  }
  const double gap = dot(gapVector, unitNormal) / unitNormalNorm;
  const double distanceSquared = dot(gapVector, gapVector);
  if (!std::isfinite(gap) || !std::isfinite(distanceSquared)) {
    return invalidProjection();
  }
  return std::make_tuple(true, contactPointLocalCoords, gap, distanceSquared, unitNormal, t1, t2);
}






  std::tuple <bool, FloatArrayF<1>, double, double, FloatArrayF<2>,FloatArrayF<2>>
FEContactSurface :: computeContactPointLocalCoordinates_2d(ContactPoint *cp, ContactElement *contactElement, TimeStep *tStep)
{
  FEInterpolation2d *interpolation = dynamic_cast<FEInterpolation2d *>(contactElement->giveInterpolation());
  FloatArrayF<1> contactPointLocalCoords;
  const auto invalidProjection = [&contactPointLocalCoords]() {
    const double infinity = std::numeric_limits<double>::infinity();
    return std::make_tuple(false, contactPointLocalCoords, infinity, infinity,
                           FloatArrayF<2> {}, FloatArrayF<2> {});
  };
  if (interpolation == nullptr) {
    return invalidProjection();
  }

  if (contactElement->giveGeometryType() != EGT_line_1) {
    return invalidProjection();
  }

  // Apply Newton's method to the unconstrained smooth-curve CPP. An exterior
  // root is not clamped to an endpoint, which would be a separate point-contact
  // problem with different kinematics.
  FEIElementDeformedGeometryWrapper cellgeo(contactElement, tStep);
  Coordinates r, ro;
  FloatArrayF<2> gapVector;
  double dFdXi;
  FloatArrayF<2> t1, normal;

  cp->giveUpdatedCoordinates(r, tStep);
  double error = std::numeric_limits<double>::infinity();
  int iter = 0;
  constexpr int maxIter = 100;
  bool converged = false;

  auto evaluateHessian = [&](const FloatArrayF<2> &currentGap,
                             double tangentSquared, double &hessian) {
    FloatMatrix d2Ndxi2;
    interpolation->surfaceEvald2Ndxi2(d2Ndxi2, contactPointLocalCoords);

    FloatArray dRodXidXi;
    for (int i = 1; i <= contactElement->giveNumberOfNodes(); ++i) {
      dRodXidXi.add(d2Ndxi2.at(i, 1), FloatArrayF<2>(cellgeo.giveVertexCoordinates(i).at(1), cellgeo.giveVertexCoordinates(i).at(2)));
    }
    hessian = tangentSquared - dRodXidXi.dotProduct(currentGap);
    const double scale = std::max(tangentSquared, std::abs(hessian));
    return std::isfinite(hessian) && std::isfinite(scale) && scale > 0.0
        && hessian > cppAlgebraTolerance * scale;
  };

  while (iter < maxIter) {
    t1 = interpolation->surfaceEvalBaseVectorsAt(1, contactPointLocalCoords, cellgeo);
    interpolation->local2global(ro, contactPointLocalCoords, cellgeo);
    gapVector = FloatArrayF<2>(r[0]-ro[0], r[1]-ro[1]);
    if (!isFinite(contactPointLocalCoords) || !isFinite(t1) || !isFinite(gapVector)) {
      return invalidProjection();
    }

    const double tangentSquared = dot(t1, t1);
    if (!std::isfinite(tangentSquared) || !(tangentSquared > 0.0)) {
      return invalidProjection();
    }

    dFdXi = -dot(t1, gapVector);
    error = fabs(dFdXi);
    const double gapNorm = norm(gapVector);
    const double errorScale = stationarityScale2d(tangentSquared, gapNorm);
    if (!std::isfinite(error) || !std::isfinite(errorScale) || !(errorScale > 0.0)) {
      return invalidProjection();
    }

    normal = {+t1.at(2), -t1.at(1)};
    if (!isFinite(normal)) {
      return invalidProjection();
    }
    if (error <= cppRelativeTolerance * errorScale) {
      converged = true;
      break;
    }

    double hessian;
    if (!evaluateHessian(gapVector, tangentSquared, hessian)) {
      return invalidProjection();
    }
    const double increment = dFdXi / hessian;
    if (!std::isfinite(increment)) {
      return invalidProjection();
    }
    contactPointLocalCoords.at(1) -= increment;
    iter++;
  }

  if (!converged) {
    return invalidProjection();
  }

  // A stationary endpoint root is a valid one-sided curve limit. A genuinely
  // exterior root is rejected instead of being converted into point contact.
  const double unconstrainedCoordinate = contactPointLocalCoords.at(1);
  if (unconstrainedCoordinate < -1.0 - domainTolerance
      || unconstrainedCoordinate > 1.0 + domainTolerance) {
    return invalidProjection();
  }
  contactPointLocalCoords.at(1) = std::clamp(unconstrainedCoordinate, -1.0, 1.0);

  t1 = interpolation->surfaceEvalBaseVectorsAt(1, contactPointLocalCoords, cellgeo);
  interpolation->local2global(ro, contactPointLocalCoords, cellgeo);
  gapVector = FloatArrayF<2>(r[0]-ro[0], r[1]-ro[1]);
  if (!isFinite(contactPointLocalCoords) || !isFinite(t1) || !isFinite(gapVector)) {
    return invalidProjection();
  }

  const double tangentSquared = dot(t1, t1);
  if (!std::isfinite(tangentSquared) || !(tangentSquared > 0.0)) {
    return invalidProjection();
  }
  dFdXi = -dot(t1, gapVector);
  const double gapNorm = norm(gapVector);
  const double finalErrorScale = stationarityScale2d(tangentSquared, gapNorm);
  if (!std::isfinite(finalErrorScale) || !(finalErrorScale > 0.0)
      || std::abs(dFdXi) > cppRelativeTolerance * finalErrorScale) {
    return invalidProjection();
  }

  double hessian;
  if (!evaluateHessian(gapVector, tangentSquared, hessian)) {
    return invalidProjection();
  }

  normal = {+t1.at(2), -t1.at(1)};
  const double normalNorm = norm(normal);
  if (!isFinite(normal) || !std::isfinite(normalNorm) || !(normalNorm > 0.0)) {
    return invalidProjection();
  }
  const double gap = dot(gapVector, normal) / normalNorm;
  const double distanceSquared = dot(gapVector, gapVector);
  if (!std::isfinite(gap) || !std::isfinite(distanceSquared)) {
    return invalidProjection();
  }
  return std::make_tuple(true, contactPointLocalCoords, gap, distanceSquared, normal, t1);
}

int
FEContactSurface :: giveNumberOfEdges(ContactElement *element) const
{
  const auto geometryType = element->giveGeometryType();
  if (geometryType == EGT_triangle_1) {
    return 3;
  }
  if (geometryType == EGT_quad_1_interface) {
    return 4;
  }
  return 0;
}

int
FEContactSurface :: giveNumberOfVertices(ContactElement *element) const
{
  const auto geometryType = element->giveGeometryType();
  if (geometryType == EGT_line_1) {
    return 2;
  }
  return this->giveNumberOfEdges(element);
}

FloatArray
FEContactSurface :: computeIncidentPseudoNormal3d(
    ContactElement *representative, int firstNode, int secondNode,
    TimeStep *tStep) const
{
  const auto evaluateNormal = [tStep](ContactElement *element) {
    FloatArray localCoordinates;
    if (element->giveGeometryType() == EGT_triangle_1) {
      localCoordinates = Vec2(1.0 / 3.0, 1.0 / 3.0);
    } else {
      localCoordinates = Vec2(0.0, 0.0);
    }
    auto *interpolation = dynamic_cast<FEInterpolation3d *>(element->giveInterpolation());
    if (interpolation == nullptr) {
      return FloatArrayF<3> {};
    }
    auto [tangent1, tangent2] = interpolation->surfaceEvalBaseVectorsAt(
      0, localCoordinates, FEIElementDeformedGeometryWrapper(element, tStep));
    FloatArrayF<3> normal = cross(tangent1, tangent2);
    const double normalNorm = norm(normal);
    if (!isFinite(normal) || !(normalNorm > 0.0)) {
      return FloatArrayF<3> {};
    }
    return normal / normalNorm;
  };

  FloatArrayF<3> referenceNormal = evaluateNormal(representative);
  if (!(norm(referenceNormal) > 0.0)) {
    return {};
  }

  FloatArrayF<3> pseudoNormal;
  for (int i = 1; i <= this->giveNumberOfContactElements(); ++i) {
    ContactElement *element = const_cast<FEContactSurface *>(this)->giveContactElement_InSet(i);
    bool containsFirst = false;
    bool containsSecond = secondNode < 0;
    for (int node = 1; node <= element->giveNumberOfNodes(); ++node) {
      const int nodeNumber = element->giveNode(node)->giveNumber();
      containsFirst = containsFirst || nodeNumber == firstNode;
      containsSecond = containsSecond || nodeNumber == secondNode;
    }
    if (!containsFirst || !containsSecond) {
      continue;
    }

    FloatArrayF<3> incidentNormal = evaluateNormal(element);
    if (!(norm(incidentNormal) > 0.0)) {
      continue;
    }
    if (dot(incidentNormal, referenceNormal) < 0.0) {
      incidentNormal *= -1.0;
    }
    pseudoNormal += incidentNormal;
  }

  const double pseudoNorm = norm(pseudoNormal);
  if (!(pseudoNorm > 0.0)) {
    pseudoNormal = referenceNormal;
  } else {
    pseudoNormal /= pseudoNorm;
  }
  return pseudoNormal;
}

FloatArray
FEContactSurface :: computeIncidentPseudoNormal2d(
    ContactElement *representative, int vertexNode, TimeStep *tStep) const
{
  const auto evaluateNormal = [tStep](ContactElement *element) {
    auto *interpolation = dynamic_cast<FEInterpolation2d *>(element->giveInterpolation());
    if (interpolation == nullptr) {
      return FloatArrayF<2> {};
    }
    const FloatArray localCoordinates = Vec1(0.0);
    FloatArrayF<2> tangent = interpolation->surfaceEvalBaseVectorsAt(
      1, localCoordinates, FEIElementDeformedGeometryWrapper(element, tStep));
    FloatArrayF<2> normal = { tangent.at(2), -tangent.at(1) };
    const double normalNorm = norm(normal);
    if (!isFinite(normal) || !(normalNorm > 0.0)) {
      return FloatArrayF<2> {};
    }
    return normal / normalNorm;
  };

  FloatArrayF<2> referenceNormal = evaluateNormal(representative);
  if (!(norm(referenceNormal) > 0.0)) {
    return {};
  }

  FloatArrayF<2> pseudoNormal;
  for (int i = 1; i <= this->giveNumberOfContactElements(); ++i) {
    ContactElement *element = const_cast<FEContactSurface *>(this)->giveContactElement_InSet(i);
    bool containsVertex = false;
    for (int node = 1; node <= element->giveNumberOfNodes(); ++node) {
      containsVertex = containsVertex
        || element->giveNode(node)->giveNumber() == vertexNode;
    }
    if (!containsVertex) {
      continue;
    }

    FloatArrayF<2> incidentNormal = evaluateNormal(element);
    if (!(norm(incidentNormal) > 0.0)) {
      continue;
    }
    if (dot(incidentNormal, referenceNormal) < 0.0) {
      incidentNormal *= -1.0;
    }
    pseudoNormal += incidentNormal;
  }

  const double pseudoNorm = norm(pseudoNormal);
  if (!(pseudoNorm > 0.0)) {
    pseudoNormal = referenceNormal;
  } else {
    pseudoNormal /= pseudoNorm;
  }
  return pseudoNormal;
}

ContactProjection
FEContactSurface :: findContactPointOnEdge_3d(
    ContactPoint *contactPoint, ContactElement *element, int edgeIndex,
    TimeStep *tStep)
{
  ContactProjection projection;
  projection.featureType = ContactFeatureType::Edge;
  projection.elementId = element->giveNumber();
  projection.featureIndex = edgeIndex;

  const auto localVertices = giveFeatureLocalVertices(element);
  const auto [firstVertex, secondVertex] =
    giveLocalEdgeVertices(element, edgeIndex);
  if (firstVertex == 0 || secondVertex == 0) {
    return projection;
  }

  const FloatArrayF<3> first = giveUpdatedNodeCoordinates3d(
    element->giveNode(firstVertex), tStep);
  const FloatArrayF<3> second = giveUpdatedNodeCoordinates3d(
    element->giveNode(secondVertex), tStep);
  const FloatArrayF<3> edge = second - first;
  const double edgeSquared = dot(edge, edge);
  if (!(edgeSquared > 0.0) || !std::isfinite(edgeSquared)) {
    return projection;
  }

  Coordinates slaveCoordinates;
  contactPoint->giveUpdatedCoordinates(slaveCoordinates, tStep);
  const FloatArrayF<3> slave(slaveCoordinates);
  double edgeCoordinate = dot(slave - first, edge) / edgeSquared;
  if (!std::isfinite(edgeCoordinate)
      || edgeCoordinate < -domainTolerance
      || edgeCoordinate > 1.0 + domainTolerance) {
    return projection;
  }
  edgeCoordinate = std::clamp(edgeCoordinate, 0.0, 1.0);
  // The open edge interior is the one-dimensional feature domain. Endpoints
  // belong to the vertex projection so a branch-frozen edge Jacobian never
  // straddles an edge/vertex kink.
  if (edgeCoordinate <= domainTolerance
      || edgeCoordinate >= 1.0 - domainTolerance) {
    return projection;
  }

  const FloatArrayF<3> closest = first + edgeCoordinate * edge;
  const FloatArrayF<3> gapVector = slave - closest;
  const double distanceSquared = dot(gapVector, gapVector);
  if (!std::isfinite(distanceSquared)) {
    return projection;
  }

  const FloatArray pseudoNormalArray = this->computeIncidentPseudoNormal3d(
    element, element->giveNode(firstVertex)->giveNumber(),
    element->giveNode(secondVertex)->giveNumber(), tStep);
  if (pseudoNormalArray.giveSize() != 3) {
    return projection;
  }
  const FloatArrayF<3> pseudoNormal(pseudoNormalArray);
  const double orientation =
    dot(gapVector, pseudoNormal) < 0.0 ? -1.0 : 1.0;
  const double distance = std::sqrt(std::max(0.0, distanceSquared));
  FloatArrayF<3> normal = pseudoNormal;
  if (distance > cppAlgebraTolerance) {
    normal = orientation * gapVector / distance;
  }
  FloatArrayF<3> secondTangent = cross(normal, edge);
  if (!(norm(secondTangent) > 0.0)) {
    return projection;
  }

  projection.localCoordinates =
    (1.0 - edgeCoordinate) * localVertices[firstVertex - 1]
    + edgeCoordinate * localVertices[secondVertex - 1];
  projection.gap = orientation * distance;
  projection.distanceSquared = distanceSquared;
  projection.normal = normal;
  projection.tangent1 = edge;
  projection.tangent2 = secondTangent;
  projection.valid = true;
  return projection;
}

ContactProjection
FEContactSurface :: findContactPointOnVertex_3d(
    ContactPoint *contactPoint, ContactElement *element, int vertexIndex,
    TimeStep *tStep)
{
  ContactProjection projection;
  projection.featureType = ContactFeatureType::Vertex;
  projection.elementId = element->giveNumber();
  projection.featureIndex = vertexIndex;

  const auto localVertices = giveFeatureLocalVertices(element);
  if (vertexIndex < 1 || vertexIndex > static_cast<int>(localVertices.size())) {
    return projection;
  }

  const FloatArrayF<3> vertex = giveUpdatedNodeCoordinates3d(
    element->giveNode(vertexIndex), tStep);
  Coordinates slaveCoordinates;
  contactPoint->giveUpdatedCoordinates(slaveCoordinates, tStep);
  const FloatArrayF<3> slave(slaveCoordinates);
  const FloatArrayF<3> gapVector = slave - vertex;
  const double distanceSquared = dot(gapVector, gapVector);
  if (!std::isfinite(distanceSquared)) {
    return projection;
  }

  const FloatArray pseudoNormalArray = this->computeIncidentPseudoNormal3d(
    element, element->giveNode(vertexIndex)->giveNumber(), -1, tStep);
  if (pseudoNormalArray.giveSize() != 3) {
    return projection;
  }
  const FloatArrayF<3> pseudoNormal(pseudoNormalArray);
  const double orientation =
    dot(gapVector, pseudoNormal) < 0.0 ? -1.0 : 1.0;
  const double distance = std::sqrt(std::max(0.0, distanceSquared));
  FloatArrayF<3> normal = pseudoNormal;
  if (distance > cppAlgebraTolerance) {
    normal = orientation * gapVector / distance;
  }

  const FloatArrayF<3> axis =
    std::abs(normal.at(1)) < 0.9
      ? FloatArrayF<3> { 1.0, 0.0, 0.0 }
      : FloatArrayF<3> { 0.0, 1.0, 0.0 };
  FloatArrayF<3> tangent1 = cross(axis, normal);
  const double tangentNorm = norm(tangent1);
  if (!(tangentNorm > 0.0)) {
    return projection;
  }
  tangent1 /= tangentNorm;
  const FloatArrayF<3> tangent2 = cross(normal, tangent1);

  projection.localCoordinates = localVertices[vertexIndex - 1];
  projection.gap = orientation * distance;
  projection.distanceSquared = distanceSquared;
  projection.normal = normal;
  projection.tangent1 = tangent1;
  projection.tangent2 = tangent2;
  projection.valid = true;
  return projection;
}

ContactProjection
FEContactSurface :: findContactPointOnVertex_2d(
    ContactPoint *contactPoint, ContactElement *element, int vertexIndex,
    TimeStep *tStep)
{
  ContactProjection projection;
  projection.featureType = ContactFeatureType::Vertex;
  projection.elementId = element->giveNumber();
  projection.featureIndex = vertexIndex;

  const auto localVertices = giveFeatureLocalVertices(element);
  if (vertexIndex < 1 || vertexIndex > static_cast<int>(localVertices.size())) {
    return projection;
  }

  const FloatArrayF<2> vertex = giveUpdatedNodeCoordinates2d(
    element->giveNode(vertexIndex), tStep);
  Coordinates slaveCoordinates;
  contactPoint->giveUpdatedCoordinates(slaveCoordinates, tStep);
  const FloatArrayF<2> slave(slaveCoordinates);
  const FloatArrayF<2> gapVector = slave - vertex;
  const double distanceSquared = dot(gapVector, gapVector);
  if (!std::isfinite(distanceSquared)) {
    return projection;
  }

  const FloatArray pseudoNormalArray = this->computeIncidentPseudoNormal2d(
    element, element->giveNode(vertexIndex)->giveNumber(), tStep);
  if (pseudoNormalArray.giveSize() != 2) {
    return projection;
  }
  const FloatArrayF<2> pseudoNormal(pseudoNormalArray);
  const double orientation =
    dot(gapVector, pseudoNormal) < 0.0 ? -1.0 : 1.0;
  const double distance = std::sqrt(std::max(0.0, distanceSquared));
  FloatArrayF<2> normal = pseudoNormal;
  if (distance > cppAlgebraTolerance) {
    normal = orientation * gapVector / distance;
  }
  const FloatArrayF<2> tangent = { -normal.at(2), normal.at(1) };

  projection.localCoordinates = localVertices[vertexIndex - 1];
  projection.gap = orientation * distance;
  projection.distanceSquared = distanceSquared;
  projection.normal = normal;
  projection.tangent1 = tangent;
  projection.valid = true;
  return projection;
}





}
