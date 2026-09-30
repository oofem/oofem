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

#include "contactsearch.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "floatmatrix.h"
#include "node.h"
#include "domain.h"
#include "fecontactsurface.h"
#include "contactpoint.h"
#include "contactpair.h"
#include "contactelement.h"
#include "datastream.h"
#include "contextioerr.h"
#include <limits>
#include <algorithm>
#include <cmath>


namespace oofem {

namespace {

bool isBetterContactProjection(
    double distanceSquared, double gap, int elementId,
    ContactFeatureType featureType, int featureIndex,
    double bestDistanceSquared, double bestGap, int bestElementId,
    ContactFeatureType bestFeatureType, int bestFeatureIndex,
    int currentElementId, ContactFeatureType currentFeatureType,
    int currentFeatureIndex, double facetOwnershipHysteresis = 0.0)
{
  if (bestElementId < 0) {
    return true;
  }

  const bool isCurrent = elementId == currentElementId
      && featureType == currentFeatureType
      && featureIndex == currentFeatureIndex;
  const bool bestIsCurrent = bestElementId == currentElementId
      && bestFeatureType == currentFeatureType
      && bestFeatureIndex == currentFeatureIndex;

  const double scale = std::max({1.0, std::abs(distanceSquared), std::abs(bestDistanceSquared)});
  double effectiveTolerance = 64.0 * std::numeric_limits<double>::epsilon() * scale;
  // At a near-degenerate shared edge/ridge between adjacent facets, distinct
  // candidates can have distances that differ only by noise-level amounts
  // that swap sign call to call, causing the owned facet to flip on every
  // residual/line-search evaluation and starving Newton of a stable
  // linearization. Require a genuinely, not just infinitesimally, better
  // distance to unseat the pair's currently owned facet.
  if (isCurrent != bestIsCurrent) {
    effectiveTolerance = std::max(effectiveTolerance, facetOwnershipHysteresis * scale);
  }

  if (distanceSquared < bestDistanceSquared - effectiveTolerance) {
    return true;
  }
  if (distanceSquared > bestDistanceSquared + effectiveTolerance) {
    return false;
  }

  const bool penetrating = gap <= 0.0;
  const bool bestPenetrating = bestGap <= 0.0;
  if (penetrating != bestPenetrating) {
    return penetrating;
  }

  if (isCurrent != bestIsCurrent) {
    return isCurrent;
  }

  if (featureType != bestFeatureType) {
    return static_cast<int>(featureType) > static_cast<int>(bestFeatureType);
  }
  if (featureIndex != bestFeatureIndex) {
    return featureIndex < bestFeatureIndex;
  }
  return elementId < bestElementId;
}

}


void
ContactSearchAlgorithm :: commitSearchState(TimeStep *tStep)
{
  for (auto &contactPair : contactPairs) {
    contactPair->commitSearchState(tStep);
  }
}

void
ContactSearchAlgorithm :: saveContext(DataStream &stream) const
{
  if (!stream.write(static_cast<int>(contactPairs.size()))) {
    THROW_CIOERR(CIO_IOERR);
  }
  for (const auto &contactPair : contactPairs) {
    contactPair->saveContext(stream);
  }
}

void
ContactSearchAlgorithm :: restoreContext(DataStream &stream)
{
  int pairCount = 0;
  if (!stream.read(pairCount)) {
    THROW_CIOERR(CIO_IOERR);
  }
  if (pairCount != static_cast<int>(contactPairs.size())) {
    OOFEM_ERROR("ContactSearchAlgorithm: stored pair count %d does not match reconstructed count %zu",
                pairCount, contactPairs.size());
  }
  for (auto &contactPair : contactPairs) {
    contactPair->restoreContext(stream);
  }
}


  ContactSearchAlgorithm_Surface2FESurface :: ContactSearchAlgorithm_Surface2FESurface(FEContactSurface *scs, FEContactSurface *mcs, Domain *d, int sd) : ContactSearchAlgorithm(d)
{
  this->slaveContactSurface = scs;
  this->masterContactSurface = mcs;
  this->surface_dimension  = sd;
}


void
ContactSearchAlgorithm_Surface2FESurface :: createContactPairs()
{
  contactPairs.clear();
  for ( int slave_ce = 1; slave_ce <= this->slaveContactSurface->giveNumberOfContactElements(); slave_ce++ ) {
    for ( GaussPoint *gp : * this->slaveContactSurface->giveContactElement_InSet(slave_ce)->giveDefaultIntegrationRulePtr() ) {
      auto cp_slave = std::make_unique<FEContactPoint_Slave>(this->slaveContactSurface, this->slaveContactSurface->giveContactElement_InSet(slave_ce)->giveNumber(),  surface_dimension, gp);
      auto contactPair = std::make_unique<ContactPair>(std::move(cp_slave));
      contactPairs.emplace_back(std::move(contactPair));
    }
  }
  // Seed slave bounds from the reference configuration before any prescribed
  // initial guess can move them across a master surface.
  ContactSearchAlgorithm :: commitSearchState(nullptr);

  // Seed swept broad-phase bounds with the undeformed master configuration so
  // the very first prescribed increment cannot cross contact undetected.
  committedMasterSearchAABBs.resize(this->masterContactSurface->giveNumberOfContactElements());
  for (int masterIndex = 1; masterIndex <= this->masterContactSurface->giveNumberOfContactElements(); ++masterIndex) {
    committedMasterSearchAABBs[masterIndex - 1] = this->masterContactSurface->giveContactElement_InSet(masterIndex)->computeAABB();
  }
}

void
ContactSearchAlgorithm_Surface2FESurface :: commitSearchState(TimeStep *tStep)
{
  ContactSearchAlgorithm :: commitSearchState(tStep);

  const int nMasterElements = this->masterContactSurface->giveNumberOfContactElements();
  committedMasterSearchAABBs.resize(nMasterElements);
  for (int i = 1; i <= nMasterElements; ++i) {
    committedMasterSearchAABBs[i - 1] =
        this->masterContactSurface->giveContactElement_InSet(i)->computeUpdatedAABB(tStep);
  }
}

void
ContactSearchAlgorithm_Surface2FESurface :: saveContext(DataStream &stream) const
{
  ContactSearchAlgorithm :: saveContext(stream);
  if (!stream.write(static_cast<int>(committedMasterSearchAABBs.size()))) {
    THROW_CIOERR(CIO_IOERR);
  }
  for (const AABB &aabb : committedMasterSearchAABBs) {
    const double bounds[6] = {
      aabb.min.x, aabb.min.y, aabb.min.z,
      aabb.max.x, aabb.max.y, aabb.max.z
    };
    if (!stream.write(bounds, 6)) {
      THROW_CIOERR(CIO_IOERR);
    }
  }
}

void
ContactSearchAlgorithm_Surface2FESurface :: restoreContext(DataStream &stream)
{
  ContactSearchAlgorithm :: restoreContext(stream);
  int masterCount = 0;
  if (!stream.read(masterCount) || masterCount < 0) {
    THROW_CIOERR(CIO_IOERR);
  }
  if (masterCount != masterContactSurface->giveNumberOfContactElements()) {
    OOFEM_ERROR("ContactSearchAlgorithm: stored master count %d does not match reconstructed count %d",
                masterCount, masterContactSurface->giveNumberOfContactElements());
  }
  committedMasterSearchAABBs.resize(masterCount);
  for (AABB &aabb : committedMasterSearchAABBs) {
    double bounds[6];
    if (!stream.read(bounds, 6)) {
      THROW_CIOERR(CIO_IOERR);
    }
    aabb = AABB(Vector(bounds[0], bounds[1], bounds[2]),
                Vector(bounds[3], bounds[4], bounds[5]));
  }
}

double
ContactSearchAlgorithm_Surface2FESurface :: computeSearchPadding(TimeStep *tStep)
{
  if (configuredSearchPadding >= 0.0) {
    return configuredSearchPadding;
  }

  double characteristicLength = 0.0;
  for (int i = 1; i <= this->masterContactSurface->giveNumberOfContactElements(); i++) {
    auto contactElement = this->masterContactSurface->giveContactElement_InSet(i);
    characteristicLength = std::max(characteristicLength, contactElement->computeUpdatedAABB(tStep).diagonalLength());
  }
  // Swept bounds cover motion already present in the current iterate.  A finite
  // proximity envelope is still required because the active-set search can run
  // before a prescribed displacement/Newton correction closes the current gap;
  // without it, the unconstrained problem may converge without another search.
  return std::max(1.e-12, 0.1 * characteristicLength);
}

AABB
ContactSearchAlgorithm_Surface2FESurface :: computeSweptMasterSearchAABB(int masterIndex, TimeStep *tStep) const
{
  ContactElement *element = this->masterContactSurface->giveContactElement_InSet(masterIndex);
  AABB currentAABB = element->computeUpdatedAABB(tStep);
  AABB sweptAABB = currentAABB;
  if (committedMasterSearchAABBs.size() == static_cast<std::size_t>(this->masterContactSurface->giveNumberOfContactElements())) {
    sweptAABB.merge(committedMasterSearchAABBs[masterIndex - 1]);
  }
  return sweptAABB;
}

void
ContactSearchAlgorithm_Surface2FESurface :: updateMasterSearchAABBs(TimeStep *tStep)
{
  const int nMasterElements = this->masterContactSurface->giveNumberOfContactElements();
  masterSearchAABBs.resize(nMasterElements);
  searchPadding = this->computeSearchPadding(tStep);

  for (int i = 1; i <= nMasterElements; i++) {
    masterSearchAABBs[i - 1] = this->computeSweptMasterSearchAABB(i, tStep);
  }
  for (auto &masterAABB : masterSearchAABBs) {
    masterAABB.expandBy(searchPadding);
  }
}

std::vector<int>
ContactSearchAlgorithm_Surface2FESurface :: findCandidateMasterContactElements(ContactPair *cp, TimeStep *tStep)
{
  std::vector<int> candidates;
  AABB slaveAABB = cp->computeSweptSlaveAABB(tStep);
  slaveAABB.expandBy(searchPadding);
  for (size_t i = 0; i < masterSearchAABBs.size(); i++) {
    if (slaveAABB.intersects(masterSearchAABBs[i])) {
      candidates.push_back(static_cast<int>(i) + 1);
    }
  }

  this->appendCurrentActiveMasterCandidate(cp, candidates);
  return candidates;
}

void
ContactSearchAlgorithm_Surface2FESurface :: appendCurrentActiveMasterCandidate(
    ContactPair *cp, std::vector<int> &candidates) const
{
  if (!cp->hasActiveContact()) {
    return;
  }

  auto *currentMaster = dynamic_cast<FEContactPoint *>(cp->giveMasterContactPoint());
  if (currentMaster == nullptr) {
    return;
  }

  const int currentElementId = currentMaster->giveContactElementId();
  const int nMasterElements = this->masterContactSurface->giveNumberOfContactElements();
  for (int i = 1; i <= nMasterElements; ++i) {
    if (this->masterContactSurface->giveContactElement_InSet(i)->giveNumber() == currentElementId) {
      if (std::find(candidates.begin(), candidates.end(), i) == candidates.end()) {
        candidates.push_back(i);
      }
      return;
    }
  }
}

void
ContactSearchAlgorithm_Surface2FESurface :: updatePairFrom3dCandidates(ContactPair *cp, const std::vector<int> &masterCandidateIndices, TimeStep *tStep)
{
  auto slavePoint = dynamic_cast<FEContactPoint_Slave*> (cp->giveSlaveContactPoint());
  int closestContactElementId = -1;
  int currentContactElementId = -1;
  ContactFeatureType currentFeatureType = ContactFeatureType::Surface;
  int currentFeatureIndex = 0;
  if (auto *currentMaster = dynamic_cast<FEContactPoint *>(cp->giveMasterContactPoint())) {
    currentContactElementId = currentMaster->giveContactElementId();
    currentFeatureType = currentMaster->giveContactFeatureType();
    currentFeatureIndex = currentMaster->giveContactFeatureIndex();
  }
  ContactFeatureType closestFeatureType = ContactFeatureType::Surface;
  int closestFeatureIndex = 0;
  FloatArray contactPointLocalCoordinates;
  double gap = std::numeric_limits<double>::infinity();
  double distanceSquared = std::numeric_limits<double>::infinity();
  FloatArrayF<3> normalVector, tangentVector1, tangentVector2;
  FloatArrayF<3> projectionDirection;
  if (directionalProjection) {
    projectionDirection = this->slaveContactSurface->computeContactPointNormal_3d(
      slavePoint, tStep);
    if (!(norm(projectionDirection) > 0.0)) {
      cp->clearContactState();
      return;
    }
  }
  const auto projectOnSurface = [&](ContactElement *element) {
    if (directionalProjection) {
      return this->masterContactSurface->findContactPointAlongDirectionInElement_3d(
        slavePoint, element, projectionDirection, tStep);
    }
    return this->masterContactSurface->findContactPointInElement_3d(
      slavePoint, element, tStep);
  };

  // Keep the current smooth master facet while its unconstrained projection
  // remains inside that facet.  This gives each integration point a persistent
  // facet owner and prevents a global closest-point search from alternating
  // between adjacent facets at a common edge during Newton iterations.  Once
  // the projection exits the current facet, the global candidate search below
  // performs the facet transition.  This is also the update order used by
  // other codes' sliding-elastic interfaces (old facet first, global search
  // second).
  bool currentSurfaceStillValid = false;
  if (currentContactElementId >= 0
      && currentFeatureType == ContactFeatureType::Surface) {
    ContactElement *currentElement =
      this->masterContactSurface->giveContactElement(currentContactElementId);
    if (currentElement != nullptr
        && slavePoint->giveContactElementId() != currentElement->giveNumber()) {
      auto [inElement, localCoords, newGap, newDistanceSquared, normal, t1, t2] =
        projectOnSurface(currentElement);
      if (inElement) {
        currentSurfaceStillValid = true;
        gap = newGap;
        distanceSquared = newDistanceSquared;
        normalVector = normal;
        tangentVector1 = t1;
        tangentVector2 = t2;
        contactPointLocalCoordinates = localCoords;
        closestContactElementId = currentElement->giveNumber();
        closestFeatureType = ContactFeatureType::Surface;
        closestFeatureIndex = 0;
      }
    }
  }

  const std::vector<int> noCandidates;
  const std::vector<int> &candidatesToSearch =
    currentSurfaceStillValid ? noCandidates : masterCandidateIndices;
  for (int i : candidatesToSearch) {
    auto contactElement = this->masterContactSurface->giveContactElement_InSet(i);
    if (slavePoint->giveContactElementId() == contactElement->giveNumber()) {
      continue;
    }
    auto [inElement, localCoords, newGap, newDistanceSquared, normal, t1, t2] =
      projectOnSurface(contactElement);
    if (inElement && isBetterContactProjection(
          newDistanceSquared, newGap, contactElement->giveNumber(),
          ContactFeatureType::Surface, 0,
          distanceSquared, gap, closestContactElementId,
          closestFeatureType, closestFeatureIndex,
          currentContactElementId, currentFeatureType, currentFeatureIndex,
          facetOwnershipHysteresis)) {
      gap = newGap;
      distanceSquared = newDistanceSquared;
      normalVector = normal;
      tangentVector1 = t1;
      tangentVector2 = t2;
      contactPointLocalCoordinates = localCoords;
      closestContactElementId = contactElement->giveNumber();
      closestFeatureType = ContactFeatureType::Surface;
      closestFeatureIndex = 0;
    }

    if (generalizedFeatures) {
      for (int edge = 1;
           edge <= this->masterContactSurface->giveNumberOfEdges(contactElement);
           ++edge) {
        ContactProjection projection =
          this->masterContactSurface->findContactPointOnEdge_3d(
            slavePoint, contactElement, edge, tStep);
        if (projection.valid && isBetterContactProjection(
              projection.distanceSquared, projection.gap, projection.elementId,
              projection.featureType, projection.featureIndex,
              distanceSquared, gap, closestContactElementId,
              closestFeatureType, closestFeatureIndex,
              currentContactElementId, currentFeatureType, currentFeatureIndex,
              facetOwnershipHysteresis)) {
          gap = projection.gap;
          distanceSquared = projection.distanceSquared;
          normalVector = projection.normal;
          tangentVector1 = projection.tangent1;
          tangentVector2 = projection.tangent2;
          contactPointLocalCoordinates = projection.localCoordinates;
          closestContactElementId = projection.elementId;
          closestFeatureType = projection.featureType;
          closestFeatureIndex = projection.featureIndex;
        }
      }

      for (int vertex = 1;
           vertex <= this->masterContactSurface->giveNumberOfVertices(contactElement);
           ++vertex) {
        ContactProjection projection =
          this->masterContactSurface->findContactPointOnVertex_3d(
            slavePoint, contactElement, vertex, tStep);
        if (projection.valid && isBetterContactProjection(
              projection.distanceSquared, projection.gap, projection.elementId,
              projection.featureType, projection.featureIndex,
              distanceSquared, gap, closestContactElementId,
              closestFeatureType, closestFeatureIndex,
              currentContactElementId, currentFeatureType, currentFeatureIndex,
              facetOwnershipHysteresis)) {
          gap = projection.gap;
          distanceSquared = projection.distanceSquared;
          normalVector = projection.normal;
          tangentVector1 = projection.tangent1;
          tangentVector2 = projection.tangent2;
          contactPointLocalCoordinates = projection.localCoordinates;
          closestContactElementId = projection.elementId;
          closestFeatureType = projection.featureType;
          closestFeatureIndex = projection.featureIndex;
        }
      }
    }
  }

  if (closestContactElementId != -1) {
    auto master_point = std::make_unique<FEContactPoint_Master>(
      this->masterContactSurface, closestContactElementId, 2,
      contactPointLocalCoordinates, closestFeatureType, closestFeatureIndex);
    cp->setMasterContactPoint(std::move(master_point));
    cp->setNormalGap(gap);
    cp->setNormalVector(normalVector);
    cp->setTangentVector1(tangentVector1);
    cp->setTangentVector2(tangentVector2);
  } else {
    cp->clearContactState();
  }
}

void
ContactSearchAlgorithm_Surface2FESurface :: updatePairFrom2dCandidates(ContactPair *cp, const std::vector<int> &masterCandidateIndices, TimeStep *tStep)
{
  auto slavePoint = dynamic_cast<FEContactPoint_Slave*> (cp->giveSlaveContactPoint());
  int closestContactElementId = -1;
  int currentContactElementId = -1;
  ContactFeatureType currentFeatureType = ContactFeatureType::Surface;
  int currentFeatureIndex = 0;
  if (auto *currentMaster = dynamic_cast<FEContactPoint *>(cp->giveMasterContactPoint())) {
    currentContactElementId = currentMaster->giveContactElementId();
    currentFeatureType = currentMaster->giveContactFeatureType();
    currentFeatureIndex = currentMaster->giveContactFeatureIndex();
  }
  ContactFeatureType closestFeatureType = ContactFeatureType::Surface;
  int closestFeatureIndex = 0;
  FloatArray contactPointLocalCoordinates;
  double gap = std::numeric_limits<double>::infinity();
  double distanceSquared = std::numeric_limits<double>::infinity();
  FloatArrayF<2> normalVector, tangentVector1;

  // Preserve the current smooth segment until its unconstrained projection
  // leaves the segment.  The global search is then responsible for selecting
  // the next segment.  Besides matching the 3-D facet policy above, this avoids
  // artificial contact-map switching at shared line endpoints.
  bool currentSurfaceStillValid = false;
  if (currentContactElementId >= 0
      && currentFeatureType == ContactFeatureType::Surface) {
    ContactElement *currentElement =
      this->masterContactSurface->giveContactElement(currentContactElementId);
    if (currentElement != nullptr
        && slavePoint->giveContactElementId() != currentElement->giveNumber()) {
      auto [inElement, localCoords, newGap, newDistanceSquared, normal, t1] =
        this->masterContactSurface->findContactPointInElement_2d(
          slavePoint, currentElement, tStep);
      if (inElement) {
        currentSurfaceStillValid = true;
        gap = newGap;
        distanceSquared = newDistanceSquared;
        normalVector = normal;
        tangentVector1 = t1;
        contactPointLocalCoordinates = localCoords;
        closestContactElementId = currentElement->giveNumber();
        closestFeatureType = ContactFeatureType::Surface;
        closestFeatureIndex = 0;
      }
    }
  }

  const std::vector<int> noCandidates;
  const std::vector<int> &candidatesToSearch =
    currentSurfaceStillValid ? noCandidates : masterCandidateIndices;
  for (int i : candidatesToSearch) {
    auto masterContactElement = this->masterContactSurface->giveContactElement_InSet(i);
    if (slavePoint->giveContactElementId() == masterContactElement->giveNumber()) {
      continue;
    }
    auto [inElement, localCoords, newGap, newDistanceSquared, normal, t1] =
      this->masterContactSurface->findContactPointInElement_2d(slavePoint, masterContactElement, tStep);
    if (inElement && isBetterContactProjection(
          newDistanceSquared, newGap, masterContactElement->giveNumber(),
          ContactFeatureType::Surface, 0,
          distanceSquared, gap, closestContactElementId,
          closestFeatureType, closestFeatureIndex,
          currentContactElementId, currentFeatureType, currentFeatureIndex,
          facetOwnershipHysteresis)) {
      gap = newGap;
      distanceSquared = newDistanceSquared;
      normalVector = normal;
      tangentVector1 = t1;
      contactPointLocalCoordinates = localCoords;
      closestContactElementId = masterContactElement->giveNumber();
      closestFeatureType = ContactFeatureType::Surface;
      closestFeatureIndex = 0;
    }

    if (generalizedFeatures) {
      for (int vertex = 1;
           vertex <= this->masterContactSurface->giveNumberOfVertices(masterContactElement);
           ++vertex) {
        ContactProjection projection =
          this->masterContactSurface->findContactPointOnVertex_2d(
            slavePoint, masterContactElement, vertex, tStep);
        if (projection.valid && isBetterContactProjection(
              projection.distanceSquared, projection.gap, projection.elementId,
              projection.featureType, projection.featureIndex,
              distanceSquared, gap, closestContactElementId,
              closestFeatureType, closestFeatureIndex,
              currentContactElementId, currentFeatureType, currentFeatureIndex,
              facetOwnershipHysteresis)) {
          gap = projection.gap;
          distanceSquared = projection.distanceSquared;
          normalVector = projection.normal;
          tangentVector1 = projection.tangent1;
          contactPointLocalCoordinates = projection.localCoordinates;
          closestContactElementId = projection.elementId;
          closestFeatureType = projection.featureType;
          closestFeatureIndex = projection.featureIndex;
        }
      }
    }
  }

  if (closestContactElementId != -1) {
    auto masterContactPoint = std::make_unique<FEContactPoint_Master>(
      this->masterContactSurface, closestContactElementId, 1,
      contactPointLocalCoordinates, closestFeatureType, closestFeatureIndex);
    cp->setMasterContactPoint(std::move(masterContactPoint));
    cp->setNormalGap(gap);
    cp->setNormalVector(normalVector);
    cp->setTangentVector1(tangentVector1);
  } else {
    cp->clearContactState();
  }
}

  ContactSearchAlgorithm_Surface2FESurface_3d ::   ContactSearchAlgorithm_Surface2FESurface_3d(FEContactSurface *scs, FEContactSurface *mcs, Domain *d) :   ContactSearchAlgorithm_Surface2FESurface(scs, mcs, d, 2)
{
}


void
ContactSearchAlgorithm_Surface2FESurface_3d :: updateContactPairs(TimeStep *tStep)
{
  this->updateMasterSearchAABBs(tStep);
  for (auto &cp : contactPairs) {
    auto candidates = this->findCandidateMasterContactElements(cp.get(), tStep);
    this->updatePairFrom3dCandidates(cp.get(), candidates, tStep);
  }

}

  

  ContactSearchAlgorithm_Surface2FESurface_2d ::   ContactSearchAlgorithm_Surface2FESurface_2d(FEContactSurface *scs, FEContactSurface *mcs, Domain *d) :   ContactSearchAlgorithm_Surface2FESurface(scs, mcs, d, 1)
{
}
  

void
ContactSearchAlgorithm_Surface2FESurface_2d :: updateContactPairs(TimeStep *tStep)
{
  this->updateMasterSearchAABBs(tStep);
  for (auto &cp : contactPairs) {
    auto candidates = this->findCandidateMasterContactElements(cp.get(), tStep);
    this->updatePairFrom2dCandidates(cp.get(), candidates, tStep);
  }
  
}


  
} // end namespace oofem
