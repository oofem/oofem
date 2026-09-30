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


#include "contactpair.h"
#include "contactelement.h"
#include "node.h"
#include "datastream.h"
#include "contextioerr.h"
#include <algorithm>
#include <limits>
#include <cmath>

namespace oofem {

ContactPair ::  ContactPair(std::unique_ptr<ContactPoint> s) : slave(std::move(s)), tractionVector()
{
  const int tangentSize = slave ? slave->giveSurfaceDimension() : 0;
  tractionVector.resize(tangentSize);
  tractionVector.zero();
  tempTractionVector.resize(tangentSize);
  tempTractionVector.zero();
  dXi = 0.0;
  temp_dXi = 0.0;
}

void
ContactPair :: setNormalGap(double ng)
{
  normal_gap = ng;
  auto *slavePoint = dynamic_cast<FEContactPoint_Slave *>(slave.get());
  if (slavePoint != nullptr) {
    slavePoint->giveContactElement()->setContactOutputState(
      slavePoint->giveIntegrationPoint(), ng, 0.0, 0);
  }
}

void
ContactPair :: setOutputContactState(double pressure, int status)
{
  auto *slavePoint = dynamic_cast<FEContactPoint_Slave *>(slave.get());
  if (slavePoint != nullptr) {
    const double outputGap = std::isfinite(normal_gap) ? normal_gap : 0.0;
    slavePoint->giveContactElement()->setContactOutputState(
      slavePoint->giveIntegrationPoint(), outputGap, pressure, status);
  }
}


const FloatArray &
ContactPair :: giveTangentVector(int i) const
{
  if(i == 1) {
    return tangentVector1;
  } else if(i == 2) {
    return tangentVector2;
  } else {
    OOFEM_ERROR("ContactPair:Wrong number of tangent vector");
  }
}

const FloatArray &
ContactPair :: givePreviousTangentVector(int i) const
{
  if(i == 1) {
    // @todo
    if (previousTangentVector1.giveSize() == 0) {
      return tangentVector1;
    }
    return previousTangentVector1;
  } else if(i == 2) {
    // @todo
    if (previousTangentVector2.giveSize() == 0) {
      return tangentVector2;
    }
    return previousTangentVector2;
  } else {
    OOFEM_ERROR("ContactPair:Wrong number of tangent vector");
  }
}

std::vector<FloatArray>
ContactPair :: giveTangentVectors() const {
  std::vector<FloatArray> ret;
  const int tangentSize = slave ? slave->giveSurfaceDimension() : 0;
  for (int i = 1; i <= tangentSize; i++) {
    ret.emplace_back(giveTangentVector(i));
  }
  return ret;
}

std::vector<FloatArray>
ContactPair :: givePreviousTangentVectors() const {
  std::vector<FloatArray> ret;
  const int tangentSize = slave ? slave->giveSurfaceDimension() : 0;
  for (int i = 1; i <= tangentSize; i++) {
    ret.emplace_back(givePreviousTangentVector(i));
  }
  return ret;
}

void
ContactPair :: computeNmatrix(FloatMatrix &answer)
{
  FloatMatrix answer_slave, answer_master;
  this->master->computeNmatrix(answer_master);
  this->slave->computeNmatrix(answer_slave);
  answer_master.times(-1);
  //
  auto master_ncols = answer_master.giveNumberOfColumns();
  auto master_nrows  = answer_master.giveNumberOfRows();
  auto slave_ncols  = answer_slave.giveNumberOfColumns();
  //
  answer.resize(master_nrows, master_ncols+slave_ncols);
  //
  answer.setSubMatrix(answer_master, 1, 1);
  answer.setSubMatrix(answer_slave, 1, master_ncols+1);
}



  

void
ContactPair :: compute_dNdxi_matrices(std::vector<FloatMatrix> &dNdxi)
{
  FloatMatrix dNs, dNm;
  auto sd = this->slave->giveSurfaceDimension();
  FloatMatrix dN(sd+1, 2 * (sd+1));
  for(int i = 1;  i <= master->giveSurfaceDimension(); i++) {
    dN.zero();
    master->compute_dNdxi_matrix(dNm, i);
    slave->compute_dNdxi_matrix(dNs, i);
    //
    auto master_ncols = dNm.giveNumberOfColumns();
    auto master_nrows = dNm.giveNumberOfRows();
    auto slave_ncols  = dNs.giveNumberOfColumns();
    //
    dN.resize(master_nrows, slave_ncols + master_ncols);
    dNm.times(-1.);
    //
    dN.setSubMatrix(dNm, 1, 1);
    //@todo: node-2-surface vs surface-2-surface!?!
    //dN.setSubMatrix(dNs, 1, master_ncols+1);
    dNdxi.emplace_back(dN);
  }
}

 

void
ContactPair :: computeCurvature(FloatMatrix &G, TimeStep *tStep)
{
  FloatArray unitNormal = this->normalVector;
  const double normalNorm = unitNormal.computeNorm();
  if (normalNorm <= 0.0) {
    OOFEM_ERROR("ContactPair: zero normal vector in curvature computation");
  }
  unitNormal /= normalNorm;
  this->master->computeCurvature(G, unitNormal, tStep);
}

void
ContactPair :: computeSecondBaseVectors(std::vector<std::vector<FloatArray>> &answer, TimeStep *tStep)
{
  this->master->computeSecondBaseVectors(answer, tStep);
}
 

double
ContactPair :: giveContactArea(TimeStep *)
{
  const double contactArea = slave->giveIntegrationWeight() * slave->giveReferenceSurfaceMeasure();
  if (contactArea <= 0.0) {
    OOFEM_ERROR("ContactPair: non-positive contact integration area");
  }
  return contactArea;
}

  
void
ContactPair :: initContactPoint()
{
  if (referenceContactPointInit == false) {
    referenceContactPointCoords = slave->giveLocalCoordinates();
    contactPointCoords = master->giveGlobalCoordinates();
    referenceContactPointInit = true;
  }
}



  

void
ContactPair :: giveLocationArray(const IntArray &dofs, IntArray &loc, const UnknownNumberingScheme &ns) const
{
  this->giveRowLocationArray(dofs, loc, ns);
}

void
ContactPair :: giveRowLocationArray(const IntArray &dofs, IntArray &loc,
                                    const UnknownNumberingScheme &ns) const
{
  IntArray loc_slave;
  this->master->giveLocationArray(loc, dofs, ns);
  this->slave->giveLocationArray(loc_slave, dofs, ns);
  loc.followedBy(loc_slave);
}

void
ContactPair :: giveColumnLocationArray(const IntArray &dofs, IntArray &loc,
                                       const UnknownNumberingScheme &ns) const
{
  loc.clear();
  for (Node *node : this->giveLinearizationNodes()) {
    IntArray nodeLoc;
    node->giveLocationArray(dofs, nodeLoc, ns);
    loc.followedBy(nodeLoc);
  }
}

std::vector<Node *>
ContactPair :: giveResidualNodes() const
{
  std::vector<Node *> nodes;
  const auto appendElementNodes = [&nodes](ContactElement *element) {
    if (element == nullptr) {
      return;
    }
    for (int i = 1; i <= element->giveNumberOfNodes(); ++i) {
      nodes.push_back(element->giveNode(i));
    }
  };

  appendElementNodes(master ? master->giveContactElement() : nullptr);
  appendElementNodes(slave ? slave->giveContactElement() : nullptr);
  return nodes;
}

ContactElement *
ContactPair :: givePreviousMasterContactElement() const
{
  if (!this->hasProjectionHistory()) {
    return nullptr;
  }
  auto *masterPoint = dynamic_cast<FEContactPoint *>(master.get());
  if (masterPoint == nullptr) {
    return nullptr;
  }
  return masterPoint->giveContactElementOnSurface(previousMasterElementId);
}

bool
ContactPair :: hasMasterFacetTransition() const
{
  auto *masterPoint = dynamic_cast<FEContactPoint *>(master.get());
  return masterPoint != nullptr && this->hasProjectionHistory()
      && masterPoint->giveContactElementId() != previousMasterElementId;
}

ContactFeatureType
ContactPair :: giveCurrentMasterFeatureType() const
{
  auto *masterPoint = dynamic_cast<FEContactPoint *>(master.get());
  return masterPoint ? masterPoint->giveContactFeatureType()
                     : ContactFeatureType::Surface;
}

int
ContactPair :: giveCurrentMasterFeatureIndex() const
{
  auto *masterPoint = dynamic_cast<FEContactPoint *>(master.get());
  return masterPoint ? masterPoint->giveContactFeatureIndex() : 0;
}

bool
ContactPair :: hasMasterFeatureTransition() const
{
  auto *masterPoint = dynamic_cast<FEContactPoint *>(master.get());
  return masterPoint != nullptr && this->hasProjectionHistory()
      && (masterPoint->giveContactElementId() != previousMasterElementId
          || masterPoint->giveContactFeatureType() != previousMasterFeatureType
          || masterPoint->giveContactFeatureIndex() != previousMasterFeatureIndex);
}

std::vector<Node *>
ContactPair :: giveLinearizationNodes() const
{
  std::vector<Node *> nodes = this->giveResidualNodes();
  ContactElement *previousMaster = this->givePreviousMasterContactElement();
  if (previousMaster == nullptr) {
    return nodes;
  }

  for (int i = 1; i <= previousMaster->giveNumberOfNodes(); ++i) {
    Node *node = previousMaster->giveNode(i);
    if (std::find(nodes.begin(), nodes.end(), node) == nodes.end()) {
      nodes.push_back(node);
    }
  }
  return nodes;
}


void
ContactPair :: clearContactState()
{
  this->setOutputContactState(0.0, 0);
  master.reset();
  normal_gap = std::numeric_limits<double>::infinity();
  normalVector.clear();
  tangentVector1.clear();
  tangentVector2.clear();

  // Search misses are iteration-local; committed friction history is cleared only after convergence.
  referenceContactPointCoords.clear();
  tempReferenceContactPointCoords.clear();
  contactPointCoords.clear();
  referenceContactPointInit = false;

  const int tangentSize = slave ? slave->giveSurfaceDimension() : 0;
  tempTractionVector.resize(tangentSize);
  tempTractionVector.zero();
  tempAccumulatedPlasticSlip = 0.0;
  temp_dXi = 0.0;
}

void
ContactPair :: updateYourself(TimeStep *tStep)
{
  const bool activeContact = this->hasActiveContact();
  if(this->hasMasterContact()) {
    referenceContactPointInit = true;
    Coordinates updatedMasterCoords;
    master->giveUpdatedCoordinates(updatedMasterCoords, tStep);
    previousContactPointCoords = updatedMasterCoords;
    previousMasterLocalCoordinates = master->giveLocalCoordinates();
    auto *masterPoint = static_cast<FEContactPoint *>(master.get());
    previousMasterElementId = masterPoint->giveContactElementId();
    previousMasterFeatureType = masterPoint->giveContactFeatureType();
    previousMasterFeatureIndex = masterPoint->giveContactFeatureIndex();
    previousContactActive = activeContact;
    previousNormalVector = normalVector;
    previousTangentVector1 = tangentVector1;
    previousTangentVector2 = tangentVector2;
    if (activeContact) {
      tractionVector = tempTractionVector;
      accumulatedPlasticSlip = tempAccumulatedPlasticSlip;
    } else {
      tractionVector.zero();
      tempTractionVector.zero();
      accumulatedPlasticSlip = 0.0;
      tempAccumulatedPlasticSlip = 0.0;
    }
  } else {
    referenceContactPointInit = false;
    previousContactPointCoords.clear();
    previousMasterLocalCoordinates.clear();
    previousMasterElementId = -1;
    previousMasterFeatureType = ContactFeatureType::Surface;
    previousMasterFeatureIndex = 0;
    previousContactActive = false;
    previousNormalVector.clear();
    previousTangentVector1.clear();
    previousTangentVector2.clear();
    tractionVector.zero();
    tempTractionVector.zero();
    accumulatedPlasticSlip = 0.0;
    tempAccumulatedPlasticSlip = 0.0;
  }


}





void
ContactPair :: computeVectorOf(ValueModeType u, TimeStep *tStep, FloatArray &answer)
{
  FloatArray s_vec;
  this->master->computeVectorOf(u, tStep, answer); 
  this->slave->computeVectorOf(u, tStep, s_vec);
  //
  int offset = answer.giveSize();
  answer.copySubVector(s_vec,offset+1);

}
  
FloatArray
ContactPair :: computeContactPointDisplacement(TimeStep *tStep) const
{
  const int spatialDimension = (slave ? slave->giveSurfaceDimension() : 2) + 1;
  if (previousMasterElementId < 0 || previousMasterLocalCoordinates.giveSize() == 0) {
    return FloatArray(spatialDimension);
  }
  auto *masterPoint = dynamic_cast<FEContactPoint *>(master.get());
  if (masterPoint == nullptr) {
    OOFEM_ERROR("ContactPair: objective friction update requires an FE master contact point");
  }
  Coordinates currentProjection, previousProjectionInCurrentConfiguration;
  master->giveUpdatedCoordinates(currentProjection, tStep);
  masterPoint->giveUpdatedCoordinatesOnElement(previousProjectionInCurrentConfiguration,
                                                previousMasterElementId,
                                                previousMasterLocalCoordinates, tStep);
  const Coordinates fullDisplacement = currentProjection - previousProjectionInCurrentConfiguration;
  FloatArray displacement(spatialDimension);
  for (int k = 1; k <= spatialDimension; ++k) {
    displacement.at(k) = fullDisplacement.at(k);
  }
  return displacement;
}

AABB
ContactPair :: computeSlaveAABB(TimeStep *tStep) const
{
  Coordinates coords;
  slave->giveUpdatedCoordinates(coords, tStep);
  double x = coords.at(1);
  double y = coords.at(2);
  double z = (coords.giveSize() > 2) ? coords.at(3) : 0.0;
  AABB aabb(Vector(x, y, z), Vector(x, y, z));
  return aabb;
}

AABB
ContactPair :: computeSweptSlaveAABB(TimeStep *tStep) const
{
  AABB currentAABB = this->computeSlaveAABB(tStep);
  AABB sweptAABB = currentAABB;
  if (committedSlaveSearchAABBInitialized) {
    sweptAABB.merge(committedSlaveSearchAABB);
  }
  return sweptAABB;
}

void
ContactPair :: commitSearchState(TimeStep *tStep)
{
  committedSlaveSearchAABB = this->computeSlaveAABB(tStep);
  committedSlaveSearchAABBInitialized = true;
}

void
ContactPair :: saveContext(DataStream &stream) const
{
  contextIOResultType result = CIO_OK;
  if ((result = previousContactPointCoords.storeYourself(stream)) != CIO_OK
      || (result = previousMasterLocalCoordinates.storeYourself(stream)) != CIO_OK
      || !stream.write(previousMasterElementId)
      || !stream.write(static_cast<int>(previousMasterFeatureType))
      || !stream.write(previousMasterFeatureIndex)
      || !stream.write(previousContactActive)
      || (result = previousNormalVector.storeYourself(stream)) != CIO_OK
      || (result = previousTangentVector1.storeYourself(stream)) != CIO_OK
      || (result = previousTangentVector2.storeYourself(stream)) != CIO_OK
      || (result = tractionVector.storeYourself(stream)) != CIO_OK
      || !stream.write(accumulatedPlasticSlip)
      || !stream.write(committedSlaveSearchAABBInitialized)) {
    THROW_CIOERR(result == CIO_OK ? CIO_IOERR : result);
  }

  if (committedSlaveSearchAABBInitialized) {
    const double bounds[6] = {
      committedSlaveSearchAABB.min.x, committedSlaveSearchAABB.min.y, committedSlaveSearchAABB.min.z,
      committedSlaveSearchAABB.max.x, committedSlaveSearchAABB.max.y, committedSlaveSearchAABB.max.z
    };
    if (!stream.write(bounds, 6)) {
      THROW_CIOERR(CIO_IOERR);
    }
  }
}

void
ContactPair :: restoreContext(DataStream &stream)
{
  contextIOResultType result = CIO_OK;
  int storedFeatureType = static_cast<int>(ContactFeatureType::Surface);
  if ((result = previousContactPointCoords.restoreYourself(stream)) != CIO_OK
      || (result = previousMasterLocalCoordinates.restoreYourself(stream)) != CIO_OK
      || !stream.read(previousMasterElementId)
      || !stream.read(storedFeatureType)
      || !stream.read(previousMasterFeatureIndex)
      || !stream.read(previousContactActive)
      || (result = previousNormalVector.restoreYourself(stream)) != CIO_OK
      || (result = previousTangentVector1.restoreYourself(stream)) != CIO_OK
      || (result = previousTangentVector2.restoreYourself(stream)) != CIO_OK
      || (result = tractionVector.restoreYourself(stream)) != CIO_OK
      || !stream.read(accumulatedPlasticSlip)
      || !stream.read(committedSlaveSearchAABBInitialized)) {
    THROW_CIOERR(result == CIO_OK ? CIO_IOERR : result);
  }
  if (storedFeatureType < static_cast<int>(ContactFeatureType::Vertex)
      || storedFeatureType > static_cast<int>(ContactFeatureType::Surface)) {
    OOFEM_ERROR("ContactPair: invalid stored master feature type");
  }
  previousMasterFeatureType = static_cast<ContactFeatureType>(storedFeatureType);

  if (committedSlaveSearchAABBInitialized) {
    double bounds[6];
    if (!stream.read(bounds, 6)) {
      THROW_CIOERR(CIO_IOERR);
    }
    committedSlaveSearchAABB =
      AABB(Vector(bounds[0], bounds[1], bounds[2]), Vector(bounds[3], bounds[4], bounds[5]));
  }

  tempTractionVector = tractionVector;
  tempAccumulatedPlasticSlip = accumulatedPlasticSlip;
  previousContactActive = previousContactActive && this->hasProjectionHistory();
}

};
