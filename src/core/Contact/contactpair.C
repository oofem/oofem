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

namespace oofem {

ContactPair ::  ContactPair(std::unique_ptr<ContactPoint> s) : slave(std::move(s)), tractionVector()
{
  tractionVector.resize(2);
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
  for (int i = 1; i <= 2; i++) {
    ret.emplace_back(giveTangentVector(i));
  }
  return ret;
}

std::vector<FloatArray>
ContactPair :: givePreviousTangentVectors() const {
  std::vector<FloatArray> ret;
  for (int i = 1; i <= 2; i++) {
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
  this->master->computeCurvature(G, this->normalVector, tStep);
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
  IntArray loc_slave;
  this->master->giveLocationArray(loc, dofs, ns);
  this->slave->giveLocationArray(loc_slave, dofs, ns);
  //
  loc.followedBy(loc_slave);
  
}


void
ContactPair :: updateYourself(TimeStep *tStep)
{
  if(this->giveNormalGap() <= 0 && this->inContact()) {
    referenceContactPointInit = true;
    previousContactPointCoords = master->giveGlobalCoordinates();
    tractionVector = tempTractionVector;
    previousNormalVector = normalVector;
    previousTangentVector1 = tangentVector1;
    previousTangentVector2 = tangentVector2;
  } else {
    referenceContactPointInit = false;
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
ContactPair :: computeContactPointDisplacement() const
{
  // @todo
  if (previousContactPointCoords.giveSize() == 0) {
    FloatArray ret = master->giveGlobalCoordinates();
    ret.times(0);
    return ret;
  }
  FloatArray contactPointCoords = master->giveGlobalCoordinates();
  return contactPointCoords - previousContactPointCoords;
}

AABB
ContactPair :: computeSlaveAABB()
{
  // TODO
  auto coords = slave->giveGlobalCoordinates();
  double x = coords.at(1);
  double y = coords.at(2);
  double z = coords.at(3);
  AABB aabb(Vector(x, y, z), Vector(x, y, z));
  //
  aabb.min.x -= 0.5;
  aabb.min.y -= 0.5;
  aabb.min.z -= 0.5;
  //
  aabb.max.x += 0.5;
  aabb.max.y += 0.5;
  aabb.max.z += 0.5;
  //
  return aabb;
}

};

