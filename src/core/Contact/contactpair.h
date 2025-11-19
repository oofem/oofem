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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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

#ifndef contactpair_h
#define contactpair_h

#include "contactpoint.h"

namespace oofem {
class IntArray;
class FloatArray;
class FloatMatrix;

/**
 * Abstract base class for all contact finite elements. Derived classes should be  base
 * classes for specific analysis type (for example base class for structural analysis,
 * thermal analysis or magnetostatics one). These derived classes then declare
 * analysis-specific part of interface and they provide default implementation
 * for these methods.
 * This abstract class declares (and possibly implements) general data and methods
 * common to all element types. General methods for obtaining characteristic vectors,
 * matrices and values are introduced and should be used instead of calling directly
 * specific member functions (these must be overloaded by derived analysis-specific
 * classes in order to invoke proper method according to type of component requested).
 */

class ContactPair
{
protected:
  std::unique_ptr<ContactPoint> master;
  std::unique_ptr<ContactPoint> slave;
  double normal_gap = 0;
  FloatArray normalVector;
  FloatArray previousNormalVector;
  FloatArray tangentVector1;
  FloatArray previousTangentVector1;
  FloatArray tangentVector2;
  FloatArray previousTangentVector2;
  //
  FloatArray referenceContactPointCoords;
  FloatArray tempReferenceContactPointCoords;
  FloatArray contactPointCoords;
  //
  FloatArray previousContactPointCoords;
  bool referenceContactPointInit = false;
  //traction vectors
  FloatArray tractionVector;
  FloatArray tempTractionVector;
  // dxi;
  double dXi;
  double temp_dXi;
  

public:
  ContactPair(std::unique_ptr<ContactPoint> slave);
  ~ContactPair(){;}
  //
  virtual bool inContact() {
    if (!master) return false;
    return master->inContact();
  }
  //
  virtual void computeNmatrix(FloatMatrix &answer);
  virtual void compute_dNdxi_matrices(std::vector<FloatMatrix> &answer);
  virtual void computeCurvature(FloatMatrix &G, TimeStep *tStep);
  //
  ContactPoint *giveMasterContactPoint() {return master.get();}
  ContactPoint *giveSlaveContactPoint()  {return slave.get();}
  void giveLocationArray(const IntArray &dofs, IntArray &loc, const UnknownNumberingScheme &ns) const;
  double giveNormalGap() { return normal_gap;}
  void setNormalGap(double ng){this->normal_gap = ng;}
  //
  const FloatArray &giveNormalVector() const {return normalVector;} 
  const FloatArray &givePreviousNormalVector() const {return previousNormalVector;}
  const FloatArray &giveTangentVector(int i) const;
  const FloatArray &givePreviousTangentVector(int i) const;
  std::vector<FloatArray> giveTangentVectors() const;
  std::vector<FloatArray> givePreviousTangentVectors() const;
  //``
  const FloatArray &giveLocalCoordinates() const {return slave->giveLocalCoordinates();}
  void initContactPoint();

  FloatArray computeContactPointDisplacement() const;
  //
  void setNormalVector(const FloatArray &nv)   {this->normalVector = nv;}
  void setTangentVector1(const FloatArray &tv1)   {this->tangentVector1 = tv1;}
  void setTangentVector2(const FloatArray &tv2)   {this->tangentVector2 = tv2;}

  //
  void setTempTractionVector(const FloatArray &tv) {this->tempTractionVector = tv;}
  const FloatArray &giveTractionVector() const {return tractionVector;}
  //  
  void setMasterContactPoint(std::unique_ptr<ContactPoint> m) {master = std::move(m);}
  void setSlaveContactPoint(std::unique_ptr<ContactPoint> s)  {slave = std::move(s);}
  //
  virtual void updateYourself(TimeStep *tStep);
  void computeVectorOf(ValueModeType u, TimeStep *tStep, FloatArray &answer);
  //
  AABB computeSlaveAABB();
};
 
} // end namespace oofem
#endif //contactpair_h
