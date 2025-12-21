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

#ifndef contactpair_h
#define contactpair_h

#include "contactpoint.h"

namespace oofem {
class IntArray;
class FloatArray;
class FloatMatrix;



 /**
 * @brief Represents a contact interaction between a master and a slave contact point.
 *
 * The ContactPair class encapsulates all data and operations associated with a single
 * master–slave contact interaction in a contact mechanics formulation. It stores pointers
 * to the corresponding master and slave ContactPoint objects and maintains the quantities
 * required to evaluate contact constraints, such as normal gap,
 * contact forces (tractions), and their temporary or incremental values.
 *
 * The class provides functionality for:
 * - Updating the contact state in time,
 * - Computing contact-related vectors 
 * - Storing and accessing current and temporary traction vectors,
 * - Computing an axis-aligned bounding box (AABB) for the slave side, typically used for
 *   contact search and detection algorithms.
 *
 * A ContactPair does not own the global contact algorithm itself; instead, it serves as a
 * lightweight container and evaluator for a single detected contact interaction, used by
 * higher-level classes.
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
  /**
 * @brief Constructs a contact pair with an initialized slave contact point.
 *
 * The master contact point is typically assigned later by a contact search/projection
 * procedure once a suitable master entity has been found.
 */
  ContactPair(std::unique_ptr<ContactPoint> slave);
  ~ContactPair(){;}
  //
  /**
   * @brief Returns true if the pair is currently considered in contact.
   *
   * The default implementation delegates the decision to the master contact point
   * (if available). If no master is assigned, the pair cannot be in contact.
   */
  virtual bool inContact() {
    if (!master) return false;
    return master->inContact();
  }
  //
  /**
   * @brief Computes the contact interpolation matrix (N-matrix) for this pair.
   *
   * Used to map nodal/DOF quantities to the contact point (e.g., for assembling
   * contact contributions). The exact definition depends on the underlying
   * contact point types.
   *
   * @param answer Output matrix.
   */
  virtual void computeNmatrix(FloatMatrix &answer);
  /**
   * @brief Computes derivatives of shape functions w.r.t. local surface parameters.
   *
   * Typically used for tangent construction, curvature evaluation, and consistent
   * linearization of contact constraints.
   *
   * @param answer Output list of derivative matrices (one per local parameter).
   */
  virtual void compute_dNdxi_matrices(std::vector<FloatMatrix> &answer);
  /**
   * @brief Computes curvature-related quantities of the contact surface at the contact point.
   *
   * @param G     Output curvature matrix/tensor representation.
   * @param tStep Current time step.
   */
  virtual void computeCurvature(FloatMatrix &G, TimeStep *tStep);
  //
  /**
   * @brief Returns the master contact point (non-owning pointer).
   */
  ContactPoint *giveMasterContactPoint() {return master.get();}
  /**
 * @brief Returns the slave contact point (non-owning pointer).
 */
  ContactPoint *giveSlaveContactPoint()  {return slave.get();}
  /**
   * @brief Builds a location array for the pair based on selected DOFs.
   *
   * Used during assembly of contact contributions into the global system.
   *
   * @param dofs DOF mask/IDs to include.
   * @param loc  Output location array.
   * @param ns   Numbering scheme.
   */
  void giveLocationArray(const IntArray &dofs, IntArray &loc, const UnknownNumberingScheme &ns) const;
  /**
 * @brief Returns the current normal gap (signed separation) of the pair.
 */
  double giveNormalGap() { return normal_gap;}
  /**
 * @brief Sets the current normal gap (signed separation) of the pair.
 */
  void setNormalGap(double ng){this->normal_gap = ng;}
  //
  /**
   * @brief Returns the current unit normal vector associated with the contact configuration.
   */
  const FloatArray &giveNormalVector() const {return normalVector;}
  /**
   * @brief Returns the unit normal vector from the previous stored state.
   */
  const FloatArray &givePreviousNormalVector() const {return previousNormalVector;}
  /**
   * @brief Returns the i-th current tangent vector at the contact point.
   *
   * @param i Tangent index.
   */
  const FloatArray &giveTangentVector(int i) const;
  /**
   * @brief Returns the i-th tangent vector from the previous stored state.
   *
   * @param i Tangent index.
   */
  const FloatArray &givePreviousTangentVector(int i) const;
  /**
   * @brief Returns all current tangent vectors (typically one or two, depending on surface dimension).
   */
  std::vector<FloatArray> giveTangentVectors() const;
  /**
   * @brief Returns the slave local coordinates used to parameterize the contact point.
   */
  std::vector<FloatArray> givePreviousTangentVectors() const;
  /**
   * @brief Returns the slave local coordinates used to parameterize the contact point.
   */
  const FloatArray &giveLocalCoordinates() const {return slave->giveLocalCoordinates();}
  /**
   * @brief Initializes contact-point related quantities for the pair.
   *
   * Typically called after assigning master/slave points and before first evaluation.
   */
  void initContactPoint();
  /**
   * @brief Computes displacement at the contact point
   *
   * Useful for gap/penetration evaluation and traction update procedures.
   */
  FloatArray computeContactPointDisplacement() const;
  /**
   * @brief Sets the current normal vector.
   */
  void setNormalVector(const FloatArray &nv)   {this->normalVector = nv;}
  /**
   * @brief Sets the first tangent vector.
   */
  void setTangentVector1(const FloatArray &tv1)   {this->tangentVector1 = tv1;}
  /**
   * @brief Sets the second tangent vector.
   */
  void setTangentVector2(const FloatArray &tv2)   {this->tangentVector2 = tv2;}
  /**
   * @brief Sets a temporary traction vector (e.g., predictor or iteration-local value).
   */  
  void setTempTractionVector(const FloatArray &tv) {this->tempTractionVector = tv;}
  /**
   * @brief Returns the current traction vector associated with the contact constraint.
   */
  const FloatArray &giveTractionVector() const {return tractionVector;}
  /**
   * @brief Assigns a master contact point for this pair.
   */
  void setMasterContactPoint(std::unique_ptr<ContactPoint> m) {master = std::move(m);}
  /**
   * @brief Assigns a slave contact point for this pair.
   */
  void setSlaveContactPoint(std::unique_ptr<ContactPoint> s)  {slave = std::move(s);}
  /**
   * @brief Updates internal state of the contact pair for the given time step.
   *
   * Typical responsibilities include storing previous normals/tangents, updating
   * contact kinematics, and synchronizing traction/gap-related data.
   *
   * @param tStep Current time step.
   */
  virtual void updateYourself(TimeStep *tStep);
  /**
   * @brief Computes a vector quantity for the pair in a given value mode.
   *
   * This is a convenience wrapper to retrieve quantities such as displacements
   * or other state vectors at the involved contact points.
   *
   * @param u      Requested value mode/quantity selector.
   * @param tStep  Current time step.
   * @param answer Output vector.
   */
  void computeVectorOf(ValueModeType u, TimeStep *tStep, FloatArray &answer);
  /**
   * @brief Computes an axis-aligned bounding box (AABB) for the slave side.
   *
   * Typically used to accelerate contact search and broad-phase collision tests.
   *
   * @return Bounding box enclosing the slave contact entity.
   */
  AABB computeSlaveAABB();
};
 
} // end namespace oofem
#endif //contactpair_h
