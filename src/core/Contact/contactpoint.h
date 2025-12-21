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

#ifndef contactpoint_h
#define contactpoint_h

#include "gausspoint.h"
#include "fecontactsurface.h"
#include "feinterpol.h"

namespace oofem {
class IntArray;
class FloatArray;

/**
 * @brief Represents a discrete contact point used in contact mechanics formulations.
 *
 * The ContactPoint class encapsulates all information associated with a single contact
 * point on a potential contact surface or boundary. It stores references to the underlying
 * finite element, local and global coordinates of the contact point, and the associated
 * degrees of freedom. The class serves as a geometric and kinematic representation of
 * contact locations used in master–slave contact algorithms.
 *
 * ContactPoint provides functionality for:
 * - Evaluating and storing the spatial position of the contact point,
 * - Accessing displacement, velocity, or other field values in different value modes,
 * - Managing local coordinate systems and surface normals associated with the contact point,
 * - Supporting projection and update operations required during contact detection and
 *   enforcement.
 *
 * The class does not enforce contact constraints by itself; instead, it acts as a
 * low-level building block used by ContactPair and higher-level contact algorithms
 * to assemble contact contributions to the global system of equations.
 */


class ContactPoint
{
protected:
  int surface_dimension;
public:
  ContactPoint() {}
  ~ContactPoint(){;}
  /**
   * @brief Computes the interpolation matrix (N-matrix) for this contact point.
   *
   * The N-matrix maps element/nodal DOFs to values at the contact point location.
   *
   * @param answer Output matrix.
   */
  virtual void computeNmatrix(FloatMatrix &answer) = 0;
  /**
   * @brief Computes derivative of shape functions w.r.t. the i-th local surface parameter.
   *
   * @param Bs Output derivative matrix.
   * @param i  Local coordinate index.
   */
  virtual void compute_dNdxi_matrix(FloatMatrix &Bs, int i) = 0;
  /**
   * @brief Computes curvature-related surface quantities at this contact point.
   *
   * @param G      Output curvature matrix/tensor representation.
   * @param normal Current surface normal used for curvature evaluation.
   * @param tStep  Current time step.
   */  
  virtual void computeCurvature(FloatMatrix &G, const FloatArray &normal, TimeStep *tStep) = 0;
  /**
   * @brief Returns the local (parametric) coordinates of the contact point.
   *
   * The returned reference must remain valid while the ContactPoint exists.
   *
   * @return Local coordinates (e.g., ξ,η on a face/segment).
   */
  virtual const FloatArray &giveLocalCoordinates() = 0;
  /**
   * @brief Returns the current global coordinates of the contact point.
   *
   * @return Global position vector.
   */
  virtual FloatArray giveGlobalCoordinates() = 0;
  /**
   * @brief Builds a location array for assembling quantities related to this contact point.
   *
   * @param locationArray Output location array.
   * @param dofIDArry     Requested DOF IDs/mask.
   * @param s             Numbering scheme.
   * @return True if the location array could be built, false otherwise.
   */
  virtual bool giveLocationArray(IntArray &locationArray,const IntArray &dofIDArry, 
				 const UnknownNumberingScheme &s) const = 0;
  /**
   * @brief Extracts the unknown vector associated with this contact point.
   *
   * Used to obtain values (e.g. displacements) corresponding to a given DOF mask
   * and value mode, with optional padding.
   *
   * @param answer  Output vector.
   * @param dofMask DOF IDs/mask.
   * @param mode    Value mode (current/incremental/etc.).
   * @param tStep   Current time step.
   * @param padding If true, pads missing entries to match requested layout.
   */
  virtual void giveUnknownVector(FloatArray &answer, const IntArray &dofMask, ValueModeType mode, TimeStep *tStep, bool padding = false) = 0;
  /**
   * @brief Returns the surface normal vector at the contact point.
   *
   * @return Normal vector.
   */
  virtual FloatArray giveNormalVector() = 0;
  /**
   * @brief Called to update internal state of the contact point for the time step.
   *
   * Default implementation is empty; subclasses may cache geometry/kinematics.
   *
   * @param tStep Current time step.
   */
  virtual void updateYourself(TimeStep *tStep){;}
  /**
   * @brief Initializes the contact point.
   *
   * Default implementation is empty; subclasses may set initial
   * coordinates, etc.
   */
  virtual void init(){;}
  /**
   * @brief Returns whether this contact point is currently in contact.
   *
   * @return True if in contact, false otherwise.
   */
  virtual bool inContact() = 0;
  /**
   * @brief Computes a vector quantity of the contact point for a given value mode.
   *
   * Typically used to retrieve displacement/velocity/etc. at the point.
   *
   * @param u      Value mode/quantity selector.
   * @param tStep  Current time step.
   * @param answer Output vector.
   */
  virtual void computeVectorOf(ValueModeType u, TimeStep *tStep, FloatArray &answer) = 0;
  /**
   * @brief Returns updated coordinates of the contact point for the given time step.
   *
   * @param coords Output coordinates.
   * @param tStep  Current time step.
   */
  virtual void giveUpdatedCoordinates(FloatArray &coords, TimeStep* tStep) = 0;
  /**
   * @brief Returns the surface dimension associated with this contact point.
   *
   * For example: 1 for a segment in 2D, 2 for a face in 3D.
   */
  int giveSurfaceDimension(){return surface_dimension;}

};


class FEContactPoint : public ContactPoint
{
protected:
  int contactElementId;
  std::unique_ptr<FEContactSurface> contactSurface;

public:
  FEContactPoint(FEContactSurface *cs, int ceId, int sd) : ContactPoint(), contactElementId(ceId), contactSurface(std::move(cs)){this->surface_dimension  = sd;}
  ~FEContactPoint(){;}
  //
  void computeNmatrix(FloatMatrix &answer) override;
  //  
  void compute_dNdxi_matrix(FloatMatrix &Bs, int i) override;
  //
  FloatArray giveNormalVector() override;
  //
  void computeVectorOf(ValueModeType mode, TimeStep *tStep, FloatArray &answer) override;
  //  
  void computeCurvature(FloatMatrix &G, const FloatArray &normal, TimeStep *tStep) override;
  //
  bool giveLocationArray(IntArray &locationArray,const IntArray &dofIDArry, 
			 const UnknownNumberingScheme &s) const override;
  //
  void giveUnknownVector(FloatArray &answer, const IntArray &dofMask, ValueModeType mode, TimeStep *tStep, bool padding = false) override;
  void giveUpdatedCoordinates(FloatArray &coords, TimeStep* tStep) override;   
  ///////////////////////////////////////////////////////////////////////////////////////
  bool inContact() override {return(contactElementId < 0 ? false : true);}
  //
  const FloatArray &giveLocalCoordinates() override = 0;
  FloatArray giveGlobalCoordinates() override = 0;
  //
   FEInterpolation* giveInterpolation();
  int giveContactElementId(){return contactElementId;}
  void setContactElementId(int ceId){contactElementId = ceId;}
  /////////////////////////////////////////////////////////////////////////////////////////
};

  
class FEContactPoint_Slave : public FEContactPoint
{
protected:
  std::unique_ptr<GaussPoint> slave_point;
  
public:
  FEContactPoint_Slave(FEContactSurface *cs, int ceId, int sd,GaussPoint *gp)  : FEContactPoint(cs, ceId, sd),slave_point(std::move(gp)){;}
  ~FEContactPoint_Slave(){;}
  //
  const FloatArray &giveLocalCoordinates() override {return slave_point->giveNaturalCoordinates();}
  FloatArray giveGlobalCoordinates() override {return slave_point->giveGlobalCoordinates();}
};

  
  
class FEContactPoint_Master : public FEContactPoint
{
protected:
  FloatArray localCoordinates;  
  
public:
  FEContactPoint_Master(FEContactSurface *cs, int ceId, int sd, FloatArray lc)  : FEContactPoint(cs, ceId, sd),localCoordinates(lc){;}
  ~FEContactPoint_Master(){;}
  //
  const FloatArray &giveLocalCoordinates() override {return this->localCoordinates;}
  FloatArray giveGlobalCoordinates() override;
};

  
} // end namespace oofem
#endif //contactpoint_h
