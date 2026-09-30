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
#include "contactprojection.h"
#include <vector>

namespace oofem {
class IntArray;
class FloatArray;
class ContactElement;

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
  virtual ~ContactPoint() = default;
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
  virtual void computeSecondBaseVectors(std::vector<std::vector<FloatArray>> &answer, TimeStep *tStep) = 0;
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
  virtual Coordinates giveGlobalCoordinates() = 0;
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
   * @brief Returns the quadrature weight associated with this contact point.
   *
   * Master points created by closest-point projection are not integration points,
   * so they use unit weight. Slave points override this with their Gauss weight.
   */
  virtual double giveIntegrationWeight() { return 1.0; }
  /**
   * @brief Returns the current surface/line measure at this contact point.
   *
   * The returned value excludes the quadrature weight.
   */
  virtual double giveSurfaceMeasure(TimeStep *tStep) = 0;
  /**
   * @brief Returns the reference surface/line measure at this contact point.
   *
   * The returned value excludes the quadrature weight. Penalty contact uses this
   * measure so the area factor is fixed with respect to displacement increments.
   */
  virtual double giveReferenceSurfaceMeasure() = 0;
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
  virtual void giveUpdatedCoordinates(Coordinates &coords, TimeStep* tStep) = 0;
  virtual ContactElement *giveContactElement() { return nullptr; }
  /**
   * @brief Returns the surface dimension associated with this contact point.
   *
   * For example: 1 for a segment in 2D, 2 for a face in 3D.
   */
  int giveSurfaceDimension() const {return surface_dimension;}

};


class FEContactPoint : public ContactPoint
{
protected:
  int contactElementId;
  FEContactSurface *contactSurface = nullptr;
  ContactFeatureType featureType = ContactFeatureType::Surface;
  int featureIndex = 0;

public:
  FEContactPoint(FEContactSurface *cs, int ceId, int sd,
                 ContactFeatureType ft = ContactFeatureType::Surface,
                 int fi = 0)
    : ContactPoint(), contactElementId(ceId), contactSurface(cs),
      featureType(ft), featureIndex(fi)
  {
    this->surface_dimension = sd;
  }
  ~FEContactPoint() override = default;
  //
  void computeNmatrix(FloatMatrix &answer) override;
  //  
  void compute_dNdxi_matrix(FloatMatrix &Bs, int i) override;
  //
  FloatArray giveNormalVector() override;
  double giveSurfaceMeasure(TimeStep *tStep) override;
  double giveReferenceSurfaceMeasure() override;
  //
  void computeVectorOf(ValueModeType mode, TimeStep *tStep, FloatArray &answer) override;
  //  
  void computeCurvature(FloatMatrix &G, const FloatArray &normal, TimeStep *tStep) override;
  void computeSecondBaseVectors(std::vector<std::vector<FloatArray>> &answer, TimeStep *tStep) override;
  //
  bool giveLocationArray(IntArray &locationArray,const IntArray &dofIDArry, 
			 const UnknownNumberingScheme &s) const override;
  //
  void giveUnknownVector(FloatArray &answer, const IntArray &dofMask, ValueModeType mode, TimeStep *tStep, bool padding = false) override;
  void giveUpdatedCoordinates(Coordinates &coords, TimeStep* tStep) override;
  void giveUpdatedCoordinatesOnElement(Coordinates &coords, int elementId,
                                       const FloatArray &localCoords, TimeStep *tStep) const;
  ///////////////////////////////////////////////////////////////////////////////////////
  bool inContact() override {return(contactElementId < 0 ? false : true);}
  //
  const FloatArray &giveLocalCoordinates() override = 0;
  Coordinates giveGlobalCoordinates() override = 0;
  //
  FEInterpolation *giveInterpolation();
  ContactElement *giveContactElement() override;
  ContactElement *giveContactElementOnSurface(int elementId) const;
  int giveContactElementId(){return contactElementId;}
  void setContactElementId(int ceId){contactElementId = ceId;}
  ContactFeatureType giveContactFeatureType() const { return featureType; }
  int giveContactFeatureIndex() const { return featureIndex; }
  /////////////////////////////////////////////////////////////////////////////////////////
};

  
class FEContactPoint_Slave : public FEContactPoint
{
protected:
  GaussPoint *slave_point = nullptr;
  
public:
  FEContactPoint_Slave(FEContactSurface *cs, int ceId, int sd,GaussPoint *gp)  : FEContactPoint(cs, ceId, sd),slave_point(gp){;}
  ~FEContactPoint_Slave() override = default;
  //
  const FloatArray &giveLocalCoordinates() override {return slave_point->giveNaturalCoordinates();}
  Coordinates giveGlobalCoordinates() override {return slave_point->giveGlobalCoordinates();}
  double giveIntegrationWeight() override { return slave_point->giveWeight(); }
  GaussPoint *giveIntegrationPoint() const { return slave_point; }
};

  
  
class FEContactPoint_Master : public FEContactPoint
{
protected:
  FloatArray localCoordinates;  
  
public:
  FEContactPoint_Master(FEContactSurface *cs, int ceId, int sd, FloatArray lc,
                        ContactFeatureType ft = ContactFeatureType::Surface,
                        int fi = 0)
    : FEContactPoint(cs, ceId, sd, ft, fi), localCoordinates(lc) { }
  ~FEContactPoint_Master() override = default;
  //
  const FloatArray &giveLocalCoordinates() override {return this->localCoordinates;}
  Coordinates giveGlobalCoordinates() override;
};

  
} // end namespace oofem
#endif //contactpoint_h
