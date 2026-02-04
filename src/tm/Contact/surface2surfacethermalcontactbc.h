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
 *               Copyright (C) 1993 - 2023   Borek Patzak
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
#ifndef surface2surfacethermalcontactbc_h
#define surface2surfacethermalcontactbc_h


#include "Contact/contactbc.h"
#include "Contact/contactpoint.h"
#include "thermalcontactbc.h"
#include "ContactSurface/thermalfecontactsurface.h"

///@name Input fields for _IFT_ContactElement
//@{
#define _IFT_Surface2SurfaceThermalContactBoundaryCondition_Name "s2sthermalcontact"

#define _IFT_Surface2SurfaceThermalContactBoundaryCondition_slaveSurfaceNum "slavesurface"
#define _IFT_Surface2SurfaceThermalContactBoundaryCondition_masterSurfaceNum "mastersurface"


//@}

namespace oofem {
class Domain;
class SparseMtrx;
class TimeStep;
class DofManager;
class GaussPoint;
class UnknownNumberingScheme;
class FloatMatrix;
class IntegrationRule;
class Node;


/**
 * @class Surface2SurfaceThermalContactBoundaryCondition
 * @brief Boundary condition for surface-to-surface thermal contact.
 *
 * This class implements thermal contact between two discretized surfaces
 * (master and slave), defined by sets of nodes and segments. It maintains
 * corresponding contact pairs and evaluates heat transfer across the
 * interface.
 *
 * The size of the thermal contact interface is user-defined, and all
 * operations are performed by iterating over all detected contact pairs.
 *
 * Responsibilities:
 * - Management of master and slave surface sets
 * - Setup and execution of the contact search algorithm
 * - Assembly of thermal contact contributions to:
 *   - global stiffness (tangent) matrix
 *   - global internal force vector
 * - Construction of location arrays for DOF coupling induced by contact
 *
 * This class represents a *geometric contact layer*, while the actual
 * thermal constitutive behavior (gap-dependent conductance) is handled
 * by the base class ThermalContactBoundaryCondition.
 */


class OOFEM_EXPORT Surface2SurfaceThermalContactBoundaryCondition : public ThermalContactBoundaryCondition
{
private:
  /**
   * Input identifier of the master surface.
   *
   * Refers to a surface definition in the input file.
   */
  int masterSurfaceNumber;
  /**
   * Input identifier of the slave surface.
   *
   * Refers to a surface definition in the input file.
   */
  int slaveSurfaceNumber;
  //    std::shared_ptr<ThermalFEContactSurface> contactSurface;
  /**
   * Pointer to the slave thermal contact surface object.
   *
   * Contains the discretized representation of the slave surface
   * used for contact detection and thermal interaction.
   */
  ThermalFEContactSurface *slaveContactSurface;
  /**
   * Pointer to the master thermal contact surface object.
   *
   * Contains the discretized representation of the master surface
   * used for contact detection and thermal interaction.
   */
  ThermalFEContactSurface *masterContactSurface;
  /**
   * Flag indicating whether contact pairs have been initialized.
   *
   * Used to avoid repeated initialization of contact pairs
   * during nonlinear iterations.
   */
  bool initializedPairs = false;
  
 public:

    /// Constructor.
    Surface2SurfaceThermalContactBoundaryCondition(int n, Domain *d) : ThermalContactBoundaryCondition(n, d) { }
    /// Destructor.
    virtual ~Surface2SurfaceThermalContactBoundaryCondition() {};

    void initializeFrom(const std::shared_ptr<InputRecord> &ir) override;
    virtual const char *giveClassName() const override { return "Surface2SurfaceThermalContactBoundaryCondition"; }
    virtual const char *giveInputRecordName() const override { return _IFT_Surface2SurfaceThermalContactBoundaryCondition_Name; }
    /**
     * Sets up the contact search algorithm.
     *
     * Builds geometric search structures and initializes
     * contact detection between master and slave surfaces.
     */
    void setupContactSearchAlgorithm() override;
 private:
     /**
     * Constructs location arrays for coupled DOFs.
     *
     * This method defines the row and column DOF indices that interact
     * through thermal contact, enabling correct assembly into the
     * global system matrices.
     *
     * @param rows Row index arrays
     * @param cols Column index arrays
     * @param type Type of matrix contribution (e.g. stiffness, mass, etc.)
     * @param r_s Row numbering scheme
     * @param c_s Column numbering scheme
     */
    void giveLocationArrays(std::vector< IntArray > &rows, std::vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) override;
    /**
     * Initializes data structures for a new nonlinear iteration.
     *
     * Used to update contact pairs
     * @param tStep Current time step
     * @param iter Nonlinear iteration counter
     */
    void initForNewIteration(TimeStep *tStep, int iter) override;
};
} // end namespace oofem
#endif // surface2surfacethermalcontact_h
