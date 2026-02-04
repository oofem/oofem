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
 *               Copyright (C) 1993 - 2026   Borek Patzak
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
#ifndef thermalcontactbc_h
#define thermalcontactbc_h


#include "Contact/contactbc.h"
#include "Contact/contactpoint.h"

///@name Input fields for 
//@{
//@}
#define _IFT_ThermalContactBoundaryCondition_gapConductanceFunctionNumber "gapconductivityfunction"

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
 * @file thermalcontactbc.h
 * @brief Declaration of thermal contact boundary condition for OOFEM.
 *
 * This file defines a boundary condition class handling heat transfer
 * across contacting surfaces. The formulation is based on a gap-dependent
 * thermal conductance, which contributes to both the global matrix and rhs vector.
 */

class OOFEM_EXPORT ThermalContactBoundaryCondition : public ContactBoundaryCondition
{
private:
   /**
     * Function number defining the gap-dependent thermal conductance.
     *
     * The corresponding function is expected to be defined in the
     * input file and queried during the evaluation of the thermal
     * contact response.
     */
  int gapConductanceFunctionNumber;
 public:

    /// Constructor.
    ThermalContactBoundaryCondition(int n, Domain *d) : ContactBoundaryCondition(n, d) { }
    /// Destructor.
    virtual ~ThermalContactBoundaryCondition() {};
    void initializeFrom(const std::shared_ptr<InputRecord> &ir) override;
private:
    /**
     * Computes the matrix contribution due to thermal contact.
     *
     * This method assembles the linearization of the thermal contact
     * contribution for a given contact pair and time step.
     *
     * @param answer Resulting tangent matrix contribution
     * @param cp Pointer to the active contact pair
     * @param tStep Current time step
     */
    void  computeTangentFromContact(FloatMatrix &answer, ContactPair *cp, TimeStep *tStep) override;
     /**
     * Computes the rhs vector contribution due to thermal contact.
     *
     * This method evaluates the heat flux across the contact interface
     * and assembles its contribution to the internal force vector.
     *
     * @param answer Resulting internal force contribution
     * @param cp Pointer to the active contact pair
     * @param tStep Current time step
     */
    void computeInternalForcesFromContact(FloatArray &answer, ContactPair *cp, TimeStep *tStep) override;
     /**
     * Evaluates the thermal conductance as a function of the normal gap.
     *
     * This function provides the constitutive relation for thermal
     * contact, typically decreasing with increasing separation.
     *
     * @param gap Normal gap between contacting surfaces
     * @return Thermal conductance value
     */
    virtual double evaluateGapConductance(double gap);

};
} // end namespace oofem
#endif // thermalcontactbc_h
