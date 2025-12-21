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

#ifndef contactbc_h
#define contactbc_h

#include "activebc.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "feinterpol.h"
#include "contactpair.h"
#include "contactsearch.h"
#include "vtkbaseexportmodule.h"

///@name Input fields for general element.
//@{
//@}

namespace oofem {
class TimeStep;
class GaussPoint;
class IntArray;
class UnknownNumberingScheme;
class SparseMtrx;
class ContactPair;
class ContactSearchAlgorithm;

/**
 * @brief Base class for contact boundary conditions.
 *
 * The ContactBC class represents boundary conditions associated with contact
 * interactions in a finite element analysis. It serves as a common base for
 * concrete contact boundary condition implementations responsible for
 * enforcing contact constraints and assembling their contributions to the
 * global system of equations.
 *
 * Contact boundary conditions typically operate on contact pairs identified
 * by a contact search algorithm and may enforce contact using different
 * formulations, such as penalty methods, Lagrange multipliers, or augmented
 * Lagrangian approaches.
 *
 * This class does not define a specific contact enforcement strategy; instead,
 * it provides a unified interface for contact-related boundary conditions used
 * by the engineering model.
 */

  class OOFEM_EXPORT ContactBoundaryCondition: public ActiveBoundaryCondition
{
protected:
    /// contactSearchAlgorithm
  std::unique_ptr<ContactSearchAlgorithm> contactSearchAlgorithm;
    int updateEachNthIter = 1;
public:
    /**
     * Constructor. Creates an element with number n belonging to domain aDomain.
     * @param n Element's number
     * @param aDomain Pointer to the domain to which element belongs.
     */
    ContactBoundaryCondition(int n, Domain * d) : ActiveBoundaryCondition(n, d){;}
    /// Virtual destructor.
    virtual ~ContactBoundaryCondition(){;}
    
    /**
     * @brief Returns the current list of detected contact pairs (modifiable).
     *
     * The pairs are provided by the configured ContactSearchAlgorithm instance.
     */
    std::vector<std::unique_ptr<ContactPair>>& getContactPairs() {
      return contactSearchAlgorithm->getContactPairs();
    };
    const std::vector<std::unique_ptr<ContactPair>>& getContactPairs() const {
        return contactSearchAlgorithm->getContactPairs();
    };
    
    /**
     * @brief Called at the beginning of a new nonlinear iteration.
     *
     * Implementations typically decide whether contact pairs must be updated
     * (e.g., every N-th iteration) and refresh contact kinematics if needed.
     *
     * @param tStep Current time step.
     * @param iter  Iteration counter.
     */
    virtual void initForNewIteration(TimeStep *tStep, int iter);
    /**
     * @brief Assembles contact contributions to the global matrix.
     *
     * @param answer Global sparse matrix to assemble into.
     * @param tStep  Current time step.
     * @param type   Assembly type selector (as used by the framework).
     * @param r_s    Row numbering scheme.
     * @param c_s    Column numbering scheme.
     * @param scale  Scaling factor.
     * @param lock   Optional lock for thread-safe assembly.
     */
    void assemble(SparseMtrx &answer, TimeStep *tStep, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale, void *lock) override;
    /**
     * @brief Assembles contact contributions to the global right-hand side vector.
     *
     * @param answer Output global vector.
     * @param tStep  Current time step.
     * @param type   Assembly type selector.
     * @param mode   Value mode.
     * @param s      Numbering scheme.
     * @param eNorms Optional element norms.
     * @param lock   Optional lock for thread-safe assembly.
     */
    virtual void assembleVector(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode, const UnknownNumberingScheme &s, FloatArray *eNorms, void *lock) override;
    /**
     * @brief Post-initialization hook called after construction and input parsing.
     *
     * Typical responsibilities include creating/configuring the contact search algorithm
     * and initializing contact data structures.
     */
    void postInitialize() override;
    /**
     * @brief Assembles extrapolated (predictor) contact forces.
     *
     * Used e.g. for output or predictor-corrector schemes.
     */
    void assembleExtrapolatedForces(FloatArray &answer, TimeStep *tStep);
 private:
    /**
     * @brief Builds a location array for a given contact pair and numbering scheme.
     *
     * @param loc Output location array.
     * @param ns  Numbering scheme.
     * @param cp  Contact pair.
     */
    void giveLocationArray(IntArray &loc, const UnknownNumberingScheme &ns, const ContactPair *cp) const;
    /**
     * @brief Builds row/column location arrays for contact assembly (pure virtual).
     *
     * This defines the DOF layout used by the concrete contact boundary condition.
     *
     * @param rows Output row location arrays.
     * @param cols Output column location arrays.
     * @param r_s  Row numbering scheme.
     * @param c_s  Column numbering scheme.
     */
    void giveLocationArrays(std::vector< IntArray > &rows, std::vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) override = 0;
    /**
     * @brief Computes a consistent tangent contribution for one contact pair (pure virtual).
     *
     * Concrete contact formulations (penalty, AL, LM, friction) provide their own
     * linearization and fill @p answer accordingly.
     *
     * @param answer Output tangent matrix contribution.
     * @param cp     Contact pair.
     * @param tStep  Current time step.
     */
    virtual void  computeTangentFromContact(FloatMatrix &answer, ContactPair *cp, TimeStep *tStep) = 0;
    /**
     * @brief Computes internal force contribution for one contact pair (pure virtual).
     *
     * Concrete contact formulations compute and assemble the corresponding force vector.
     *
     * @param answer Output internal force contribution.
     * @param cp     Contact pair.
     * @param tStep  Current time step.
     */
    virtual void computeInternalForcesFromContact(FloatArray &answer, ContactPair *cp, TimeStep *tStep) = 0;
    /**
     * @brief Returns the configured contact search algorithm
     */
    virtual  ContactSearchAlgorithm* giveContactSearchAlgorithm(){return contactSearchAlgorithm.get();}
    /**
     * @brief Creates and configures the contact search algorithm (pure virtual).
     *
     * Concrete boundary conditions decide which search algorithm to use and how
     * it is parameterized.
     */
    virtual  void setupContactSearchAlgorithm() = 0;
    /**
     * @brief Updates internal BC state for the given time step.
     *
     * Typically triggers contact search updates and stores time-step dependent data.
     */
    void updateYourself(TimeStep *tStep) override;
 public:
    /**
     * @brief Exports contact-related data (e.g., for VTK output).
     *
     * @param vtkPieces Output regions/pieces.
     * @param shift     Coordinate shift (for periodicity or visualization).
     * @param tStep     Current time step.
     */
    void giveExportData(std::vector< ExportRegion> &vtkPieces, FloatArray shift, TimeStep *tStep );



};


} // end namespace oofem
#endif //contactbc_h
