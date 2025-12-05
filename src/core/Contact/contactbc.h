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
    std::vector<std::unique_ptr<ContactPair>>& getContactPairs() {
        return contactSearchAlgorithm->getContactPairs();
    };
    const std::vector<std::unique_ptr<ContactPair>>& getContactPairs() const {
        return contactSearchAlgorithm->getContactPairs();
    };

    ContactBoundaryCondition(int n, Domain * d) : ActiveBoundaryCondition(n, d){;}
    /// Virtual destructor.
    virtual ~ContactBoundaryCondition(){;}
    //
    virtual void initForNewIteration(TimeStep *tStep, int iter);

    void assemble(SparseMtrx &answer, TimeStep *tStep, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale, void *lock) override;
    //
    virtual void assembleVector(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode, const UnknownNumberingScheme &s, FloatArray *eNorms, void *lock) override;
    //
    void postInitialize() override;
    
    void assembleExtrapolatedForces(FloatArray &answer, TimeStep *tStep);
 private:
    void giveLocationArray(IntArray &loc, const UnknownNumberingScheme &ns, const ContactPair *cp) const;
    //
    void giveLocationArrays(std::vector< IntArray > &rows, std::vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) override = 0;

    virtual void  computeTangentFromContact(FloatMatrix &answer, ContactPair *cp, TimeStep *tStep) = 0;
    virtual void computeInternalForcesFromContact(FloatArray &answer, ContactPair *cp, TimeStep *tStep) = 0;

    virtual  ContactSearchAlgorithm* giveContactSearchAlgorithm(){return contactSearchAlgorithm.get();}
    virtual  void setupContactSearchAlgorithm() = 0;
    void updateYourself(TimeStep *tStep) override;
 public:
    void giveExportData(std::vector< ExportRegion> &vtkPieces, FloatArray shift, TimeStep *tStep );



};


} // end namespace oofem
#endif //contactbc_h
