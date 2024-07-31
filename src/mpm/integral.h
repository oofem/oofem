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
#ifndef integral_h
#define integral_h

#include "mpm.h"
#include "set.h"
#include "intarray.h"
#include "domain.h"
#include "sparsemtrx.h"

namespace oofem {
    /**
    * @brief Class representing weak form integral
    * 
    */
   class Integral {
    public:

        Set &set;
        const Term &term;
        Domain *domain;
        /// @brief Constructor, creates an integral of given term over entities in given set 
        /// @param d
        /// @param s 
        /// @param t 
        Integral (Domain* d, Set& s, const Term& t ) : set(s), term(t) {
            this->domain = d;
        }   
        /// @brief Initialize the integal domain 
        void initialize() {
            for (auto i: this->set.giveElementList()) { // loop over elements
                // introduce necessary dofs and set-up integration rules
                this->term.initializeCell(*(domain->giveElement(i)));
            }
        } 
        // evaluate term contribution to weak form on given cell at given point 
        void assemble_lhs (SparseMtrx& dest, const UnknownNumberingScheme &s, TimeStep* tStep) const {
            FloatMatrix contrib;
            IntArray locr, locc;

            for (auto i: this->set.giveElementList()) { // loop over elements
                MPElement *e = dynamic_cast<MPElement*>(this->domain->giveElement(i));
                if (e) {
                    this->getElementTermCodeNumbers(locr, locc, e, this->term, s);
                    // determine integration rule (this has to be set up on element, as we need to track history there)
                    // the IR is created by term.initialize, we just need mechanism to get this IR
                    // specific terms can have specific integration requirements (reduced integration, etc)
                    // at the same time same rules should be shared between terms->
                    // ->need to querry/store/identify element IR based on NIP.
                    IntegrationRule* ir =  this->term.giveElementIntegrationRule(e);
                    e->integrateTerm_dw(contrib, this->term, ir, tStep); // @todo IR 
                    // assemble
                    dest.assemble (locr, locc, contrib);
                }
            }
        }
        // evaluate contribution (all vars known) on given cell
        void assemble_rhs (FloatArray& dest, const UnknownNumberingScheme &s, TimeStep* tstep) {
            FloatArray contrib;
            IntArray locr, locc;

            for (auto i: this->set.giveElementList()) { // loop over elements
                MPElement *e = dynamic_cast<MPElement*>(this->domain->giveElement(i));
                if (e) {
                    this->getElementTermCodeNumbers(locr, locc, e, this->term, s);
                    // determine integration rule (this has to be set up on element, as we need to track history there)
                    // the IR is created by term.initialize, we just need mechanism to get this IR
                    // specific terms can have specific integration requirements (reduced integration, etc)
                    // at the same time same rules should be shared between terms->
                    // ->need to querry/store/identify element IR based on NIP.
                    IntegrationRule* ir =  this->term.giveElementIntegrationRule(e);
                    e->integrateTerm_c(contrib, this->term, ir, tstep); // @todo IR 
                    // assemble
                    dest.assemble (contrib, locr);
                }
            }
           

        }
        void getElementTermCodeNumbers (IntArray &locr, IntArray &locc, Element* e, const Term& t, const UnknownNumberingScheme &s) const {
            IntArray nodes, internalDofMans, loc;
            locr.resize(0);
            locc.resize(0);
            // term.field and its interpolation determines row code numbers
            // from interpolation get node and internalDofMan lists
            t.field.interpolation.giveCellDofMans(nodes, internalDofMans, e);
            // loop over dof managers to get code numbers
            for (int i: nodes) {
                // get mode dofID mask
                e->giveDofManager(i)->giveLocationArray(t.field.getDofManDofIDs(), loc, s);
                locr.followedBy(loc);
            }
            for (int i: internalDofMans) {
                e->giveInternalDofManager(i)->giveLocationArray(t.field.getDofManDofIDs(), loc, s);
                locr.followedBy(loc);
            }

            // term.testField and its interpolation determines column code numbers
            // from interpolation get node and internalDofMan lists
            t.testField.interpolation.giveCellDofMans(nodes, internalDofMans, e);
            // loop over dof managers to get code numbers
            for (int i: nodes) {
                // get mode dofID mask
                e->giveDofManager(i)->giveLocationArray(t.field.getDofManDofIDs(), loc, s);
                locc.followedBy(loc);
            }
            for (int i: internalDofMans) {


                e->giveInternalDofManager(i)->giveLocationArray(t.field.getDofManDofIDs(), loc, s);
                locc.followedBy(loc);
            }
        }
   };

} // namespace oofem

#endif // integral_h