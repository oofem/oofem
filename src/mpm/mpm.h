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
#ifndef mpm_h
#define mpm_h

/**
 * Multiphysics module 
 * Classes:
 * - ElementBase(Element) defining geometry
 * - Variable class representing unknown field (or test feld) in a weak psolution. The variable has its interpolation, type (scalar, vector), size.
    When test field, it keeps reference to its primary (unknown) variable. The history parameter dermines how many time steps to remember. 
 * - Term class represnting a term to evaluate on element. Paramaters element(geometry), variables
 * - Element - responsible for defining and performing integration (of terms), assembly of term contributions. 
 */

#include "element.h"
#include "gausspoint.h"
#include "feinterpol.h"
#include "intarray.h"
#include "classfactory.h"


namespace oofem {

class Variable {
    public:
    enum VariableType {
        scalar,
        vector
    };

    enum VariableQuantity {
        Displacement,
        Temperature
    };

    FEInterpolation& interpolation;  
    Variable* dualVar; //? or just bool?
    VariableType type;
    VariableQuantity q;
    int size;
    IntArray dofIDs;

    Variable (FEInterpolation& i, Variable::VariableQuantity q, Variable::VariableType t, int size, Variable* dual = NULL) : 
        interpolation(i), 
        dualVar(dual), 
        q(q) {
        this->type = t;
        this->size = size;
    }

    /// Returns DodIF mask in node
    const IntArray& getDofManDofIDs () {return this->dofIDs;}
};

class Term {
    public:
    Variable& field;
    Variable& testField;

    public:
    Term (Variable& unknownField, Variable &testField) : field(unknownField), testField(testField) {}
    
    // evaluate term contribution to weak form on given cell at given point 
    virtual void evaluate_dw (FloatMatrix& , Element& cell, GaussPoint* gp, TimeStep* tStep) =0;
    // evaluate contribution (all vars known) on given cell
    virtual void evaluate_c (FloatArray&, Element& cell, GaussPoint* gp, TimeStep* tStep)=0;
    virtual void getDimensions_dw(Element& cell) =0;
    virtual void initializeCell(Element& cell) =0;
};


/**
* Element code sample:
* term=Poison(Variable(interpolation, Temperature, 1), Variable(Interpolation, Temperature,1));
* this-> Assemble(integrationRule, term, destination); // where to integrate (volume, surface,edge?)
* 
*/
class MPElement : public Element {
   public:

    MPElement (int n, Domain * aDomain) : 
        Element(n, aDomain)
    {}

    void initialize () {
        // loop over variables and allocate nodal dofs (for unknownFields)
    }

    void integrateTerm_dw (FloatMatrix& answer, Term& term, IntegrationRule* iRule, TimeStep* tstep) {
        // need for integration domain and rule.
        // who should determine integration domain? Element or term? Term is just integrand, not integral
        // so integral type (surface, volume, etc) defined by element ---
        FloatMatrix dw;
        for ( GaussPoint *igp : * iRule ) {
            term.evaluate_dw(dw, *this, igp, tstep);
            dw.times(this->computeVolumeAround(igp));
            answer.add(dw);
        }
    } 

    /// @brief  returns local code number on element to assemble variable/term contribution
    /// @param answer 
    /// @param v 
    virtual void getLocalCodeNumbers (IntArray& answer, Variable& v) =0;
    

    void assembleTermContribution (FloatMatrix& answer, FloatMatrix& contrib, Term& t) {
        IntArray uloc, tloc;
        this->getLocalCodeNumbers(uloc, t.field);
        this->getLocalCodeNumbers(tloc, t.testField);
        answer.assemble(contrib, uloc, tloc);
    }
};




} // end namespace oofem
#endif // mpm_h