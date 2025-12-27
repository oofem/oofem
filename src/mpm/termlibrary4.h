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
#ifndef termlibrary4_h
#define termlibrary4_h

#include "mpm.h"
#include "boundaryload.h"
#include "termlibrary.h"
#include "field.h"

// file containing various term definitions 

namespace oofem {

 /**
 * @brief Advection problem mass matrix $S=(N)^T N$. 
 */
class NTN : public Term {
    protected:
    public:
    NTN (const Variable *testField, const Variable* unknownField) ;

    /**
     * @brief Evaluates the linearization of term (the lhs contribution)
     * 
     * @param answer 
     * @param e 
     * @param coords 
     */
    void evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override;
    /**
     * @brief Evaluates Internal forces vector, i.e. $w^T(\grad N)^T f(p)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override;
    void getDimensions(Element& cell) const override;
    void initializeCell(Element& cell) const override;
    
};

 /**
 * @brief A advection equation term $T=(\grad N)^T N$. 
 */
class dnTaN : public Term {
    protected:
    FieldPtr velocity;
    public:
    dnTaN (const Variable *testField, const Variable* unknownField, FieldPtr velocity);

    /**
     * @brief Evaluates the linearization of term (the lhs contribution)
     * 
     * @param answer 
     * @param e 
     * @param coords 
     */
    void evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override;
    /**
     * @brief Evaluates Internal forces vector, i.e. $w^T(\grad N)^T f(p)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override;
    void getDimensions(Element& cell) const override;
    void initializeCell(Element& cell) const override;
    
};



} // end namespace oofem
#endif // termlibrary_h
