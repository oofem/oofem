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
#ifndef termlibrary_h
#define termlibrary_h

#include "mpm.h"

namespace oofem {
/**
 * @brief A Linear momentum balance equation term ($w^TB^T\sigma(u)$)
 * 
 */
class wTBTSigTerm : public Term {
    protected:
    public:
    wTBTSigTerm (Variable& unknownField, Variable &testField) ;

    /**
     * @brief Evaluates the linearization of $B^T\sigma(u)$, i.e. $B^TDBu$
     * 
     * @param answer 
     * @param e 
     * @param coords 
     */
    void evaluate_dw (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) override;
    /**
     * @brief Evaluates Internal forces vector, i.e. $b^T\sigma(u)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate_c (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) override;
    void getDimensions_dw(Element& cell) override;
    void initializeCell(Element& cell) override;

    protected:
    /**
     * @brief Evaluates B matrix; i.e. $LN$ where $L$ is operator matrix and $N$ is interpolation matrix of unknowns
     * 
     * @param answer B matrix
     * @param v 
     * @param interpol 
     * @param cell 
     * @param coords 
     */
    void grad(FloatMatrix& answer, Variable &v, FEInterpolation& interpol, Element& cell, const FloatArray& coords) ;
    
};

} // end namespace oofem
#endif // termlibrary_h