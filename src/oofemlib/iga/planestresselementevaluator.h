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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef planestresselementevaluator_h
#define planestresselementevaluator_h

#include "structuralelementevaluator.h"

namespace oofem {
/**
 * General purpose Plane stress structural element evaluator.
 */
class PlaneStressStructuralElementEvaluator : public StructuralElementEvaluator
{
public:
    PlaneStressStructuralElementEvaluator() : StructuralElementEvaluator() { }

protected:
    /**
     * Assemble interpolation matrix at given IP.
     * In case of IGAElements, N is assumed to contain only nonzero interpolation functions.
     */
    void computeNMatrixAt(FloatMatrix &answer, GaussPoint *gp);
    /**
     * Assembles the strain-displacement matrix of the receiver at given integration point.
     * In case of IGAElements, B is assumed to contain only contribution from nonzero interpolation functions.
     */
    void computeBMatrixAt(FloatMatrix &answer, GaussPoint *gp);
    double computeVolumeAround(GaussPoint *gp);
    void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const {
        answer.setValues(2, D_u, D_v);
    }
}; // end of PlaneStressStructuralElementEvaluator definition
} // end namespace oofem
#endif //planestresselementevaluator_h
