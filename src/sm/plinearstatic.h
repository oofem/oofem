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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//
// Class PLinearStatic
//
#ifndef plinearstatic_h
#define plinearstatic_h

#ifdef __PARALLEL_MODE

 #ifndef __MAKEDEPEND
  #include <stdio.h>
 #endif
 #include "linearstatic.h"
 #include "sparsemtrx.h"

namespace oofem {
///@todo Remove this class
///@deprecated Using normal linear static should work fine even in parallell.
class PLinearStatic : public LinearStatic
{
    /*
     * This class implements Parallel LinearStatic Engineering problem.
     * Multiple loading works only if linear elastic material (such as isoLE)  is used.
     * (Other non-linear materials kepp load history, so such multiple loading
     * will cause that next step will be assumed as new load increment,
     * not the total new load). Because they always copute real stresses acording
     * to reached strain state, they are not able to respond to linear analysis.
     *
     * DESCRIPTION:
     * Solution of this problem is series of loading cases, maintained as sequence of
     * time-steps. This solution is in form of linear equation system Ax=b
     * TASK:
     * Creating Numerical method for solving Ax=b
     * Interfacing Numerical method to Elements
     * Managing time  steps
     */

protected:

public:
    PLinearStatic(int i, EngngModel *_master = NULL) : LinearStatic(i, _master)
    { }
    ~PLinearStatic()
    { }
    NumericalMethod *giveNumericalMethod(TimeStep *);
    IRResultType initializeFrom(InputRecord *ir);

    /**
     * Assembles characteristic vector of required type into given vector.
     * Overloaded in order to properly handle nodal loading of shared DofManagers.
     * According to general rules, all shared nodes records on all partitions must
     * contain loading. It is therefore necessary to localize loading only on on partition or
     * localize scaled loading on all partitions to guarantee the proper value.
     * The last is used.
     * @param answer assembled vector
     * @param tStep time step, when answer is assembled.
     * @param type characterisctic components of type type are requsted
     * from dofManagers/elements and assembled.
     */
    virtual double assembleVectorFromDofManagers(FloatArray &, TimeStep *, EquationID ut,
                                                 CharType type, ValueModeType mode,
                                                 const UnknownNumberingScheme &s, Domain *domain);
};
} // end namespace oofem
#endif
#endif // plinearstatic_h
