/* $Header: /home/cvs/bp/oofem/oofemlib/src/stationaryflow.h,v 1.14 2003/04/06 14:08:26 bp Exp $ */
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
// Class StationaryFlow
//

#ifndef stationaryflow_h
#define stationaryflow_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "engngm.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"

namespace oofem {

class StationaryFlow : public EngngModel
{
    /*
     * This class implements StationaryFlow Engineering problem.
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
    SparseMtrx *conductivityMatrix;
    FloatArray loadVector;
    FloatArray fluxVector;
    /// Numerical method used to solve the problem
    SparseLinearSystemNM *nMethod;


public:
    StationaryFlow(int i, EngngModel *_master = NULL) : EngngModel(i, _master), loadVector(), fluxVector()
    { conductivityMatrix = NULL;
      ndomains = 1; }
    ~StationaryFlow()
    { delete  conductivityMatrix; }
    // solving
    //void solveYourself ();
    void solveYourselfAt(TimeStep *);
    //int requiresNewLhs () {return 0;}

    double giveUnknownComponent(EquationID, ValueModeType, TimeStep *, Domain *, Dof *);
    /** Stores receiver state to output stream */
    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /** Restores the receiver state previously written in stream */
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    TimeStep *giveNextStep();
    NumericalMethod *giveNumericalMethod(TimeStep *);
    void               terminate(TimeStep *);
    /** DOF printing routine. Called by DofManagers to print Dof specific part.
     * Dof class provides component printing routines, but emodel is responsible
     * for what will be printed at DOF level.
     * @param stream output stream
     * @param iDof dof to be processed
     * @param atTime solution step
     */
    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);

    // identification
    const char *giveClassName() const { return "StationaryFlow"; }
    classType giveClassID() const { return StationaryFlowClass; }
    int isFlow() { return 1; }
    fMode giveFormulation() { return TL; }
};

} // end namespace oofem
#endif // stationaryflow_h
