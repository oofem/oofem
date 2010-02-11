/* $Header: /home/cvs/bp/oofem/sm/src/plinearstatic.C,v 1.6.4.1 2004/04/05 15:19:47 bp Exp $ */
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
#ifdef __PARALLEL_MODE

#include "plinearstatic.h"
#include "fetisolver.h"
#include "usrdefsub.h"

namespace oofem {

NumericalMethod *PLinearStatic :: giveNumericalMethod(TimeStep *atTime)
// only one has reason for LinearStatic
//     - SolutionOfLinearEquations

{
    if ( ( solverType == ST_Petsc ) || ( solverType == ST_Feti ) ) {
        if ( nMethod ) {
            return nMethod;
        }

        nMethod = CreateUsrDefSparseLinSolver(solverType, 1, this->giveDomain(1), this);
    }

    if ( nMethod == NULL ) {
        _error("giveNumericalMethod: linear solver creation failed (unknown type or no parallel support)");
    }

    return nMethod;
}

IRResultType
PLinearStatic :: initializeFrom(InputRecord *ir)
{
    LinearStatic :: initializeFrom(ir);
    //this->giveNumericalMethod (giveCurrentStep())-> instanciateFrom (ir);

    return IRRT_OK;
}


void PLinearStatic :: assembleVectorFromDofManagers(FloatArray &answer, TimeStep *tStep, EquationID ut,
                                                    CharType type, ValueModeType mode, 
						    const UnknownNumberingScheme& s,
						    Domain *domain)
{
    /*
     * Assembles characteristic vector of required type into given vector.
     * Overloaded in order to properly handle nodal loading of shared DofManagers.
     * According to general rules, all shared nodes records on all partitions must
     * contain loading. It is therefore necessary to localize loading only on on partition or
     * localize scaled loading on all partitions to guarantee the proper value.
     * The last is used.
     */
    int i;
    IntArray loc;
    FloatArray charVec;
    double scale;

    int nnode = domain->giveNumberOfDofManagers();

    if ( type == NodalLoadVector ) {
        DofManager *node;
        for ( i = 1; i <= nnode; i++ ) {
            node = domain->giveDofManager(i);
            node->giveCompleteLocationArray(loc,s);
            node->computeLoadVectorAt(charVec, tStep, mode);
            if ( node->giveParallelMode() == DofManager_shared ) {
                scale = 1. / ( node->givePartitionList()->giveSize() + 1 );
                charVec.times(scale);
            }

            if ( charVec.giveSize() ) {
                answer.assemble(charVec, loc);
            }
        }
    } else {
      EngngModel :: assembleVectorFromDofManagers(answer, tStep, EID_MomentumBalance, type, mode, s, domain);
    }
}

} // end namespace oofem
#endif
