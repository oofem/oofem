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

#include "dynamicrelaxationsolver.h"
#include "timestep.h"
#include "classfactory.h"
#include "exportmodulemanager.h"
#include "engngm.h"
#include "domain.h"
#include "dofmanager.h"
#include "element.h"
#include "unknownnumberingscheme.h"

namespace oofem {

#define ERROR_NORM_SMALL_NUM 1.e-6
#define MAX_REL_ERROR_BOUND 1.e20

REGISTER_SparseNonLinearSystemNM(DynamicRelaxationSolver)

DynamicRelaxationSolver :: DynamicRelaxationSolver(Domain *d, EngngModel *m) : NRSolver(d, m)
{
}


IRResultType
DynamicRelaxationSolver :: initializeFrom(InputRecord *ir)
{
    return NRSolver :: initializeFrom(ir);
}


NM_Status
DynamicRelaxationSolver :: solve(SparseMtrx &k, FloatArray &R, FloatArray *R0, FloatArray *iR,
                  FloatArray &X, FloatArray &dX, FloatArray &F,
                  const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                  int &nite, TimeStep *tStep)

{
    // residual, iteration increment of solution, total external force
    FloatArray rhs, ddX, RT, X_0, X_n, X_n1, M;

    double RRT;
    int neq = X.giveSize();
    NM_Status status = NM_None;
    bool converged, errorOutOfRangeFlag;
    ParallelContext *parallel_context = engngModel->giveParallelContext( this->domain->giveNumber() );
    if ( engngModel->giveProblemScale() == macroScale ) {
        OOFEM_LOG_INFO("DRSolver: Iteration");
        if ( rtolf.at(1) > 0.0 ) {
            OOFEM_LOG_INFO(" ForceError");
        }
        if ( rtold.at(1) > 0.0 ) {
            OOFEM_LOG_INFO(" DisplError");
        }
        OOFEM_LOG_INFO("\n----------------------------------------------------------------------------\n");
    }

    // compute total load R = R+R0
    l = 1.0;
    RT = R;
    if ( R0 ) {
        RT.add(* R0);
    }

    RRT = parallel_context->localNorm(RT);
    RRT *= RRT;

    ddX.resize(neq);
    ddX.zero();

    X_0 = X;
    X_n = X_0;
    X_n1 = X_0;

    // Compute the mass "matrix" (lumped, only storing the diagonal)
    M.resize(neq);
    M.zero();
    engngModel->assembleVector(M, tStep, LumpedMassVectorAssembler(), VM_Total, EModelDefaultEquationNumbering(), domain);

    double Le = -1.0;
    for ( auto &elem : domain->giveElements() ) {
        double size = elem->computeMeanSize();
        if ( Le < 0 || Le >= size ) {
            Le = size;
        }
    }

    for ( nite = 0; ; ++nite ) {
        // Compute the residual
        engngModel->updateComponent(tStep, InternalRhs, domain);
        rhs.beDifferenceOf(RT, F);

        converged = this->checkConvergence(RT, F, rhs, ddX, X, RRT, internalForcesEBENorm, nite, errorOutOfRangeFlag);
        if ( errorOutOfRangeFlag ) {
            status = NM_NoSuccess;
            dX.zero();
            X.zero();
            OOFEM_WARNING("Divergence reached after %d iterations", nite);
            break;
        } else if ( converged && ( nite >= minIterations ) ) {
            break;
        } else if ( nite >= nsmax ) {
            OOFEM_LOG_DEBUG("Maximum number of iterations reached\n");
            break;
        }

        double density = 1.;
        double lambda = 210e9;
        double mu = 210e9;
        double c = sqrt((lambda + 2*mu) / density);
        double dt = 0.25 * Le / c;
        double alpha = 0.1 / dt;
        printf("dt = %e\n", dt);
        for ( int j = 0; j < neq; ++j ) {
            //M * x'' + C*x' * dt = rhs * dt*dt
            X[j] = rhs[j] * dt * dt / M[j] - ( -2*X_n1[j] + X_n[j] ) - alpha * (X_n1[j] - X_n[j]) * dt;
        }
        X_n = X_n1;
        X_n1 = X;

        dX.beDifferenceOf(X, X_0);
        tStep->incrementStateCounter(); // update solution state counter
        tStep->incrementSubStepNumber();
    }

    status |= NM_Success;
    solved = 1;

    return status;
}

} // end namespace oofem

