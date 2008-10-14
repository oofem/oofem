/* $Header: /home/cvs/bp/oofem/sm/src/adaptnlinearstatic.h,v 1.7 2003/04/06 14:08:30 bp Exp $ */
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
// Class AdaptiveNonLinearStatic
//

#ifndef adaptnlinearstatic_h
#define adaptnlinearstatic_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "nlinearstatic.h"
#include "sparsemtrx.h"
#include "errorestimator.h"
#include "meshpackagetype.h"

class AdaptiveNonLinearStatic : public NonLinearStatic
{
    /*
     * This class implements Adaptive Non - LinearStatic Engineering problem.
     * DESCRIPTION:
     * Solution of this problem is performed  as a series of increments (loading or displacement).
     * At start of Each increment we assemble new tangent stiffness, and iteratively trying
     * to fullfill balance of external and real internal forces
     * at end of load step (see numerical method ).
     *
     * The error is estimated at the end of each load increment, and based on reached error,
     * the computation continues, or remeshing happens. The solution is then mapped to
     * new mesh. The solution may be taken from current state or if error is too high then is
     * taken from previous equilibrium state.
     *
     * TASK:
     * Creating Numerical method for solving nonlinear problem.
     * Assembling tangent stiffness matrix
     * Interfacing Numerical method to Elements
     * Managing time  steps
     */

protected:
    FloatArray d2_totalDisplacement, d2_incrementOfDisplacement;
    MeshPackageType meshPackage;
    int equilibrateMappedConfigurationFlag;
    /// Error estimator
    ErrorEstimator *ee;
    /**
     * Array storing the load levels reached in correpsonding timesteps.
     * It is necessary to keep track of this load level history,
     * because after adaptive restart one has to assemble
     * the initial and total load vectors on new mesh ->
     * and for this it is necessary to know the history of loading.
     * The size of this array is equal to numberOfSolutionSteps and
     * should be stored/restored in every context file.
     */
    FloatArray timeStepLoadLevels;



public:
    AdaptiveNonLinearStatic(int i, EngngModel *_master = NULL);
    ~AdaptiveNonLinearStatic();
    // solving
    // void solveYourself ();
    void solveYourselfAt(TimeStep *);
    virtual void               updateYourself(TimeStep *);
    //void terminate (TimeStep *);
    IRResultType initializeFrom(InputRecord *ir);
    /** Service for accessing ErrorEstimator corresponding to particular domain */
    ErrorEstimator *giveDomainErrorEstimator(int n) { return ee; }

    double giveUnknownComponent(EquationID, ValueModeType, TimeStep *, Domain *, Dof *);

    /**
     * Returns the loaad level cooresponding to given solution step number
     */
    double giveTimeStepLoadLevel(int istep);
    /**
     * Initializes the newly generated discretization state acording to previous solution.
     * This process should typically include restoring old solution, instanciating newly
     * generated domain(s) and by mapping procedure.
     */
    virtual int                initializeAdaptive(int stepNumber);
    /**
     * Initializes the receiver state acording to state of given source problem.
     * This process should typically include mapping of source solution, internal variable mapping procedures and
     * optionally restoring new global equilibrium.
     */
    virtual int                initializeAdaptiveFrom(EngngModel *sourceProblem);
    /**
     * Remaps the solution state to newly given domain. This includes mapping of source solution, 
     * internal variable mapping procedures and optionally restoring new global equilibrium.
     * Given domain becomes new domain of receiver.
     */
    int                        adaptiveRemap (Domain *dNew);


    /**
     * Restores the  state of model from output stream. Restores not only the receiver state,
     * but also same function is invoked for all DofManagers and Elements in associated
     * domain. Note that by restoring element  context also contexts of all associated
     * integration points (and material statuses) are restored.
     */
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    void updateDomainLinks();

    /// Returns class name of the receiver.
    virtual const char *giveClassName() const { return "AdaptiveNonLinearStatic"; }
    /// Returns classType id of receiver.
    virtual classType giveClassID() const { return AdaptiveNonLinearStaticClass; }

protected:
    void assembleInitialLoadVector(FloatArray &loadVector, FloatArray &loadVectorOfPrescribed,
                                   AdaptiveNonLinearStatic *sourceProblem, int domainIndx, TimeStep *atTime);
    /*
     * void assembleCurrentTotalLoadVector (FloatArray& loadVector, FloatArray& loadVectorOfPrescribed,
     *                 AdaptiveNonLinearStatic* sourceProblem, int domainIndx, TimeStep* atTime);
     */
};

#endif // adaptnlinearstatic_h
