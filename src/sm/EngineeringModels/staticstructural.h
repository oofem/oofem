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

#ifndef staticstructural_h
#define staticstructural_h

#include "../sm/EngineeringModels/structengngmodel.h"
#include "../sm/EngineeringModels/xfemsolverinterface.h"
#include "sparsenonlinsystemnm.h"
#include "sparsemtrxtype.h"
#include "staggeredsolver.h"

#define _IFT_StaticStructural_Name "staticstructural"
#define _IFT_StaticStructural_deltat "deltat"

namespace oofem {
class SparseMtrx;

/**
 * Solves a static structural problem.
 * @author Mikael Öhman
 */
class StaticStructural : public StructuralEngngModel, public XfemSolverInterface
{
protected:
    FloatArray internalForces;
    FloatArray *solution;
    FloatArray eNorm;
    SparseMtrx *stiffnessMatrix;

    PrimaryField *field;

    SparseMtrxType sparseMtrxType;

    SparseNonLinearSystemNM *nMethod;

    double deltaT;
    
public:
    StaticStructural(int i, EngngModel * _master = NULL);
    virtual ~StaticStructural();
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void solveYourself();
    virtual void solveYourselfAt(TimeStep *tStep);

    virtual void terminate(TimeStep *tStep);

    virtual void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d);
    
    virtual double giveUnknownComponent(ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);

    virtual void updateDomainLinks();

    virtual int forceEquationNumbering();

    virtual TimeStep *giveNextStep();
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);

    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep);
    virtual fMode giveFormulation() { return TL; }

    virtual void updatePrimaryField(ValueModeType mode, TimeStep *tStep, const FloatArray &vectorToStore);

    // identification
    virtual const char *giveInputRecordName() const { return _IFT_StaticStructural_Name; }
    virtual const char *giveClassName() const { return "StaticStructural"; }

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    
#ifdef __PARALLEL_MODE
    int estimateMaxPackSize(IntArray &commMap, CommunicationBuffer &buff, int packUnpackType);
#endif
};
} // end namespace oofem
#endif // staticstructural_h
