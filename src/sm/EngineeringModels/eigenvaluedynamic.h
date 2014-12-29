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

#ifndef eigenvaluedynamic_h
#define eigenvaluedynamic_h

#include "engngm.h"
#include "sparsegeneigenvalsystemnm.h"
#include "sparsemtrx.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "geneigvalsolvertype.h"

///@name Input fields for EigenValueDynamic
//@{
#define _IFT_EigenValueDynamic_Name "eigenvaluedynamic"
#define _IFT_EigenValueDynamic_nroot "nroot"
#define _IFT_EigenValueDynamic_rtolv "rtolv"
#define _IFT_EigenValueDynamic_stype "stype"
//@}

namespace oofem {
/**
 * This class implements way for examining eigenvalues and eigenvectors in
 * dynamic problems.
 *
 * Solution of this problem is base on equation in the form of: @f$ K\cdot y=w M\cdot y @f$
 * Currently eigenvalue problem is solved using subspace iteration.
 * Tasks:
 * - Assembling the governing equation in the form @f$ K\cdot y=wM\cdot y@f$.
 * - Creating Numerical method for @f$ K\cdot y=wM\cdot y@f$.
 * - Interfacing Numerical method to Elements.
 */
class EigenValueDynamic : public EngngModel
{
private:
    std :: unique_ptr< SparseMtrx > stiffnessMatrix;
    std :: unique_ptr< SparseMtrx > massMatrix;
    SparseMtrxType sparseMtrxType;
    FloatMatrix eigVec;
    FloatArray eigVal;
    int numberOfRequiredEigenValues;
    int activeVector;
    int restoreFlag;
    /// Relative tolerance.
    double rtolv;
    /// Numerical method used to solve the problem.
    std :: unique_ptr< SparseGeneralEigenValueSystemNM > nMethod;
    GenEigvalSolverType solverType;

public:
    EigenValueDynamic(int i, EngngModel * _master = NULL) : EngngModel(i, _master)
    {
        numberOfSteps = 1;
        ndomains = 1;
    }
    virtual ~EigenValueDynamic() { }

    virtual void solveYourselfAt(TimeStep *tStep);
    virtual void terminate(TimeStep *tStep);
    virtual void updateYourself(TimeStep *tStep);

    virtual double giveUnknownComponent(ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual TimeStep *giveNextStep();
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);
    virtual void setActiveVector(int i) { activeVector = i; }
    virtual int resolveCorrespondingEigenStepNumber(void *obj);

    virtual double giveEigenValue(int eigNum) { return eigVal.at(eigNum); }

    // identification
    virtual const char *giveClassName() const { return "EigenValueDynamic"; }
};
} // end namespace oofem
#endif // eigenvaluedynamic_h
