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

#ifndef staggeredsolver_h
#define staggeredsolver_h

#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "convergedreason.h"
#include "sparsemtrx.h"
#include "floatarray.h"
#include "linesearch.h"
#include "nrsolver.h"
#include "unknownnumberingscheme.h"

#include <set>
#include <vector>


///@name Input fields for StaggeredSolver
//@{
#define _IFT_StaggeredSolver_Name "staggeredsolver"
#define _IFT_StaggeredSolver_DofIdList "dofidlist"
#define _IFT_StaggeredSolver_DofIdListPositions "idpos"
//@}

namespace oofem {
class Domain;
class EngngModel;


class CustomEquationNumbering : public UnknownNumberingScheme
{
protected:
    bool prescribed;
    int numEqs;
    int numPresEqs;

public:
    CustomEquationNumbering();

    IntArray dofIdArray; // should be private
    bool isDefault() const override { return true; }
    void setDofIdArray(IntArray array) { this->dofIdArray = std::move(array); }

    int giveDofEquationNumber(Dof *dof) const override;

    int giveRequiredNumberOfDomainEquation() const override { return numEqs; }

    int giveNewEquationNumber() { return ++numEqs; }
    int giveNewPrescribedEquationNumber() { return ++numPresEqs; }
    int giveNumEquations() { return this->numEqs; }
    int giveNumPresEquations() { return this->numPresEqs; }
};

/**
 * Support struct to handle all the split up variables used during the solving step.
 */
class DofGrouping
{
public:
    std :: vector< std :: unique_ptr< SparseMtrx > > stiffnessMatrixList;
    std :: vector< FloatArray > fIntList;
    std :: vector< FloatArray > fExtList;

    std :: vector< IntArray > locArrayList;
    std :: vector< FloatArray > X;
    std :: vector< FloatArray > dX;
    std :: vector< FloatArray > ddX;

    DofGrouping(const std :: vector< CustomEquationNumbering > &numberings, Domain * m);
    void giveTotalLocationArray(IntArray &locationArray, const UnknownNumberingScheme &s, Domain *d);
};

/**
 * The staggered solver will perform Newton iterations on subsets of DofIDs, in a staggered manner.
 * @author Jim Brouzoulis
 */
class OOFEM_EXPORT StaggeredSolver : public NRSolver
{
private:
    IntArray totalIdList;
    IntArray idPos;
    std :: vector< CustomEquationNumbering > UnknownNumberingSchemeList;

    bool checkConvergenceDofIdArray(FloatArray &RT, FloatArray &F, FloatArray &rhs, FloatArray &ddX, FloatArray &X,
                          double RRT, const FloatArray &internalForcesEBENorm, int nite, bool &errorOutOfRange, TimeStep *tStep, IntArray &dofIdArray);

public:
    StaggeredSolver(Domain * d, EngngModel * m);
    virtual ~StaggeredSolver() {}

    // Overloaded methods:
    ConvergedReason solve(SparseMtrx &k, FloatArray &R, FloatArray *R0,
                    FloatArray &X, FloatArray &dX, FloatArray &F,
                    const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                    int &nite, TimeStep *) override;

    void initializeFrom(InputRecord &ir) override;

    const char *giveClassName() const override { return "StaggeredSolver"; }
    const char *giveInputRecordName() const override { return _IFT_StaggeredSolver_Name; }
};
} // end namespace oofem
#endif // staggeredsolver_h
