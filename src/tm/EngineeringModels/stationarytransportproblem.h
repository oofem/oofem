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

#ifndef stationarytransportproblem_h
#define stationarytransportproblem_h

#include "engngm.h"
#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "sparsemtrx.h"
#include "primaryfield.h"

///@name Input fields for StationaryTransportProblem
//@{
#define _IFT_StationaryTransportProblem_Name "stationaryproblem"
#define _IFT_StationaryTransportProblem_exportfields "exportfields"
#define _IFT_StationaryTransportProblem_keepTangent "keeptangent"
//@}

namespace oofem {

/**
 * This class represents stationary transport problem.
 * @author Mikael Öhman (among others)
 */
class StationaryTransportProblem : public EngngModel
{
protected:
    SparseMtrxType sparseMtrxType = SMT_Skyline;
    /// This field stores solution vector. For fixed size of problem, the PrimaryField is used, for growing/decreasing size, DofDistributedPrimaryField applies.
    std :: unique_ptr< PrimaryField > UnknownsField;

    std :: unique_ptr< SparseMtrx > conductivityMatrix;
    FloatArray internalForces;
    FloatArray eNorm;

    /// Numerical method used to solve the problem
    std::unique_ptr<SparseNonLinearSystemNM> nMethod;

    bool keepTangent = false;

public:
    StationaryTransportProblem(int i, EngngModel * _master);

    void solveYourselfAt(TimeStep *tStep) override;
    void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d) override;
    void updateSolution(FloatArray &solutionVector, TimeStep *tStep, Domain *d) override;
    void updateInternalRHS(FloatArray &answer, TimeStep *tStep, Domain *d, FloatArray *eNorm) override;
    void updateMatrix(SparseMtrx &mat, TimeStep *tStep, Domain *d) override;
    double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof) override;
    FieldPtr giveField (FieldType key, TimeStep *) override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    void updateDomainLinks() override;

    TimeStep *giveNextStep() override;
    NumericalMethod *giveNumericalMethod(MetaStep *mStep) override;

    void initializeFrom(InputRecord &ir) override;

    int checkConsistency() override;

    // identification
    const char *giveInputRecordName() const { return _IFT_StationaryTransportProblem_Name; }
    const char *giveClassName() const override { return "StationaryTransportProblem"; }
    fMode giveFormulation() override { return TL; }
};
} // end namespace oofem
#endif // stationarytransportproblem_h
