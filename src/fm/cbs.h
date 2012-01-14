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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef cbs_h
#define cbs_h

#include "engngm.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "primaryfield.h"
//<RESTRICTED_SECTION>
#include "materialinterface.h"
//</RESTRICTED_SECTION>

namespace oofem {
/**
 * This class represents CBS algorithm for solving incompressible Navier-Stokes equations
 */
class CBS : public EngngModel
{
protected:
    /// Numerical method used to solve the problem.
    SparseLinearSystemNM *nMethod;

    LinSystSolverType solverType;
    SparseMtrxType sparseMtrxType;

    SparseMtrx *lhs;
    /// Pressure field.
    PrimaryField PressureField;
    /// Velocity field.
    PrimaryField VelocityField;
    FloatArray deltaAuxVelocity;
    FloatArray prescribedTractionPressure;
    FloatArray nodalPrescribedTractionPressureConnectivity;

    /// Lumped mass matrix.
    FloatArray mm;
    /// Sparse consistent mass.
    SparseMtrx *mss;

    /// Time step and its minimal value.
    double deltaT, minDeltaT;
    /// Integration constants.
    double theta [ 2 ];

    int initFlag;
    /// Consistent mass flag.
    int consistentMassFlag;

    int numberOfMomentumEqs, numberOfConservationEqs;
    int numberOfPrescribedMomentumEqs, numberOfPrescribedConservationEqs;

    bool equationScalingFlag;
    /// Length scale.
    double lscale;
    /// Velocity scale.
    double uscale;
    /// Density scale.
    double dscale;
    /// Reynolds number
    double Re;

    //<RESTRICTED_SECTION>
    // material interface representation for multicomponent flows
    MaterialInterface *materialInterface;
    //</RESTRICTED_SECTION>
public:
    CBS(int i, EngngModel *_master = NULL) : EngngModel(i, _master), PressureField(this, 1, FT_Pressure, EID_ConservationEquation, 1),
        VelocityField(this, 1, FT_Velocity, EID_MomentumBalance, 1) {
        initFlag = 1;
        lhs = NULL;
        ndomains = 1;
        nMethod = NULL;
        numberOfMomentumEqs = numberOfConservationEqs = numberOfPrescribedMomentumEqs = numberOfPrescribedConservationEqs = 0;
        consistentMassFlag = 0;
        equationScalingFlag = false;
        lscale = uscale = dscale = 1.0;
        //<RESTRICTED_SECTION>
        materialInterface = NULL;
        //</RESTRICTED_SECTION>
    }
    ~CBS() {
        //<RESTRICTED_SECTION>
        if ( materialInterface ) { delete materialInterface; }
        //</RESTRICTED_SECTION>
    }

    void solveYourselfAt(TimeStep *tStep);

    virtual void updateYourself(TimeStep *tStep);

    double giveUnknownComponent(EquationID eid, ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);
    double giveUnknownComponent(UnknownType ut, ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);
    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    void   updateDomainLinks();

    TimeStep *giveNextStep();
    TimeStep *giveSolutionStepWhenIcApply();
    NumericalMethod *giveNumericalMethod(TimeStep *tStep);

    IRResultType initializeFrom(InputRecord *ir);

    virtual int checkConsistency();
    // identification
    const char *giveClassName() const { return "CBS"; }
    classType giveClassID() const { return CBSClass; }
    fMode giveFormulation() { return TL; }

    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep);

    virtual int giveNumberOfEquations(EquationID eid);
    virtual int giveNumberOfPrescribedEquations(EquationID eid);
    virtual int giveNumberOfDomainEquations(int, EquationID eid);
    virtual int giveNumberOfPrescribedDomainEquations(int, EquationID eid);

    virtual int giveNewEquationNumber(int domain, DofIDItem);
    virtual int giveNewPrescribedEquationNumber(int domain, DofIDItem);

    virtual bool giveEquationScalingFlag() { return equationScalingFlag; }
    virtual double giveVariableScale(VarScaleType varId);

protected:
    /**
     * Updates nodal values
     * (calls also this->updateDofUnknownsDictionary for updating dofs unknowns dictionaries
     * if model supports changes of static system). The element internal state update is also forced using
     * updateInternalState service.
     */
    void updateInternalState(TimeStep *tStep);
    void applyIC(TimeStep *tStep);
    void assembleAlgorithmicPartOfRhs(FloatArray &rhs, EquationID ut, TimeStep *tStep, int nite);
};
} // end namespace oofem
#endif // cbs_h
