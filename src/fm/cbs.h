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

///@name Input fields for CBS
//@{
#define _IFT_CBS_deltat "deltat"
#define _IFT_CBS_mindeltat "mindeltat"
#define _IFT_CBS_cmflag "cmflag"
#define _IFT_CBS_theta1 "theta1"
#define _IFT_CBS_theta2 "theta2"
#define _IFT_CBS_scaleflag "scaleflag"
#define _IFT_CBS_lscale "lscale"
#define _IFT_CBS_uscale "uscale"
#define _IFT_CBS_dscale "dscale"
#define _IFT_CBS_miflag "miflag"
//@}

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
    virtual ~CBS() {
        //<RESTRICTED_SECTION>
        if ( materialInterface ) { delete materialInterface; }
        //</RESTRICTED_SECTION>
    }

    virtual void solveYourselfAt(TimeStep *tStep);

    virtual void updateYourself(TimeStep *tStep);

    virtual double giveUnknownComponent(EquationID eid, ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);
    virtual double giveUnknownComponent(UnknownType ut, ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual void updateDomainLinks();

    virtual TimeStep *giveNextStep();
    virtual TimeStep *giveSolutionStepWhenIcApply();
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int checkConsistency();
    // identification
    virtual const char *giveClassName() const { return "CBS"; }
    virtual classType giveClassID() const { return CBSClass; }
    virtual fMode giveFormulation() { return TL; }

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
