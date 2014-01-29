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

#ifndef pfem_h
#define pfem_h

#include "engngm.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "primaryfield.h"
#include "dofdistributedprimaryfield.h"
//<RESTRICTED_SECTION>
#include "materialinterface.h"
//</RESTRICTED_SECTION>
#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

#include "pfemnumberingschemes.h"

///@name Input fields for PFEM
//@{
#define _IFT_PFEM_Name "pfem"
#define _IFT_CBS_deltat "deltat"
#define _IFT_CBS_mindeltat "mindeltat"
#define _IFT_CBS_cmflag "cmflag"
#define _IFT_PFEM_alphashapecoef "alphashapecoef"
#define _IFT_PFEM_particalRemovalRatio "removalratio"

//@}

namespace oofem {
/**
 * This class represents PFEM method for solving incompressible Navier-Stokes equations
 */
class PFEM : public EngngModel
{
protected:
    /// Numerical method used to solve the problem
    SparseLinearSystemNM *nMethod;

    LinSystSolverType solverType;
    SparseMtrxType sparseMtrxType;

    SparseMtrx *lhs;
    SparseMtrx *avLhs, *pLhs, *vLhs;

    /// Pressure field
    PrimaryField PressureField;
    //DofDistributedPrimaryField PressureField;
    /// Velocity field
    PrimaryField VelocityField;
    //DofDistributedPrimaryField VelocityField;
    FloatArray deltaAuxVelocity;
    FloatArray AuxVelocity;
    FloatArray prescribedTractionPressure;
    FloatArray nodalPrescribedTractionPressureConnectivity;

    /** lumped mass matrix */
    FloatArray mm;
    /** Sparse consistent mass */
    SparseMtrx *mss;

    /// time step and its minimal value
    double deltaT, minDeltaT;

    double alphaShapeCoef;

	double particleRemovalRatio;

    int initFlag;


    int numberOfMomentumEqs, numberOfConservationEqs;
    int numberOfPrescribedMomentumEqs, numberOfPrescribedConservationEqs;

    /// Pressure numbering
    PressureNumberingScheme pns;
    /// Auxiliary Velocity numbering
    AuxVelocityNumberingScheme avns;
    /// Velocity numbering
    VelocityNumberingScheme vns;
    /// Prescribed velocity numbering
    VelocityNumberingScheme prescribedVns;

public:
    PFEM(int i, EngngModel * _master = NULL) :
        EngngModel(i, _master)
        , avLhs(NULL)
        , pLhs(NULL)
        , vLhs(NULL)
        , PressureField(this, 1, FT_Pressure, EID_ConservationEquation, 1)
        , VelocityField(this, 1, FT_Velocity, EID_MomentumBalance, 1)
        , pns()
        , vns(false)
        , prescribedVns(true)

    {
        initFlag = 1;
        lhs = NULL;
        ndomains = 1;
        nMethod = NULL;
        numberOfMomentumEqs = numberOfConservationEqs = numberOfPrescribedMomentumEqs = numberOfPrescribedConservationEqs = 0;
    }
    ~PFEM() { }

    void solveYourselfAt(TimeStep *);
    /**
     * Updates nodal values
     * (calls also this->updateDofUnknownsDictionary for updating dofs unknowns dictionaries
     * if model supports changes of static system). The element internal state update is also forced using
     * updateInternalState service.giv
     */
    virtual void               updateYourself(TimeStep *);

    double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof);

    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);


    TimeStep *giveNextStep();
    TimeStep *giveSolutionStepWhenIcApply();
    NumericalMethod *giveNumericalMethod(MetaStep *);

    //equation numbering using PressureNumberingScheme
    virtual int forceEquationNumbering(int id);
    /**
     * Forces equation renumbering on all domains associated to engng model.
     * All equation numbers in all domains for all dofManagers are invalidated,
     * and new equation numbers are generated starting from 1 on each domain.
     * It will update numberOfEquations variable accordingly.
     * Should be used at startup to force equation numbering and therefore sets numberOfEquations.
     * Must be used if model supports changes of static system to assign  new valid equation numbers
     * to dofManagers.
     */
    virtual int       forceEquationNumbering() { return EngngModel :: forceEquationNumbering(); }

    virtual int requiresUnknownsDictionaryUpdate() { return true; }


    /// Initialization from given input record
    IRResultType initializeFrom(InputRecord *ir);

    // consistency check
    virtual int checkConsistency(); // returns nonzero if o.k.
    // identification
    const char *giveClassName() const { return "PFEM"; }
    classType giveClassID()      const { return PFEMClass; }
    fMode giveFormulation() { return AL; }
    //    fMode giveFormulation() { return TL; }

    /** DOF printing routine. Called by DofManagers to print Dof specific part.
     *  Dof class provides component printing routines, but emodel is responsible
     *  for what will be printed at DOF level.
     *  @param stream output stream
     *  @param iDof dof to be processed
     *  @param atTime solution step
     */
    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);

    virtual int giveNumberOfDomainEquations(int, const UnknownNumberingScheme &num);

    virtual int giveNewEquationNumber(int domain, DofIDItem);
    virtual int giveNewPrescribedEquationNumber(int domain, DofIDItem);

    virtual bool giveBCEnforcementFlag(EquationID eid);
    virtual bool giveBCEnforcementFlag(PrimaryField &field) { return true; }


    virtual int giveUnknownDictHashIndx(ValueModeType mode, TimeStep *stepN); // { return ( int ) mode; }
    virtual void updateDofUnknownsDictionary(DofManager *inode, TimeStep *tStep);

    /// Writes pressures into the dof unknown dictionaries
    void updateDofUnknownsDictionaryPressure(DofManager *inode, TimeStep *tStep);
    /// Writes velocities into the dof unknown dictionaries
    void updateDofUnknownsDictionaryVelocities(DofManager *inode, TimeStep *tStep);
    void resetEquationNumberings();
protected:
    /**
     * Updates nodal values
     * (calls also this->updateDofUnknownsDictionary for updating dofs unknowns dictionaries
     * if model supports changes of static system). The element internal state update is also forced using
     * updateInternalState service.
     */
    void updateInternalState(TimeStep *);
    void applyIC(TimeStep *);
};
} // end namespace oofem
#endif // pfem_h
