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

// Implementation according paper:
// S.R. Idelsohn, E. Onate, F. Del Pin
// The particle finite element method: a powerful tool to solve
// incompressible flows with free-surface and breaking waves

#include "pfem.h"
#include "nummet.h"
#include "ldltfact.h"
//#include "imlsolver.h"
#include "timestep.h"
#include "metastep.h"
#include "element.h"
#include "dofmanager.h"
#include "elementside.h"
#include "dof.h"
#include "initialcondition.h"

#include "verbose.h"
#include "connectivitytable.h"
#include "cbselement.h"
#include "pfemelement.h"
#include "tr1_2d_pfem.h"
#include "delaunaytriangulator.h"
#include "load.h"
#include "pfemparticle.h"
#include "octreelocalizert.h" //changed from "octreelocalizer.h"
#include "spatiallocalizer.h"
#include "classfactory.h"
//#include "usrdefsub.h"
#include "mathfem.h"
#include "datastream.h"
//<RESTRICTED_SECTION>
#include "leplic.h"
//</RESTRICTED_SECTION>
#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif
#ifdef TIME_REPORT
 #ifndef __MAKEDEPEND
  #include <time.h>
 #endif
//#include "clock.h"
#endif
#include "contextioerr.h"

namespace oofem {
REGISTER_EngngModel(PFEM);

NumericalMethod *PFEM :: giveNumericalMethod(MetaStep *mStep)
{
    if ( nMethod ) {
        return nMethod;
    }

    nMethod = classFactory.createSparseLinSolver(solverType, this->giveDomain(1), this);
    if ( nMethod == NULL ) {
        _error("giveNumericalMethod: linear solver creation failed");
    }

    return nMethod;
}


int
PFEM :: forceEquationNumbering(int id)
{
    resetEquationNumberings();
    // forces equation renumbering for current time step
    // intended mainly for problems with changes of static system
    // during solution
    // OUTPUT:
    // sets this->numberOfEquations and this->numberOfPrescribedEquations and returns this value

    // first initialize default numbering (for velocity unknowns only)

    int i, j, ndofs, nnodes;
    DofManager *inode;
    Domain *domain = this->giveDomain(id);
    TimeStep *currStep = this->giveCurrentStep();
    IntArray loc;
    Dof *jDof;
    DofIDItem type;

    this->domainNeqs.at(id) = 0;
    this->domainPrescribedNeqs.at(id) = 0;


    nnodes = domain->giveNumberOfDofManagers();

    // number first velocities
    for ( i = 1; i <= nnodes; i++ ) {
        inode = domain->giveDofManager(i);
        ndofs = inode->giveNumberOfDofs();
        for ( j = 1; j <= ndofs; j++ ) {
            jDof  =  inode->giveDof(j);
            type  =  jDof->giveDofID();
            if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                jDof->askNewEquationNumber(currStep);
            }
        }
    }

    // invalidate element local copies of location arrays
    //    for ( i = 1; i <= nelem; i++ ) {
    //        domain->giveElement(i)->invalidateLocationArray();
    //    }

    // initialize the separate pressure equation numbering
    pns.init(domain, currStep);
    numberOfConservationEqs =  pns.giveTotalNumberOfEquations();
    numberOfPrescribedConservationEqs = pns.giveTotalNumberOfPrescribedEquations();

    avns.init(domain);

    return domainNeqs.at(id);
}



IRResultType
PFEM :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";            // Required by IR_GIVE_FIELD macro
    IRResultType result;                            // Required by IR_GIVE_FIELD macro

    EngngModel :: initializeFrom(ir);
    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_lstype);
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    sparseMtrxType = ( SparseMtrxType ) val;

    IR_GIVE_FIELD(ir, deltaT, _IFT_CBS_deltat);
    minDeltaT = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, minDeltaT, _IFT_CBS_mindeltat);

    alphaShapeCoef = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaShapeCoef, _IFT_PFEM_alphashapecoef);


    return IRRT_OK;
}



double
PFEM :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of dof
{
    //#ifdef DEBUG
    //  if (dof->giveDofID() != P_f && mode != VM_Intermediate) // uses other equation numbering
    //  {
    //  if ( dof->__giveEquationNumber() == 0 ) {
    //        _error("giveUnknownComponent: invalid equation number");
    //    }
    //  }
    //#endif
    if ( mode == VM_Intermediate ) {
        if ( AuxVelocity.isNotEmpty() ) {
            return AuxVelocity.at( avns.giveDofEquationNumber(dof) );
        } else {
            return 0.;
        }
    } else {
        if ( this->requiresUnknownsDictionaryUpdate() ) {
            int hash = this->giveUnknownDictHashIndx(mode, tStep);
            if ( dof->giveUnknowns()->includes(hash) ) {
                return dof->giveUnknowns()->at(hash);
            } else {
                OOFEM_ERROR2( "giveUnknown:  Dof unknowns dictionary does not contain unknown of value mode (%s)", __ValueModeTypeToString(mode) );
                return 0.; // to make compiler happy
            }
        } else {
            if ( dof->giveDofID() == P_f ) {                     // pressures
                return PressureField.giveUnknownValue(dof, mode, tStep);
            } else {                     // velocities
                return VelocityField.giveUnknownValue(dof, mode, tStep);
            }
        }
    }
}


TimeStep *
PFEM :: giveSolutionStepWhenIcApply()
{
    if ( stepWhenIcApply == NULL ) {
        stepWhenIcApply = new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0, 0.0, deltaT, 0);
    }

    return stepWhenIcApply;
}

TimeStep *
PFEM :: giveNextStep()
{
    int istep = this->giveNumberOfFirstStep();
    int i, nelem, nnodes;
    double totalTime = 0;
    StateCounterType counter = 1;
    delete previousStep;
    DofManager *dman = NULL;
    Domain *domain = this->giveDomain(1);
    Dof *jDof = NULL;
    DofIDItem type;


    if ( currentStep == NULL ) {
        // first step -> generate initial step
        currentStep = new TimeStep( *giveSolutionStepWhenIcApply() );
    } else {
        istep =  currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
        domain->clearElements();
    }

    DelaunayTriangulator myMesher(domain, alphaShapeCoef);
    myMesher.generateMesh();

    nnodes = domain->giveNumberOfDofManagers();
    for ( i = 1; i <= nnodes; i++ ) {
        PFEMParticle* particle = dynamic_cast< PFEMParticle * >( domain->giveDofManager(i) );
		particle->setFree();
		//particle->storeCoordinatesTimeStepBegin();
    }

    nelem = domain->giveNumberOfElements();
    if ( nelem ) {
        for ( i = 1; i <= nelem; i++ ) {
            Element *element = domain->giveElement(i);
            element->checkConsistency();

            for ( int j = 1; j <= element->giveNumberOfDofManagers(); j++ ) {
                dynamic_cast< PFEMParticle * >( element->giveDofManager(j) )->setFree(false);
            }
        }
    } else {
        VERBOSE_PRINTS("Mesh generation failed", "0 elements created");
    }


    nnodes = domain->giveNumberOfDofManagers();
    for ( i = 1; i <= nnodes; i++ ) {
        dman = domain->giveDofManager(i);
        for ( int j = 1; j <= dman->giveNumberOfDofs(); j++ ) {
            jDof  =  dman->giveDof(j);
            type  =  jDof->giveDofID();
            if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                if ( jDof->giveBcId() ) {
                    dynamic_cast< PFEMParticle * >(dman)->setFree(false);
                    break;
                }
            }
        }
    }

    previousStep = currentStep;

    // FORCE EQUATION NUMBERING
    this->forceEquationNumbering();

    double ndt = ( ( PFEMElement * ) domain->giveElement(1) )->computeCriticalTimeStep(previousStep);
    // check for critical time step
	TR1_2D_PFEM* ielem = NULL;
    for ( i = 2; i <= domain->giveNumberOfElements(); i++ ) {
		ielem = dynamic_cast<TR1_2D_PFEM*>(domain->giveElement(i));
		double idt = ielem->computeCriticalTimeStep(previousStep);
		if ( idt < ndt)
		{
			ndt = idt;
			printf("Reducing time step due to element #%i \n", i);
		}
		
    }

    // UP error computeCriticalTimeStep

    //COdt = max(dt*0.8, minDeltaT);
    ndt = min (ndt, deltaT);
	ndt = max (ndt, minDeltaT);


    if ( currentStep != NULL ) {
        totalTime = currentStep->giveTargetTime() + ndt;
    }

    currentStep = new TimeStep(istep, this, 1, totalTime, ndt, counter);
    // time and dt variables are set eq to 0 for staics - has no meaning

    OOFEM_LOG_INFO( "SolutionStep %d : t = %e, dt = %e\n", istep, totalTime * this->giveVariableScale(VST_Time), ndt * this->giveVariableScale(VST_Time) );

    return currentStep;
}

void
PFEM :: solveYourselfAt(TimeStep *tStep)
{
    int auxmomneq = this->giveNumberOfDomainEquations(1, avns);
    int momneq =  this->giveNumberOfDomainEquations(1, vns);
    int presneq =  this->giveNumberOfDomainEquations(1, pns);

    double deltaT = tStep->giveTimeIncrement();

    FloatArray rhs(auxmomneq);


    double d_pnorm = 1.0;
    double d_vnorm = 1.0;

    //  if ( initFlag ) {
    deltaAuxVelocity.resize(auxmomneq);
    AuxVelocity.resize(auxmomneq);

    if ( avLhs ) {
        delete avLhs;
    }

    avLhs = classFactory.createSparseMtrx(sparseMtrxType);
    if ( avLhs == NULL ) {
        _error("solveYourselfAt: sparse matrix creation failed");
    }

    avLhs->buildInternalStructure(this, 1, EID_AuxMomentumBalance, avns);

    // ??? stepWhenIcApply
    this->assemble( avLhs, tStep, EID_AuxMomentumBalance, LumpedMassMatrix, avns, this->giveDomain(1) );

//	double ulala = avLhs->at(444,444);
    if ( pLhs ) {
        delete pLhs;
    }

    pLhs = classFactory.createSparseMtrx(sparseMtrxType);
    if ( pLhs == NULL ) {
        _error("solveYourselfAt: sparse matrix creation failed");
    }

    pLhs->buildInternalStructure(this, 1, EID_ConservationEquation, pns);

    this->assemble( pLhs, tStep, EID_ConservationEquation, PressureLaplacianMatrix, pns, this->giveDomain(1) );

    pLhs->times(deltaT);

    //////////////////////////////////////////////////////////////////////////

    if ( vLhs ) {
        delete vLhs;
    }

    vLhs = classFactory.createSparseMtrx(sparseMtrxType);
    if ( vLhs == NULL ) {
        _error("solveYourselfAt: sparse matrix creation failed");
    }

    vLhs->buildInternalStructure(this, 1, EID_MomentumBalance, vns);
    this->assemble( vLhs, tStep, EID_MomentumBalance, LumpedMassMatrix, vns, this->giveDomain(1) );

    initFlag = 0;

    if ( tStep->giveNumber() == giveNumberOfFirstStep() ) {
        TimeStep *stepWhenIcApply = tStep->givePreviousStep();
        this->applyIC(stepWhenIcApply);
    }

    VelocityField.advanceSolution(tStep);
    PressureField.advanceSolution(tStep);

    FloatArray *velocityVector = VelocityField.giveSolutionVector(tStep);
    FloatArray *pressureVector = PressureField.giveSolutionVector(tStep);

    FloatArray velocityVectorThisStep(*velocityVector);
    FloatArray velocityVectorLastStep(*velocityVector);

    FloatArray pressureVectorThisStep(*pressureVector);
    FloatArray pressureVectorLastStep(*pressureVector);

    if ( tStep->isTheFirstStep() ) {
        for ( int i = 1; i <= this->giveDomain(1)->giveNumberOfDofManagers(); i++ ) {
            this->updateDofUnknownsDictionary(this->giveDomain(1)->giveDofManager(i), tStep);
        }
    }

    int iteration = 0;
    do {
        iteration++;

        velocityVectorLastStep = velocityVectorThisStep;
        pressureVectorLastStep = pressureVectorThisStep;

        /************************** STEP 1 - calculates auxiliary velocities *************************/

        rhs.resize(auxmomneq);
        rhs.zero();

        bool tStepNumberTricked = false;
        //in order not to get bc-value in getUnknown
        //will be reset

        if ( tStep->givePreviousStep()->isIcApply() ) {
            tStep->givePreviousStep()->setNumber(-1);
            tStepNumberTricked = true;
        }

        if ( iteration > 1 ) {
            this->assembleVectorFromElements( rhs, tStep, EID_AuxMomentumBalance, LaplaceVelocityVector, VM_Total, avns, this->giveDomain(1) );
        }


        rhs.times(-1.0);

        this->assembleVectorFromElements( rhs, tStep, EID_AuxMomentumBalance, LoadVector, VM_Total, avns, this->giveDomain(1) );

        rhs.times(deltaT);

        this->assembleVectorFromElements( rhs, tStep->givePreviousStep(), EID_AuxMomentumBalance, MassVelocityVector, VM_Total, avns, this->giveDomain(1) );

        if ( tStepNumberTricked ) {
            tStep->givePreviousStep()->setNumber(0);
        }

        this->giveNumericalMethod(this->giveMetaStep( tStep->giveMetaStepNumber()));
        nMethod->solve(avLhs, & rhs, & AuxVelocity);

        /************************* STEP 2 - calculates pressure ************************************/

        rhs.resize(presneq);
        rhs.zero();

        this->assembleVectorFromElements( rhs, tStep, EID_ConservationEquation, PrescribedRhsVector, VM_Total, pns, this->giveDomain(1) );

        rhs.times(-1.0);

        this->assembleVectorFromElements( rhs, tStep, EID_ConservationEquation, DivergenceVelocityVector, VM_Total, pns, this->giveDomain(1) );

        this->giveNumericalMethod(this->giveMetaStep( tStep->giveMetaStepNumber() ));
        pressureVector->resize(presneq);
        nMethod->solve(pLhs, & rhs, pressureVector);

        for ( int i = 1; i <= this->giveDomain(1)->giveNumberOfDofManagers(); i++ ) {
            this->updateDofUnknownsDictionaryPressure(this->giveDomain(1)->giveDofManager(i), tStep);
        }

        /**************************** STEP 3 - velocity correction step ********************************/



        rhs.resize(momneq);
        rhs.zero();
        velocityVector->resize(momneq);

        this->assembleVectorFromElements( rhs, tStep, EID_MomentumBalance, PressureGradientVector, VM_Total, vns, this->giveDomain(1) );

        rhs.times(-1.0);

        rhs.times(deltaT);

        this->assembleVectorFromElements( rhs, tStep, EID_MomentumBalance, MassAuxVelocityVector, VM_Total, vns, this->giveDomain(1) );

        // deactivated, problem of prescribed pressure solved by improvement in alpha shape free surface definition
        // this->assembleVectorFromElements(rhs, tStep, EID_MomentumBalance, PrescribedPressureRhsVector, VM_Total, vns, this->giveDomain(1) );

        this->giveNumericalMethod(this->giveMetaStep( tStep->giveMetaStepNumber() ));

        nMethod->solve(vLhs, & rhs, velocityVector);

        for (int i = 1; i <= this->giveDomain(1)->giveNumberOfDofManagers(); i++)
        {
			PFEMParticle* particle = dynamic_cast<PFEMParticle*>(this->giveDomain(1)->giveDofManager(i));
			this->updateDofUnknownsDictionaryVelocities(particle, tStep);

			//particle->resetNodalCoordinates();
			//particle->updateNodalCoordinates(tStep);
        }

        velocityVectorThisStep = * velocityVector;
        pressureVectorThisStep = * pressureVector;

        if ( iteration == 1 ) {
            d_vnorm = velocityVectorThisStep.computeNorm();
            d_pnorm = pressureVectorThisStep.computeNorm();
        } else {
            FloatArray diffVelocity = velocityVectorThisStep;
            diffVelocity.subtract(velocityVectorLastStep);
			/*for (int i = 1; i <= diffVelocity.giveSize(); i++)
			{
				if(fabs(diffVelocity.at(i)) > 1.e-2)
				{
					
					Dof* dof = vns.giveDofToEquationNumber(this->giveDomain(1), i);
					if (dof)
					{
						DofManager* dofman = dof->giveDofManager();
						printf("Too large difference at time step %i, iteration %i \n", tStep->giveNumber(), iteration);
						printf("Velocity thisStep %f - lastStep %f = difference %f \n", velocityVectorThisStep.at(i), velocityVectorLastStep.at(i), diffVelocity.at(i));
						printf("DofManager #%i at (%f, %f), DofID #%s \n", dofman->giveNumber(), dofman->giveCoordinate(1), dofman->giveCoordinate(2), __DofIDItemToString(dof->giveDofID()));
						char str[80];
						scanf ("%s",str);
					}
					else
						OOFEM_ERROR("Finding Dof to equation Number failed");
				}
			}*/

			d_vnorm = 0.0;
			for (int i = 1; i <= diffVelocity.giveSize(); i++)
			{
				d_vnorm = max(d_vnorm, abs(diffVelocity.at(i)/velocityVectorLastStep.at(i)));
			}
            //d_vnorm = diffVelocity.computeNorm() / velocityVectorLastStep.computeNorm();

            FloatArray diffPressure = pressureVectorThisStep;
            diffPressure.subtract(pressureVectorLastStep);
			/*double normDiffPressure = diffPressure.computeNorm();
			if (normDiffPressure < 1.e-9)
				d_pnorm = 0.0;
			else
				d_pnorm = normDiffPressure / pressureVectorLastStep.computeNorm();*/
			d_pnorm = 0.0;
			for (int i = 1; i <= diffPressure.giveSize(); i++)
			{
				d_pnorm = max(d_pnorm, abs(diffPressure.at(i)/pressureVectorLastStep.at(i)));
			}
        }
    } while ( ( d_vnorm > 1.e-8 || d_pnorm > 1.e-8 ) && iteration < 50 );

    if ( iteration > 49 ) {
        OOFEM_ERROR("Maximal iteration count exceded");
    } else {
        printf("\n %i iterations performed.\n", iteration);
    }

    Domain *d = this->giveDomain(1);
    for ( int i = 1; i <= this->giveDomain(1)->giveNumberOfDofManagers(); i++ ) {
        PFEMParticle *particle = dynamic_cast< PFEMParticle * >( d->giveDofManager(i) );
		
		//resetting coordinates of all particles to avoid doing it again in updateYourself()
        //particle->resetNodalCoordinates();

		if ( particle->isFree() ) {
            for ( int j = 1; j <= particle->giveNumberOfDofs(); j++ ) {
                Dof *dof = particle->giveDof(j);
                DofIDItem type  =  dof->giveDofID();
                if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                    int eqnum = dof->giveEquationNumber(vns);
                    if ( eqnum ) {
                        double previousValue = dof->giveUnknown( VM_Total, tStep->givePreviousStep() );
                        Load *load = d->giveLoad(3);
                        FloatArray gVector;
                        load->computeComponentArrayAt(gVector, tStep, VM_Total);
                        previousValue += gVector.at(j) * deltaT;
                        velocityVector->at(eqnum) = previousValue;
                    }
                }
            }
        }
    }


    tStep->incrementStateCounter();
	this->updateYourself( this->giveCurrentStep() );
}


void
PFEM :: updateYourself(TimeStep *stepN)
{
    this->updateInternalState(stepN);
    EngngModel :: updateYourself(stepN);
}



void
PFEM :: updateInternalState(TimeStep *stepN)
{
    int j, nnodes;
    Domain *domain;

    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);

        nnodes = domain->giveNumberOfDofManagers();
        if ( requiresUnknownsDictionaryUpdate() ) {
            for ( j = 1; j <= nnodes; j++ ) {
                this->updateDofUnknownsDictionary(domain->giveDofManager(j), stepN);
            }
        }

        int nelem = domain->giveNumberOfElements();
        for ( j = 1; j <= nelem; j++ ) {
            domain->giveElement(j)->updateInternalState(stepN);
        }
    }
}

int
PFEM :: giveUnknownDictHashIndx(ValueModeType mode, TimeStep *stepN)
{
    if ( ( stepN == this->giveCurrentStep() ) || ( stepN == this->givePreviousStep() ) ) {
        return ( stepN->giveNumber() % 2 ) * 100 + mode;
    } else {
        _error("giveUnknownDictHashIndx: unsupported solution step");
    }

    return 0;
}

void
PFEM :: updateDofUnknownsDictionary(DofManager *inode, TimeStep *tStep)
{
    int ndofs = inode->giveNumberOfDofs();

    for ( int i = 1; i <= ndofs; i++ ) {
        Dof *iDof = inode->giveDof(i);
        int eqNum = 0;
        double val = 0;
        if ( iDof->hasBc(tStep) ) {                 // boundary condition
            val = iDof->giveBcValue(VM_Total, tStep);
        } else {
            if ( iDof->giveDofID() == P_f ) {
                eqNum = pns.giveDofEquationNumber(iDof);
                FloatArray *vect = PressureField.giveSolutionVector(tStep);
                if ( vect->giveSize() > 0 ) {                     // in the first step -> zero will be set
                    val = vect->at(eqNum);
                }
            } else {                 // velocities
                eqNum = iDof->__giveEquationNumber();
                FloatArray *vect = VelocityField.giveSolutionVector(tStep);
                if ( vect->giveSize() > 0 ) {
                    val = vect->at(eqNum);
                }
            }
        }

        iDof->updateUnknownsDictionary(tStep, VM_Total, val);
    }
}

void
PFEM :: updateDofUnknownsDictionaryPressure(DofManager *inode, TimeStep *tStep)
{
    Dof *iDof = inode->giveDofWithID(P_f);
    double val = 0;
    if ( iDof->hasBc(tStep) ) {
        val = iDof->giveBcValue(VM_Total, tStep);
    } else {
        int eqNum = pns.giveDofEquationNumber(iDof);
        FloatArray *vect = PressureField.giveSolutionVector(tStep);
        val = vect->at(eqNum);
        //iDof->updateUnknownsDictionary(tStep, VM_Total, val);
    }

    iDof->updateUnknownsDictionary(tStep, VM_Total, val);
}

void
PFEM :: updateDofUnknownsDictionaryVelocities(DofManager *inode, TimeStep *tStep)
{
    int ndofs = inode->giveNumberOfDofs();

    for ( int i = 1; i <= ndofs; i++ ) {
        Dof *iDof = inode->giveDof(i);
        int eqNum = 0;
        double val = 0;
        if ( iDof->hasBc(tStep) ) { // boundary condition
            val = iDof->giveBcValue(VM_Total, tStep);
        } else {
            DofIDItem type = iDof->giveDofID();
            if ( type == V_u || type == V_v || type == V_w ) {
                eqNum = iDof->__giveEquationNumber();
                FloatArray *vect = VelocityField.giveSolutionVector(tStep);
                if ( vect->giveSize() > 0 ) {
                    val = vect->at(eqNum);
                }
            }
        }

        iDof->updateUnknownsDictionary(tStep, VM_Total, val);
    }
}


// NOT ACTIVE
contextIOResultType
PFEM :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves state variable - displacement vector
//
{
    contextIOResultType iores;
    int closeFlag = 0;
    FILE *file = NULL;

    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, this->giveCurrentStep()->giveNumber(),
                                    this->giveCurrentStep()->giveVersion(), contextMode_write) ) {
            THROW_CIOERR(CIO_IOERR);                             // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = EngngModel :: saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = PressureField.saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = VelocityField.saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = prescribedTractionPressure.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }                                                                   // ensure consistent records

    return CIO_OK;
}


// NOT ACTIVE
contextIOResultType
PFEM :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restore state variable - displacement vector
//
{
    contextIOResultType iores;
    int closeFlag = 0;
    int istep, iversion;
    FILE *file = NULL;

    this->resolveCorrespondingStepNumber(istep, iversion, obj);

    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, istep, iversion, contextMode_read) ) {
            THROW_CIOERR(CIO_IOERR);                     // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = EngngModel :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = PressureField.restoreContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = VelocityField.restoreContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = prescribedTractionPressure.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }                                                                   // ensure consistent records

    return CIO_OK;
}


int
PFEM :: checkConsistency()
{
    // check internal consistency
    // if success returns nonzero
    int i, nelem;
    Element *ePtr;
    PFEMElement *sePtr;
    Domain *domain = this->giveDomain(1);

    nelem = domain->giveNumberOfElements();
    // check for proper element type

    for ( i = 1; i <= nelem; i++ ) {
        ePtr = domain->giveElement(i);
        sePtr = dynamic_cast< PFEMElement * >(ePtr);
        if ( sePtr == NULL ) {
            _warning2("Element %d has no PFEM base", i);
            return 0;
        }
    }

    EngngModel :: checkConsistency();

    return 1;
}


void
PFEM :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
{
    DofIDItem type = iDof->giveDofID();
    if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
        iDof->printSingleOutputAt(stream, atTime, 'x', VM_Intermediate);
        iDof->printSingleOutputAt(stream, atTime, 'v', VM_Total);
    } else if ( ( type == P_f ) ) {
        iDof->printSingleOutputAt(stream, atTime, 'p', VM_Total);
    } else {
        _error("printDofOutputAt: unsupported dof type");
    }
}

void
PFEM :: applyIC(TimeStep *stepWhenIcApply)
{
    Domain *domain = this->giveDomain(1);
    int mbneq =  this->giveNumberOfDomainEquations(1, vns);
    int pdneq =  this->giveNumberOfDomainEquations(1, pns);
    FloatArray *velocityVector, *pressureVector;

#ifdef VERBOSE
    OOFEM_LOG_INFO("Applying initial conditions\n");
#endif
    int nDofs, j, k, jj;
    int nman  = domain->giveNumberOfDofManagers();
    DofManager *node;
    Dof *iDof;
    DofIDItem type;

    VelocityField.advanceSolution(stepWhenIcApply);
    velocityVector = VelocityField.giveSolutionVector(stepWhenIcApply);
    velocityVector->resize(mbneq);
    velocityVector->zero();

    PressureField.advanceSolution(stepWhenIcApply);
    pressureVector = PressureField.giveSolutionVector(stepWhenIcApply);
    pressureVector->resize(pdneq);
    pressureVector->zero();


    for ( j = 1; j <= nman; j++ ) {
        node = domain->giveDofManager(j);
        nDofs = node->giveNumberOfDofs();

        for ( k = 1; k <= nDofs; k++ ) {
            // ask for initial values obtained from
            // bc (boundary conditions) and ic (initial conditions)
            iDof  =  node->giveDof(k);
            if ( !iDof->isPrimaryDof() ) {
                continue;
            }

            jj = iDof->__giveEquationNumber();
            type = iDof->giveDofID();

            if ( jj ) {
                if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                    velocityVector->at(jj) = iDof->giveUnknown(VM_Total, stepWhenIcApply);
                } else {
                    pressureVector->at(jj) = iDof->giveUnknown(VM_Total, stepWhenIcApply);
                }
            }
        }
    }

    // update element state according to given ic
    int nelem = domain->giveNumberOfElements();
    PFEMElement *element;

    for ( j = 1; j <= nelem; j++ ) {
        element = ( PFEMElement * ) domain->giveElement(j);
        element->updateInternalState(stepWhenIcApply);
        element->updateYourself(stepWhenIcApply);
    }
}

int
PFEM :: giveNewEquationNumber(int domain, DofIDItem id)
{
    if ( ( id == V_u ) || ( id == V_v ) || ( id == V_w ) ) {
        return this->vns.askNewEquationNumber();
    } else {
        _error("giveNewEquationNumber:: Unknown DofIDItem");
    }

    return 0;
}

int
PFEM :: giveNewPrescribedEquationNumber(int domain, DofIDItem id)
{
    if ( ( id == V_u ) || ( id == V_v ) || ( id == V_w ) ) {
        return prescribedVns.askNewEquationNumber();
    } else {
        _error("giveNewPrescribedEquationNumber:: Unknown DofIDItem");
    }

    return 0;
}

int PFEM :: giveNumberOfDomainEquations(int d, const UnknownNumberingScheme &numberingScheme) {            //
    // returns number of equations of current problem
    // this method is implemented here, because some method may add some
    // conditions in to system and this may results into increased number of
    // equations.
    //

    if ( !equationNumberingCompleted ) {
        this->forceEquationNumbering();
    }

    return numberingScheme.giveRequiredNumberOfDomainEquation();
}

// NOT ACTIVE called in masterdof - outcommented at the moment
bool
PFEM ::  giveBCEnforcementFlag(EquationID eid)
{
    if ( eid == EID_AuxMomentumBalance ) {
        return false;
    } else {
        return true;
    }
}

void
PFEM :: resetEquationNumberings()
{
    vns.reset();
    prescribedVns.reset();
    pns.reset();
    avns.reset();
}
} // end namespace oofem
