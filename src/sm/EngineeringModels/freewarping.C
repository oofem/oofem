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

#include "freewarping.h"
#include "crosssection.h"
#include "Elements/trwarp.h"
#include "nummet.h"
#include "timestep.h"
#include "element.h"
#include "dof.h"
#include "sparsemtrx.h"
#include "verbose.h"
#include "Elements/structuralelement.h"
#include "unknownnumberingscheme.h"
#include "Elements/structuralelementevaluator.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"

//#define THROW_CIOERR(e) throw ContextIOERR(e, __FILE__, __LINE__); // km???

#ifdef __PARALLEL_MODE
 #include "problemcomm.h"
 #include "communicator.h"
#endif

#include <typeinfo>

namespace oofem {
REGISTER_EngngModel(FreeWarping);

FreeWarping :: FreeWarping(int i, EngngModel *_master) : StructuralEngngModel(i, _master), loadVector(), displacementVector()
{
    ndomains = 1;
    initFlag = 1;
    solverType = ST_Direct;
}


FreeWarping :: ~FreeWarping()
{
}


NumericalMethod *FreeWarping :: giveNumericalMethod(MetaStep *mStep)
{
    if ( isParallel() ) {
        if ( ( solverType == ST_Petsc ) || ( solverType == ST_Feti ) ) {
            nMethod.reset( classFactory.createSparseLinSolver(solverType, this->giveDomain(1), this) );
        }
    } else {
        nMethod.reset( classFactory.createSparseLinSolver(solverType, this->giveDomain(1), this) );
    }

    if ( !nMethod ) {
        OOFEM_ERROR("linear solver creation failed for lstype %d", solverType);
    }

    return nMethod.get();
}

IRResultType
FreeWarping :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = StructuralEngngModel :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_lstype);
    solverType = ( LinSystSolverType ) val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    sparseMtrxType = ( SparseMtrxType ) val;




#ifdef __PARALLEL_MODE
    if ( isParallel() ) {
        commBuff = new CommunicatorBuff( this->giveNumberOfProcesses() );
        communicator = new NodeCommunicator(this, commBuff, this->giveRank(),
                                            this->giveNumberOfProcesses());
    }

#endif


    return IRRT_OK;
}


void
FreeWarping :: computeCenterOfGravity()
{
    int noCS = this->giveDomain(1)->giveNumberOfCrossSectionModels(); //number of warping Crosssections
    CG.resize(noCS, 2);
    //	cg.resize(2);
    FloatArray Sx, Sy, moments, domainArea;
    Sx.resize(noCS);
    Sy.resize(noCS);
    domainArea.resize(noCS);
    moments.resize(2);

    for ( int i = 1; i <= this->giveDomain(1)->giveNumberOfElements(); ++i ) {
        int j = this->giveDomain(1)->giveElement(i)->giveCrossSection()->giveNumber();
        Tr_Warp *trwarp = dynamic_cast< Tr_Warp * >( this->giveDomain(1)->giveElement(i) );
        if ( trwarp ) {
            trwarp->computeFirstMomentOfArea(moments);
            Sx.at(j) += moments.at(1);
            Sy.at(j) += moments.at(2);

            domainArea.at(j) += fabs( this->giveDomain(1)->giveElement(i)->computeArea() );
        } else                   {
            OOFEM_ERROR("Error during dynamic_cast");
        }
    }
    for ( int j = 1; j <= noCS; ++j ) {
        double A = domainArea.at(j);
        if ( A != 0 ) {
            CG.at(j, 1) = Sx.at(j) / A;
            CG.at(j, 2) = Sy.at(j) / A;
        } else   {
            OOFEM_ERROR("Zero crosssection area");
        }
    }
    fprintf(outputStream, "\n Center of gravity:");
    for ( int j = 1; j <= noCS; ++j ) {
        fprintf( outputStream, "\n  CrossSection %d  x = %f     y = %f", j, CG.at(j, 1), CG.at(j, 2) );
    }
}

void
FreeWarping :: computeResultAtCenterOfGravity(TimeStep *tStep)
{
    int noCS = this->giveDomain(1)->giveNumberOfCrossSectionModels(); //number of warping Crosssections
    SolutionAtCG.resize(noCS);
    Element *closestElement;
    FloatArray lcoords,  closest, lcg;
    SpatialLocalizer *sp = this->giveDomain(1)->giveSpatialLocalizer();
    sp->init();
    lcoords.resize(2);
    closest.resize(2);
    lcg.resize(2);

    for ( int j = 1; j <= noCS; ++j ) {
        lcg.at(1) = CG.at(j, 1);
        lcg.at(2) = CG.at(j, 2);
        closestElement = sp->giveElementClosestToPoint(lcoords, closest, lcg, 0);

        StructuralElement *sE = dynamic_cast< StructuralElement * >(closestElement);
        FloatArray u, r, b;
        FloatMatrix N;
        sE->computeNmatrixAt(lcoords, N);
        sE->computeVectorOf(VM_Total, tStep, u);
        u.resizeWithValues(3);
        r.beProductOf(N, u);

        SolutionAtCG.at(j) = r.at(1);
    }
}


double FreeWarping :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of equation eq
// This function translates this request to numerical method language
{
    int eq = dof->__giveEquationNumber();
#ifdef DEBUG
    if ( eq == 0 ) {
        OOFEM_ERROR("invalid equation number");
    }
#endif

    if ( tStep != this->giveCurrentStep() ) {
        OOFEM_ERROR("unknown time step encountered");
        return 0.;
    }

    switch ( mode ) {
    case VM_Total:
    case VM_Incremental:
        if ( displacementVector.isNotEmpty() ) {
            return displacementVector.at(eq);
        } else {
            return 0.;
        }

    default:
        OOFEM_ERROR("Unknown is of undefined type for this problem");
    }

    return 0.;
}


TimeStep *FreeWarping :: giveNextStep()
{
    int istep = this->giveNumberOfFirstStep();
    StateCounterType counter = 1;

    if ( currentStep ) {
        istep = currentStep->giveNumber() + 1;
        counter = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = std :: move(currentStep);
    currentStep.reset( new TimeStep(istep, this, 1, ( double ) istep, 0., counter) );
    return currentStep.get();
}


void FreeWarping :: solveYourself()
{
    this->computeCenterOfGravity();

    if ( this->isParallel() ) {
 #ifdef __VERBOSE_PARALLEL
        // force equation numbering before setting up comm maps
      int neq = this->giveNumberOfDomainEquations(1, EModelDefaultEquationNumbering());
        OOFEM_LOG_INFO("[process rank %d] neq is %d\n", this->giveRank(), neq);
 #endif
	this->initializeCommMaps();
    }

    StructuralEngngModel :: solveYourself();
}



void FreeWarping :: solveYourselfAt(TimeStep *tStep)
{
    //
    // creates system of governing eq's and solves them at given time step
    //
    // first assemble problem at current time step

    if ( initFlag ) {
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Assembling stiffness matrix\n");
#endif

        //
        // first step  assemble stiffness Matrix
        //
        stiffnessMatrix.reset( classFactory.createSparseMtrx(sparseMtrxType) );
        if ( !stiffnessMatrix ) {
            OOFEM_ERROR("sparse matrix creation failed");
        }

        stiffnessMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );

        this->assemble( *stiffnessMatrix, tStep, TangentAssembler(TangentStiffness),
                        EModelDefaultEquationNumbering(), this->giveDomain(1) );
        //
        // original stiffnes Matrix is singular (no Dirichlet b.c's exist for free warping problem)
        // thus one diagonal element is made stiffer (for each warping crosssection)
        this->updateStiffnessMatrix(stiffnessMatrix.get());

        initFlag = 0;
    }

#ifdef VERBOSE
    OOFEM_LOG_DEBUG("Assembling load\n");
#endif

    //
    // allocate space for displacementVector
    //
    displacementVector.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) );
    displacementVector.zero();

    //
    // assembling the load vector
    //
    loadVector.resize( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) );
    loadVector.zero();
    this->assembleVector( loadVector, tStep, ExternalForceAssembler(), VM_Total,
                          EModelDefaultEquationNumbering(), this->giveDomain(1) );

    //
    // internal forces (from Dirichlet b.c's, or thermal expansion, etc.)
    //
    // no such forces exist for the free warping problem
    /*
     * FloatArray internalForces( this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() ) );
     * internalForces.zero();
     * this->assembleVector( internalForces, tStep, InternalForceAssembler(), VM_Total,
     *                   EModelDefaultEquationNumbering(), this->giveDomain(1) );
     *
     * loadVector.subtract(internalForces);
     */

    this->updateSharedDofManagers(loadVector, EModelDefaultEquationNumbering(), ReactionExchangeTag);

    //
    // set-up numerical model
    //
    this->giveNumericalMethod( this->giveMetaStep( tStep->giveMetaStepNumber() ) );

    //
    // call numerical model to solve arose problem
    //
#ifdef VERBOSE
    OOFEM_LOG_INFO("\n\nSolving ...\n\n");
#endif
    NM_Status s = nMethod->solve(*stiffnessMatrix, loadVector, displacementVector);
    if ( !( s & NM_Success ) ) {
        OOFEM_ERROR("No success in solving system.");
    }

    //
    // update computed "displacementVector" with the respect to the center of gravity
    //
    this->computeResultAtCenterOfGravity(tStep);
    //
    // results are shifted to be zero at the center of gravity (for each warping crosssection)
    //
    this->updateComputedResults(displacementVector, tStep);


    // update solution state counter
    tStep->incrementStateCounter();
}


void
FreeWarping :: terminate(TimeStep *tStep)
{
    fflush( this->giveOutputStream() );
    this->printReactionForces(tStep, 1); // computed reaction forces have the mening of torsional stiffness
    StructuralEngngModel :: terminate(tStep);
}


void
FreeWarping :: updateDomainLinks()
{
    EngngModel :: updateDomainLinks();
    this->giveNumericalMethod( this->giveCurrentMetaStep() )->setDomain( this->giveDomain(1) );
}


int
FreeWarping :: estimateMaxPackSize(IntArray &commMap, DataStream &buff, int packUnpackType)
{
    int count = 0, pcount = 0;
    Domain *domain = this->giveDomain(1);

    if ( packUnpackType == 0 ) { ///@todo Fix this old ProblemCommMode__NODE_CUT value
        for ( int map: commMap ) {
            for ( Dof *jdof: *domain->giveDofManager( map ) ) {
                if ( jdof->isPrimaryDof() && ( jdof->__giveEquationNumber() ) ) {
                    count++;
                } else {
                    pcount++;
                }
            }
        }

        // --------------------------------------------------------------------------------
        // only pcount is relevant here, since only prescribed components are exchanged !!!!
        // --------------------------------------------------------------------------------

        return ( buff.givePackSizeOfDouble(1) * pcount );
    } else if ( packUnpackType == 1 ) {
        for ( int map: commMap ) {
            count += domain->giveElement( map )->estimatePackSize(buff);
        }

        return count;
    }

    return 0;
}


void
FreeWarping :: updateStiffnessMatrix(SparseMtrx *answer)
{
    // increase diagonal stiffness (coresponding to the 1st node of 1st element ) for each crosssection
    for ( int j = 1; j <= this->giveDomain(1)->giveNumberOfCrossSectionModels(); j++ ) {
        for ( auto &elem : this->giveDomain(1)->giveElements() ) {
            int CSnumber = elem->giveCrossSection()->giveNumber();
            if ( CSnumber == j ) {
                IntArray locationArray;
                EModelDefaultEquationNumbering s;
                elem->giveLocationArray(locationArray, s);
                if ( locationArray.at(1) != 0 ) {
                    int cn = locationArray.at(1);
                    if ( answer->at(cn, cn) != 0 ) {
                        answer->at(cn, cn) *= 2;
                    } else {
                        answer->at(cn, cn) = 1000;
                    }
                }
                break;
            }
        }
    }
}




void
FreeWarping :: updateComputedResults(FloatArray &answer, TimeStep *tStep)
{
    // value of result in the center of gravity is interpolated
    // and substracted from the original solution vector
    // (individualy for each crosssection)


    // set up vector of crosssections numbers
    // through all nodes
    int length = answer.giveSize();
    IntArray CSindicator(length);
    EModelDefaultEquationNumbering s;
    for ( int i = 1; i <= this->giveDomain(1)->giveNumberOfElements(); i++ ) {
        int CSnumber = this->giveDomain(1)->giveElement(i)->giveCrossSection()->giveNumber();
        IntArray locationArray;
        this->giveDomain(1)->giveElement(i)->giveLocationArray(locationArray, s);
        for ( int j = 1; j <= 3; j++ ) {
            int locationNumber = locationArray.at(j);
            if ( locationNumber != 0 ) {
                CSindicator.at(locationNumber) =  CSnumber;
            }
        }
    }

    // substracr answer for corresponding warping crosssection
    for ( int i = 1; i <= length; i++ ) {
        answer.at(i) -= SolutionAtCG.at( CSindicator.at(i) );
    }
}
} // end namespace oofem
